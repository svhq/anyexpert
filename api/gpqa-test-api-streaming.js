const express = require('express');
const router = express.Router();
const workflowEngine = require('../src/workflow-engine');
const OpenRouterService = require('../src/openrouter-service');
const mathAgent = require('../src/math-agent');
const winston = require('winston');
const { v4: uuidv4 } = require('uuid');
const { EventEmitter } = require('events');

// Configure logger for API
const apiLogger = winston.createLogger({
    level: 'info',
    format: winston.format.combine(
        winston.format.timestamp(),
        winston.format.json()
    ),
    transports: [
        new winston.transports.Console({
            format: winston.format.combine(
                winston.format.colorize(),
                winston.format.simple()
            )
        })
    ]
});

// Initialize services
const openRouterService = new OpenRouterService();

// Event emitter for streaming updates
const streamEmitter = new EventEmitter();

// Streaming endpoint for real-time updates
router.get('/test-gpqa-stream/:requestId', (req, res) => {
    const { requestId } = req.params;
    
    // Set headers for SSE
    res.writeHead(200, {
        'Content-Type': 'text/event-stream',
        'Cache-Control': 'no-cache',
        'Connection': 'keep-alive',
        'Access-Control-Allow-Origin': '*'
    });
    
    // Send initial connection message
    res.write(`data: ${JSON.stringify({ type: 'connected', requestId })}\n\n`);
    
    // Listen for events for this request
    const sendUpdate = (data) => {
        res.write(`data: ${JSON.stringify(data)}\n\n`);
    };
    
    streamEmitter.on(requestId, sendUpdate);
    
    // Clean up on disconnect
    req.on('close', () => {
        streamEmitter.removeListener(requestId, sendUpdate);
    });
});

// Direct model testing endpoint with streaming
router.post('/test-gpqa', async (req, res) => {
    const requestId = `gpqa-test-${uuidv4()}`;
    const startTime = Date.now();
    
    const { question, query, config } = req.body;
    
    // Send requestId immediately so client can connect to stream
    res.json({ requestId, status: 'started' });
    
    // Initialize response tracking
    const responseData = {
        questionId: question.id,
        requestId,
        timestamp: new Date().toISOString(),
        config,
        systemLogs: [],
        thinking: [],
        response: null,
        modelAnswer: null,
        correct: null,
        responseTime: null,
        tokensUsed: null,
        error: null
    };

    // Add system log and emit to stream
    const addLog = (level, type, message) => {
        const log = {
            timestamp: new Date().toISOString(),
            level,
            type,
            message
        };
        responseData.systemLogs.push(log);
        console.log(`[${requestId}] ${level.toUpperCase()} - ${type}: ${message}`);
        
        // Emit log update
        streamEmitter.emit(requestId, {
            type: 'log',
            data: log
        });
    };

    // Emit thinking update
    const addThinking = (thought) => {
        responseData.thinking.push(thought);
        streamEmitter.emit(requestId, {
            type: 'thinking',
            data: thought
        });
    };

    try {
        addLog('info', 'start', `Testing question ${question.id} with ${config.model}`);
        
        // Set model for OpenRouter service
        openRouterService.setModel(config.model);
        
        let finalResponse;
        
        switch (config.routingMode) {
            case 'direct':
                // Direct model call without routing
                addLog('info', 'routing', 'Using direct mode - no workflow routing');
                
                const systemPrompt = `You are an expert in ${question.category}. Answer the following multiple choice question by analyzing it carefully and selecting the correct option. Provide detailed reasoning for your answer. At the end of your response, clearly state your final answer as "The answer is (X)" where X is the letter of the correct option.`;
                
                const messages = [
                    { role: 'system', content: systemPrompt },
                    { role: 'user', content: query }
                ];
                
                addLog('info', 'api_call', 'Sending request to OpenRouter');
                
                // Stream partial responses if available
                streamEmitter.emit(requestId, {
                    type: 'status',
                    data: 'Waiting for model response...'
                });
                
                const apiResponse = await openRouterService.sendRequest(messages, {
                    max_tokens: config.maxTokens,
                    temperature: 0.3
                });
                
                finalResponse = apiResponse.content;
                responseData.tokensUsed = apiResponse.usage?.total_tokens || 0;
                
                // Stream the full response
                streamEmitter.emit(requestId, {
                    type: 'response',
                    data: finalResponse
                });
                
                // Extract thinking process (for direct mode, use the full response)
                addThinking(finalResponse);
                
                break;
                
            case 'workflow':
                // Use workflow engine with full routing
                addLog('info', 'routing', 'Using workflow engine with full routing');
                
                // Hook into workflow engine logs
                const originalLog = console.log;
                console.log = (...args) => {
                    originalLog(...args);
                    const message = args.join(' ');
                    if (message.includes(requestId) || message.includes('agent') || message.includes('routing')) {
                        streamEmitter.emit(requestId, {
                            type: 'workflow_log',
                            data: message
                        });
                    }
                };
                
                try {
                    streamEmitter.emit(requestId, {
                        type: 'status',
                        data: 'Analyzing query and routing to appropriate agent...'
                    });
                    
                    const workflowResult = await workflowEngine.answer(query, {});
                    finalResponse = workflowResult.content;
                    
                    // Stream the response
                    streamEmitter.emit(requestId, {
                        type: 'response',
                        data: finalResponse
                    });
                    
                    // Extract metadata
                    responseData.tokensUsed = workflowResult.metadata?.tokensUsed || 0;
                    
                    // Add reasoning steps
                    if (workflowResult.reasoning && Array.isArray(workflowResult.reasoning)) {
                        workflowResult.reasoning.forEach(reason => addThinking(reason));
                    }
                    
                    // Log routing decisions
                    if (workflowResult.metadata?.agent) {
                        addLog('info', 'routing_decision', `Routed to: ${workflowResult.metadata.agent}`);
                    }
                    
                } finally {
                    // Restore original console.log
                    console.log = originalLog;
                }
                
                break;
                
            case 'math':
                // Direct math agent call
                addLog('info', 'routing', 'Using math agent directly');
                
                try {
                    streamEmitter.emit(requestId, {
                        type: 'status',
                        data: 'Processing mathematical query...'
                    });
                    
                    const mathResult = await mathAgent.handle(
                        query, 
                        { isMath: true, needsExecution: true, complexity: 'complex' },
                        {},
                        { requestId }
                    );
                    finalResponse = mathResult.content;
                    
                    // Stream the response
                    streamEmitter.emit(requestId, {
                        type: 'response',
                        data: finalResponse
                    });
                    
                    // Extract execution results
                    if (mathResult.executionResults) {
                        mathResult.executionResults.forEach((exec, idx) => {
                            addThinking(`Step ${idx + 1} Code:\n${exec.code}`);
                            addThinking(`Step ${idx + 1} Output:\n${exec.result.output}`);
                        });
                    }
                    
                    responseData.tokensUsed = mathResult.metadata?.tokensUsed || 0;
                    
                } catch (mathError) {
                    addLog('error', 'math_agent', mathError.message);
                    throw mathError;
                }
                
                break;
        }
        
        // Store response
        responseData.response = finalResponse;
        responseData.responseTime = Date.now() - startTime;
        
        // Extract model's answer
        const answerMatch = finalResponse.match(/[Tt]he (?:answer|final answer) is \(?([A-D])\)?/);
        if (answerMatch) {
            responseData.modelAnswer = answerMatch[1];
            responseData.correct = responseData.modelAnswer === question.correct_answer;
            addLog('info', 'answer_extraction', `Model answered: ${responseData.modelAnswer}, Expected: ${question.correct_answer}`);
        } else {
            addLog('warning', 'answer_extraction', 'Could not extract answer from response');
        }
        
        addLog('info', 'complete', `Test completed in ${responseData.responseTime}ms`);
        
        // Send final complete event
        streamEmitter.emit(requestId, {
            type: 'complete',
            data: responseData
        });
        
    } catch (error) {
        addLog('error', 'test_error', error.message);
        responseData.error = error.message;
        responseData.responseTime = Date.now() - startTime;
        
        streamEmitter.emit(requestId, {
            type: 'error',
            data: error.message
        });
    }
});

// Health check
router.get('/health', (req, res) => {
    res.json({
        status: 'ok',
        services: {
            workflowEngine: 'active',
            openRouter: 'active',
            mathAgent: 'active',
            e2b: workflowEngine.config.E2B_API_KEY ? 'configured' : 'not configured'
        }
    });
});

module.exports = router;