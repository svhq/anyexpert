const express = require('express');
const router = express.Router();
const workflowEngine = require('../src/workflow-engine');
const OpenRouterService = require('../src/openrouter-service');
const mathAgent = require('../src/math-agent');
const winston = require('winston');
const { v4: uuidv4 } = require('uuid');

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

// Direct model testing endpoint
router.post('/test-gpqa', async (req, res) => {
    const requestId = `gpqa-test-${uuidv4()}`;
    const startTime = Date.now();
    
    const { question, query, config } = req.body;
    
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

    // Add system log
    const addLog = (level, type, message) => {
        const log = {
            timestamp: new Date().toISOString(),
            level,
            type,
            message
        };
        responseData.systemLogs.push(log);
        apiLogger[level](`[${requestId}] ${type}: ${message}`);
    };

    try {
        addLog('info', 'start', `Testing question ${question.id} with ${config.model}`);
        
        // Set model for OpenRouter service
        openRouterService.model = config.model;
        
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
                
                const apiResponse = await openRouterService.sendRequest(messages, {
                    max_tokens: config.maxTokens,
                    temperature: 0.3
                });
                
                finalResponse = apiResponse.content;
                responseData.tokensUsed = apiResponse.usage?.total_tokens || 0;
                
                // Extract thinking process (for direct mode, use the full response)
                responseData.thinking = [finalResponse];
                
                break;
                
            case 'workflow':
                // Use workflow engine with full routing
                addLog('info', 'routing', 'Using workflow engine with full routing');
                
                // Temporarily override timeout
                const originalTimeout = workflowEngine.config.E2B_TIMEOUT;
                workflowEngine.config.E2B_TIMEOUT = config.timeout;
                
                try {
                    const workflowResult = await workflowEngine.answer(query, {});
                    finalResponse = workflowResult.content;
                    
                    // Extract metadata
                    responseData.tokensUsed = workflowResult.metadata?.tokensUsed || 0;
                    responseData.thinking = workflowResult.reasoning || [];
                    
                    // Log routing decisions
                    if (workflowResult.metadata?.agent) {
                        addLog('info', 'routing_decision', `Routed to: ${workflowResult.metadata.agent}`);
                    }
                    
                } finally {
                    // Restore original timeout
                    workflowEngine.config.E2B_TIMEOUT = originalTimeout;
                }
                
                break;
                
            case 'math':
                // Direct math agent call
                addLog('info', 'routing', 'Using math agent directly');
                
                try {
                    const mathResult = await mathAgent.handle(
                        query, 
                        { isMath: true, needsExecution: true, complexity: 'complex' },
                        {},
                        { requestId }
                    );
                    finalResponse = mathResult.content;
                    
                    // Extract execution results
                    if (mathResult.executionResults) {
                        mathResult.executionResults.forEach((exec, idx) => {
                            responseData.thinking.push(`Step ${idx + 1} Code:\n${exec.code}`);
                            responseData.thinking.push(`Step ${idx + 1} Output:\n${exec.result.output}`);
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
        
        res.json(responseData);
        
    } catch (error) {
        addLog('error', 'test_error', error.message);
        responseData.error = error.message;
        responseData.responseTime = Date.now() - startTime;
        
        res.status(500).json(responseData);
    }
});

// Get test history
router.get('/test-history/:questionId', (req, res) => {
    // This could be expanded to store and retrieve test history from a database
    res.json({ message: 'Test history not yet implemented' });
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