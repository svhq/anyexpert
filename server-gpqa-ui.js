const express = require('express');
const cors = require('cors');
const path = require('path');
const winston = require('winston');
const config = require('./config');

// Import API routes
const gpqaTestRoutes = require('./api/gpqa-test-api-streaming');

const app = express();
const PORT = process.env.PORT || 3456;

// Configure logger
const logger = winston.createLogger({
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
        }),
        new winston.transports.File({ 
            filename: path.join(__dirname, 'logs', 'gpqa-ui-server.log') 
        })
    ]
});

// Middleware
app.use(cors());
app.use(express.json({ limit: '50mb' }));
app.use(express.urlencoded({ extended: true, limit: '50mb' }));

// Request logging
app.use((req, res, next) => {
    logger.info(`${req.method} ${req.path}`, {
        query: req.query,
        body: req.body ? Object.keys(req.body) : undefined
    });
    next();
});

// API routes
app.use('/api', gpqaTestRoutes);

// Serve the testing UI
app.get('/', (req, res) => {
    res.sendFile(path.join(__dirname, 'gpqa-testing-ui.html'));
});

// Serve GPQA questions data
app.get('/data/gpqa-questions.json', (req, res) => {
    res.sendFile(path.join(__dirname, 'gpqa-10-questions-clean.json'));
});

// Error handling middleware
app.use((err, req, res, next) => {
    logger.error('Server error:', err);
    res.status(500).json({
        error: 'Internal server error',
        message: err.message
    });
});

// Start server
app.listen(PORT, () => {
    logger.info(`GPQA Testing UI server running on http://localhost:${PORT}`);
    logger.info('Configuration:', {
        model: config.MODEL,
        e2b: config.E2B_API_KEY ? 'configured' : 'not configured',
        openRouterKey: config.OPENROUTER_API_KEY ? 'configured' : 'not configured'
    });
    
    // Check E2B service
    if (!config.E2B_API_KEY) {
        logger.warn('E2B_API_KEY not configured - code execution will not work');
    }
    
    // Log available routes
    logger.info('Available endpoints:');
    logger.info('  GET  / - Testing UI');
    logger.info('  POST /api/test-gpqa - Run GPQA test');
    logger.info('  GET  /api/health - Health check');
});