const express = require('express');
const cors = require('cors');
const path = require('path');
const config = require('./config');

// Import API routes
const gpqaTestRoutes = require('./api/gpqa-test-api-simple');

const app = express();
const PORT = process.env.PORT || 3456;

// Middleware
app.use(cors());
app.use(express.json({ limit: '50mb' }));
app.use(express.urlencoded({ extended: true, limit: '50mb' }));

// Request logging
app.use((req, res, next) => {
    console.log(`[${new Date().toISOString()}] ${req.method} ${req.path}`);
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
    console.error('Server error:', err);
    res.status(500).json({
        error: 'Internal server error',
        message: err.message
    });
});

// Start server
app.listen(PORT, () => {
    console.log(`\n🚀 GPQA Testing UI Server Running`);
    console.log(`================================`);
    console.log(`📍 URL: http://localhost:${PORT}`);
    console.log(`📊 API: http://localhost:${PORT}/api/test-gpqa`);
    console.log(`🔧 E2B Service: ${config.codeExecution?.serviceUrl || 'http://localhost:8001'}`);
    console.log(`🤖 Default Model: ${config.MODEL || 'google/gemini-2.5-flash-lite'}`);
    console.log(`\n✅ Ready for testing!`);
});