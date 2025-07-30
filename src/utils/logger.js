const pino = require('pino');
const path = require('path');
const fs = require('fs');

// Ensure logs directory exists
const logsDir = path.join(__dirname, '../../logs');
if (!fs.existsSync(logsDir)) {
  fs.mkdirSync(logsDir, { recursive: true });
}

// Configure pino with pretty printing for development
const isDevelopment = process.env.NODE_ENV !== 'production';

const logger = pino({
  level: process.env.LOG_LEVEL || 'info',
  formatters: {
    level: (label) => {
      return { level: label };
    }
  },
  timestamp: pino.stdTimeFunctions.isoTime,
  base: {
    pid: process.pid,
    hostname: require('os').hostname()
  },
  // Pretty print in development
  transport: isDevelopment ? {
    target: 'pino-pretty',
    options: {
      colorize: true,
      translateTime: 'yyyy-mm-dd HH:MM:ss',
      ignore: 'pid,hostname'
    }
  } : undefined
});

// Production file logging
if (!isDevelopment) {
  // Create file stream for production
  const fileStream = pino.destination({
    dest: path.join(logsDir, 'app.log'),
    sync: false
  });

  // Override logger
  module.exports = pino({
    level: process.env.LOG_LEVEL || 'info',
    timestamp: pino.stdTimeFunctions.isoTime
  }, fileStream);
} else {
  module.exports = logger;
}

// Helper functions for structured logging
module.exports.logRequest = (requestId, phase, data) => {
  logger.info({
    requestId,
    phase,
    ...data
  });
};

module.exports.logError = (requestId, error, context) => {
  logger.error({
    requestId,
    error: {
      message: error.message,
      stack: error.stack,
      code: error.code
    },
    context
  });
};

module.exports.logMetrics = (metrics) => {
  logger.info({
    type: 'metrics',
    ...metrics
  });
};

// Enhanced logging functions for system monitoring
module.exports.logJsonParsingFailure = (requestId, service, content, error) => {
  logger.warn({
    requestId,
    service,
    type: 'json_parsing_failure',
    contentLength: content?.length || 0,
    contentPreview: content?.substring(0, 200) || '',
    error: error.message
  });
};

module.exports.logSearchDecision = (requestId, query, needsSearch, rationale, method) => {
  logger.info({
    requestId,
    type: 'search_decision',
    query,
    needsSearch,
    rationale,
    analysisMethod: method // 'json', 'fallback', 'heuristic'
  });
};

module.exports.logSerperResponse = (requestId, queryCount, resultCount, hasUndefinedUrls) => {
  logger.info({
    requestId,
    type: 'serper_response',
    queryCount,
    resultCount,
    hasUndefinedUrls
  });
};

module.exports.logResponseQuality = (requestId, responseLength, wasTruncated, citationCount) => {
  logger.info({
    requestId,
    type: 'response_quality',
    responseLength,
    wasTruncated,
    citationCount
  });
};

module.exports.logCodeExecution = (requestId, language, codeLength, success) => {
  logger.info({
    requestId,
    type: 'code_execution',
    language,
    codeLength,
    success
  });
};