/**
 * Central Tool Configuration System
 * Defines available tools and capabilities for different API modes
 */

// API Modes
const API_MODES = {
  NONE: 'none',        // No tools - pure reasoning only
  SEARCH: 'search',    // Web search only
  FULL: 'full'         // All tools (search + code execution)
};

// Get current API mode from environment
const getCurrentApiMode = () => {
  const mode = process.env.API_MODE || 'full';
  if (!Object.values(API_MODES).includes(mode)) {
    console.warn(`Invalid API_MODE "${mode}", defaulting to "full"`);
    return API_MODES.FULL;
  }
  return mode;
};

// Tool Definitions
const TOOL_DEFINITIONS = {
  search_web: {
    type: 'function',
    function: {
      name: 'search_web',
      description: 'Search the web for current information, facts, and recent events. Returns relevant snippets and sources.',
      parameters: {
        type: 'object',
        properties: {
          query: { 
            type: 'string',
            description: 'The search query'
          },
          count: { 
            type: 'integer',
            description: 'Number of results to return',
            default: 5
          }
        },
        required: ['query']
      }
    }
  },
  
  run_code: {
    type: 'function',
    function: {
      name: 'run_code',
      description: 'Execute Python code in an isolated sandbox. Available libraries: numpy, scipy, numpy_financial, sympy, mpmath, pandas, matplotlib, seaborn, scikit-learn, biopython, splicekit, deeptools, beautifulsoup4, requests.',
      parameters: {
        type: 'object',
        properties: {
          code: { 
            type: 'string',
            description: 'Python code to execute'
          },
          timeout: { 
            type: 'integer',
            description: 'Execution timeout in milliseconds',
            default: 30000
          }
        },
        required: ['code']
      }
    }
  }
};

// Available Python Libraries (for documentation/prompts)
const PYTHON_LIBRARIES = {
  mathematical: ['math', 'decimal', 'fractions', 'mpmath'],
  numerical: ['numpy', 'scipy', 'numpy_financial'],
  symbolic: ['sympy'],
  data_analysis: ['pandas', 'statistics', 'statsmodels'],
  visualization: ['matplotlib', 'seaborn', 'plotly'],
  machine_learning: ['scikit-learn'],
  bioinformatics: ['biopython', 'splicekit', 'deeptools'],
  web_utilities: ['beautifulsoup4', 'requests', 'openpyxl', 'xlrd'],
  date_time: ['python-dateutil']
};

// API Mode Configurations
const API_MODE_CONFIGS = {
  [API_MODES.NONE]: {
    name: 'No Tools',
    description: 'Pure reasoning without external tools',
    tools: [],
    systemPromptAddition: '',
    capabilities: {
      webSearch: false,
      codeExecution: false,
      pythonLibraries: []
    }
  },
  
  [API_MODES.SEARCH]: {
    name: 'Search Only',
    description: 'Web search capability only',
    tools: [TOOL_DEFINITIONS.search_web],
    systemPromptAddition: '\n\n**Tools Available:** You have access to web search for current information and fact-checking. Use search strategically when you need current/external information that may not be in your training data.',
    capabilities: {
      webSearch: true,
      codeExecution: false,
      pythonLibraries: []
    }
  },
  
  [API_MODES.FULL]: {
    name: 'Full Tools',
    description: 'All tools available (search + code execution)',
    tools: [TOOL_DEFINITIONS.search_web, TOOL_DEFINITIONS.run_code],
    systemPromptAddition: '\n\n**Tools Available:** You have access to: (1) Web search for current information and fact-checking, (2) Python code execution sandbox for calculations, analysis, and verification. When uncertainty exists or calculations are needed, use tools strategically: search for current/external info, code for calculations/verification. You may use multiple tools in parallel when beneficial.',
    capabilities: {
      webSearch: true,
      codeExecution: true,
      pythonLibraries: Object.values(PYTHON_LIBRARIES).flat()
    }
  }
};

// Get configuration for current API mode
const getCurrentConfig = () => {
  const mode = getCurrentApiMode();
  return {
    mode,
    ...API_MODE_CONFIGS[mode]
  };
};

// Generate library documentation text
const generateLibraryDocumentation = (libraries) => {
  if (libraries.length === 0) return '';
  
  const grouped = {
    mathematical: libraries.filter(lib => PYTHON_LIBRARIES.mathematical.includes(lib)),
    numerical: libraries.filter(lib => PYTHON_LIBRARIES.numerical.includes(lib)),
    symbolic: libraries.filter(lib => PYTHON_LIBRARIES.symbolic.includes(lib)),
    data_analysis: libraries.filter(lib => PYTHON_LIBRARIES.data_analysis.includes(lib)),
    visualization: libraries.filter(lib => PYTHON_LIBRARIES.visualization.includes(lib)),
    machine_learning: libraries.filter(lib => PYTHON_LIBRARIES.machine_learning.includes(lib)),
    bioinformatics: libraries.filter(lib => PYTHON_LIBRARIES.bioinformatics.includes(lib)),
    web_utilities: libraries.filter(lib => PYTHON_LIBRARIES.web_utilities.includes(lib)),
    date_time: libraries.filter(lib => PYTHON_LIBRARIES.date_time.includes(lib))
  };
  
  let docs = '\n\nAvailable Python libraries:\n';
  
  if (grouped.mathematical.length > 0) {
    docs += `- ${grouped.mathematical.join(', ')} (mathematical operations)\n`;
  }
  if (grouped.numerical.length > 0) {
    docs += `- ${grouped.numerical.join(', ')} (numerical computing)\n`;
  }
  if (grouped.symbolic.length > 0) {
    docs += `- ${grouped.symbolic.join(', ')} (symbolic math)\n`;
  }
  if (grouped.data_analysis.length > 0) {
    docs += `- ${grouped.data_analysis.join(', ')} (data analysis)\n`;
  }
  if (grouped.visualization.length > 0) {
    docs += `- ${grouped.visualization.join(', ')} (visualization)\n`;
  }
  if (grouped.machine_learning.length > 0) {
    docs += `- ${grouped.machine_learning.join(', ')} (machine learning)\n`;
  }
  if (grouped.bioinformatics.length > 0) {
    docs += `- ${grouped.bioinformatics.join(', ')} (bioinformatics)\n`;
  }
  if (grouped.web_utilities.length > 0) {
    docs += `- ${grouped.web_utilities.join(', ')} (web utilities)\n`;
  }
  if (grouped.date_time.length > 0) {
    docs += `- ${grouped.date_time.join(', ')} (date/time handling)\n`;
  }
  
  return docs.trim();
};

module.exports = {
  API_MODES,
  getCurrentApiMode,
  getCurrentConfig,
  TOOL_DEFINITIONS,
  PYTHON_LIBRARIES,
  generateLibraryDocumentation
};