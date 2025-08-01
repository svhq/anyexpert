const fs = require('fs');
const path = require('path');

/**
 * Model configuration loader and manager
 */
class ModelLoader {
  constructor() {
    this.configPath = path.join(__dirname, 'config.json');
    this.loadConfig();
  }

  /**
   * Load model configurations from config.json
   */
  loadConfig() {
    try {
      const configData = fs.readFileSync(this.configPath, 'utf8');
      this.config = JSON.parse(configData);
      this.models = this.config.models;
    } catch (error) {
      console.error('Failed to load model config:', error.message);
      this.models = {};
    }
  }

  /**
   * Get all available models
   * @returns {Object} All model configurations
   */
  getAllModels() {
    return this.models;
  }

  /**
   * Get a specific model configuration
   * @param {string} modelKey - The model key (e.g., 'glm-4.5-air')
   * @returns {Object} Model configuration
   */
  getModel(modelKey) {
    if (!this.models[modelKey]) {
      throw new Error(`Model '${modelKey}' not found. Available models: ${Object.keys(this.models).join(', ')}`);
    }
    return this.models[modelKey];
  }

  /**
   * Get model ID for OpenRouter API
   * @param {string} modelKey - The model key
   * @returns {string} OpenRouter model ID
   */
  getModelId(modelKey) {
    const model = this.getModel(modelKey);
    return model.id;
  }

  /**
   * Get current model from environment or default
   * @returns {Object} Current model configuration
   */
  getCurrentModel() {
    const envModel = process.env.OPENROUTER_MODEL;
    
    // Find model by ID if env var is set
    if (envModel) {
      const modelEntry = Object.entries(this.models).find(([key, model]) => model.id === envModel);
      if (modelEntry) {
        return { key: modelEntry[0], ...modelEntry[1] };
      }
    }
    
    // Return default model
    const defaultKey = this.config.defaultModel || 'glm-4.5-air';
    return { key: defaultKey, ...this.models[defaultKey] };
  }

  /**
   * Calculate cost for token usage
   * @param {string} modelKey - The model key
   * @param {number} inputTokens - Number of input tokens
   * @param {number} outputTokens - Number of output tokens
   * @returns {Object} Cost breakdown
   */
  calculateCost(modelKey, inputTokens, outputTokens) {
    const model = this.getModel(modelKey);
    const pricing = model.pricing;
    
    // Convert to millions for pricing calculation
    const inputMillions = inputTokens / 1_000_000;
    const outputMillions = outputTokens / 1_000_000;
    
    const inputCost = inputMillions * pricing.input;
    const outputCost = outputMillions * pricing.output;
    const totalCost = inputCost + outputCost;
    
    return {
      inputTokens,
      outputTokens,
      totalTokens: inputTokens + outputTokens,
      inputCost: parseFloat(inputCost.toFixed(6)),
      outputCost: parseFloat(outputCost.toFixed(6)),
      totalCost: parseFloat(totalCost.toFixed(6)),
      currency: pricing.currency,
      model: model.name
    };
  }

  /**
   * Get models by feature
   * @param {string} feature - Feature to filter by (e.g., 'tools', 'streaming')
   * @returns {Array} Models that support the feature
   */
  getModelsByFeature(feature) {
    return Object.entries(this.models)
      .filter(([key, model]) => model.features[feature])
      .map(([key, model]) => ({ key, ...model }));
  }

  /**
   * Get models by price range
   * @param {number} maxInputPrice - Maximum input price per million tokens
   * @param {number} maxOutputPrice - Maximum output price per million tokens
   * @returns {Array} Models within price range
   */
  getModelsByPriceRange(maxInputPrice, maxOutputPrice) {
    return Object.entries(this.models)
      .filter(([key, model]) => 
        model.pricing.input <= maxInputPrice && 
        model.pricing.output <= maxOutputPrice
      )
      .map(([key, model]) => ({ key, ...model }));
  }

  /**
   * Get testing profile models
   * @param {string} profile - Profile name (quick, comprehensive, budget, premium)
   * @returns {Array} Model keys for the profile
   */
  getTestingProfile(profile) {
    return this.config.testingProfiles[profile] || [];
  }

  /**
   * List all models with basic info
   */
  listModels() {
    console.log('\n📋 Available Models:\n');
    Object.entries(this.models).forEach(([key, model]) => {
      console.log(`${key}:`);
      console.log(`  Name: ${model.name}`);
      console.log(`  Provider: ${model.provider}`);
      console.log(`  Pricing: $${model.pricing.input}/$${model.pricing.output} per 1M tokens`);
      console.log(`  Context: ${model.limits.contextWindow.toLocaleString()} tokens`);
      console.log(`  Features: ${Object.entries(model.features)
        .filter(([k, v]) => v)
        .map(([k]) => k)
        .join(', ')}`);
      console.log('');
    });
  }
}

// Singleton instance
let instance = null;

module.exports = {
  ModelLoader,
  
  // Get singleton instance
  getInstance() {
    if (!instance) {
      instance = new ModelLoader();
    }
    return instance;
  }
};

// If run directly, list all models
if (require.main === module) {
  const loader = new ModelLoader();
  loader.listModels();
  
  console.log('\n💰 Cost Examples:');
  const exampleTokens = { input: 1000, output: 500 };
  ['glm-4.5-air', 'o4-mini-high', 'claude-3.5-sonnet'].forEach(modelKey => {
    const cost = loader.calculateCost(modelKey, exampleTokens.input, exampleTokens.output);
    console.log(`${modelKey}: $${cost.totalCost} (${exampleTokens.input} in, ${exampleTokens.output} out)`);
  });
}