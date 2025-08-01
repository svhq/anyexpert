const fs = require('fs');
const path = require('path');
const { ModelLoader } = require('./models/model-loader');

/**
 * Model Swap Utility - Easy way to switch between models
 */
class ModelSwapper {
  constructor() {
    this.modelLoader = new ModelLoader();
    this.envPath = path.join(__dirname, '.env');
  }

  /**
   * Get current model from environment
   */
  getCurrentModel() {
    const current = this.modelLoader.getCurrentModel();
    return current;
  }

  /**
   * Swap to a new model
   */
  swapModel(modelKey) {
    // Validate model exists
    const model = this.modelLoader.getModel(modelKey);
    
    // Read current .env file
    let envContent = '';
    try {
      envContent = fs.readFileSync(this.envPath, 'utf8');
    } catch (error) {
      console.error('Warning: .env file not found, creating new one');
    }
    
    // Update or add OPENROUTER_MODEL line
    const modelLine = `OPENROUTER_MODEL=${model.id}`;
    const lines = envContent.split('\n');
    let modelLineFound = false;
    
    for (let i = 0; i < lines.length; i++) {
      if (lines[i].startsWith('OPENROUTER_MODEL=')) {
        lines[i] = modelLine;
        modelLineFound = true;
        break;
      }
    }
    
    if (!modelLineFound) {
      // Add after the API key line
      const apiKeyIndex = lines.findIndex(line => line.startsWith('OPENROUTER_API_KEY='));
      if (apiKeyIndex !== -1) {
        lines.splice(apiKeyIndex + 1, 0, modelLine);
      } else {
        lines.push(modelLine);
      }
    }
    
    // Write back to file
    fs.writeFileSync(this.envPath, lines.join('\n'));
    
    // Update current process env
    process.env.OPENROUTER_MODEL = model.id;
    
    return model;
  }

  /**
   * Display swap information
   */
  displaySwapInfo(fromModel, toModel) {
    console.log('\n🔄 MODEL SWAP SUCCESSFUL\n');
    
    console.log('📤 Previous Model:');
    console.log(`   Name: ${fromModel.name}`);
    console.log(`   ID: ${fromModel.id}`);
    console.log(`   Cost: $${fromModel.pricing.input}/$${fromModel.pricing.output} per 1M tokens`);
    
    console.log('\n📥 New Model:');
    console.log(`   Name: ${toModel.name}`);
    console.log(`   ID: ${toModel.id}`);
    console.log(`   Cost: $${toModel.pricing.input}/$${toModel.pricing.output} per 1M tokens`);
    console.log(`   Features: ${Object.entries(toModel.features)
      .filter(([k, v]) => v)
      .map(([k]) => k)
      .join(', ')}`);
    
    // Cost comparison
    const costRatio = {
      input: toModel.pricing.input / (fromModel.pricing.input || 0.01),
      output: toModel.pricing.output / (fromModel.pricing.output || 0.01)
    };
    
    if (fromModel.pricing.input === 0) {
      console.log(`\n💰 Cost Impact: Moving from FREE to PAID model`);
    } else {
      console.log(`\n💰 Cost Impact:`);
      console.log(`   Input: ${costRatio.input > 1 ? '🔺' : '🔻'} ${Math.abs(costRatio.input - 1) * 100}% ${costRatio.input > 1 ? 'more' : 'less'} expensive`);
      console.log(`   Output: ${costRatio.output > 1 ? '🔺' : '🔻'} ${Math.abs(costRatio.output - 1) * 100}% ${costRatio.output > 1 ? 'more' : 'less'} expensive`);
    }
    
    console.log('\n✅ .env file updated successfully');
    console.log('📝 Restart your application to use the new model\n');
  }

  /**
   * Interactive model selection
   */
  async interactiveSwap() {
    const readline = require('readline').createInterface({
      input: process.stdin,
      output: process.stdout
    });
    
    const question = (query) => new Promise(resolve => readline.question(query, resolve));
    
    console.log('\n🤖 INTERACTIVE MODEL SWAP\n');
    
    // Show current model
    const current = this.getCurrentModel();
    console.log(`Current Model: ${current.name} (${current.key})\n`);
    
    // List all models
    console.log('Available Models:');
    const models = this.modelLoader.getAllModels();
    const modelKeys = Object.keys(models);
    
    modelKeys.forEach((key, index) => {
      const model = models[key];
      console.log(`${index + 1}. ${key}`);
      console.log(`   ${model.name} - $${model.pricing.input}/$${model.pricing.output} per 1M`);
    });
    
    // Get user choice
    const choice = await question('\nEnter model number or key: ');
    readline.close();
    
    let selectedKey;
    if (isNaN(choice)) {
      selectedKey = choice;
    } else {
      const index = parseInt(choice) - 1;
      if (index >= 0 && index < modelKeys.length) {
        selectedKey = modelKeys[index];
      }
    }
    
    if (!selectedKey || !models[selectedKey]) {
      console.error('❌ Invalid selection');
      return;
    }
    
    // Perform swap
    const newModel = this.swapModel(selectedKey);
    this.displaySwapInfo(current, { ...newModel, key: selectedKey });
  }
}

// Command line interface
if (require.main === module) {
  const swapper = new ModelSwapper();
  const modelKey = process.argv[2];
  
  if (!modelKey) {
    // Interactive mode
    swapper.interactiveSwap().catch(console.error);
  } else if (modelKey === 'current') {
    // Show current model
    const current = swapper.getCurrentModel();
    console.log('\n📍 Current Model:');
    console.log(`   Key: ${current.key}`);
    console.log(`   Name: ${current.name}`);
    console.log(`   ID: ${current.id}`);
    console.log(`   Provider: ${current.provider}`);
  } else if (modelKey === 'list') {
    // List all models
    const loader = swapper.modelLoader;
    loader.listModels();
  } else {
    // Direct swap
    try {
      const current = swapper.getCurrentModel();
      const newModel = swapper.swapModel(modelKey);
      swapper.displaySwapInfo(current, { ...newModel, key: modelKey });
    } catch (error) {
      console.error(`❌ Error: ${error.message}`);
      console.log('\nUsage:');
      console.log('  node swap-model.js              # Interactive mode');
      console.log('  node swap-model.js <model-key>  # Direct swap');
      console.log('  node swap-model.js current      # Show current model');
      console.log('  node swap-model.js list         # List all models');
      process.exit(1);
    }
  }
}

module.exports = ModelSwapper;