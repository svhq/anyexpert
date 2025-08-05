/**
 * Flexible JSON parser that handles various model output formats
 * Designed to work with multiple LLM providers
 */

class FlexibleJSONParser {
  /**
   * Extract JSON from various text formats
   * @param {string} text - The text containing JSON
   * @returns {Object|null} - Parsed JSON object or null
   */
  static extract(text) {
    if (!text || typeof text !== 'string') {
      return null;
    }

    // Try these patterns in order of specificity
    const patterns = [
      // 1. Markdown code blocks with json label
      /```json\s*\n?([\s\S]*?)\n?```/,
      
      // 2. Markdown code blocks without label
      /```\s*\n?([\s\S]*?)\n?```/,
      
      // 3. JSON with curly braces (greedy to get complete object)
      /(\{[\s\S]*\})/,
      
      // 4. JSON array format
      /(\[[\s\S]*\])/
    ];

    for (const pattern of patterns) {
      const match = text.match(pattern);
      if (match) {
        try {
          // Get the captured group (1) or the full match (0)
          const jsonStr = match[1] || match[0];
          
          // Clean up common issues
          const cleaned = jsonStr
            .replace(/^[^{[]+/, '') // Remove text before JSON
            .replace(/[^}\]]+$/, '') // Remove text after JSON
            .trim();
          
          // Try to parse
          const parsed = JSON.parse(cleaned);
          return parsed;
        } catch (e) {
          // Continue to next pattern
          continue;
        }
      }
    }

    // Last resort: try to find JSON-like structure
    try {
      // Look for object-like structure with quotes
      const objectMatch = text.match(/\{[^{}]*"[^"]+"\s*:[^{}]*\}/);
      if (objectMatch) {
        return JSON.parse(objectMatch[0]);
      }
    } catch (e) {
      // Failed
    }

    return null;
  }

  /**
   * Parse JSON with detailed error info
   * @param {string} text - The text containing JSON
   * @returns {Object} - { success: boolean, data: Object|null, error: string|null }
   */
  static parse(text) {
    try {
      const extracted = this.extract(text);
      if (extracted) {
        return {
          success: true,
          data: extracted,
          error: null
        };
      }
      return {
        success: false,
        data: null,
        error: 'No valid JSON found in text'
      };
    } catch (error) {
      return {
        success: false,
        data: null,
        error: error.message
      };
    }
  }

  /**
   * Test parser with various formats
   */
  static test() {
    const testCases = [
      '{"key": "value"}',
      '```json\n{"key": "value"}\n```',
      'Here is the JSON: {"key": "value"}',
      'The result is:\n```\n{"key": "value"}\n```\nThat\'s it.',
      'Response: {"needsSearch": "Y", "rationale": "test"}',
      '```json{"key":"value"}```'
    ];

    console.log('Testing FlexibleJSONParser:\n');
    testCases.forEach((test, i) => {
      const result = this.parse(test);
      console.log(`Test ${i + 1}: ${result.success ? '✅' : '❌'}`);
      if (result.success) {
        console.log('  Parsed:', JSON.stringify(result.data));
      } else {
        console.log('  Error:', result.error);
      }
    });
  }
}

module.exports = FlexibleJSONParser;