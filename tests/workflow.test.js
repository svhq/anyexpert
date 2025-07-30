const workflowEngine = require('../src/workflow-engine');
const queryAnalyzer = require('../src/query-analyzer');
const searchPlanner = require('../src/search-planner');
const webSearch = require('../src/web-search');

// Mock dependencies
jest.mock('../src/query-analyzer');
jest.mock('../src/search-planner');
jest.mock('../src/web-search');
jest.mock('../src/utils/logger');

describe('Workflow Engine', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  describe('Direct Answer Flow', () => {
    test('returns direct answer when search is not needed', async () => {
      // Mock query analyzer to say no search needed
      queryAnalyzer.assess.mockResolvedValue({
        needsSearch: false,
        rationale: 'Basic mathematical calculation'
      });

      const query = 'What is 2 + 2?';
      const result = await workflowEngine.answer(query, {});

      expect(queryAnalyzer.assess).toHaveBeenCalledWith(query, {});
      expect(searchPlanner.generate).not.toHaveBeenCalled();
      expect(webSearch.runBatch).not.toHaveBeenCalled();
      expect(result.searchPerformed).toBe(false);
    });
  });

  describe('Search Flow', () => {
    test('performs search when needed and respects confidence threshold', async () => {
      // Mock query analyzer to require search
      queryAnalyzer.assess.mockResolvedValue({
        needsSearch: true,
        rationale: 'Current information needed'
      });

      // Mock search planner
      searchPlanner.generate
        .mockResolvedValueOnce(['query 1', 'query 2'])
        .mockResolvedValueOnce(['refined query 1', 'refined query 2']);

      // Mock web search results
      const mockResults = [
        { title: 'Result 1', url: 'http://example.com/1', snippet: 'Test snippet 1' },
        { title: 'Result 2', url: 'http://example.edu/2', snippet: 'Test snippet 2' }
      ];
      webSearch.runBatch.mockResolvedValue(mockResults);

      const query = 'Latest AI developments 2025';
      const result = await workflowEngine.answer(query, {});

      expect(queryAnalyzer.assess).toHaveBeenCalled();
      expect(searchPlanner.generate).toHaveBeenCalled();
      expect(webSearch.runBatch).toHaveBeenCalled();
      expect(result.metadata.rounds).toBeGreaterThanOrEqual(1);
    });

    test('stops searching when confidence threshold is reached', async () => {
      queryAnalyzer.assess.mockResolvedValue({
        needsSearch: true,
        rationale: 'Need current data'
      });

      searchPlanner.generate.mockResolvedValue(['test query']);
      
      // Mock results with .edu domain for high confidence
      webSearch.runBatch.mockResolvedValue([
        { title: 'Academic Result', url: 'http://stanford.edu/ai', snippet: 'Academic content' },
        { title: 'Gov Result', url: 'http://nsf.gov/research', snippet: 'Government content' }
      ]);

      const result = await workflowEngine.answer('AI research funding', {});
      
      // Should stop after first round due to authoritative sources
      expect(result.metadata.rounds).toBe(1);
      expect(result.metadata.finalConfidence).toBeGreaterThanOrEqual(0.8);
    });

    test('performs maximum 3 rounds of search', async () => {
      queryAnalyzer.assess.mockResolvedValue({
        needsSearch: true,
        rationale: 'Complex topic'
      });

      searchPlanner.generate.mockResolvedValue(['query']);
      
      // Mock results with low confidence (no authoritative sources)
      webSearch.runBatch.mockResolvedValue([
        { title: 'Blog Post', url: 'http://blog.com/post', snippet: 'Blog content' }
      ]);

      const result = await workflowEngine.answer('Complex technical question', {});
      
      // Should perform exactly 3 rounds (max)
      expect(result.metadata.rounds).toBe(3);
      expect(searchPlanner.generate).toHaveBeenCalledTimes(3);
    });
  });

  describe('Error Handling', () => {
    test('handles search errors gracefully', async () => {
      queryAnalyzer.assess.mockResolvedValue({
        needsSearch: true,
        rationale: 'Need search'
      });

      searchPlanner.generate.mockResolvedValue(['query']);
      webSearch.runBatch.mockRejectedValue(new Error('API Error'));

      await expect(workflowEngine.answer('Test query', {}))
        .rejects.toThrow('API Error');
    });
  });

  describe('Performance', () => {
    test('completes within reasonable time', async () => {
      queryAnalyzer.assess.mockResolvedValue({
        needsSearch: false,
        rationale: 'Direct answer'
      });

      const start = Date.now();
      await workflowEngine.answer('Quick question', {});
      const duration = Date.now() - start;

      // Should complete quickly for direct answers
      expect(duration).toBeLessThan(1000); // 1 second for direct answer
    });
  });
});

describe('Integration Tests', () => {
  test('end-to-end workflow with mocked external APIs', async () => {
    // This would be a more comprehensive integration test
    // For now, marking as todo
    expect(true).toBe(true);
  });
});