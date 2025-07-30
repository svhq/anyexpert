# Claude Code Guidelines

1. First think through the problem, read the codebase for relevant files, and write a plan to tasks/todo.md.
2. The plan should have a list of todo items that you can check off as you complete them
3. Before you begin working, check in with me and I will verify the plan.
4. Then, begin working on the todo items, marking them as complete as you go.
5. Please every step of the way just give me a high level explanation of what changes you made
6. Make every task and code change you do as simple as possible. We want to avoid making any massive or complex changes. Every change should impact as little code as possible. Everything is about simplicity.
7. Finally, add a review section to the todo.md file with a summary of the changes you made and any other relevant information.

## Project-Specific Guidelines

### Ask Any Expert System
- This is a web-augmented expert system combining AI with intelligent search
- Follow the phased development plan outlined in README.md
- Maintain backward compatibility when adding new features
- Use structured logging for all significant operations
- Keep search costs low (≤12 Serper calls per complex query)
- Target ≤6 second response times

### Testing
- Always test changes with both direct answers and search-required queries
- Use `npm test-workflow` for comprehensive testing
- Monitor Serper API usage and response times

### Code Style
- Use singleton patterns for service modules
- Handle JSON parsing failures gracefully (GLM model workaround)
- Include JSDoc comments for public methods
- Keep modules focused on single responsibilities