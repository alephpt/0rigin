# Contributing to Kuramoto Model Implementation

## Development Standards

This project maintains strict professional development standards to ensure code quality, maintainability, and scientific rigor.

### Code Quality Requirements

#### File and Function Size Limits
- **Files**: Maximum 500 lines (split into modules if exceeded)
- **Functions**: Target < 50 lines (complex numerical methods may be 50-55 lines)
- **Nesting**: Maximum depth 3 levels

These limits are MANDATORY and enforced during code review.

#### Style and Documentation
- **PEP 8 compliance**: All code must follow PEP 8 style guidelines
- **Type hints**: Required for all function signatures
- **Docstrings**: NumPy-style docstrings for all public APIs
- **Comments**: Explain why, not what (code should be self-documenting)

#### Code Structure
- Clean, readable, self-documenting code
- Descriptive naming (no cryptic abbreviations)
- Comprehensive error handling (no silent failures)
- Proper logging at appropriate levels
- Zero compiler/linter warnings tolerated
- No TODO comments in production code

### Architecture Principles

**Core Principle**: Know it early → lock it down. Know it late → keep it flexible.

**Applied**:
- **State Management**: Compile-time knowable → `const`/`immutable`, runtime-contingent → dynamic
- **Type Safety**: Static knowledge → strict typing, dynamic knowledge → flexible with validation
- **Separation of Concerns**: Clear boundaries between compute logic and runtime context
- **Performance**: Prefer O(n) or O(log n) over O(n²) where possible

### Testing Requirements

All new functionality must include comprehensive tests:

#### Test Coverage
- Unit tests for all new functions and classes
- Integration tests for end-to-end workflows
- Validation tests against known analytical solutions
- Target: 100% coverage for critical paths

#### Running Tests
```bash
PYTHONPATH=src:$PYTHONPATH pytest tests/
```

All tests must pass before code can be merged.

#### Test Structure
- Use pytest fixtures for common setup
- Test edge cases and error conditions
- Validate numerical accuracy where applicable
- Include docstrings explaining what is being tested

### Zero Tolerance Policy

The following are NEVER acceptable in production code:

- ❌ Duplicate files or components
- ❌ Stub implementations or mock data
- ❌ Hardcoded credentials or secrets
- ❌ Orphaned code or commented-out blocks
- ❌ Unhandled edge cases
- ❌ Copy-paste code duplication
- ❌ Functions or files exceeding size limits

### Scientific Validation

For numerical algorithms:
- Compare against analytical solutions where available
- Verify convergence properties
- Check conservation laws (energy, etc.)
- Validate against published results
- Document all assumptions and limitations

### Git Workflow

#### Branch Strategy
1. Create feature branch from main: `git checkout -b feature/your-feature-name`
2. Make changes with clear, atomic commits
3. Run full test suite locally
4. Submit pull request with description

#### Commit Messages
- Use clear, descriptive commit messages
- Format: `<type>: <description>`
- Types: `feat`, `fix`, `refactor`, `test`, `docs`, `perf`
- Example: `feat: Add Klein-Gordon field solver with CFL stability`

#### Pull Request Requirements
- Description of changes and motivation
- All tests pass
- Code review approval
- No merge conflicts with main

### Code Review Checklist

Before submitting a pull request, verify:

- [ ] All tests pass locally
- [ ] No code duplication
- [ ] Files < 500 lines, functions < 50 lines
- [ ] Type hints present for all function signatures
- [ ] NumPy-style docstrings complete
- [ ] No TODOs or FIXMEs in production code
- [ ] No hardcoded values or secrets
- [ ] Error handling comprehensive
- [ ] Scientific accuracy validated (if applicable)
- [ ] Examples updated (if public API changed)
- [ ] Documentation updated (if needed)

### Adding New Features

When adding new functionality:

1. **Plan**: Review existing architecture, identify extension points
2. **Design**: Follow established patterns and principles
3. **Implement**: Write clean, tested, documented code
4. **Validate**: Test thoroughly, including edge cases
5. **Document**: Update docs, add examples if appropriate
6. **Review**: Submit PR with clear description

### Module Organization

Follow the existing module structure:

```
src/kuramoto/
├── core/               # Core Kuramoto model
├── distributions/      # Frequency distributions
├── solvers/           # ODE solvers
├── analysis/          # Analysis tools
├── visualization/     # Plotting functions
└── field_theory/      # MSFT extensions
```

Place new code in the appropriate module. If creating a new module, discuss with maintainers first.

### Dependencies

- **Core dependencies**: Keep minimal (NumPy, SciPy, Matplotlib)
- **Optional dependencies**: Mark clearly as optional (e.g., py-pde)
- **New dependencies**: Justify necessity, prefer well-maintained packages
- **Version constraints**: Specify minimum versions in requirements.txt

### Documentation

#### Code Documentation
- All public functions, classes, and modules must have docstrings
- Use NumPy docstring format
- Include mathematical formulas where relevant
- Provide usage examples in docstrings

#### User Documentation
- Update README.md for user-facing changes
- Update ARCHITECTURE.md for structural changes
- Add examples for new features
- Keep docs concise and accurate

### Performance Considerations

- Use NumPy vectorization for array operations
- Avoid Python loops over large arrays
- Profile before optimizing
- Document performance characteristics in docstrings
- Consider memory usage for large-scale simulations

### Getting Help

- **Questions**: Open a GitHub issue with the `question` label
- **Bugs**: Open an issue with the `bug` label and include reproduction steps
- **Feature requests**: Open an issue with the `enhancement` label
- **Discussions**: Use GitHub Discussions for design questions

## Recognition

Contributors will be acknowledged in commit history and may be listed in project documentation for significant contributions.

## License

By contributing, you agree that your contributions will be licensed under the same terms as the project (Research use - 0rigin Project).

---

**Thank you for contributing to the Kuramoto Model implementation!**
