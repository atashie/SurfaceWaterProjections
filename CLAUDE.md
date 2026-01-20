# Surface Water Projections

## Overview

Research repository for surface water modeling, streamflow analysis, and flood prediction.

## Project Structure

```
SurfaceWaterProjections/
├── streamflowSignatures/       # Hydrological signature extraction (R)
├── EA-LSTM_globalStreamflowModeling/  # Deep learning streamflow models (Python)
├── floodModeling/              # Flood prediction models
├── surfaceWaterModeling_*.r    # Legacy R modeling scripts
└── CLAUDE.md                   # This file
```

## Subproject Documentation

Each subproject has its own `CLAUDE.md` with specific conventions:
- `streamflowSignatures/CLAUDE.md` - R-based signature extraction system

## Development Practices

### Test-Driven Development (TDD)

Use the `superpowers:test-driven-development` skill for all new features and bug fixes.

**Core Principle:** No production code without a failing test first.

**Red-Green-Refactor Cycle:**
1. **RED** - Write one minimal failing test
2. **Verify RED** - Watch it fail for the expected reason
3. **GREEN** - Write minimal code to pass
4. **Verify GREEN** - Confirm all tests pass
5. **REFACTOR** - Clean up while staying green

**Iron Law:**
```
NO PRODUCTION CODE WITHOUT A FAILING TEST FIRST
```

Write code before the test? Delete it. Start over. No exceptions.

**When to Use TDD:**
- New features
- Bug fixes
- Refactoring
- Behavior changes

**Verification Checklist:**
- [ ] Every new function/method has a test
- [ ] Watched each test fail before implementing
- [ ] Each test failed for expected reason (feature missing, not typo)
- [ ] Wrote minimal code to pass each test
- [ ] All tests pass
- [ ] Output pristine (no errors, warnings)
- [ ] Edge cases and errors covered

**Common Rationalizations to Avoid:**
| Excuse | Reality |
|--------|---------|
| "Too simple to test" | Simple code breaks. Test takes 30 seconds. |
| "I'll test after" | Tests passing immediately prove nothing. |
| "Already manually tested" | Ad-hoc ≠ systematic. No record, can't re-run. |
| "TDD will slow me down" | TDD faster than debugging. |

### Code Quality

- Use structured logging (see `streamflowSignatures/config.R`)
- Validate inputs at function entry points
- Validate outputs before returning
- Handle errors with meaningful messages

### Independent Code Review

After implementing substantial code changes (>20 lines added/modified/deleted), invoke the `you-review` skill before committing. This spawns a separate Claude Sonnet agent to review changes with fresh eyes.

**Why a separate agent?** The reviewer has no context from the implementation conversation, providing unbiased feedback on bugs, security issues, and code quality.

**Triggers:**
- Before commits with substantial changes
- After completing plan-mode implementations
- On explicit request (`/you-review`)

## Languages & Tools

- **R**: Hydrological analysis, signature extraction
- **Python**: Deep learning models (PyTorch/TensorFlow)
- **Git**: Version control

## Available Claude Skills

- `superpowers:test-driven-development` - TDD workflow guidance
- `you-review` - Independent code review via Claude Sonnet agent (invoke after substantial changes)
- `scientific-skills:matplotlib` - Visualization
- `scientific-skills:geopandas` - Geospatial analysis
- `scientific-skills:aeon` - Time series ML
- `scientific-skills:networkx` - Network analysis
