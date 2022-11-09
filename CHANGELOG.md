# Master (becomes Release v0.3)
- Added linear solver capability to accept matrix-valued rhs #76
- Added `EnhancedAssumedStrains` that can be used to decorated a linear elasticity element with several EAS method #74
- Added codespell workflow #70 (CI checks now for grammar and typos in comments and variable names)
- Added Path following technique to using, e.g. `Arc length method` the scalar subsidiary equation can be implemented indepentently #80

# Release v0.2 (Apollodorus)
- Moved from Catch2 testing environment to `Dune::TestSuite`
- Generalized local function test and added `linearStrains` expression

