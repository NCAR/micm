# Debug Test

Debug a failing MICM test with tracing support.

## Usage

- `/debug-test <test_name>` - Debug a specific test

## Instructions

When the user invokes this skill with a test name:

1. **Run the test first** to confirm failure:
   ```bash
   ctest --test-dir /Users/fillmore/EarthSystem/MICM/build --output-on-failure -R "<test_name>" -V
   ```

2. **If the test passes**, report success and exit.

3. **If the test fails**, analyze the output:
   - Look for NaN, Inf, or exploding values
   - Check for assertion failures
   - Identify which solver component failed

4. **For numerical failures** (NaN/Inf/explosion):
   - Offer to add debug tracing to `rosenbrock.inl`
   - Key trace points:
     - After `AlphaMinusJacobian`: print diagonal values
     - After linear solve: print K values
     - In forcing: print residual values

5. **Add tracing code** (with user confirmation):
   ```cpp
   // In rosenbrock.inl Solve() loop, after LinearFactor call:
   std::cerr << "=== Step " << result.stats_.number_of_steps_ << " ===" << std::endl;
   std::cerr << "H = " << H << ", alpha = " << (1.0 / (H * parameters.gamma_[0])) << std::endl;
   for (std::size_t i = 0; i < state.jacobian_.NumColumns(); ++i)
     std::cerr << "Diag[" << i << "] = " << state.jacobian_.DiagonalElement(0, i) << std::endl;
   ```

6. **Rebuild and rerun** with tracing enabled.

7. **Analyze trace output**:
   - Compare diagonal magnitudes (should be similar)
   - Check K value scaling (should scale with H)
   - Look for sign errors

8. **Clean up** tracing code when done.

## Common Issues

| Symptom | Likely Cause | Check |
|---------|--------------|-------|
| Values explode to 10^16+ | Ill-conditioned matrix | Diagonal magnitudes |
| NaN after linear solve | Singular matrix | Zero diagonals |
| Wrong steady state | Sign error in Jacobian | SubtractJacobianTerms signs |
| Constraint not satisfied | Missing projection step | Post-solve enforcement |
