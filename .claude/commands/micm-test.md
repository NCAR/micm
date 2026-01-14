# MICM Test Runner

Build and run MICM tests efficiently.

## Usage

- `/micm-test` - Build and run all tests
- `/micm-test build` - Build only (no tests)
- `/micm-test <pattern>` - Run tests matching pattern (e.g., `constraint`, `rosenbrock`)

## Instructions

When the user invokes this skill:

1. **Parse the argument** (if any):
   - No argument or empty: run all tests
   - `build`: build only
   - Any other string: use as test filter pattern

2. **Build the project**:
   ```bash
   cmake --build /Users/fillmore/EarthSystem/MICM/build -j$(sysctl -n hw.ncpu)
   ```

3. **Run tests** (unless build-only):
   - All tests: `ctest --test-dir /Users/fillmore/EarthSystem/MICM/build --output-on-failure`
   - Filtered: `ctest --test-dir /Users/fillmore/EarthSystem/MICM/build --output-on-failure -R "<pattern>"`

4. **Report results**:
   - Show pass/fail count
   - If failures, show which tests failed
   - Suggest `/debug-test <name>` for failed tests

## Example Output

```
Building MICM...
Build successful.

Running tests matching "constraint"...
3/3 tests passed:
  - constraint
  - constraint_set
  - equilibrium (contains constraint tests)
```
