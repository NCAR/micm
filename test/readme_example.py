import re
import subprocess
import sys
import pandas as pd

README_FILE = "README.md"
CODE_FILE = "foo_chem.cpp"
EXECUTABLE = "./foo_chem"

def extract_code():
    """Extract C++ code from README.md and save it to foo_chem.cpp."""
    with open(README_FILE, "r") as f:
        readme = f.read()
    
    match = re.search(r"```c\+\+\n([\s\S]+?)\n```", readme)
    if not match:
        print("Error: Could not find C++ code block in README.")
        sys.exit(1)
    
    code = match.group(1).strip()
    
    with open(CODE_FILE, "w") as f:
        f.write(code)

def extract_expected_output():
    """Extract expected output from README.md and return as a pandas DataFrame."""
    with open(README_FILE, "r") as f:
        readme = f.read()
    
    match = re.search(r"Output:\n```([\s\S]+?)```", readme)
    if not match:
        print("Error: Could not find expected output in README.")
        sys.exit(1)
    
    output_lines = match.group(1).strip().split("\n")
    return parse_output_to_dataframe(output_lines)

def parse_output_to_dataframe(output_lines):
    """Convert output lines to a pandas DataFrame."""
    # Extract column headers
    headers = output_lines[0].split(",")
    headers = [h.strip() for h in headers]

    # Extract numerical data
    data = []
    for line in output_lines[1:]:
        values = [float(v.strip()) for v in line.split(",")]
        data.append(values)

    return pd.DataFrame(data, columns=headers)

def run_test():
    """Run the test to verify the README example builds, runs, and produces the correct output."""
    extract_code()
    expected_df = extract_expected_output()

    print("Compiling...")
    compile_command = "g++ -o foo_chem foo_chem.cpp -Iinclude -std=c++20"
    subprocess.run(compile_command, shell=True, check=True)

    print("Running executable...")
    result = subprocess.run(EXECUTABLE, shell=True, capture_output=True, text=True)
    
    actual_lines = result.stdout.strip().split("\n")
    actual_df = parse_output_to_dataframe(actual_lines)

    print("Comparing outputs...")
    # Compare numerical values
    differences = (actual_df - expected_df).abs()


    assert(expected_df.shape == actual_df.shape)
    assert(not (differences > 1e-6).any().any())

    print("Test Passed: Output matches expected numerical values.")

if __name__ == "__main__":
    run_test()
