import re
import subprocess
import pandas as pd
import pytest

README_FILE = "README.md"
CODE_FILE = "foo_chem.cpp"
EXECUTABLE = "./foo_chem"


def extract_code(readme_content):
    """Extract C++ code from README content."""
    match = re.search(r"```c\+\+\n([\s\S]+?)\n```", readme_content)
    assert match, "Could not find C++ code block in README."
    return match.group(1).strip()


def extract_expected_output(readme_content):
    """Extract expected output from README content and return as a pandas DataFrame."""
    match = re.search(r"Output:\n```([\s\S]+?)```", readme_content)
    assert match, "Could not find expected output in README."
    output_lines = match.group(1).strip().split("\n")
    return parse_output_to_dataframe(output_lines)


def parse_output_to_dataframe(output_lines):
    """Convert output lines to a pandas DataFrame."""
    headers = [h.strip() for h in output_lines[0].split(",")]
    data = [[float(v.strip()) for v in line.split(",")]
            for line in output_lines[1:]]
    return pd.DataFrame(data, columns=headers)


@pytest.fixture
def readme_content():
    with open(README_FILE, "r") as f:
        return f.read()


@pytest.fixture
def test_execution(tmp_path, readme_content):
    """Compile and run the extracted C++ code, then compare output."""
    extracted_code = extract_code(readme_content)
    expected_output_df = extract_expected_output(readme_content)

    code_file = tmp_path / CODE_FILE
    executable = tmp_path / "foo_chem"

    with open(code_file, "w") as f:
        f.write(extracted_code)

    compile_command = f"g++ -o {executable} {code_file} -Iinclude -std=c++20"
    subprocess.run(compile_command, shell=True, check=True)

    result = subprocess.run(str(executable), shell=True,
                            capture_output=True, text=True)
    actual_lines = result.stdout.strip().split("\n")
    actual_output_df = parse_output_to_dataframe(actual_lines)

    assert actual_output_df.shape == expected_output_df.shape, "Output shape mismatch."
    assert actual_output_df.columns.equals(expected_output_df.columns), "Output columns do not match expected columns."
    assert actual_output_df.equals(expected_output_df), "Output values do not match expected values."


def test_readme_example(test_execution):
    """Run the test_execution fixture to validate the README example."""
    pass  # The actual test logic is handled in the fixture.
