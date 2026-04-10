# RepOrtR Test Style Guide

## Overview
This document defines the coding standards, formatting conventions, and best practices for all test files in the RepOrtR pipeline.

## Core Principles

### 1. **Consistency First**
- All tests use the same formatting symbols and style
- Consistent error message structure
- Uniform section headers and output format

### 2. **User-Friendly Error Messages**
- Every error includes actionable "Install with" commands
- **All Python library installations use conda** (not pip)
- Clear, specific installation instructions
- Helpful recommendations for troubleshooting

### 3. **Comprehensive Coverage**
- Test each component exactly once
- No redundant functionality
- Essential functionality prioritized over edge cases

### 4. **Clear Organization**
- Logical test progression
- Consistent section structure
- Easy-to-read output format

### 5. **Conda-First Installation**
- **All Python libraries must be installed via conda**
- Use `conda install -c conda-forge` for Python packages
- Use `conda install -c bioconda` for bioinformatics tools
- Never recommend `pip install` for Python libraries
- This ensures consistent environments and dependency management

## Formatting Standards

### Symbols and Icons
```python
# Consistent formatting constants
SUCCESS = "✅"
ERROR = "❌"
WARNING = "⚠️"
INFO = "ℹ️"
```

### Header Formatting
```python
def print_header(title):
    """Print consistent header"""
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")

def print_section(title):
    """Print consistent section header"""
    print(f"\n--- {title} ---")
```

### Error Message Structure
```python
def check_dependency(name, check_func, install_cmd=None):
    """Check dependency with consistent formatting"""
    try:
        result = check_func()
        if result:
            print(f"{SUCCESS} {name} is available")
            return True
        else:
            print(f"{ERROR} {name} not working properly")
            if install_cmd:
                print(f"   Install with: {install_cmd}")
            return False
    except Exception as e:
        print(f"{ERROR} {name} not available: {e}")
        if install_cmd:
            print(f"   Install with: {install_cmd}")
        return False
```

## Test Structure Guidelines

### 1. **File Organization**
```python
#!/usr/bin/env python3
"""
Brief description of test purpose
Reference: See TEST_STYLE_GUIDE.md for formatting standards
"""

import os
import sys
import subprocess
import time
import yaml
import shutil
from pathlib import Path

# Consistent formatting constants
SUCCESS = "✅"
ERROR = "❌"
WARNING = "⚠️"
INFO = "ℹ️"

# Helper functions
def print_header(title):
    """Print consistent header"""
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")

def print_section(title):
    """Print consistent section header"""
    print(f"\n--- {title} ---")

# Test functions
def test_component_name():
    """Test description"""
    print_section("Component Name")
    # Test logic here
    return True/False

# Main function
def main():
    """Run all tests"""
    print_header("Test Suite Name")
    # Test execution logic
    return success_status
```

### 2. **Test Function Naming**
- Use descriptive names: `test_file_existence()`, `test_config_loading()`
- Follow snake_case convention
- Include component name in function name

### 3. **Section Organization**
Each test should follow this structure:
1. **Section Header** - Clear description of what's being tested
2. **Prerequisites** - Check if required files/tools exist
3. **Test Logic** - Actual testing code
4. **Results** - Clear pass/fail with details
5. **Cleanup** - Remove temporary files/directories

## Error Handling Standards

### 1. **Dependency Checks**
```python
# Always provide installation commands
# All Python libraries use conda (not pip)
dependencies = [
    ("Biopython", lambda: __import__("Bio"), "conda install -c conda-forge biopython"),
    ("PyYAML", lambda: __import__("yaml"), "conda install -c conda-forge pyyaml"),
    ("FastQC", lambda: subprocess.run(["fastqc", "--version"], 
                                     capture_output=True, timeout=10).returncode == 0,
     "conda install -c bioconda fastqc"),
]
```

### 2. **File Existence Checks**
```python
required_files = [
    ("Snakefile", "Snakefile"),
    ("config.yaml", "config.yaml"),
    ("prepareReadsRE.py", "prepareReadsRE.py"),
]

missing_files = []
for name, path in required_files:
    if os.path.exists(path):
        print(f"{SUCCESS} {name} exists")
    else:
        print(f"{ERROR} {name} not found: {path}")
        missing_files.append(name)
```

### 3. **Subprocess Error Handling**
```python
try:
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    
    if result.returncode != 0:
        print(f"{ERROR} Command failed")
        if "No module named 'Bio'" in result.stderr:
            print(f"   Install with: conda install -c conda-forge biopython")
        elif "No module named 'yaml'" in result.stderr:
            print(f"   Install with: conda install -c conda-forge pyyaml")
        else:
            print(f"   Error: {result.stderr}")
        return False
    
    print(f"{SUCCESS} Command completed successfully")
    return True
    
except subprocess.TimeoutExpired:
    print(f"{ERROR} Command timed out")
    return False
except Exception as e:
    print(f"{ERROR} Command error: {e}")
    return False
```

## Output Formatting Standards

### 1. **Success Messages**
- Use ✅ symbol
- Include relevant details (file sizes, counts, etc.)
- Be specific about what passed

### 2. **Error Messages**
- Use ❌ symbol
- Include specific error details
- Always provide "Install with" command when applicable
- Suggest alternative solutions when possible

### 3. **Warning Messages**
- Use ⚠️ symbol
- For non-critical issues
- Provide context about why it's a warning

### 4. **Information Messages**
- Use ℹ️ symbol
- For general information
- Provide helpful context

## Summary and Reporting

### 1. **Test Summary Format**
```python
print_header("Test Summary")
print(f"Tests passed: {total_passed}/{len(tests)}")
print(f"Success rate: {(total_passed/len(tests))*100:.1f}%")

print(f"\nDetailed Results:")
for test_name, passed in results.items():
    status = f"{SUCCESS} PASS" if passed else f"{ERROR} FAIL"
    print(f"  {test_name:<25} {status}")
```

### 2. **Dependency Status Report**
```python
print(f"\nDependency Status:")
for dep, available in deps.items():
    status = f"{SUCCESS} Available" if available else f"{ERROR} Missing"
    print(f"  {dep:<15} {status}")
```

### 3. **Recommendations Section**
```python
print(f"\nRecommendations:")
if not deps.get("Biopython", False):
    print(f"  - Install Biopython: conda install -c conda-forge biopython")
if not deps.get("PyYAML", False):
    print(f"  - Install PyYAML: conda install -c conda-forge pyyaml")
if not deps.get("Snakemake", False):
    print(f"  - Install Snakemake: conda install -c conda-forge snakemake")
```

## Best Practices

### 1. **Test Independence**
- Each test should be independent
- Clean up after each test
- Don't rely on previous test results

### 2. **Timeout Handling**
- Always set reasonable timeouts for subprocess calls
- Handle timeout exceptions gracefully
- Provide clear timeout error messages

### 3. **Resource Management**
- Use temporary directories for test files
- Clean up created files and directories
- Handle file permissions appropriately

### 4. **Error Recovery**
- Provide fallback options when possible
- Don't fail completely if optional components are missing
- Give users clear next steps

### 5. **Documentation**
- Include docstrings for all functions
- Document test purpose and expected behavior
- Reference this style guide in file headers

## Example Implementation

See `test_consolidated.py` for a complete implementation of these standards.

## Maintenance

- Update this guide when adding new formatting standards
- Ensure all test files follow these guidelines
- Review and update as the codebase evolves 