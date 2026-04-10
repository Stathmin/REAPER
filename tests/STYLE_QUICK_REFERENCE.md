# Test Style Quick Reference

## Symbols
```python
SUCCESS = "✅"  # Pass/success
ERROR = "❌"    # Fail/error
WARNING = "⚠️"  # Warning
INFO = "ℹ️"     # Information
```

## Header Functions
```python
def print_header(title):
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")

def print_section(title):
    print(f"\n--- {title} ---")
```

## Error Message Template
```python
print(f"{ERROR} Component not available: {error_details}")
print(f"   Install with: {install_command}")
```

## Conda Installation Commands
```python
# Python libraries (use conda-forge)
"conda install -c conda-forge biopython"
"conda install -c conda-forge pyyaml"
"conda install -c conda-forge snakemake"

# Bioinformatics tools (use bioconda)
"conda install -c bioconda fastqc"
"conda install -c bioconda bbmap"
```

## Test Structure
```python
def test_component():
    """Test description"""
    print_section("Component Name")
    
    # Prerequisites
    if not os.path.exists(required_file):
        print(f"{ERROR} Required file not found")
        return False
    
    # Test logic
    try:
        # Test code here
        print(f"{SUCCESS} Component working")
        return True
    except Exception as e:
        print(f"{ERROR} Test failed: {e}")
        print(f"   Install with: {install_cmd}")
        return False
```

## File Header Template
```python
#!/usr/bin/env python3
"""
Brief description of test purpose
Reference: See TEST_STYLE_GUIDE.md for formatting standards
"""
```

## Summary Template
```python
print_header("Test Summary")
print(f"Tests passed: {total_passed}/{len(tests)}")
print(f"Success rate: {(total_passed/len(tests))*100:.1f}%")

print(f"\nDetailed Results:")
for test_name, passed in results.items():
    status = f"{SUCCESS} PASS" if passed else f"{ERROR} FAIL"
    print(f"  {test_name:<25} {status}")
```

## Key Principles
1. **Consistency** - Use same symbols and format everywhere
2. **User-Friendly** - Always provide "Install with" commands
3. **Conda-First** - All Python libraries use conda (not pip)
4. **Comprehensive** - Test each component once, no redundancy
5. **Clear** - Easy-to-read output with logical structure

See `TEST_STYLE_GUIDE.md` for complete documentation. 