# RepOrtR Testing Guide

This guide covers the comprehensive testing framework for the RepOrtR pipeline.

## Overview

The testing framework validates both the main pipeline (`reportr` environment) and the RepeatExplorer integration (`repeatexplorer` environment).

## Test Architecture

```
Test Runner (reportr env)
    ↓
Unit Tests (reportr env)
    ↓
Integration Tests (reportr env)
    ↓
RepeatExplorer Tests (repeatexplorer env via wrapper)
```

## Test Organization (P0-aligned)

### Root Level
- `run_tests.py` - Main test runner; runs consolidated suite and prints dependency status.

### tests/ Directory
- `conftest.py` - Pytest setup (repo root on sys.path, cwd reset before each test).
- `test_consolidated.py` - File existence, config load, dependencies, prepareReadsRE, seqclust; recommended CLI entrypoint.
- `test_pipeline_integration_modular.py` - Snakefile_modular, workflows/*.smk, config validation rules.
- `test_post_tarean_pipeline.py` - Post-TAREAN pipeline and quality-gating step (used in Docker full test).
- `test_multi_project.py` - Multi-project and ProjectManager (project creation, sample handling, caching).
- `run_snakemake_workflow_modular.py` - Optional modular workflow runner (e.g. `--validate-only`).

## Running Tests

### Recommended
```bash
conda activate reportr
python3 run_tests.py
```

### Pytest (full suite, CI / container)
```bash
conda activate reportr
pytest -q tests
```

### Specific Suites
```bash
# Consolidated suite (CLI)
python3 tests/test_consolidated.py

# Modular workflow integration
pytest -q tests/test_pipeline_integration_modular.py

# Post-TAREAN pipeline
pytest -q tests/test_post_tarean_pipeline.py

# Multi-project
pytest -q tests/test_multi_project.py
```

## Test Categories

### 1. Consolidated (`test_consolidated.py`)
- **Purpose**: One-shot health check (files, config, deps, prepareReadsRE, seqclust).
- **Environment**: reportr.
- **Run**: `python3 run_tests.py` or `python3 tests/test_consolidated.py`.

### 2. Modular workflow (`test_pipeline_integration_modular.py`)
- **Purpose**: Validate Snakefile_modular, config validation rules, dry-run.
- **Environment**: reportr.

### 3. Post-TAREAN pipeline (`test_post_tarean_pipeline.py`)
- **Purpose**: Orchestrator, steps (blast, summary, quality_gating), JSON output.
- **Environment**: reportr; used in Docker full test.

### 4. Multi-project (`test_multi_project.py`)
- **Purpose**: ProjectManager, workflow config, directory layout, caching.
- **Environment**: reportr.

## Test Data

### Input Files
- `data/reads/test_R1.fq` - Test forward reads (10,000 reads)
- `data/reads/test_R2.fq` - Test reverse reads (10,000 reads)

### Generated Files
- `tests/cleaned_reads/` - Simulated cleaned reads
- `tests/interlaced/` - Mock prepared reads
- `tests/interlaced/test_sample_RepExRES/` - Mock clustering output

## Test Results

### Expected Output
```
Tests passed: 4/6
Success rate: 66.7%

Detailed Results:
  File Existence            ✅ PASS
  Configuration Loading     ✅ PASS
  Dependencies              ✅ PASS
  Test Data Quality         ❌ FAIL (Biopython missing)
  Prepare Reads Script      ❌ FAIL (Biopython missing)
  RepeatExplorer Clustering ✅ PASS
```

### Dependency Status
- **Biopython**: Required for read validation
- **PyYAML**: Required for config parsing
- **FastQC**: Optional for quality control
- **BBDuk**: Optional for read cleaning
- **Snakemake**: Optional for workflow execution

## Environment-Specific Testing

### reportr Environment Tests
```bash
# Activate environment
conda activate reportr

# Run tests
python3 run_tests.py
```

### repeatexplorer Environment Tests
```bash
# Test seqclust inside the repeatexplorer env (in-repo binary)
conda run -n repeatexplorer ./repex_tarean/seqclust --help
```

### Cross-Environment Tests
```bash
# From reportr environment
conda activate reportr

# Test seqclust via conda-run
conda run -n repeatexplorer ./repex_tarean/seqclust --help
```

## Troubleshooting

### Common Test Failures

1. **"Biopython not available"**
   ```bash
   conda activate reportr
   conda install -c conda-forge biopython
   ```

2. **"seqclust not found"**
   ```bash
   # Check if repeatexplorer environment exists
   conda env list | grep repeatexplorer
   
   # Reinstall if missing
   python3 install_reportr.py
   ```

3. **"Test data not found"**
   ```bash
   # Check test data exists
   ls data/reads/test_R*.fq
   
   # Create test data if missing
   python3 tests/create_test_data.py
   ```

### Test Environment Issues

1. **Wrong environment active**
   ```bash
   echo $CONDA_DEFAULT_ENV  # Should show "reportr"
   conda activate reportr
   ```

2. **Missing dependencies**
   ```bash
   conda list | grep -E "(biopython|pyyaml|snakemake)"
   ```

3. **Path issues**
   ```bash
   # Check in-repo binary exists (built from repex_tarean)
   ls -la ./repex_tarean/seqclust
   ```

## TAREAN Runs: Logs and Stable Launch

### Logging layers

- **Sample-level logs**: `logs/seqclust_{project}_{sample}.log`
  - Primary place to inspect TAREAN/seqclust behaviour.
  - Captures seqclust stdout/stderr and R/pyRserve stack traces.
  - Includes structured markers such as `TAREAN_SUCCESS`, `TAREAN_WARNING`, and `TAREAN_ERROR` with exit codes.
- **Workflow-level logs**: `.snakemake/log/<timestamp>.snakemake.log`
  - Snakemake’s own log of rule execution and `RuleException` summaries.
  - Points back to the relevant `logs/seqclust_{project}_{sample}.log` for details.
- **Optional progress tracking**: when `global.progress_tracking.enabled: true` in
  `[projects/global_config.yaml](/root/development/RepOrtR/projects/global_config.yaml)`,
  additional progress artifacts are written under `logs/progress_{project}_{sample}.*`.

### Stable launch for many assemblies

- **Recommended command** for large projects:

  ```bash
  snakemake -s Snakefile_modular \
    --configfile projects/global_config.yaml \
    --cores N \
    --keep-going \
    all
  ```

  - `--keep-going` ensures that one failing sample does not stop independent ones.
  - The `all` target (see `workflows/report_rules.smk`) drives project-level completion.

- **Resuming after a failure**:
  - Fix configuration or input for the failing sample.
  - Rerun the same Snakemake command; completed jobs are skipped based on existing outputs and tokens.

### TAREAN core outputs vs late reporting

- RepeatExplorer/TAREAN may fail late in reporting (e.g. supercluster reports) even when core outputs are already present.
- The `run_tarean` rule in
  `[workflows/core_rules.smk](/root/development/RepOrtR/workflows/core_rules.smk)`:
  - Treats the presence of core TAREAN outputs (key DB files and `TAREAN_consensus_rank_*.fasta`)
    as a **success with warnings** when the final reporting step fails.
  - Writes `TAREAN_COMPLETE_WITH_WARNINGS` to the per-sample `tarean.done` token and logs a
    `TAREAN_WARNING` line with the non-zero exit code.
  - Allows downstream RepOrtR post-TAREAN steps to proceed without rerunning the unstable
    reporting part, while still surfacing the underlying error in the seqclust log.

## Test Development

### Adding New Tests

1. **Create test file** in `tests/` with naming `test_*.py` for pytest discovery.
2. Use repo root in path (conftest sets cwd and sys.path).
3. For CLI-style suites, add to the consolidated runner or keep as pytest-only.

### Test Style Guide

- Use descriptive test names
- Include setup and teardown
- Handle missing dependencies gracefully
- Provide clear error messages
- Test both success and failure cases

### Example Test Structure

```python
def test_example():
    """Test description"""
    # Setup
    setup_test_environment()
    
    # Test
    result = run_test_function()
    
    # Assert
    assert result is True
    
    # Cleanup
    cleanup_test_files()
```

## Performance Testing

### Benchmark Tests
- File processing speed
- Memory usage
- CPU utilization
- Disk I/O performance

### Load Testing
- Large file processing
- Multiple concurrent operations
- Resource exhaustion scenarios

## Continuous Integration

### Automated Testing
- Run tests on every commit
- Validate both environments
- Check cross-environment compatibility
- Report test coverage

### Test Reports
- Generate HTML reports
- Track test history
- Identify performance regressions
- Monitor dependency changes

## Support

For test-related issues:
1. Check troubleshooting section
2. Run individual test files
3. Verify environment setup
4. Check test data integrity
5. Review test logs for details 