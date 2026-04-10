# Modular Workflow Structure

This directory contains the modularized Snakemake workflow components for better maintainability and stability.

## Structure

### `core_rules.smk` - STABLE COMPONENTS
**Purpose**: Core pipeline functionality that rarely changes
**Contains**:
- Read cleaning and preparation rules
- TAREAN execution rules
- Project validation rules

**Maintenance**: These rules are well-tested and should only be modified for critical bug fixes.

### `test_rules.smk` - SYSTEM VALIDATION
**Purpose**: System testing and validation before running full analyses
**Contains**:
- System resource testing
- Read preparation testing
- TAREAN testing with small datasets
- Comparative analysis testing

**Maintenance**: Test rules can be enhanced for better validation but core logic should remain stable.

### `analysis_rules.smk` - ANALYSIS PIPELINES
**Purpose**: Post-TAREAN analysis and comparative analysis
**Contains**:
- Post-TAREAN BLAST analysis
- Comparative read preparation
- Comparative TAREAN analysis
- Quality control (FastQC)

**Maintenance**: This module can be extended with new analysis types while keeping existing functionality stable.

### `report_rules.smk` - REPORTING AND AGGREGATION
**Purpose**: Final reporting and project aggregation
**Contains**:
- Project summary generation
- Final report aggregation
- Utility rules (cleaning, validation)

**Maintenance**: Report formats can be enhanced, but aggregation logic should remain stable.

## Benefits of Modular Structure

### 1. **Stability Protection**
- Core rules are isolated and protected from experimental changes
- Test rules ensure system stability before running expensive analyses
- Each module can be version-controlled independently

### 2. **Maintainability**
- Clear separation of concerns
- Easier to locate and fix issues
- Reduced risk of breaking stable components

### 3. **Extensibility**
- New analysis types can be added to `analysis_rules.smk`
- New test types can be added to `test_rules.smk`
- Report formats can be enhanced in `report_rules.smk`

### 4. **Development Workflow**
- Experimental changes can be made in specific modules
- Stable components remain untouched during development
- Easy rollback of problematic changes

### 5. **Configuration Compliance**
- All values are dynamically loaded from configuration files
- No hardcoded sample names, file paths, or parameters
- Automatic validation of configuration completeness
- Tool path validation ensures FastQC and BBDuk write to correct directories

## Usage

### Running the Modular Workflow
```bash
# Use the modular Snakefile
snakemake -s Snakefile_modular --cores 4

# Or include specific modules in a custom workflow
include: "workflows/core_rules.smk"
include: "workflows/analysis_rules.smk"
```

### Development Guidelines

1. **Core Rules (`core_rules.smk`)**
   - Only modify for critical bug fixes
   - Extensive testing required for any changes
   - Document all changes thoroughly

2. **Test Rules (`test_rules.smk`)**
   - Can be enhanced for better validation
   - Add new test types as needed
   - Ensure tests are comprehensive

3. **Analysis Rules (`analysis_rules.smk`)**
   - Can be extended with new analysis types
   - Keep existing functionality stable
   - Add new analysis methods here

4. **Report Rules (`report_rules.smk`)**
   - Can enhance report formats
   - Add new aggregation methods
   - Keep core aggregation logic stable

## Migration from Monolithic Snakefile

The original `Snakefile` contains all rules in one file. The modular structure provides:

- **Better organization**: Related rules are grouped together
- **Easier maintenance**: Changes are isolated to specific modules
- **Improved testing**: Each module can be tested independently
- **Enhanced stability**: Core components are protected from experimental changes

## Stability Matrix

| Component | Stability Level | Change Frequency | Testing Required |
|-----------|----------------|------------------|------------------|
| Core Rules | High | Rare | Extensive |
| Test Rules | Medium | Occasional | Moderate |
| Analysis Rules | Medium | Regular | Moderate |
| Report Rules | Low | Frequent | Light |
| Config Validation | High | Rare | Extensive |

This structure ensures that the most critical components (core pipeline functionality) remain stable while allowing flexibility for enhancements and new features.

## Configuration Compliance

### **No Hardcoded Values**
- All sample names are dynamically loaded from `config["projects"][project]["samples"]`
- All file paths use wildcards: `{project}`, `{sample}`
- All parameters are loaded from configuration files
- Comparative species are dynamically determined from `config["projects"][project]["comparative_species"]`

### **Tool Validation**
- FastQC writes to `projects/{project}/samples/{sample}/fastqc/`
- BBDuk writes to `projects/{project}/samples/{sample}/filtered_reads/`
- All tool paths are validated before workflow execution
- Tool availability is checked during configuration validation

### **Configuration Validation**
- Required configuration keys are validated
- Project-specific settings are checked
- TAREAN parameters are verified
- Global settings are confirmed 