# RepOrtR Post-TAREAN Pipeline

A consolidated and streamlined pipeline for repeat analysis after TAREAN processing. This pipeline combines functionality from multiple redundant scripts into a coherent, modular structure.

## Overview

The original `post_tarean` folder contained many redundant scripts with overlapping functionality. This consolidated pipeline removes redundancy and provides:

- **Unified BLAST analysis** (`blast_consolidated.py`)
- **Main pipeline** (`pipeline.py`) 
- **Utility functions** (`utils.py`)
- **Configuration management** (`config.py`)
- **Comprehensive documentation**

## Files Structure

### Core Pipeline Files

- **`pipeline.py`** - Main pipeline entry point
- **`blast_consolidated.py`** - Consolidated BLAST analysis module
- **`utils.py`** - Common utility functions
- **`config.py`** - Configuration management

### Legacy Files (Redundant)

The following files have been consolidated and are no longer needed:

#### BLAST Variants (consolidated into `blast_consolidated.py`)
- `blast.py`
- `blast_safe.py` 
- `blast_Kroupin.py`
- `blastold.py`
- `super_blast.py`

#### Representation Variants (consolidated into `pipeline.py`)
- `represenation.py`
- `represenation_safe.py`
- `represenation_Kroupin.py`
- `represenation_Kroupin_abomination.py`
- `represenation (copy).py`

#### Word Processing (consolidated into `pipeline.py`)
- `worder.py`
- `worder_new.py`

#### Analysis Scripts (consolidated into `pipeline.py`)
- `annotate_repeatome_structure.py`
- `repeat_sums_by_group.py`
- `types.py`

#### Utility Scripts (consolidated into `utils.py`)
- `rm_homology.py`
- `screenshot_repex.py`
- `katya_csv.py`
- `kate_rep.py`

## Installation

### Prerequisites

```bash
# Required packages
pip install pandas numpy matplotlib

# Optional packages (for enhanced functionality)
pip install xlsxwriter python-docx plotly rapidfuzz pyyaml
```

### BLAST Installation

Ensure BLAST+ is installed and accessible in your PATH:

```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# CentOS/RHEL
sudo yum install blast

# macOS
brew install blast

# Verify installation
blastn -version
```

## Usage

### Basic Usage

```bash
# Run full pipeline for a sample
python pipeline.py SAMPLE_NAME

# Run with custom output directory
python pipeline.py SAMPLE_NAME --output-dir custom_reports

# Run only BLAST analysis
python pipeline.py SAMPLE_NAME --blast-only

# Generate reports only (skip BLAST)
python pipeline.py SAMPLE_NAME --report-only

# Verbose output
python pipeline.py SAMPLE_NAME --verbose
```

### Advanced Usage

```python
from pipeline import RepeatAnalyzer
from config import get_config

# Create custom configuration
config = get_config()
config.update_config(
    blast_only=True,
    output_dir='my_reports'
)

# Initialize analyzer
analyzer = RepeatAnalyzer(config)

# Run analysis
success = analyzer.run_full_pipeline('SAMPLE_NAME')
```

### Configuration

Create a custom configuration file `config.yaml`:

```yaml
blast:
  default_tasks: ['megablast', 'dc-megablast']
  default_threads: 16
  evalue_threshold: 0.05

analysis:
  output_dir: 'custom_reports'
  create_excel: true
  create_word: true

filtering:
  levenshtein_threshold: 85.0
  remove_self_blasts: true
```

Load custom configuration:

```python
from config import create_config_from_file

config = create_config_from_file('config.yaml')
analyzer = RepeatAnalyzer(config)
```

## Pipeline Components

### 1. BLAST Analysis (`blast_consolidated.py`)

**Features:**
- Automatic database discovery
- Multiple BLAST tasks (megablast, dc-megablast, blastn)
- X3 sequence variants support
- Coverage calculation and classification
- Levenshtein distance filtering
- Priority-based result sorting

**Usage:**
```python
from blast_consolidated import BLASTAnalyzer

analyzer = BLASTAnalyzer()
results = analyzer.run_blast(
    subjects=['SAMPLE1', 'SAMPLE2'],
    dbs=['local', 'ncbi'],
    tasks=['megablast', 'dc-megablast'],
    num_threads=20
)
```

### 2. Data Parsing (`utils.py`)

**Functions:**
- `parse_tarean_data()` - Parse TAREAN output
- `parse_copy_data()` - Parse copy number data
- `parse_comparative_data()` - Parse comparative analysis
- `get_cyphers()` - Load cypher mappings
- `get_ncbi_naming()` - Load NCBI naming data

### 3. Report Generation (`pipeline.py`)

**Report Types:**
- Excel reports with multiple sheets
- Word documents with formatted text
- CSV data exports
- Summary statistics

### 4. Configuration Management (`config.py`)

**Features:**
- YAML configuration files
- Validation of settings
- Default configurations
- Runtime configuration updates

## Database Requirements

The pipeline expects the following database structure:

```
.
├── ncbi_repeats_db/
│   ├── ncbi_repeats.fasta
│   ├── ncbi_repeats_x3.fasta
│   └── *.nhr, *.nin, *.nsq (BLAST index files)
├── local_db_solo/
│   ├── multifasta.fasta
│   ├── multifasta_x3.fasta
│   └── *.nhr, *.nin, *.nsq
├── comparatives_db/
│   ├── COMPBASE.fasta
│   ├── COMPBASE_x3.fasta
│   └── *.nhr, *.nin, *.nsq
└── important_db/
    ├── reference.fasta
    └── *.nhr, *.nin, *.nsq
```

## Input Data Structure

For each sample, the pipeline expects:

```
SAMPLE_NAME/
├── *TAREAN* (FASTA files with repeat sequences)
├── CLUSTER_TABLE.csv
├── COMPARATIVE_ANALYSIS_COUNTS.csv (optional)
└── seqclust/clustering/clusters/dir_*/graph_layout.png
```

## Output Structure

```
reports/
├── SAMPLE_NAME_repeat_analysis.xlsx
├── SAMPLE_NAME_repeat_analysis.docx
├── SAMPLE_NAME_repeat_analysis.csv
└── SAMPLE_NAME_blast_results.csv (if --blast-only)
```

## Error Handling

The pipeline includes comprehensive error handling:

- **File validation** - Checks for required input files
- **Database validation** - Verifies BLAST database availability
- **Graceful degradation** - Continues processing when optional components fail
- **Detailed logging** - Provides informative error messages
- **Retry mechanisms** - Automatically retries failed operations

## Performance Optimization

- **Parallel processing** - Multi-threaded BLAST analysis
- **Memory management** - Efficient data handling for large datasets
- **Chunked processing** - Processes data in manageable chunks
- **Temporary file cleanup** - Automatically removes temporary files

## Migration from Legacy Scripts

### From `blast.py` variants:
```python
# Old way
from blast import blaster
results = blaster(subjects=['SAMPLE'], dbs=['local'])

# New way
from blast_consolidated import BLASTAnalyzer
analyzer = BLASTAnalyzer()
results = analyzer.run_blast(subjects=['SAMPLE'], dbs=['local'])
```

### From `represenation.py` variants:
```python
# Old way
python represenation.py SAMPLE_NAME

# New way
python pipeline.py SAMPLE_NAME
```

### From `worder.py` variants:
```python
# Old way
python worder.py SAMPLE_NAME local_or_comp

# New way
python pipeline.py SAMPLE_NAME --output-dir reports
```

## Troubleshooting

### Common Issues

1. **BLAST not found**
   ```bash
   # Check BLAST installation
   which blastn
   blastn -version
   ```

2. **Database not found**
   ```bash
   # Check database files
   ls -la ncbi_repeats_db/
   ls -la local_db_solo/
   ```

3. **Permission errors**
   ```bash
   # Check file permissions
   chmod +x pipeline.py
   chmod -R 755 reports/
   ```

4. **Memory issues**
   ```python
   # Reduce thread count
   config.update_config(blast={'default_threads': 8})
   ```

### Debug Mode

```bash
# Enable debug logging
python pipeline.py SAMPLE_NAME --verbose

# Check configuration
python -c "from config import get_config; print(get_config()._config_to_dict())"
```

## Contributing

When adding new functionality:

1. **Use the modular structure** - Add functions to appropriate modules
2. **Update configuration** - Add new settings to `config.py`
3. **Add documentation** - Update this README
4. **Include tests** - Add test cases for new functionality
5. **Follow naming conventions** - Use descriptive function and variable names

## License

This pipeline is part of the RepOrtR project. Please refer to the main project license.

## Support

For issues and questions:

1. Check the troubleshooting section
2. Review the configuration options
3. Enable verbose logging for detailed error messages
4. Check the input data structure requirements 