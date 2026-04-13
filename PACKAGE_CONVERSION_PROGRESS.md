# HybPhaser Package Conversion - Progress Report

## Completed Tasks

### 1. Package Infrastructure ✓
- Created `DESCRIPTION` file with proper metadata
- Set up `NAMESPACE` file for roxygen2
- Created `.Rbuildignore` to exclude original scripts
- Set up test infrastructure with testthat

### 2. Directory Structure ✓
```
HybPhaser/
├── R/                          # Package functions
│   ├── HybPhaser-package.R    # Package documentation
│   ├── utils.R                # Utility functions
│   ├── count_snps.R           # SNP counting (from 1a_count_snps.R)
│   └── data_access.R          # Access to scripts and examples
├── inst/
│   ├── scripts/               # Bash scripts
│   │   ├── 1_generate_consensus_sequences.sh
│   │   ├── 1_generate_consensus_sequences_single-core.sh
│   │   └── 2_extract_mapped_reads.sh
│   └── extdata/               # Example data
│       ├── clade_references.csv
│       └── phasing_prep.csv
├── tests/
│   └── testthat/             # Unit tests
│       ├── test-utils.R
│       ├── test-count_snps.R
│       └── test-data_access.R
├── man/                       # Documentation (generated)
└── vignettes/                # Tutorials (to be created)
```

### 3. Functions Created ✓

#### Utility Functions (`R/utils.R`)
- `read_config()` - Read and validate configuration files
- `validate_paths()` - Check file/directory existence and type
- `create_output_dirs()` - Create output directory structure

#### SNP Analysis (`R/count_snps.R`)
- `count_snps()` - Main function to count SNPs in consensus sequences
  - Replaces `1a_count_snps.R` script
  - Returns both SNP proportions and sequence lengths
  - Saves results as .Rds files
- `.seq_stats()` - Internal helper for sequence statistics

#### Data Access (`R/data_access.R`)
- `hybphaser_scripts()` - Get path to bash scripts
- `hybphaser_example()` - Get path to example data
- `run_hybphaser_script()` - Execute bash scripts from R

### 4. Tests Written ✓
All functions have comprehensive unit tests using testthat:
- `test-utils.R` - Tests for utility functions
- `test-count_snps.R` - Tests for SNP counting
- `test-data_access.R` - Tests for data access functions

## Next Steps

### Remaining R Scripts to Convert
1. `1b_assess_dataset.R` - Dataset optimization and assessment
2. `1c_generate_sequence_lists.R` - Generate sequence lists
3. `2a_prepare_bbsplit_script.R` - BBSplit clade association prep
4. `2b_collate_bbsplit_results.R` - Collate BBSplit results
5. `3a_prepare_phasing_script.R` - Phasing script preparation
6. `3b_collate_phasing_stats.R` - Collate phasing statistics
7. `4_merge_sequence_lists.R` - Merge phased with normal sequences

### Documentation
- Generate roxygen2 documentation
- Update README for package usage
- Create vignette showing complete workflow

### Validation
- Run `devtools::check()`
- Run all tests
- Verify bash scripts are accessible

## Testing the Package

To test the current state:

```r
# Load development package
devtools::load_all()

# Generate documentation
devtools::document()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Example Usage (Current Functionality)

```r
library(HybPhaser)

# Count SNPs in consensus sequences
results <- count_snps(
  path_to_output_folder = "path/to/hybphaser_output",
  fasta_file_with_targets = "path/to/targets.fasta",
  targets_file_format = "DNA",
  path_to_namelist = "path/to/namelist.txt",
  intronerated_contig = FALSE
)

# Access bash scripts
script_path <- hybphaser_scripts("1_generate_consensus_sequences.sh")

# Run bash script
run_hybphaser_script(
  "1_generate_consensus_sequences.sh",
  args = c("-n", "namelist.txt", "-p", "hybpiper_output"),
  dry_run = TRUE  # Preview command
)

# Access example data
clade_refs <- read.csv(hybphaser_example("clade_references.csv"))
```
