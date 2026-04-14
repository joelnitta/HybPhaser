# Docker Integration for HybPhaser

## Summary

Added Docker wrapper functions to enable running complete HybPhaser 
workflows from R without manually managing command-line tools.

## New Functions

### Docker Utilities (`R/docker.R`)

1. **`check_docker()`** - Verify Docker is installed and running
2. **`check_docker_image()`** - Check if HybPhaser image exists, 
   optionally pull it
3. **`.normalize_docker_path()`** - Internal function to handle path 
   translation for volume mounting (works on Unix, macOS, and Windows)
4. **`.run_docker()`** - Internal function to execute commands in Docker 
   containers

### Workflow Functions (`R/docker_workflows.R`)

1. **`run_generate_consensus_sequences()`** - Generate consensus 
   sequences from HybPiper output
   - Runs `1_generate_consensus_sequences.sh` in Docker
   - Automatically mounts directories
   - Configurable parameters for variant calling
   - Supports parallel processing

2. **`run_extract_mapped_reads()`** - Extract reads that mapped to 
   targets
   - Runs `2_extract_mapped_reads.sh` in Docker
   - Prepares data for clade association
   - Optional duplicate removal

## Features

- **Cross-platform**: Works on macOS, Linux, and Windows
- **Automatic setup**: Creates directories and handles path mapping
- **Error handling**: Validates inputs and provides clear error messages
- **Flexible**: Can use Docker or run scripts directly
- **Well-tested**: Comprehensive test suite (53 tests passing)
- **Well-documented**: Full roxygen2 documentation for all functions

## Usage Example

```r
library(HybPhaser)

# Check Docker setup
check_docker()
check_docker_image("joelnitta/hybphaser:latest", pull = TRUE)

# Run complete workflow
run_generate_consensus_sequences(
  hybpiper_dir = "path/to/hybpiper_output",
  output_dir = "path/to/hybphaser_output",
  namelist = "samples.txt",
  threads = 8
)

# Count SNPs in consensus sequences
results <- count_snps(
  path_to_output_folder = "path/to/hybphaser_output",
  fasta_file_with_targets = "targets.fasta",
  targets_file_format = "DNA",
  path_to_namelist = "samples.txt"
)

# Extract mapped reads
run_extract_mapped_reads(
  base_dir = "path/to/hybphaser_output",
  output_dir = "path/to/mapped_reads"
)
```

## Benefits

1. **Single interface**: Everything runs from R
2. **No manual installation**: Docker image includes all dependencies
3. **Reproducible**: Same environment across all systems
4. **Integrated workflow**: Seamlessly combine Docker execution with R 
   analysis
5. **User-friendly**: Automatic directory creation and path handling

## Docker Image

Uses `joelnitta/hybphaser:latest` which includes:
- BWA (read mapping)
- SAMtools (alignment processing)
- Bcftools (variant calling)
- BBMap/BBSplit (clade association and phasing)
- R and required packages
- All HybPhaser bash scripts

## Testing

All Docker functions tested with:
- Input validation
- Path handling
- Directory creation
- Docker image checking
- Actual Docker execution (when image available)

Test coverage: 11 dedicated Docker tests + integration tests

## Vignette Updates

The vignette has been updated to include:
- Docker setup instructions
- Complete end-to-end workflow example
- Comparison of Docker vs. local execution
- Docker troubleshooting guidance

## Files Created/Modified

### New Files
- `R/docker.R` - Docker utility functions
- `R/docker_workflows.R` - Workflow wrapper functions
- `tests/testthat/test-docker.R` - Docker function tests

### Modified Files
- `vignettes/HybPhaser.Rmd` - Added Docker workflow examples
- `NAMESPACE` - Exported new functions
- `man/` - Generated documentation for new functions

## Next Steps

Additional workflow functions could be added for:
- BBSplit clade association
- Phasing operations
- Sequence list merging

These would follow the same pattern as the existing Docker wrappers.
