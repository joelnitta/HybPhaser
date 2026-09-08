# Get Path to HybPhaser Bash Scripts

Returns the full path to bash scripts included with the package

## Usage

``` r
hybphaser_scripts(script_name = NULL)
```

## Arguments

- script_name:

  Name of the script. If NULL, returns the scripts directory. One of:

  - "1_generate_consensus_sequences.sh"

  - "1_generate_consensus_sequences_single-core.sh"

  - "2_extract_mapped_reads.sh"

## Value

Character path to the script or scripts directory

## Examples

``` r
# Get path to scripts directory
script_dir <- hybphaser_scripts()

# Get path to specific script
consensus_script <- hybphaser_scripts(
  "1_generate_consensus_sequences.sh"
)
```
