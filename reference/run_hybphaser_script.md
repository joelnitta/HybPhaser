# Run HybPhaser Bash Script

Execute one of the HybPhaser bash scripts with specified arguments

## Usage

``` r
run_hybphaser_script(script_name, args = character(0), dry_run = FALSE)
```

## Arguments

- script_name:

  Name of the script to run

- args:

  Character vector of arguments to pass to the script

- dry_run:

  Logical; if TRUE, prints the command instead of running it

## Value

If dry_run = FALSE, returns the exit status invisibly. If dry_run =
TRUE, returns the command string.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate consensus sequences
run_hybphaser_script(
  "1_generate_consensus_sequences.sh",
  args = c("-n", "namelist.txt", "-p", "hybpiper_output")
)

# Preview command without running
run_hybphaser_script(
  "1_generate_consensus_sequences.sh",
  args = c("-n", "namelist.txt"),
  dry_run = TRUE
)
} # }
```
