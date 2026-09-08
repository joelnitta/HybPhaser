# rhybphaser Docker Image

This image bundles the external command-line tools and the bash scripts
that the `rhybphaser` R package runs in a container
([`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md),
[`run_extract_mapped_reads()`](https://joelnitta.github.io/rhybphaser/reference/run_extract_mapped_reads.md)).
Normally you do not use it directly: the R functions call `docker run`
for you. This guide covers building/publishing the image and running it
by hand.

The R functions default to the image `joelnitta/rhybphaser:latest`;
override with their `docker_image` argument.

## Prerequisites

- Docker installed on your system
- Docker Compose (optional, but recommended)

## Included Software

All dependencies are installed via conda/bioconda:

- R (≥4.0) with packages: ape, seqinr, stringr
- BWA (≥0.7.17)
- SAMtools (≥1.9)
- BCFtools (≥1.9)
- BBMap (≥38.87) including BBSplit
- HybPiper (≥2.0)
- the bundled `rhybphaser` bash scripts

The specific versions installed are defined in `environment.yml`.

## Quick Start

### Option 1: Using Docker Compose (Recommended)

1.  **Build the image:**

    ``` bash
    docker compose build
    ```

2.  **Run the container:**

    ``` bash
    docker compose run --rm rhybphaser
    ```

3.  **Access your data:**

    - Place input files in `./data/` directory
    - Outputs will be saved to `./output/` directory

### Option 2: Using Docker directly

1.  **Build the image:**

    ``` bash
    docker build -t joelnitta/rhybphaser:latest .
    ```

2.  **Run the container:**

    ``` bash
    docker run -it --rm \
      -v $(pwd)/data:/data/input \
      -v $(pwd)/output:/data/output \
      joelnitta/rhybphaser:latest
    ```

## Usage Examples

### Running the Bundled Scripts

The container’s entry points are the two bash scripts on `PATH`:

``` bash
# Generate consensus sequences
1_generate_consensus_sequences.sh -n namelist.txt \
  -p /data/input/hybpiper_output \
  -o /data/output/hybphaser_output

# Extract mapped reads
2_extract_mapped_reads.sh -b /data/input/hybphaser_output \
  -o /data/output/mapped_reads
```

From R, call
[`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md)
/
[`run_extract_mapped_reads()`](https://joelnitta.github.io/rhybphaser/reference/run_extract_mapped_reads.md)
instead — they mount the right directories and invoke these scripts for
you.

### Interactive Session

For an interactive session where you can run multiple commands:

``` bash
docker run -it --rm \
  -v $(pwd)/data:/data/input \
  -v $(pwd)/output:/data/output \
  joelnitta/rhybphaser:latest /bin/bash
```

## Directory Structure

Inside the container: - `/opt/rhybphaser/` - bundled bash scripts and
configuration - `/data/input/` - Mounted input data directory -
`/data/output/` - Mounted output directory - `/opt/conda/` - Conda
environment with all dependencies

## Resource Management

### Memory Limits

If you encounter memory issues with Java (BBMap/BBSplit), you can:

1.  **Set memory in docker-compose.yml:**

    ``` yaml
    mem_limit: 8g
    ```

2.  **Set memory when running Docker directly:**

    ``` bash
    docker run --memory=8g -it joelnitta/rhybphaser:latest
    ```

3.  **Adjust Java memory in config.txt:**

        java_memory_usage_clade_association = "6G"
        java_memory_usage_phasing = "6G"

### CPU Limits

Control the number of CPUs used:

``` bash
docker run --cpus=4 -it joelnitta/rhybphaser:latest
```

Or in docker-compose.yml:

``` yaml
cpus: 4
```

## Complete Workflow Example

``` bash
# 1. Prepare your data structure
mkdir -p data/hybpiper_output data/namelist
mkdir -p output

# 2. Copy your HybPiper results to data/hybpiper_output/
# 3. Create namelist.txt in data/namelist/

# 4. Build and run
docker compose build
docker compose run --rm rhybphaser

# Inside container:
# 5. Generate consensus sequences
1_generate_consensus_sequences.sh \
  -n /data/input/namelist/namelist.txt \
  -p /data/input/hybpiper_output \
  -o /data/output/hybphaser_output
```

Then run the downstream analysis from R on the host with the
`rhybphaser` package
([`count_snps()`](https://joelnitta.github.io/rhybphaser/reference/count_snps.md),
[`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md),
…), pointing at `./output/hybphaser_output`.

## Troubleshooting

### Permission Issues

If you encounter permission issues with output files:

``` bash
# Linux/Mac: Run container with your user ID
docker run -it --rm --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data/input \
  -v $(pwd)/output:/data/output \
  joelnitta/rhybphaser:latest
```

### Out of Memory Errors

1.  Increase Docker’s memory allocation in Docker Desktop settings
2.  Set appropriate memory limits in config.txt
3.  Reduce the number of parallel threads in config.txt

### BBMap/BBSplit Issues

If BBMap commands fail:

``` bash
# Check if BBMap is in PATH
which bbsplit.sh

# BBMap commands should be available via conda
bbsplit.sh [options]
```

### Modifying Dependencies

To add or update dependencies:

1.  Edit `environment.yml`

2.  Rebuild the Docker image:

    ``` bash
    docker compose build --no-cache
    ```

## Publishing

The image is built and pushed to Docker Hub by
`.github/workflows/docker.yml` on pushes to `main` that touch the
`Dockerfile`, `environment.yml`, or the bundled scripts (requires the
`DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` repository secrets).

## Additional Resources

- [rhybphaser GitHub](https://github.com/joelnitta/rhybphaser)
- [Original HybPhaser method and
  code](https://github.com/LarsNauheimer/HybPhaser)
- [HybPhaser
  Manuscript](https://www.biorxiv.org/content/10.1101/%5C%202020.10.27.354589v2)
- [HybPiper Documentation](https://github.com/mossmatters/HybPiper/wiki)

## Support

For issues related to: - Docker image or `rhybphaser`: Open an issue in
this repository - The HybPhaser method itself: See the original
HybPhaser repository - Dependencies: Check respective software
documentation
