# HybPhaser Docker Usage Guide

This Docker image provides a complete environment for running HybPhaser with 
all required dependencies pre-installed.

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
- HybPhaser scripts

The specific versions installed are defined in `environment.yml`.

## Quick Start

### Option 1: Using Docker Compose (Recommended)

1. **Build the image:**
   ```bash
   docker-compose build
   ```

2. **Run the container:**
   ```bash
   docker-compose run --rm hybphaser
   ```

3. **Access your data:**
   - Place input files in `./data/` directory
   - Outputs will be saved to `./output/` directory

### Option 2: Using Docker directly

1. **Build the image:**
   ```bash
   docker build -t hybphaser:latest .
   ```

2. **Run the container:**
   ```bash
   docker run -it --rm \
     -v $(pwd)/data:/data/input \
     -v $(pwd)/output:/data/output \
     hybphaser:latest
   ```

## Usage Examples

### Running HybPhaser Scripts

Once inside the container, you can run HybPhaser scripts:

```bash
# Generate consensus sequences
1_generate_consensus_sequences.sh -n namelist.txt \
  -p /data/input/hybpiper_output \
  -o /data/output/hybphaser_output

# Run R scripts
cd /opt/hybphaser
R --no-save < 1a_count_snps.R
```

### Running with Custom Configuration

1. **Copy and edit the config file:**
   ```bash
   cp /opt/hybphaser/config.txt /data/input/config.txt
   # Edit config.txt with your settings
   ```

2. **Run the main script:**
   ```bash
   docker run -it --rm \
     -v $(pwd)/data:/data/input \
     -v $(pwd)/output:/data/output \
     hybphaser:latest \
     R -e "config_file='/data/input/config.txt'; \
           source('/opt/hybphaser/HybPhaser_Main_script.R')"
   ```

### Interactive Session

For an interactive session where you can run multiple commands:

```bash
docker run -it --rm \
  -v $(pwd)/data:/data/input \
  -v $(pwd)/output:/data/output \
  hybphaser:latest /bin/bash
```

## Directory Structure

Inside the container:
- `/opt/hybphaser/` - HybPhaser scripts and configuration
- `/data/input/` - Mounted input data directory
- `/data/output/` - Mounted output directory
- `/opt/conda/` - Conda environment with all dependencies

## Resource Management

### Memory Limits

If you encounter memory issues with Java (BBMap/BBSplit), you can:

1. **Set memory in docker-compose.yml:**
   ```yaml
   mem_limit: 8g
   ```

2. **Set memory when running Docker directly:**
   ```bash
   docker run --memory=8g -it hybphaser:latest
   ```

3. **Adjust Java memory in config.txt:**
   ```
   java_memory_usage_clade_association = "6G"
   java_memory_usage_phasing = "6G"
   ```

### CPU Limits

Control the number of CPUs used:

```bash
docker run --cpus=4 -it hybphaser:latest
```

Or in docker-compose.yml:
```yaml
cpus: 4
```

## Complete Workflow Example

```bash
# 1. Prepare your data structure
mkdir -p data/hybpiper_output data/namelist
mkdir -p output

# 2. Copy your HybPiper results to data/hybpiper_output/
# 3. Create namelist.txt in data/namelist/

# 4. Build and run
docker-compose build
docker-compose run --rm hybphaser

# Inside container:
# 5. Generate consensus sequences
1_generate_consensus_sequences.sh \
  -n /data/input/namelist/namelist.txt \
  -p /data/input/hybpiper_output \
  -o /data/output/hybphaser_output

# 6. Configure and run R scripts
cd /opt/hybphaser
# Edit config.txt as needed, then:
R --no-save < HybPhaser_Main_script.R
```

## Troubleshooting

### Permission Issues

If you encounter permission issues with output files:

```bash
# Linux/Mac: Run container with your user ID
docker run -it --rm --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data/input \
  -v $(pwd)/output:/data/output \
  hybphaser:latest
```

### Out of Memory Errors

1. Increase Docker's memory allocation in Docker Desktop settings
2. Set appropriate memory limits in config.txt
3. Reduce the number of parallel threads in config.txt

### BBMap/BBSplit Issues

If BBMap commands fail:

```bash
# Check if BBMap is in PATH
which bbsplit.sh

# BBMap commands should be available via conda
bbsplit.sh [options]
```

### Modifying Dependencies

To add or update dependencies:

1. Edit `environment.yml`
2. Rebuild the Docker image:
   ```bash
   docker-compose build --no-cache
   ```

## Additional Resources

- [HybPhaser GitHub](https://github.com/larsnauheimer/HybPhaser)
- [HybPhaser Manuscript](https://www.biorxiv.org/content/10.1101/\
2020.10.27.354589v2)
- [HybPiper Documentation](https://github.com/mossmatters/HybPiper/wiki)

## Support

For issues related to:
- Docker image: Open an issue in this repository
- HybPhaser functionality: Refer to the main HybPhaser repository
- Dependencies: Check respective software documentation
