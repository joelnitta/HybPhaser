# rhybphaser Docker Image
# Bundles the external tools and bash scripts used by the rhybphaser
# Docker workflow functions (run_generate_consensus_sequences(),
# run_extract_mapped_reads()).

FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="rhybphaser"
LABEL org.opencontainers.image.description="Toolchain for the rhybphaser \
Docker workflow: hybrid and polyploid phasing for target capture datasets"
LABEL org.opencontainers.image.source="https://github.com/joelnitta/rhybphaser"

# Switch to root for setup
USER root

# Install dependencies via conda/bioconda
# Using micromamba for faster dependency resolution
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Ensure conda environment is activated by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy the rhybphaser bash scripts and configuration
COPY --chown=$MAMBA_USER:$MAMBA_USER *.sh /opt/rhybphaser/
COPY --chown=$MAMBA_USER:$MAMBA_USER *.R /opt/rhybphaser/
COPY --chown=$MAMBA_USER:$MAMBA_USER config.txt /opt/rhybphaser/
COPY --chown=$MAMBA_USER:$MAMBA_USER example_files/ \
     /opt/rhybphaser/example_files/

# Make bash scripts executable
RUN chmod +x /opt/rhybphaser/*.sh

# Add the scripts to PATH
ENV PATH="/opt/rhybphaser:${PATH}"

# Create default directories
RUN mkdir -p /data/input /data/output && \
    chown -R $MAMBA_USER:$MAMBA_USER /data

# Switch back to non-root user
USER $MAMBA_USER

# Set the working directory for users
WORKDIR /data

# Default command
CMD ["/bin/bash"]
