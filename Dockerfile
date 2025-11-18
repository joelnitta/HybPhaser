# HybPhaser Docker Image
# Includes all dependencies for running HybPhaser v2.x

FROM mambaorg/micromamba:1.5.8

LABEL maintainer="HybPhaser"
LABEL description="Docker image for HybPhaser - hybrid and polyploid \
phasing for target capture datasets"
LABEL version="2.1"

# Switch to root for setup
USER root

# Install dependencies via conda/bioconda
# Using micromamba for faster dependency resolution
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Ensure conda environment is activated by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy HybPhaser scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER *.sh /opt/hybphaser/
COPY --chown=$MAMBA_USER:$MAMBA_USER *.R /opt/hybphaser/
COPY --chown=$MAMBA_USER:$MAMBA_USER config.txt /opt/hybphaser/
COPY --chown=$MAMBA_USER:$MAMBA_USER example_files/ \
     /opt/hybphaser/example_files/

# Make bash scripts executable
RUN chmod +x /opt/hybphaser/*.sh

# Add HybPhaser to PATH
ENV PATH="/opt/hybphaser:${PATH}"

# Create default directories
RUN mkdir -p /data/input /data/output && \
    chown -R $MAMBA_USER:$MAMBA_USER /data

# Switch back to non-root user
USER $MAMBA_USER

# Set the working directory for users
WORKDIR /data

# Default command
CMD ["/bin/bash"]
