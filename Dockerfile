# Dockerfile
FROM debian:12-slim

# Prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    samtools \
    openjdk-17-jre-headless \
    curl \
    --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

# Set up Python virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install Python dependencies in virtual environment
COPY requirements.txt /app/requirements.txt
RUN pip3 install --no-cache-dir -r /app/requirements.txt

# Create necessary directories
RUN mkdir -p /app/ref /app/samples /app/results /app/beds

# Copy files into container
COPY maf30_snps.txt /app/
COPY loh.py /app/
COPY centromeres.json /app/
COPY run_analysis.sh /app/
COPY beds/R210.bed /app/beds/

# Download VarScan from GitHub
RUN curl -L https://github.com/dkoboldt/varscan/releases/download/v2.4.6/VarScan.v2.4.6.jar -o /app/VarScan.v2.4.6.jar

# Set working directory
WORKDIR /app

# Make the run script executable
RUN chmod +x run_analysis.sh

# Default command
ENTRYPOINT ["./run_analysis.sh"]