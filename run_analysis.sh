#!/bin/bash
# run_analysis.sh

# Exit on any error
set -e

# Check if input sample name was provided
if [ -z "$1" ]; then
    echo "Error: Please provide input sample name"
    echo "Usage: $0 <sample_name>"
    exit 1
fi

SAMPLE=$1
echo "Processing sample: $SAMPLE"

# Check if required files exist
if [ ! -f "ref/ucsc_hg19.fa" ]; then
    echo "Error: Reference file not found at ref/ucsc_hg19.fa"
    exit 1
fi

if [ ! -f "samples/${SAMPLE}.bam" ]; then
    echo "Error: BAM file not found at samples/${SAMPLE}.bam"
    exit 1
fi

if [ ! -f "beds/R210.bed" ]; then
    echo "Error: BED file not found at beds/R210.bed"
    exit 1
fi

echo "Running samtools mpileup..."
samtools mpileup -l maf30_snps.txt -f ref/ucsc_hg19.fa "samples/${SAMPLE}.bam" > "results/${SAMPLE}.pileup"

echo "Running VarScan..."
java -jar VarScan.v2.4.6.jar pileup2cns "results/${SAMPLE}.pileup" \
    --min-coverage 20 \
    --min-reads2 0 \
    --min-var-freq 0 > "results/${SAMPLE}.cns"

echo "Running LOH analysis..."
python3 loh.py "results/${SAMPLE}.cns" "beds/R210.bed"

echo "Analysis complete. Results can be found in results/${SAMPLE}.loh.csv"