# LDscaff: LD-based scaffolding of de novo genome assemblies

## User Manual

*Version 1.0*  
*July 27, 2025*

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
  - [System Requirements](#system-requirements)
  - [Setting up Conda Environment](#setting-up-conda-environment)
  - [Installing LDscaff](#installing-ldscaff)
- [Usage](#usage)
  - [Quick Start](#quick-start)
  - [Full Pipeline](#full-pipeline)
  - [Step-by-Step Commands](#step-by-step-commands)
- [Command Reference](#command-reference)
  - [run](#run)
  - [extract](#extract)
  - [prepare](#prepare)
  - [match](#match)
  - [paths](#paths)
  - [assemble](#assemble)
- [Input File Formats](#input-file-formats)
- [Output Files](#output-files)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Contact Information](#contact-information)

## Introduction

LDscaff is a tool for improving genome assemblies using population-scale linkage disequilibrium (LD) information. It can connect scaffolds together based on genetic linkage detected in population samples, helping to create chromosome-scale assemblies.

The software consists of a C++ core for LD calculation and maximum weight matching, and a Python CLI (command-line interface) for executing the full pipeline. LDscaff does not require a genetic map or other reference information, making it ideal for de novo genome assembly projects.

## Installation

### System Requirements

- Linux-based operating system
- At least 8GB of RAM (16GB+ recommended for large genomes)
- 10GB+ of free disk space
- [Conda](https://docs.conda.io/en/latest/) package manager

### Setting up Conda Environment

1. Create a new conda environment with Python 3.8:

```bash
conda create -n ldscaff python=3.8
```

2. Activate the environment:

```bash
conda activate ldscaff
```

3. Install the required packages:

```bash
# Install BioPython
conda install -c conda-forge biopython

# Install PyVCF
conda install -c bioconda pyvcf
```

### Installing LDscaff

1. Clone the LDscaff repository:

```bash
git clone https://github.com/username/ldscaff.git
cd ldscaff
```

2. Compile the C++ matching program:

```bash
cd src
make
cd ..
```

3. Test the installation:

```bash
# Make sure the CLI script is executable
chmod +x pipeline/cli.py

# Run the help command
./pipeline/cli.py --help
```

## Usage

### Quick Start

To run the complete LDscaff pipeline on your data:

```bash
./pipeline/cli.py run \
  --vcf your_data.vcf \
  --input-fasta scaffolds.fasta \
  --output results_directory \
  --output-fasta improved_assembly.fasta
```

### Full Pipeline

The full pipeline consists of the following steps:

1. Extract SNP markers from the VCF file
2. Prepare sample data for the matching program
3. Run the maximum weight matching algorithm
4. Find scaffold paths based on matching results
5. Assemble the final scaffolds

Running the `run` command executes all these steps automatically. However, you can also run each step individually as described below.

### Step-by-Step Commands

For more control, you can run each step of the pipeline separately:

```bash
# 1. Extract SNP markers from VCF
./pipeline/cli.py extract \
  --vcf your_data.vcf \
  --output marker_data \
  --marker-count 100 \
  --sample-list samples.txt

# 2. Prepare sample data for matching
./pipeline/cli.py prepare \
  --marker-dir marker_data \
  --samples samples.txt \
  --output prepared_data

# 3. Run matching algorithm
./pipeline/cli.py match \
  --sample-list prepared_data/sample.list \
  --samples 100 \
  --scaffolds 500 \
  --output match_results.txt

# 4. Find scaffold paths
./pipeline/cli.py paths \
  --match-result match_results.txt \
  --scaffolds 500 \
  --threshold 0.2 \
  --output path_results.txt

# 5. Assemble final scaffolds
./pipeline/cli.py assemble \
  --paths path_results.txt \
  --input-fasta scaffolds.fasta \
  --output-fasta improved_assembly.fasta
```

## Command Reference

### run

Run the complete LDscaff pipeline.

```
./pipeline/cli.py run [options]
```

**Required arguments:**
- `--vcf FILE`: Input VCF file with population variation data
- `--output DIR`: Output directory

**Optional arguments:**
- `--samples FILE`: Sample list file (will be generated if not provided)
- `--gap-size INT`: Gap size between scaffolds in final assembly [default: 100]
- `--minor-ratio FLOAT`: SNP marker minor ratio cutoff [default: 0.2]
- `--ld-threshold FLOAT`: LD threshold for breaking weak connections [default: 0.2]
- `--marker-count INT`: Number of markers per scaffold end [default: 100]
- `--input-fasta FILE`: Input scaffold FASTA file for final assembly
- `--output-fasta FILE`: Output assembled scaffold FASTA file
- `--temp-dir DIR`: Temporary directory for matching program [default: /tmp]

### extract

Extract SNP markers from a VCF file.

```
./pipeline/cli.py extract [options]
```

**Required arguments:**
- `--vcf FILE`: Input VCF file
- `--output DIR`: Output directory for markers

**Optional arguments:**
- `--scaffolds FILE`: File containing scaffold names (will be extracted from VCF if not provided)
- `--marker-count INT`: Number of markers per scaffold end [default: 100]
- `--minor-ratio FLOAT`: SNP marker minor ratio cutoff [default: 0.2]
- `--sample-list FILE`: Output file for sample list

### prepare

Prepare sample data for the matching program.

```
./pipeline/cli.py prepare [options]
```

**Required arguments:**
- `--marker-dir DIR`: Directory containing extracted markers
- `--samples FILE`: Sample list file
- `--output DIR`: Output directory

**Optional arguments:**
- `--scaffolds FILE`: File containing scaffold names (will be extracted from marker directory if not provided)

### match

Run the maximum weight matching algorithm.

```
./pipeline/cli.py match [options]
```

**Required arguments:**
- `--sample-list FILE`: Sample list file from prepare step
- `--samples INT`: Number of samples
- `--scaffolds INT`: Number of scaffolds
- `--output FILE`: Output file for matching results

**Optional arguments:**
- `--marker-count INT`: Number of markers to use from each scaffold [default: 100]
- `--minor-ratio FLOAT`: Minor allele cutoff [default: 0.2]
- `--temp-dir DIR`: Temporary directory [default: /tmp]
- `--edge-log FILE`: Log file for edge weights [default: edge.log]

### paths

Find scaffold paths from matching results.

```
./pipeline/cli.py paths [options]
```

**Required arguments:**
- `--match-result FILE`: Match result file from match step
- `--scaffolds INT`: Number of scaffolds
- `--threshold FLOAT`: LD threshold for breaking weak connections
- `--output FILE`: Output file for path results

### assemble

Assemble scaffolds based on paths.

```
./pipeline/cli.py assemble [options]
```

**Required arguments:**
- `--paths FILE`: Path result file from paths step
- `--input-fasta FILE`: Input scaffold FASTA file
- `--output-fasta FILE`: Output assembled scaffold FASTA file

**Optional arguments:**
- `--scaffolds FILE`: File containing scaffold names (will be extracted from input FASTA if not provided)
- `--gap-size INT`: Gap size between scaffolds in output [default: 100]

## Input File Formats

### VCF File

LDscaff requires a standard VCF (Variant Call Format) file containing population variation data. The VCF file should include:

- Multiple samples (minimum 30 samples recommended)
- SNP variants (indels are not used)
- Genotype information (GT field)

The VCF file should be either sorted and indexed (for better performance) or at least organized by scaffold/chromosome.

### FASTA File

The scaffold FASTA file should contain all the scaffolds you want to assemble. Each scaffold should have a unique identifier in the header line:

```
>scaffold1
ATGCATGCATGCATGCATGC...
>scaffold2
GCTAGCTAGCTAGCTAGCTA...
```

## Output Files

### Main Output Files

- `assembled.fasta`: The final assembled scaffolds
- `paths.result`: The scaffold paths used for assembly
- `match.result`: The matching results showing connections between scaffold ends
- `edge.log`: Detailed edge weights from the matching algorithm

### Intermediate Files

- `preprocess/`: Directory containing extracted SNP markers
- `samples/`: Directory containing prepared sample data
- `sample.list`: List of sample files for the matching program

## Examples

### Example 1: Basic Assembly

Assembling a draft genome with default parameters:

```bash
./pipeline/cli.py run \
  --vcf population.vcf \
  --input-fasta draft.fasta \
  --output ldscaff_results \
  --output-fasta improved.fasta
```

### Example 2: Fine-tuning Parameters

Using custom parameters for better results:

```bash
./pipeline/cli.py run \
  --vcf population.vcf \
  --input-fasta draft.fasta \
  --output ldscaff_results \
  --output-fasta improved.fasta \
  --marker-count 150 \
  --minor-ratio 0.1 \
  --ld-threshold 0.15 \
  --gap-size 200
```

### Example 3: Step-by-Step Processing

Running the pipeline steps individually to monitor progress:

```bash
# Step 1: Extract markers
./pipeline/cli.py extract \
  --vcf population.vcf \
  --output markers \
  --marker-count 100 \
  --sample-list samples.txt

# Check marker extraction results
ls markers

# Step 2: Prepare samples
./pipeline/cli.py prepare \
  --marker-dir markers \
  --samples samples.txt \
  --output prepared

# Step 3: Run matching
./pipeline/cli.py match \
  --sample-list prepared/sample.list \
  --samples 96 \
  --scaffolds 1242 \
  --output matching.out

# Examine matching results
head matching.out

# Step 4: Find paths
./pipeline/cli.py paths \
  --match-result matching.out \
  --scaffolds 1242 \
  --threshold 0.2 \
  --output paths.out

# Look at path results
cat paths.out

# Step 5: Assemble final scaffolds
./pipeline/cli.py assemble \
  --paths paths.out \
  --input-fasta draft.fasta \
  --output-fasta improved.fasta

# Check results
grep ">" improved.fasta | wc -l
```

## Troubleshooting

### Common Issues

1. **Error: Command 'matching' not found**
   - Make sure you compiled the C++ matching program
   - Check if the matching program is in the correct location

2. **Error: Failed to open VCF file**
   - Verify the path to your VCF file
   - Check if the VCF file is readable

3. **Error: No samples found in VCF file**
   - Make sure your VCF file contains genotype information
   - Check if the VCF file is properly formatted

4. **Error: Not enough markers found**
   - Try reducing the `--minor-ratio` parameter
   - Check if your VCF file contains enough variants

5. **Error: No paths found**
   - Try reducing the `--ld-threshold` parameter
   - Check the `edge.log` file to see if any edges were detected

### Getting Help

If you encounter issues not covered in this manual, please:

1. Run the command with the `--verbose` flag to get more detailed output
2. Check the log files in the output directory
3. Contact the developers with a description of the problem and the log files

## Contact Information

For support, bug reports, or feature requests, please contact:

- Zicheng ZHAO <zachazhao2-c@my.cityu.edu.hk>

Or create an issue on the GitHub repository: https://github.com/username/ldscaff

---

*This manual was created on July 27, 2025 for LDscaff version 1.0.*
