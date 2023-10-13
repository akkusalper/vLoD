# Variant Limit of Detection Predictor (vLoD)

## Introduction

vLoD is a comprehensive open-source tool designed to predict the detectability status of alleles from variant call files (VCF) using matched sequencing data. 
vLoD calculates the likelihood of observing each variant in the context of a given sequencing error rate, true positive rate, and false positive rate. This allows users to assign a detectability score to each variant and classify variants as detectable or non-detectable.

## Features

- **Detectability Scoring**: For each variant in the VCF, vLoD computes a detectability score, which is a log-odds ratio representing the likelihood of the variant being a true positive.
- **Parallel Processing**: vLoD utilizes all available CPU cores for faster variant processing.
- **Integration with VCF**: The tool provides functionality to integrate detectability status directly into the input VCF, allowing users to quickly assess the detectability of variants in downstream analyses.

## Requirements

- Python 3.x
- Docker (for containerized execution)
- Pysam
- Pandas

## Installation

While vLoD can be run from source, we recommend using the provided Docker container for ease of use.

#### Please follow the docker page for the latest release and updates: https://hub.docker.com/r/alperakkus/vlod

## Usage Example

#### Variant Limit of Detection Prediction:
###### To evaluate the detectability of variants, use the `LOD_11_05_23_updated_14_08_23.py` script.
docker run -v $PWD:/data --rm -w /data -t --entrypoint python alperakkus/vlod:latest /usr/src/app/LOD_11_05_23_updated_14_08_23.py --input-vcf [input.vcf] --input-bam [input.bam] --input-bam-index [input.bam.bai] --output [output.xls]

#### Integrating Detectability Status into VCF:
###### To integrate detectability status into the original VCF, use the `merge_detectability.py` script.
docker run -v $PWD:/data --rm -w /data -t --entrypoint python alperakkus/vlod:latest merge_detectability.py /data/[input.vcf] /data/[input.xls] /data/[output.vcf]

### From Source

1. Clone the repository:

git clone https://github.com/akkusalper/vLoD.git
cd vLoD

2. Install required Python libraries:

### Example

#### Variant Limit of Detection Prediction:
python LOD_11_05_23_updated_14_08_23.py --input-vcf [input.vcf] --input-bam [input.bam] --input-bam-index [input.bam.bai] --output [output.xls]

#### Integrating Detectability Status into VCF:
python merge_detectability.py [input.vcf] [input.xls] [output.vcf]

## Outputs

vLoD produces an output table (in `.xls` format) with the following columns:

- **VCF_ID**: Identifier for each variant in the format `chrom_position_ref_alt`.
- **Detectability_Score**: The log-odds ratio representing the likelihood of the variant being a true positive.
- **Detectability_Condition**: A classification of the variant as either "Detectable" or "Non-detectable".
- **Coverage**: Total coverage at the variant position.
- **Variant_Reads**: Number of reads supporting the given allele.
