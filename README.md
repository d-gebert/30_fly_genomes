# 30 Fly Genomes

**Analysis of 30 Chromosome-Level *Drosophila* Genome Assemblies Reveals Dynamic Evolution of Centromeric Satellite Repeats**

## Introduction

This repository contains the scripts and code utilized in the study titled "Analysis of 30 Chromosome-Level *Drosophila* Genome Assemblies Reveals Dynamic Evolution of Centromeric Satellite Repeats" by Gebert D, Hay AD, Hoang J, Gibbon AE, Henderson IR, and Teixeira FK (doi: https://doi.org/10.1101/2024.06.17.599346). The project encompasses:

- **Hi-C Genome Scaffolding**: Utilizing Hi-C chromatin conformation data to scaffold genome assemblies to chromosome-level.
- **Hi-C-Based Muller Element Organization and Compartment Analysis**: Analyzing chromatin compartments and interactions between Muller elements.
- **Gene Annotation**: Predicting gene models using de novo transcriptome assemblies and protein homology.
- **Genome Rearrangement Analysis**: Detecting syntenic blocks and structural variations to study genome rearrangements.
- **Transposable Element Annotation**: Analysing transposable elements (TEs) across the genomes.
- **Tandem Repeat Annotation and Analysis**: Characterizing tandem repeats, including centromeric satellite repeats.

## Project Structure

The repository is organized as follows:

- `Code/HiC_scaffolding/`: Scripts for Hi-C genome scaffolding centred around YAHS.
- `Code/HiC_analysis/`: Scripts for Hi-C data processing and chromatin compartment analysis.
- `Code/Gene_annotation/`: Scripts for gene annotation based on MAKER.
- `Code/Genome_rearrangement/`: Scripts for synteny and rearrangement analysis based on GENSPACE.
- `Code/TE_annotation/`: Scripts for transposable element analysis and visualization.
- `Code/TR_analysis/`: Scripts for TRASH-based tandem repeat analysis and visualization.
- `Docs/`: Documentation and supplementary materials.

## Installation

To replicate the analyses, ensure you have the following dependencies installed:

- **Python**: Version 3.x
- **Perl**: Version 5.x
- **R**: Version 4.x

## Usage

Instructions for running each analysis are provided within each scripts usage/help messages.

## Data Availability

Due to the large size of the genomic datasets, input data is hosted externally. Please refer to the [Data Availability Statement](Docs/Data_availability.md) for information on accessing the necessary datasets.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE.md) file for details.

## Contact

For any inquiries or further information, please contact:

- **Name**: Daniel Gebert
- **Email**: dg572@cam.ac.uk
- **Institution**: University of Cambridge
