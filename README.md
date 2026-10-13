# HERA

## 1.0 Introduction

We present HERA (HRR and ROH Exploration and Annotation), a scalable, reproducible, and user-friendly Nextflow-based pipeline designed to analyze and annotate ROH and HRR. HERA identifies ROH and/or HRR in a species of interest alongside corresponding genomic islands and gene regions intersected by these segments. HERA can additionally identify unique ROH/HRR in user-defined subpopulations, allowing for more depth in downstream analyses.

## 2.0 HERA steps

<img width="861" height="496" alt="HERA workflow" src="https://github.com/user-attachments/assets/9d1bcf2d-b752-422f-b846-21e22b122aaf" />

*Summary workflow of the HERA pipeline. The HERA pipeline is summarized in both HRR mode (yellow) and ROH mode (blue). Both modes can also be run simultaneously in parallel.*

An overview of the steps involved in the HERA pipeline alongside required software is shown below.

### ROH/HRR identification and annotation

- **Characterize ROH:** using PLINK
- **Characterize HRR:** using detectRUNs
- **Find overlapping genes:** using Bedtools

### Reformat ROH/HRR and SNP islands

- **Reformat ROH/HRR:** using Python (required libraries: NumPy, Pandas, YAML)
- **Identify SNP islands:** using Python (required libraries: NumPy, Pandas, YAML)

### Generate report

- **Produce HTML report:** using Quarto
- **Embed tables and plots in Quarto report:** using R (required libraries: dplyr, stringr, tidyr, tibble, plotly, DT, kableExtra, zoo)

## 3.0 Computational resources and dependencies

HERA supports execution in both local Linux environments as well as high-performance computing (HPC) SLURM clusters. The desired environment is specified at runtime using the `-profile` option with either `local` or `slurm`. Resources can be allocated differently for either local or SLURM execution.

HERA processes can require significant computational resources, and the SLURM option is strongly recommended unless testing locally with a subsampled dataset. The total required resources will vary based on the dataset used, including SNP density, number of samples, mode of operation, among other factors. An example of the computational resources used by HERA in a test dataset is shown below.

The HERA pipeline requires installation of any Nextflow version that employs DSL2 (i.e., any version after Nextflow 20.07.1). For all other dependencies, HERA fully integrates Apptainer containerization and can run all processes using Singularity Image files (SIF). The required software for each process is shown in the `main.yaml` configuration file, and Apptainer containers can optionally be downloaded from the `containers` directory.

<img width="975" height="1013" alt="Computational resources" src="https://github.com/user-attachments/assets/52c1161b-faf2-4604-8846-0cf56e69363f" />

*Computational resources used by HERA for full ROH and HRR analysis. HERA was run in ROH mode on a SLURM cluster, and metrics were calculated on a per-process basis using Nextflow. Calculated metrics are the duration of each process, the number of CPUs used, and the total memory used.*

## 4.0 Input files

To run HERA, the working directory must contain a subdirectory named `param_files`, which contains two parameter files:

- **`main.yaml`:** This parameter file contains the information needed for pipeline execution, including the output directory, temporary file storage, and computational resource allocation for individual processes. Please see the `main.yaml` file provided for a detailed description of each parameter.

- **`ROHRR.yaml`:** This file contains specific parameters for ROH/HRR identification, genomic island detection, and gene annotation. Please see the `ROHRR.yaml` file provided for a detailed description of each parameter.

Additionally, two required files and two optional files are specified. The absolute path to these files must be provided in `main.yaml`. Please refer to the `example_files` directory for well-formatted input files:

- **VCF or PLINK binary files:** The input data must either be a Variant Call Format (VCF) file, or the output of PLINK in binary format (BED, BIM, FAM files). If binary format files are used, they must all contain the same prefix (i.e., `example.bed`, `example.bim`, `example.fam`).

- **Scaffolds:** A list containing the name of each scaffold that is to be analyzed by the pipeline.

- **Gene annotation (OPTIONAL):** A GTF or BED file containing gene annotation information for the species being studied.

- **Subpopulation file (OPTIONAL):** If the subpopulation analysis feature is set to `true` in `ROHRR.yaml`, a file classifying each sample into an appropriate subpopulation must be provided. Note that samples can belong to more than one subpopulation.

## 5.0 How to run

To run HERA, you must first prepare the directory structure as shown in Section 4.0. Next, the parameter files must be configured. Once complete, you can run the HERA pipeline using Nextflow commands. Please see the Nextflow documentation for details on Nextflow functionality.

To run HERA, use the following command:

```bash
nextflow run omicsflow98/HERA -profile <slurm or local>
```

The `slurm` profile is intended for execution on high-performance computing clusters using the SLURM workload manager, while the `local` profile is intended for execution on local systems.
