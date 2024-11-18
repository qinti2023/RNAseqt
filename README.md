# Snakemake workflow of RNAseq

This package can help you generate Snakemake workflows, enabling you to perform upstream RNA-seq analysis in one step.

**Author** :`qinti`

## Installation

Install the package using pip:

```bash
mamba env create - n RNAseq -f environment.yml
pip install RNAseqt
```

## Usage

Run the command to initialize a project:

```bash
RNAseqt /path/to/project
```

Edit the `congfig.yaml`:

```shell
vim congfig.yaml
# samples 
# reference
```

Run the workflow of RNAseq upstream

```shell
snakemake -n
snakemake --core 32
```

