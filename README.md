# ProHap Graph
Pipeline to transform the ProHap files and PSM lists into a Neo4j-based graph database. 

## Requirements

The pipeline requires Snakemake and Conda installed. You can install these following [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), _Installation via Conda/Mamba_. 

You will need a local or remote running deployment of Neo4j with the APOC plugin. We have tested the current version of this pipeline with Neo4j v.5.24.2. For a deployment of Neo4j on an Ubuntu server, [this tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-and-configure-neo4j-on-ubuntu-22-04) may be helpful.

## Input Files and Usage

* The three output files produced by [ProHap](https://github.com/ProGenNo/ProHap) - concatenated FASTA, haplotype table, and haplotype FASTA. If using one of the published ProHap databases (e.g., the [protein haplotypes of the 1000 Genomes panel](https://doi.org/10.5281/zenodo.10149277), these are the F1, F2, and F3 files. For F1, use the full rather than the simplified format.
* TSV files containing a list of PSMs from a proteomic search using these databases, annotated with [ProHap Peptide Annotator](https://github.com/ProGenNo/ProHap_PeptideAnnotator). In addition to the standard output of the Peptide Annotator, the following columns should be included:
    - `posterior_error_prob`: Estimated posterior error probability of the PSM,
    - `q-value`: Estimated Q value of the PSM,
    - `USI`: The universal spectrum identifier for this PSM (ideally also containing the peptide sequence and precursor charge),
    - `rt_Abs_error`: The difference between the observed and predicted retention time,
    - `SpectrumTitle`: 
    - `SpectrumFilename`: Name of the raw file containing the spectrum, without the suffix (e.g., `240511_S123_plasma_R1`, the PSM comes from a search on the file `240511_S123_plasma_R1.raw`). This has to match the file names in the metadata file as below.
* Metadata file (SDRF or similar) containing information on the raw files used in the search

Usage:
 1. Clone this repository: `git clone https://github.com/ProGenNo/ProHap_Graph.git; cd ProHap_Graph/;`
 2. Create a configuration file called `config.yaml` based on the instructions in `config_example.yaml`
 3. Test Snakemake with a dry-run: `snakemake --cores <# provided cores> -n -q`
 4. Run the Snakemake pipeline to create your protein database: `snakemake --cores <# provided cores> -p --use-conda`
