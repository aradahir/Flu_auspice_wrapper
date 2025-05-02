# Flu_auspice_wrapper pipeline

A Snakemake-based pipeline for automating phylogenetic analysis and visualization via the Nextstrain platform.

## About this workflow

Phylogenetic trees are created using `augur` tools by aligning sequences to provided reference files. Clade and subclade information are automatically retrieved during pipeline execution using `nextclade`, ensuring up-to-date clade assignment. Tree building is performed using MAFFT and IQ-TREE, followed by time-resolved reconstruction using `timetree`. Internal node traits are inferred using the traits command in augur, allowing annotation of ancestral nodes with clade and subclade information.

## Dependencies

Ensure the following tools and libraries are installed:
- biopython=1.85
- iqtree=2.4.0
- mafft=7.525
- nextclade=3.12.0
- nextstrain-augur=29.0.0
- phylo-treetime=0.11.4
- pandas=2.2.3


## Installation and Setup

Follow these steps to clone and install the repository (Linux-based systems):

    ```
    git clone ssh https://git@github.com:aradahir/Flu_auspice_wrapper.git
    cd .\Flu_auspice_wrapper\
    conda env create -f requirements.yml

    ```

Activate the snakemake environment for running the pipeline.
    
    ```
    conda activate nextstrian_build
    
    ```

## Usage

 - Prepare Input
    Place your .fasta and metadata files into the ./data/ directory:
    ./data/
    ├── variant_segment.fasta
    └── meta_variant_segment.tsv
     
 - Open the environment
 
        ```
        conda activate nextstrain_build
        ```
 - Run the pipeline
        -
        thread can be adjusted by changing from 1 into the specific number of thread

        ```
        snakemake -j 1
        ```

# Output Structure

├── data/
│   └── concatenated_files/
│       └── variant/
│           └── segment/
│               └── variant_segment_with_reference.fasta
│                   # Input sequences with reference added
├── nextclade/
│   └── variant_segment.db
│       # Clade and subclade assignment
├── results/
│   └── variant_segment/
│       ├── auspice/
│       │   ├── variant_segment.json : json files using for visualize the data in https://auspice.us/
│       │   └── variant_segment_auspice_root-sequence.json
│       └── tree/
│           ├── variant_segment_aa_muts.json
│           ├── variant_segment_aligned.fasta
│           ├── variant_segment_aligned.fasta.insertion.csv
│           ├── variant_segment_branch_length.json
│           ├── variant_segment_nt_muts.json
│           ├── variant_segment_traits.json
│           ├── variant_segment_traitsclade.mugration_model.txt
│           ├── variant_segment_traitssubclade.mugration_model.txt
│           ├── variant_segment_tree.nwk
│           └── variant_segment_tree_raw.nwk
├── formatted_metadata.tsv
    # Reformatted metadata used in the pipeline



[!WARNING]

This pipeline is designed for assembled influenza sequences with known subtypes.

Currently, only the HA segment is supported.

File naming is used to determine variant and segment information.

Before rerunning the pipeline, please remove the following: results/, data/, formatted_metadata.tsv


