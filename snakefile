#nextstrain build protocol
import os
import datetime
import pandas as pd
import pathlib
import glob

wildcard_constraints:
    subtype=r"h1n1pdm|h3n2|vic",
    segment=r"ha|na|pb2|pb1|pa|np|mp|ns"

os.makedirs("data", exist_ok=True)
filenames = list(pathlib.Path('./data').glob('*.fasta'))

subtypes = [filename.stem.split("_")[0] for filename in filenames]
segments = [filename.stem.split("_")[1] for filename in filenames]

reference_name = {
    "h3n2"    : "A/Darwin/6/2021_F",
    "h1n1pdm" : "A/Wisconsin/588/2019", 
    "vic"     : "B/BRISBANE/60/2008"}

header_name = {
    "h3n2"    : "A/H3N2",
    "h1n1pdm" : "A/H1N1(pdm)", 
    "vic"     : "B(Victoria)"}

def clock_rate(lineage, segment):
    # these rates are from 12y runs on 2019-10-18   
    rate = {
     ('h1n1pdm', 'ha'): 0.00329,
     ('h1n1pdm', 'na'): 0.00326,
     ('h1n1pdm', 'np'): 0.00221,
     ('h1n1pdm', 'pa'): 0.00217,
     ('h1n1pdm', 'pb1'): 0.00205,
     ('h1n1pdm', 'pb2'): 0.00277,
     ('h3n2', 'ha'): 0.00382,
     ('h3n2', 'na'): 0.00267,
     ('h3n2', 'np'): 0.00157,
     ('h3n2', 'pa'): 0.00178,
     ('h3n2', 'pb1'): 0.00139,
     ('h3n2', 'pb2'): 0.00218,
     ('vic', 'ha'): 0.00145,
     ('vic', 'na'): 0.00133,
     ('vic', 'np'): 0.00132,
     ('vic', 'pa'): 0.00178,
     ('vic', 'pb1'): 0.00114,
     ('vic', 'pb2'): 0.00106,
     ('yam', 'ha'): 0.00176,
     ('yam', 'na'): 0.00177,
     ('yam', 'np'): 0.00133,
     ('yam', 'pa'): 0.00112,
     ('yam', 'pb1'): 0.00092,
     ('yam', 'pb2'): 0.00113}
    return rate.get((lineage, segment), 0.001)

def clock_std_dev(lineage, segement):
    return 0.2*clock_rate(lineage, segement)

rule all:
    input:
        auspice_json  = expand("results/{subtype}_{segment}/auspice/{subtype}_{segment}_auspice.json", zip, subtype=subtypes, segment=segments), 
        nexus_out = expand("results/{subtype}_{segment}/tree/{subtype}_{segment}_tree.nwk", zip, subtype=subtypes, segment=segments),
        nextclade_out = expand("nextclade/{subtype}_{segment}/clade.tsv", zip, subtype=subtypes, segment=segments),
        metadata_out = expand("data/{subtype}_{segment}_metadata.tsv", zip, subtype=subtypes, segment=segments)

rule generate_concatenated_files:
    output:
        'data/concatenated_files/{subtype}/{segment}/{subtype}_{segment}_with_reference.fasta'
    shell:
        """
        mkdir -p $(dirname {output})
        cat data/{wildcards.subtype}_{wildcards.segment}.fasta \
            config/reference_fasta/reference_{wildcards.subtype}_{wildcards.segment}.fasta \
            > {output}
        """

rule clade:
    input: 
        sequence = 'data/concatenated_files/{subtype}/{segment}/{subtype}_{segment}_with_reference.fasta'
    output:
        "nextclade/{subtype}_{segment}/clade.tsv"
    shell:
        """
        nextclade3 dataset get -n flu_{wildcards.subtype}_{wildcards.segment} --output-dir nextclade/{wildcards.subtype}_{wildcards.segment} 
        nextclade3 run -j 5 -D nextclade/{wildcards.subtype}_{wildcards.segment} \
                  {input.sequence} --quiet --output-tsv {output}
        """

rule combined_metadata:
    input:
        sample_meta = "data/meta_{subtype}_{segment}.tsv",
        reference_meta = "config/metadata/reference_{subtype}_{segment}.tsv",
        nextclade = "nextclade/{subtype}_{segment}/clade.tsv"
    output:
        "data/{subtype}_{segment}_metadata.tsv"
    shell:
        """
        # First combine the metadata files
        tail -n +2 {input.sample_meta} > temp_noheader_sample.tsv
        csvtk concat -t -H {input.reference_meta} temp_noheader_sample.tsv > temp_metadata.tsv


        # Then join with nextclade data
        csvtk join -t -H --fields "strain;seqName" temp_metadata.tsv {input.nextclade} > temp_joined.tsv
        # Remove duplicate lines from the joined output
        csvtk uniq -t -H temp_joined.tsv > {output}
        
        # Clean up temp file
        #rm temp_metadata.tsv
        #rm temp_joined.tsv
        """

rule augur_index:
    input: 
        sequence = 'data/concatenated_files/{subtype}/{segment}/{subtype}_{segment}_with_reference.fasta'
    output:
        "results/{subtype}_{segment}/tree/{subtype}_{segment}_index.tsv"
    shell:
        """ 
        augur index --sequences {input.sequence} --output {output} 2>&1
        """

rule augur_align:
    input:
        sequence = 'data/concatenated_files/{subtype}/{segment}/{subtype}_{segment}_with_reference.fasta',
        reference = 'config/reference/reference_{subtype}_{segment}.gb'
    output:
        "results/{subtype}_{segment}/tree/{subtype}_{segment}_aligned.fasta"
    log:
        "logs/{subtype}_{segment}_align.log"    
    shell:
        """
        augur align \
          --sequences {input.sequence} \
          --reference-sequence {input.reference} \
          --output {output} \
          --fill-gaps \
          --remove-reference &> {log}
        """    
    
rule augur_tree:
    input: 
        'results/{subtype}_{segment}/tree/{subtype}_{segment}_aligned.fasta'
    output:
        "results/{subtype}_{segment}/tree/{subtype}_{segment}_tree_raw.nwk"
    log:
        "logs/{subtype}_{segment}_tree.log"
    shell:
        """
        augur tree --alignment {input} --output {output} &> {log}
        """


rule augur_refine:
    input:
        tree = "results/{subtype}_{segment}/tree/{subtype}_{segment}_tree_raw.nwk",
        alignment = "results/{subtype}_{segment}/tree/{subtype}_{segment}_aligned.fasta",
        metadata = "data/{subtype}_{segment}_metadata.tsv"
    output:
        tree = "results/{subtype}_{segment}/tree/{subtype}_{segment}_tree.nwk",
        node = "results/{subtype}_{segment}/tree/{subtype}_{segment}_branch_lengths.json"
    params:
        clock_rate = lambda wildcards: clock_rate(wildcards.subtype, wildcards.segment),
        clock_std_dev = lambda wildcards: clock_std_dev(wildcards.subtype, wildcards.segment),
        root = lambda wildcards: reference_name[wildcards.subtype]
    log:
        "logs/{subtype}_{segment}_refine.log"
    shell:
        """
        augur curate format-dates \
          --metadata {input.metadata} \
          --date-fields date \
          --expected-date-formats "%d/%m/%Y" --output-metadata formatted_metadata.tsv  

        augur refine \
          --tree {input.tree} \
          --alignment {input.alignment} \
          --metadata formatted_metadata.tsv \
          --output-tree {output.tree} \
          --output-node-data {output.node} \
          --timetree \
          --clock-rate {params.clock_rate} \
          --clock-std-dev {params.clock_std_dev} \
          --root {params.root} \
          --date-inference marginal \
          --keep-polytomies  &> {log}
        """
rule traits:  
    input:  
        tree =  "results/{subtype}_{segment}/tree/{subtype}_{segment}_tree.nwk",  
        metadata = "data/{subtype}_{segment}_metadata.tsv" 
    output:  
        node_data = "results/{subtype}_{segment}/tree/{subtype}_{segment}_traits.json"
    log:
        "logs/{subtype}_{segment}_traits.log"    
    shell:  
        """  
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns subclade clade \
            --confidence \
            --sampling-bias-correction 2.0 &> {log}
        """  

rule augur_ancestral:
    input:
        tree =  "results/{subtype}_{segment}/tree/{subtype}_{segment}_tree.nwk",
        alignment = "results/{subtype}_{segment}/tree/{subtype}_{segment}_aligned.fasta"
    log:
        "logs/{subtype}_{segment}_ancestral.log"
    output:
        "results/{subtype}_{segment}/tree/{subtype}_{segment}_nt_muts.json"
    shell:
        """
        augur ancestral \
          --tree {input.tree} \
          --alignment {input.alignment} \
          --output-node-data {output} \
          --inference joint &> {log}
        """

rule augur_translate:
    input:
        tree =  "results/{subtype}_{segment}/tree/{subtype}_{segment}_tree.nwk",
        ancestral_sequence = "results/{subtype}_{segment}/tree/{subtype}_{segment}_nt_muts.json",
        reference_sequence = 'config/reference/reference_{subtype}_{segment}.gb'
    log:
        "logs/{subtype}_{segment}_translate.log"    
    output:
        "results/{subtype}_{segment}/tree/{subtype}_{segment}_aa_muts.json"
    shell:
        """
        augur translate \
          --tree {input.tree} \
          --ancestral-sequences {input.ancestral_sequence} \
          --reference-sequence {input.reference_sequence} \
          --output-node-data {output} &> {log}
        """
    
rule export:
    message:
        "Exporting data files for auspice"
    input:
        tree="results/{subtype}_{segment}/tree/{subtype}_{segment}_tree.nwk",
        metadata="data/{subtype}_{segment}_metadata.tsv",
        branch_lengths = "results/{subtype}_{segment}/tree/{subtype}_{segment}_branch_lengths.json",
        aa_mut = "results/{subtype}_{segment}/tree/{subtype}_{segment}_aa_muts.json",
        nt_mut = "results/{subtype}_{segment}/tree/{subtype}_{segment}_nt_muts.json",
        traits = "results/{subtype}_{segment}/tree/{subtype}_{segment}_traits.json",
        auspice_config= "config/auspice_config.json"
    output:
        auspice_json="results/{subtype}_{segment}/auspice/{subtype}_{segment}_auspice.json",
    params:
        fields="passage type clade subclade location",
        select_fields = ','.join(["strain", "passage", "center", "type", "clade","subclade", "location"]),
        header =  lambda wildcards: header_name[wildcards.subtype]
    log:
        "logs/{subtype}_{segment}_export.log"    
    shell:
        """
        tsv-select -H -f {params.select_fields} {input.metadata} > {input.metadata}.tmp
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata}.tmp \
            --node-data {input.branch_lengths} {input.aa_mut} {input.nt_mut} {input.traits}\
            --auspice-config {input.auspice_config} \
            --color-by-metadata {params.fields} \
            --minify-json \
            --title "Influenza type: {params.header} " \
            --include-root-sequence \
            --output {output.auspice_json} &> {log};
        """
