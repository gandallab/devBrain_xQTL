from os.path import join
import os
import numpy as np
import sys

configfile: "config.yaml"

"""
rules:
- index: generate STAR index with Gencode v33, with annotation GTF as recommended

"""

rule index:
    input:
        config["GENOME_FASTA"],
        config["GTF_FILE"]
    output:
        "/u/project/gandalm/cindywen/isoform_twas/star/gencodev33_STARindex/SAindex"
    resources:
        mem_gb=4,
        time_min=120,
        num_cores=8
    shell:
        """
        {config[STAR]} --runThreadN {resources.num_cores} \
            --runMode genomeGenerate \
            --genomeDir /u/project/gandalm/cindywen/isoform_twas/star/gencodev33_STARindex \
            --genomeFastaFiles {config[GENOME_FASTA]} \
            --sjdbGTFfile {config[GTF_FILE]} \
            --sjdbOverhang 100
        """
