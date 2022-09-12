from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"


GWAS_DIC = {s["trait"]: s["file"] for s in config["GWAS_LIST"]}

# update 6-6-2022: files are deleted
rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/sLDSC/results/{GWAS_trait}_PEC_{level}.{file}",
            file=["cov", "delete", "log", "part_delete", "results"],
            level=["gene", "isoform"],
            GWAS_trait=list(GWAS_DIC.keys()),
        ),


rule ldsc:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/sLDSC/annot/PEC_{level}/PEC_{level}.{chr}.annot.gz",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/sLDSC/annot/PEC_{{level}}/PEC_{{level}}.{{chr}}.{file}",
            file=["l2.ldscore.gz", "l2.M", "l2.M_5_50", "log"],
        ),
    resources:
        mem_gb=4,
        time_min=60,
    params:
        annot_dir="/u/project/gandalm/cindywen/isoform_twas/sLDSC/annot/PEC_{level}/",
        annot_prefix="PEC_{level}",
    shell:
        """
        set +eu 
        source ~/anaconda3/bin/activate ldsc
        {config[LDSC_PATH]}ldsc.py \
            --l2 \
            --bfile {config[LDSC_PATH]}LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.{wildcards.chr} \
            --print-snps {config[LDSC_PATH]}LDSCORE/1000G_EUR_Phase3_baseline/print_snps.txt \
            --ld-wind-cm 1 \
            --annot {params.annot_dir}{params.annot_prefix}.{wildcards.chr}.annot.gz \
            --out {params.annot_dir}{params.annot_prefix}.{wildcards.chr}
        set -eu
        """


rule partition_h2:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/sLDSC/annot/PEC_{{level}}/PEC_{{level}}.{chr_num}.{file}",
            file=["l2.ldscore.gz", "l2.M", "l2.M_5_50", "log"],
            chr_num=np.arange(1, 23, 1),
        ),
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/sLDSC/results/{{GWAS_trait}}_PEC_{{level}}.{file}",
            file=["cov", "delete", "log", "part_delete", "results"],
        ),
    resources:
        mem_gb=4,
        time_min=60,
    params:
        out_dir="/u/project/gandalm/cindywen/isoform_twas/sLDSC/results/",
        annot_dir="/u/project/gandalm/cindywen/isoform_twas/sLDSC/annot/PEC_{level}/",
        annot_prefix="PEC_{level}",
        GWAS_file=lambda wildcards: GWAS_DIC[wildcards.GWAS_trait],
    shell:
        """
        set +eu 
        source ~/anaconda3/bin/activate ldsc
        {config[LDSC_PATH]}ldsc.py \
            --h2 {params.GWAS_file}.sumstats.gz \
            --ref-ld-chr {config[LDSC_PATH]}LDSCORE/1000G_EUR_Phase3_baseline/baseline.,{params.annot_dir}{params.annot_prefix}. \
            --frqfile-chr {config[LDSC_PATH]}LDSCORE/1000G_Phase3_frq/1000G.EUR.QC. \
            --w-ld-chr {config[LDSC_PATH]}LDSCORE/weights_hm3_no_hla/weights. \
            --overlap-annot \
            --print-cov \
            --print-coefficients \
            --print-delete-vals \
            --out {params.out_dir}{wildcards.GWAS_trait}_{params.annot_prefix}
        set -eu
        """
