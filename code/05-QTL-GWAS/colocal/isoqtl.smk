from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"


"""
rules:
- see eqtl.smk
"""

GWAS_DIC = {s["trait"]: s["file"] for s in config["GWAS_LIST"]}
LOCI_DIC = {s["trait"]: s["locus"] for s in config["LOCI_LIST"]}

GWAS_LOCI_TABLE = pd.DataFrame(
    [(t, l) for t in LOCI_DIC.keys() for l in LOCI_DIC[t]], columns=["trait", "locus"]
)


def get_locus_chr(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return chromosome


def get_locus_start(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    start = LOCI_TABLE.loc[int(wildcards.locus), "START"]
    return start


def get_locus_end(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    end = LOCI_TABLE.loc[int(wildcards.locus), "END"]
    return end


def get_1kg_eur_bim_chr(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
        + ".bim"
    )


def get_1kg_eur_bfile_chr(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
    )


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{trait}/locus_{locus}/extract_cis_assoc_qsub.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),


# gwas variants extracted for each locus is in results_eqtl


rule locus_isoform:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
        "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{trait}/locus_{locus}/locus_isoform.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{wildcards.trait}/locus_{wildcards.locus}/

        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_isoform.R \
            --locus {wildcards.locus} \
            --trait {wildcards.trait}
        """


# need to manually check all files are done (does snakemake's checkpoint do this?)
# also will terminate when >500 jobs are in the queue
# see sqtl_extarct_cis.sh
rule extract_cis_assoc:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/gtex.allpairs.txt.gz",
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{trait}/locus_{locus}/locus_isoform.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{trait}/locus_{locus}/extract_cis_assoc_qsub.done",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        cat {input[1]} | while read gene
        do
            qsub scripts/extract_cis_assoc.sh \
                ${{gene}} \
                {input[0]} \
                /u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{wildcards.trait}/locus_{wildcards.locus}/
        done
        touch {output[0]}
        """
