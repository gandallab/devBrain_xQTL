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


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_isoqtl/{trait}/locus_{locus}/locus_isoform.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),


# gwas variants extracted for each locus is in results_eqtl
# Ran other steps in bash


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
