from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"


GWAS_DIC = {s["trait"]: s["file"] for s in config["GWAS_LIST"]}
LOCI_DIC = {s["trait"]: s["locus"] for s in config["LOCI_LIST"]}

GWAS_LOCI_TABLE = pd.DataFrame(
    [(t, l) for t in LOCI_DIC.keys() for l in LOCI_DIC[t]], columns=["trait", "locus"]
)


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_male_eqtl/{trait}/locus_{locus}/locus_egene.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_female_eqtl/{trait}/locus_{locus}/locus_egene.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_tri1_eqtl/{trait}/locus_{locus}/locus_egene.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_tri2_eqtl/{trait}/locus_{locus}/locus_egene.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_tri1_sqtl/{trait}/locus_{locus}/locus_intron.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_tri2_sqtl/{trait}/locus_{locus}/locus_intron.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_male_sqtl/{trait}/locus_{locus}/locus_intron.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_female_sqtl/{trait}/locus_{locus}/locus_intron.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),


rule locus_egene_sex_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
        # "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{group}_eqtl/{trait}/locus_{locus}/locus_egene.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.group}_eqtl/{wildcards.trait}/locus_{wildcards.locus}/

        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_egene_sex_tri.R \
            --locus {wildcards.locus} \
            --trait {wildcards.trait} \
            --group {wildcards.group}
        """


rule locus_intron_sex_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
        # "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{group}_sqtl/{trait}/locus_{locus}/locus_intron.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.group}_sqtl/{wildcards.trait}/locus_{wildcards.locus}/

        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_intron_sex_tri.R \
            --locus {wildcards.locus} \
            --trait {wildcards.trait} \
            --group {wildcards.group}
        """
