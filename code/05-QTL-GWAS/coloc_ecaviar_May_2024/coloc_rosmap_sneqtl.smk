from os.path import join
import os
import numpy as np
import pandas as pd
import sys

configfile: "config.yaml"

GWAS_DIC = {s["gwas"]: s["file"] for s in config["GWAS_LIST"]}
LOCI_DIC = {s["gwas"]: s["locus"] for s in config["LOCI_LIST"]}
GWAS_LOCI_TABLE = pd.DataFrame(
    [(g, l) for g in LOCI_DIC.keys() for l in LOCI_DIC[g]], columns=["gwas", "locus"]
)

rule all:
    input:
        expand(
            "../out/{gwas}/locus{locus}/rosmap_Ast_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/rosmap_End_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/rosmap_Exc_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/rosmap_Inh_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/rosmap_Mic_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/rosmap_Oli_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/rosmap_OPC_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        
rule coloc_rosmap_sneqtl:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
        "/u/project/gandalm/shared/GenomicDatasets/ROSMAP_snRNA_eQTL/celltype-eqtl-sumstats.{type}.tsv.gz",
    output:
        "../out/{gwas}/locus{locus}/rosmap_{type}_coloc.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/coloc_rosmap_sneqtl.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --gwas_file {input[1]} \
            --type {wildcards.type}

        touch {output[0]}
        """
