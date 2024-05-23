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


"""
rules:
- coloc: run coloc for each GWAS locus. In R, loop over candidate features, make sure beta/SE and REF/ALT are consistent between GWAS and QTL
- concat: concatenate all significant PP4 results
- annotate: add gene names etc.
"""

rule all:
    input:
        # expand(
        #     "../out/{gwas}/coloc_sigPP4_annot.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # )
        expand(
            "../out/{gwas}/locus{locus}/MB_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/thistle_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        
        

################# run COLOC #################
rule coloc_metabrain:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
        "/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz",
    output:
        "../out/{gwas}/locus{locus}/MB_coloc.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/coloc_MB.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --gwas_file {input[1]}

        touch {output[0]}
        """

rule coloc_thistle:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
        "/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/job.out.convert",
    output:
        "../out/{gwas}/locus{locus}/thistle_coloc.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/coloc_thistle.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --gwas_file {input[1]}

        touch {output[0]}
        """

rule coloc:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
    output:
        "../out/{gwas}/locus{locus}/coloc_{annot}.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/coloc_fetal.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --annot {wildcards.annot} \
            --gwas_file {input[1]}
        touch {output[0]}
        """
# ################# concat sig PP4 results #################
rule concat:
    input:
        expand(
            "../out/{gwas}/locus{locus}/MB_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/thistle_coloc.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_eqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_isoqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_sqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_tri1_eqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_tri2_eqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_tri1_isoqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_tri2_isoqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_tri1_sqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "../out/{gwas}/locus{locus}/coloc_tri2_sqtl.done",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
    output:
        "../out/{gwas}/coloc_sigPP4.txt",
    shell:
        """
        cd ../out/{wildcards.gwas}/
        ls -1 locus*/*rds > coloc_sigPP4.txt
        """

rule annot:
    input:
         "../out/{gwas}/coloc_sigPP4.txt",
    output:
         "../out/{gwas}/coloc_sigPP4_annot.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/coloc_annot.R \
            --gwas {wildcards.gwas} 
        """

