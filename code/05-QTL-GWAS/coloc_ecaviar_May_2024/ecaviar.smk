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
In this Snakefile, generate a list of locus-gene to run eCAVIAR for each GWAS. Since this file is data-dependent, and the next steps depend on this file, run the next steps in a separate Snakefile (ecaviar_run.smk).
"""

rule all:
    input:
        # expand(
        #     "{gwas}_locus_gene_list_MB.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_probe_list.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_eqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_isoqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_sqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_tri1_eqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_tri2_eqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_tri1_isoqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_tri2_isoqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_tri1_sqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        # expand(
        #     "{gwas}_locus_fetal_tri2_sqtl.txt",
        #     gwas=GWAS_LOCI_TABLE.gwas.values,
        # ),
        expand(
            "{gwas}_1MB_locus_fetal_eqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_isoqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_sqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_tri1_eqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_tri2_eqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_tri1_isoqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_tri2_isoqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_tri1_sqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),
        expand(
            "{gwas}_1MB_locus_fetal_tri2_sqtl.txt",
            gwas=GWAS_LOCI_TABLE.gwas.values,
        ),


# MetaBrain eQTL
rule locus_egene_MB:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
        "/u/project/gandalm/shared/GenomicDatasets/MetaBrain/2021-07-23-cortex-EUR-80PCs-TopEffects.txt.gz",
    output:
        "../out/{gwas}/locus{locus}/locus_egene_MB.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/locus_egene_MB.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --gwas_file {input[1]}
        """

rule write_locus_gene_list_MB:
    input:
        expand(
            "../out/{gwas}/locus{locus}/locus_egene_MB.txt",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
    output:
        "{gwas}_locus_gene_list_MB.txt",
    shell:
        """
        cat ../out/{wildcards.gwas}/locus*/locus_egene_MB.txt > {output[0]}
        sed  -i '1i locus\tgene' {output[0]}
        """

# THISTLE sQTL
rule locus_probe:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
        "/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/job.out.convert",
    output:
        "../out/{gwas}/locus{locus}/locus_probe.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/locus_probe.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --gwas_file {input[1]}
        """

rule write_locus_probe:
    input:
        expand(
            "../out/{gwas}/locus{locus}/locus_probe.txt",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
    output:
        "{gwas}_locus_probe_list.txt",
    shell:
        """
        cat ../out/{wildcards.gwas}/locus*/locus_probe.txt > {output[0]}
        sed  -i '1i locus\tgene' {output[0]}
        """

# fetal xQTL
rule locus_feature:
    input:
        "{gwas}_loci.tsv",
        lambda wildcards: GWAS_DIC[wildcards.gwas],
    output:
        "../out_1MB/{gwas}/locus{locus}/locus_fetal_{annot}.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p ../out_1MB/{wildcards.gwas}/locus{wildcards.locus}

        Rscript scripts/locus_feature.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --gwas_file {input[1]} \
            --annot {wildcards.annot}
        """

rule write_locus_feature:
    input:
        expand(
            "../out_1MB/{gwas}/locus{locus}/locus_fetal_{{annot}}.txt",
            zip,
            gwas=GWAS_LOCI_TABLE.gwas.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
    output:
        "{gwas}_1MB_locus_fetal_{annot}.txt",
    shell:
        """
        cat ../out_1MB/{wildcards.gwas}/locus*/locus_fetal_{wildcards.annot}.txt > {output[0]}
        sed  -i '1i locus\tgene' {output[0]}
        """

