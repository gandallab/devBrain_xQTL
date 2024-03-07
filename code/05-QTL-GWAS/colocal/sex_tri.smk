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
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/colocal/results_male_eqtl/{trait}/locus_{locus}/locus_egene.txt",
        #     zip,
        #     trait=GWAS_LOCI_TABLE.trait.values,
        #     locus=GWAS_LOCI_TABLE.locus.values,
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/colocal/results_female_eqtl/{trait}/locus_{locus}/locus_egene.txt",
        #     zip,
        #     trait=GWAS_LOCI_TABLE.trait.values,
        #     locus=GWAS_LOCI_TABLE.locus.values,
        # ),
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
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_tri1_isoqtl/{trait}/locus_{locus}/locus_isoform.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_tri2_isoqtl/{trait}/locus_{locus}/locus_isoform.txt",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/colocal/results_male_sqtl/{trait}/locus_{locus}/locus_intron.txt",
        #     zip,
        #     trait=GWAS_LOCI_TABLE.trait.values,
        #     locus=GWAS_LOCI_TABLE.locus.values,
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/colocal/results_female_sqtl/{trait}/locus_{locus}/locus_intron.txt",
        #     zip,
        #     trait=GWAS_LOCI_TABLE.trait.values,
        #     locus=GWAS_LOCI_TABLE.locus.values,
        # ),


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

rule locus_isoform_sex_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
        # "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{group}_isoqtl/{trait}/locus_{locus}/locus_isoform.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.group}_isoqtl/{wildcards.trait}/locus_{wildcards.locus}/

        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_isoform_sex_tri.R \
            --locus {wildcards.locus} \
            --trait {wildcards.trait} \
            --group {wildcards.group}
        """


# # SCZ
# qsub -t 1-307 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri1_eqtl PGC3_SCZ_wave3.european.autosome.public.v3 egene 
# qsub -t 1-307 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri2_eqtl PGC3_SCZ_wave3.european.autosome.public.v3 egene 

# # ADHD
# qsub -t 1-11 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri1_isoqtl ADHD.Demontis.2019 isoform 
# qsub -t 1-11 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri2_isoqtl ADHD.Demontis.2019 isoform

# # ASD
# qsub -t 1-3 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri1_isoqtl ASD.iPSYCHPGC.2018 isoform 
# qsub -t 1-3 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri2_isoqtl ASD.iPSYCHPGC.2018 isoform

# # MDD
# qsub -t 1-101 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri1_isoqtl MDD.Howard.PGC.2019 isoform 
# qsub -t 1-101 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri2_isoqtl MDD.Howard.PGC.2019 isoform

# # BIP
# qsub -t 1-63 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri1_isoqtl pgc-bip2021-all isoform 
# qsub -t 1-63 ~/project-gandalm/isoform_twas/colocal/code/iso_s_bash/analyze.sh tri2_isoqtl pgc-bip2021-all isoform