from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"

# TEST_TABLE = pd.read_table("../results_" + str(wildcards.level) + "_ieqtl/" + str(wildcards.trait) + "/loc_mod_gene_list.txt")

ISO_TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/isoform_twas/colocal/results_isomod_ieqtl/PGC3_SCZ_wave3.european.autosome.public.v3/loc_mod_gene_list.txt")

GENE_TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/isoform_twas/colocal/results_genemod_ieqtl/PGC3_SCZ_wave3.european.autosome.public.v3/loc_mod_gene_list.txt")

def get_1kg_eur_bfile_chr(wildcards):
    LOCI_TABLE = pd.read_table("tables/" + str(wildcards.trait) + "_1Mb.txt").set_index(
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
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_genemod_ieqtl/PGC3_SCZ_wave3.european.autosome.public.v3/locus_{locus}/{module_gene}.done",
            zip,
            locus=GENE_TEST_TABLE.locus.values,
            module_gene=GENE_TEST_TABLE.module_gene.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_isomod_ieqtl/PGC3_SCZ_wave3.european.autosome.public.v3/locus_{locus}/{module_gene}.done",
            zip,
            locus=ISO_TEST_TABLE.locus.values,
            module_gene=ISO_TEST_TABLE.module_gene.values,
        ),



rule zscore:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_gwas_zscore.txt",
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_gwas_snplist.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore_mod.R \
            --level {wildcards.level} \
            --trait {wildcards.trait} \
            --locus {wildcards.locus} \
            --module_gene {wildcards.module_gene}
        """

rule ld:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_gwas_snplist.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_1kg_eur.ld",
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_eqtl.ld",
    resources:
        mem_gb=4,
        time_min=60,
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624
        
        plink --bfile {params[0]} \
              --r \
              --matrix \
              --extract {input[0]} \
              --out /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.level}_ieqtl/{wildcards.trait}/locus_{wildcards.locus}/{wildcards.module_gene}_1kg_eur
        
        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
                --r \
                --matrix \
                --extract {input[0]} \
                --out /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.level}_ieqtl/{wildcards.trait}/locus_{wildcards.locus}/{wildcards.module_gene}_eqtl
        """

rule run_ecaviar:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_eqtl.ld",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_ecaviar_col",
    resources:
        mem_gb=16,
        time_min=120,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.level}_ieqtl/{wildcards.trait}/locus_{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.module_gene}_ecaviar \
                -l {wildcards.module_gene}_1kg_eur.ld \
                -z {wildcards.module_gene}_gwas_zscore.txt \
                -l {wildcards.module_gene}_eqtl.ld \
                -z {wildcards.module_gene}.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{level}_ieqtl/{trait}/locus_{locus}/{module_gene}.done",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.level}_ieqtl/{wildcards.trait}/locus_{wildcards.locus}/

        Rscript /u/project/gandalm/cindywen/isoform_twas/colocal/code/scripts/analyze_mod.R \
                --level {wildcards.level} \
                --trait {wildcards.trait} \
                --locus {wildcards.locus} \
                --module_gene {wildcards.module_gene}

        touch {output[0]}
        """