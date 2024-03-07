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
    LOCI_TABLE = pd.read_table("tables/" + str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return chromosome


def get_locus_start(wildcards):
    LOCI_TABLE = pd.read_table("tables/" + str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    start = LOCI_TABLE.loc[int(wildcards.locus), "START"]
    return start


def get_locus_end(wildcards):
    LOCI_TABLE = pd.read_table("tables/" + str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    end = LOCI_TABLE.loc[int(wildcards.locus), "END"]
    return end


def get_1kg_eur_bim_chr(wildcards):
    LOCI_TABLE = pd.read_table("tables/" + str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
        + ".bim"
    )


def get_1kg_eur_bfile_chr(wildcards):
    LOCI_TABLE = pd.read_table("tables/" + str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
    )


# should have used ct_hcp as one wildcard


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_end_100hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_ex_90hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_in_90hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_ip_80hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_mic_80hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_opc_80hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_per_70hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_pg_80hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_rg_100hcp_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),


# gwas variants extracted for each locus is in results_eqtl
rule locus_egene:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{type}_nominal_{hcp}hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/locus_egene.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.type}_{wildcards.hcp}hcp_eqtl/{wildcards.trait}/locus_{wildcards.locus}/

        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_egene_ct.R \
            --locus {wildcards.locus} \
            --trait {wildcards.trait} \
            --type {wildcards.type} \
            --hcp {wildcards.hcp}
        """


rule extract_cis_assoc:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{type}_nominal_{hcp}hcp/all.chunks.txt.gz",
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/locus_egene.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/extract_cis_assoc.done",
    resources:
        mem_gb=6,
        time_min=240,
    shell:
        """
        cat {input[1]} | while read gene
        do
            awk -v a="${{gene}}" '$1 == a {{print}}' <(zcat {input[0]}) > /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.type}_{wildcards.hcp}hcp_eqtl/{wildcards.trait}/locus_{wildcards.locus}/${{gene}}_all_pairs.txt
        done
        touch {output[0]}
        """


rule zscore:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/extract_cis_assoc.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/zscore.done",
    resources:
        mem_gb=6,
        time_min=240,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        cat /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.type}_{wildcards.hcp}hcp_eqtl/{wildcards.trait}/locus_{wildcards.locus}/locus_egene.txt | while read gene
        do
            Rscript scripts/zscore_ct.R \
                --locus {wildcards.locus} \
                --gene ${{gene}} \
                --trait {wildcards.trait} \
                --type {wildcards.type} \
                --hcp {wildcards.hcp}
        done
        touch {output[0]}
        """


rule ld:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/zscore.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/ld.done",
    resources:
        mem_gb=4,
        time_min=60,
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.type}_{wildcards.hcp}hcp_eqtl/{wildcards.trait}/locus_{wildcards.locus}/

        cat locus_egene.txt | while read gene
        do
            plink --bfile {params[0]} \
                --r \
                --matrix \
                --extract ${{gene}}_gwas_eqtl_snp_set.txt \
                --out ${{gene}}_1kg_eur
            plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier \
                --r \
                --matrix \
                --extract ${{gene}}_gwas_eqtl_snp_set.txt \
                --out ${{gene}}_eqtl
        done
        touch {output[0]}
        """


rule run_ecaviar:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/ld.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/ecaviar.done",
    resources:
        mem_gb=8,
        num_cores=2,
        time_min=720,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.type}_{wildcards.hcp}hcp_eqtl/{wildcards.trait}/locus_{wildcards.locus}/

        cat locus_egene.txt | while read gene
        do
            /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o ${{gene}}_ecaviar \
                -l ${{gene}}_1kg_eur.ld \
                -z ${{gene}}_gwas_zscore.txt \
                -l ${{gene}}_eqtl.ld \
                -z ${{gene}}_eqtl_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        done
        touch {output[0]}
        """


rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/ecaviar.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_{type}_{hcp}hcp_eqtl/{trait}/locus_{locus}/analyze.done",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_{wildcards.type}_{wildcards.hcp}hcp_eqtl/{wildcards.trait}/locus_{wildcards.locus}

        cat locus_egene.txt | while read gene
        do
            Rscript /u/project/gandalm/cindywen/isoform_twas/colocal/code/scripts/analyze_ct.R \
                --locus {wildcards.locus} \
                --gene ${{gene}} \
                --trait {wildcards.trait} \
                --type {wildcards.type} \
                --hcp {wildcards.hcp}
        done
        touch {output[0]}
        """


# for i in end_100 ex_90 in_90 ip_80 mic_80 per_70 rg_100 pg_80 opc_80
# > do awk 'FNR==1 && NR!=1 {next;}{print}' results_${i}hcp_eqtl/MDD.Howard.PGC.2019/*/*sig.txt > results_${i}hcp_eqtl/MDD.Howard.PGC.2019/CLPP_sig.txt

# for i in end_100 ex_90 in_90 ip_80 mic_80 per_70 rg_100 pg_80 opc_80
# do cat results_${i}hcp_eqtl/MDD.Howard.PGC.2019/CLPP_sig.txt >> cell_MDD_sig.txt
# dat <- dat %>% filter(SNP_ID != "SNP_ID")
# dat <- dat %>% left_join(ref, by=c("gene"="ensg"))
# dat <- dat %>% select(c(1:6),V11,V12)
# colnames(dat)[7] <- "gene_type"
# colnames(dat)[8] <- "gene_name"
# write.table(dat, "cell_ADHD_sig.txt", col.names=T, row.names=F, quote=F, sep="\t")
