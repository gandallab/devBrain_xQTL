from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"


"""
rules:
Gene set analysis (using individual level expression and genotype)
    - plink_covar: generate covar file in plink format
    - expr_rel: remove relatatives in expression
    - geno_split: generate chr genotype files
    - effect_egene: estimate eQTL effect sizes
    - gene_set: generate egene set file
    - score_egene: estimate gene set expression scores
    - h2med_egene: estimate h2med
Overall gene expression analysis
    - score_all_gene: estimate overall gene expression score
    - h2med_all_gene: estimate h2med
Overall isoform expression
    - plink_covar_iso
    - expr_rel_iso
    - score_all_iso
    - h2med_all_iso
Overall splicing
    - plink_covar_intron
    - expr_rel_intron
    - score_all_intron
    - h2med_all_intron
"""

GWAS_DIC = {s["trait"]: s["file"] for s in config["GWAS_LIST"]}


rule all:
    input:
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/MESC/out/gene.{GWAS_trait}.all.h2med",
        #     GWAS_trait=list(GWAS_DIC.keys()),
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/MESC/out/gene.{GWAS_trait}.categories.h2med",
        #     GWAS_trait=list(GWAS_DIC.keys()),
        # ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{GWAS_trait}.all.h2med",
            GWAS_trait=list(GWAS_DIC.keys()),
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{GWAS_trait}.categories.h2med",
            GWAS_trait=list(GWAS_DIC.keys()),
        ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{GWAS_trait}.all.h2med",
        #     GWAS_trait=list(GWAS_DIC.keys()),
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{GWAS_trait}.categories.h2med",
        #     GWAS_trait=list(GWAS_DIC.keys()),
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{GWAS_trait}.all.h2med",
        #     GWAS_trait=list(GWAS_DIC.keys()),
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{GWAS_trait}.categories.h2med",
        #     GWAS_trait=list(GWAS_DIC.keys()),
        # ),


############### overall gene expression analysis ###############
rule score_all_gene:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.{{chr}}.{suffix}",
            suffix=["bed", "bim", "fam"],
        ),
        covar="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/90hcp_cov_629_plink_MESC.txt",
        expr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/gene.counts.scaled.normalized.629_MESC.tsv",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{{chr}}.{file}",
            file=["G", "ave_h2cis", "gannot.gz", "expscore.gz", "hsq", "lasso"],
        ),
    resources:
        mem_gb=6,
        time_min=1440,
        num_cores=4,
    params:
        geno="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier",
        out="/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene",
    shell:
        """
        set +eu
        source ~/anaconda3/bin/activate ldsc
        {config[MESC]}run_mesc.py \
            --compute-expscore-indiv \
            --plink-path {config[PLINK]} \
            --expression-matrix {input.expr} \
            --columns 4,1,2,5 \
            --exp-bfile {params.geno} \
            --geno-bfile {config[LDSC]}LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.{wildcards.chr} \
            --chr {wildcards.chr} \
            --out {params.out} \
            --covariates {input.covar}
        set -eu
        """


rule h2med_all_gene:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{chr}.{file}",
            file=["G", "ave_h2cis", "gannot.gz", "expscore.gz", "hsq", "lasso"],
            chr=np.arange(1, 23, 1),
        ),
    output:
        "/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{GWAS_trait}.all.h2med",
        "/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{GWAS_trait}.categories.h2med",
    resources:
        mem_gb=4,
        time_min=120,
    params:
        expscore="/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene",
        out="/u/project/gandalm/cindywen/isoform_twas/MESC/out/test.all.gene.{GWAS_trait}",
        GWAS_file=lambda wildcards: GWAS_DIC[wildcards.GWAS_trait],
    shell:
        """ 
        set +eu
        source ~/anaconda3/bin/activate ldsc
        {config[MESC]}run_mesc.py \
            --h2med {params.GWAS_file}.sumstats.gz \
            --exp-chr {params.expscore} \
            --out {params.out}
        set -eu
        """


# ############### overall isoform expression analysis ###############
# rule plink_covar_iso:
#     input:
#         "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/data/70hcp_cov_626.txt",
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/data/70hcp_cov_626_plink_MESC.txt",
#     resources:
#         num_gb=4,
#         time_min=120,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load R/3.6.0
#         Rscript scripts/plink_covar.R \
#             --input {input[0]} \
#             --out {output[0]}
#         """


# rule expr_rel_iso:
#     input:
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt",
#         "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/data/tx.counts.scaled.normalized.bed.gz",
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/data/tx.counts.scaled.normalized.626_MESC.tsv",
#     resources:
#         num_gb=4,
#         time_min=120,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load R/3.6.0
#         Rscript scripts/expr_rel.R \
#             --rel {input[0]} \
#             --input {input[1]} \
#             --out {output[0]}
#         """


# rule score_all_iso:
#     input:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.{{chr}}.{suffix}",
#             suffix=["bed", "bim", "fam", "log", "nosex"],
#         ),
#         covar="/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/data/70hcp_cov_626_plink_MESC.txt",
#         expr="/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/data/tx.counts.scaled.normalized.626_MESC.tsv",
#     output:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{{chr}}.{file}",
#             file=["G", "ave_h2cis", "gannot.gz", "expscore.gz", "hsq", "lasso"],
#         ),
#     resources:
#         mem_gb=8,
#         time_min=600,
#         num_cores=12,
#     params:
#         geno="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted",
#         out="/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso",
#     shell:
#         """
#         set +eu
#         source ~/anaconda3/bin/activate ldsc
#         {config[MESC]}run_mesc.py \
#             --compute-expscore-indiv \
#             --plink-path {config[PLINK]} \
#             --expression-matrix {input.expr} \
#             --columns 4,1,2,5 \
#             --exp-bfile {params.geno} \
#             --geno-bfile {config[LDSC]}LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.{wildcards.chr} \
#             --chr {wildcards.chr} \
#             --out {params.out} \
#             --covariates {input.covar}
#         set -eu
#         """


# rule h2med_all_iso:
#     input:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{chr}.{file}",
#             file=["G", "ave_h2cis", "gannot.gz", "expscore.gz", "hsq", "lasso"],
#             chr=np.arange(1, 23, 1),
#         ),
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{GWAS_trait}.all.h2med",
#         "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{GWAS_trait}.categories.h2med",
#     resources:
#         mem_gb=4,
#         time_min=240,
#         num_cores=4,
#     params:
#         expscore="/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso",
#         out="/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.iso.{GWAS_trait}",
#         GWAS_file=lambda wildcards: GWAS_DIC[wildcards.GWAS_trait],
#     shell:
#         """
#         set +eu
#         source ~/anaconda3/bin/activate ldsc
#         {config[MESC]}run_mesc.py \
#             --h2med {params.GWAS_file}.sumstats.gz \
#             --exp-chr {params.expscore} \
#             --out {params.out}
#         set -eu
#         """


# ############### overall splicing analysis ###############
# rule plink_covar_intron:
#     input:
#         "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/40hcp_cov_640.txt",
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/40hcp_cov_640_plink_MESC.txt",
#     resources:
#         num_gb=4,
#         time_min=120,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load R/3.6.1
#         Rscript scripts/plink_covar.R \
#             --input {input[0]} \
#             --out {output[0]}
#         """


# rule expr_rel_intron:
#     input:
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt",
#         "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/cluster/leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.bed.gz",
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/cluster/lc_MESC.tsv",
#     resources:
#         num_gb=4,
#         time_min=120,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load R/3.6.1
#         Rscript scripts/expr_rel.R \
#             --rel {input[0]} \
#             --input {input[1]} \
#             --out {output[0]}
#         """


# rule score_all_intron:
#     input:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.{{chr}}.{suffix}",
#             suffix=["bed", "bim", "fam", "log", "nosex"],
#         ),
#         covar="/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/40hcp_cov_640_plink_MESC.txt",
#         expr="/u/project/gandalm/cindywen/isoform_twas/sqtl_new/cluster/lc_MESC.tsv",
#     output:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{{chr}}.{file}",
#             file=["G", "ave_h2cis", "gannot.gz", "expscore.gz", "hsq", "lasso"],
#         ),
#     resources:
#         mem_gb=8,
#         time_min=1800,
#         num_cores=16,
#     params:
#         geno="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted",
#         out="/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron",
#     shell:
#         """
#         set +eu
#         source ~/anaconda3/bin/activate ldsc
#         {config[MESC]}run_mesc.py \
#             --compute-expscore-indiv \
#             --plink-path {config[PLINK]} \
#             --expression-matrix {input.expr} \
#             --columns 4,1,2,5 \
#             --exp-bfile {params.geno} \
#             --geno-bfile {config[LDSC]}LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.{wildcards.chr} \
#             --chr {wildcards.chr} \
#             --out {params.out} \
#             --covariates {input.covar}
#         set -eu
#         """


# rule h2med_all_intron:
#     input:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{chr}.{file}",
#             file=["G", "ave_h2cis", "gannot.gz", "expscore.gz", "hsq", "lasso"],
#             chr=np.arange(1, 23, 1),
#         ),
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{GWAS_trait}.all.h2med",
#         "/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{GWAS_trait}.categories.h2med",
#     resources:
#         mem_gb=4,
#         time_min=240,
#         num_cores=4,
#     params:
#         expscore="/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron",
#         out="/u/project/gandalm/cindywen/isoform_twas/MESC/out/all.intron.{GWAS_trait}",
#         GWAS_file=lambda wildcards: GWAS_DIC[wildcards.GWAS_trait],
#     shell:
#         """
#         set +eu
#         source ~/anaconda3/bin/activate ldsc
#         {config[MESC]}run_mesc.py \
#             --h2med {params.GWAS_file}.sumstats.gz \
#             --exp-chr {params.expscore} \
#             --out {params.out}
#         set -eu
#         """

