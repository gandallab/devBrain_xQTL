from os.path import join
import os
import numpy as np
import sys
import pandas as pd


configfile: "config.yaml"


GENES = pd.read_table(config["PAINTOR_GENE"]).set_index("gene_index", drop=False)
# POPS = ["eur", "amr", "afr"]

"""
rules:
    - all_assoc: for each gene, write the association file for each population
    - zscore: calculate zscore, generate locus file and shared variant file for each gene
    - ld: write paintor LD matrices
    - annot: write paintor annotation files
    - write_list: write paintor input list
    - run_paintor: run paintor with specified number of causal
"""


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/{gene}.results.{max_causal}causal",
            gene=GENES.gene_id,
            max_causal=[2],
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/Enrichment.Values.{max_causal}causal",
            gene=GENES.gene_id,
            max_causal=[2],
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/Log.BayesFactor.{max_causal}causal",
            gene=GENES.gene_id,
            max_causal=[2],
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/LogFile.results.{max_causal}causal",
            gene=GENES.gene_id,
            max_causal=[2],
        ),


rule all_assoc:
    input:
        eur="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_nominal_50HCP/all_assoc.txt.gz",
        amr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/amr_nominal_15HCP/all_assoc.txt.gz",
        afr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/afr_nominal_25HCP/all_assoc.txt.gz",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{pop}_assoc.txt",
            pop=["eur", "amr", "afr"],
        ),
    resources:
        num_gb=4,
        time_min=60,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{wildcards.gene}_dir
        cd /u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{wildcards.gene}_dir
        zcat {input.eur} | awk '$1 == "{wildcards.gene}" {{ print $0 }} ' > eur_assoc.txt
        zcat {input.amr} | awk '$1 == "{wildcards.gene}" {{ print $0 }} ' > amr_assoc.txt
        zcat {input.afr} | awk '$1 == "{wildcards.gene}" {{ print $0 }} ' > afr_assoc.txt
        """


rule zscore:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{pop}_assoc.txt",
            pop=["eur", "amr", "afr"],
        ),
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/{gene}",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/shared_variants.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/paintor_zscore.R \
            --gene {wildcards.gene}
        """


# PAINTOR wants signed pearson r for LD; not r2; not sqrt of r2 as it has no sign
# note plink can output NaN in LD matrix, PAINTOR will get error (all problems are in LD3, AFR)
# From https://github.com/gkichaev/PAINTOR_V3.0/issues/1,
# NA values are probably monomorphic snps. I'd recommend dropping them.
# Paintor 3 does not have faculties to deal with NA anymore.
rule ld:
    input:
        extract="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/shared_variants.txt",
        eur_vcf="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeGeneOutlier.vcf.gz",
        amr_vcf="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/amr/filtered.hg19.sorted.removeGeneOutlier.vcf.gz",
        afr_vcf="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/afr/filtered.hg19.sorted.removeGeneOutlier.vcf.gz",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{{gene}}.{file}",
            file=["LD1", "LD2", "LD3"],
        ),
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624
        cd /u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{wildcards.gene}_dir
        plink --vcf {input.eur_vcf} \
            --r \
            --matrix \
            --extract {input[0]} \
            --out eur
        plink --vcf {input.amr_vcf} \
            --r \
            --matrix \
            --extract {input[0]} \
            --out amr
        plink --vcf {input.afr_vcf} \
            --r \
            --matrix \
            --extract {input[0]} \
            --out afr
        mv eur.ld {wildcards.gene}.LD1
        mv amr.ld {wildcards.gene}.LD2
        mv afr.ld {wildcards.gene}.LD3
        """


# this is wrong and inconvenient
# rule plink_ld:
#     input:
#         extract="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/shared_variants.txt",
#         vcf="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/{population}/filtered.hg19.sorted.removeGeneOutlier.vcf.gz",
#     output:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{{population}}.{file}",
#             file=["ld.bin", "log", "nosex"],
#         ),
#     resources:
#         mem_gb=4,
#         time_min=60,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load plink/1.90b624
#         plink --r2 bin \
#             --vcf {input.vcf} \
#             --extract {input.extract} \
#             --out /u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{wildcards.gene}_dir/{wildcards.population}
#         """


# rule ld:
#     input:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{population}.{file}",
#             file=["ld.bin", "log", "nosex"],
#             population=["eur", "amr", "afr"],
#         ),
#     output:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{{gene}}.{file}",
#             file=["LD1", "LD2", "LD3"],
#         ),
#     resources:
#         mem_gb=4,
#         time_min=60,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load R/4.1.0-BIO
#         Rscript scripts/paintor_ld.R \
#             --gene {wildcards.gene}
#         """


rule annot:
    input:
        annot="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/annot/eur_variant_annot.txt.gz",
        var="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/shared_variants.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/{gene}.annotations",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/paintor_annot.R \
            --gene {wildcards.gene} 
        """


rule write_list:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/{gene}",
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{{gene}}_dir/{{gene}}.{file}",
            file=["LD1", "LD2", "LD3"],
        ),
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/{gene}.annotations",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/input.file",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        echo "{wildcards.gene}" > {output[0]}
        """


# genes in run_paintor_genes.txt have NaN in AFR LD, and 0 in AFR zscore
# variants MAC=0 in AFR genotype, because MAF was filtered before relative/outlier removal
# remove these monomorphic variants and updated input files, see paintor_fixLD.sh
# note don't run fixLD on genes w/0 NA
rule run_paintor:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/input.file",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/{gene}.results.{max_causal}causal",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/Enrichment.Values.{max_causal}causal",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/Log.BayesFactor.{max_causal}causal",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{gene}_dir/LogFile.results.{max_causal}causal",
    resources:
        mem_gb=6,
        num_cores=4,
        time_min=360,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/7.5.0
        cd /u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/{wildcards.gene}_dir/
        /u/project/gandalm/shared/apps/PAINTOR_V3.0/PAINTOR \
            -input {input[0]} \
            -LDname LD1,LD2,LD3 \
            -Zhead ZSCORE.P1,ZSCORE.P2,ZSCORE.P3 \
            -annotations TF_binding_site_d,promoter_flanking_region_d,promoter_d,open_chromatin_region_d,enhancer_d,CTCF_binding_site_d \
            -enumerate {wildcards.max_causal} \
            -in . \
            -out . \
            -Gname Enrichment.Values.{wildcards.max_causal}causal \
            -Lname Log.BayesFactor.{wildcards.max_causal}causal \
            -RESname results.{wildcards.max_causal}causal \
            -set_seed 123 
            # -max_causal {wildcards.max_causal} \
            # -mcmc

        """
