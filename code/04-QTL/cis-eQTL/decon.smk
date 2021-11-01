from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"


"""
rules:
cell type specific
    - ct_cov
    - ct_fastqtl_nominal
    - ct_merge_nominal
    - ct_call_nominal
    - ct_fastqtl_perm
    - ct_merge_perm
    - ct_call_perm
cell type/group interaction
    - make_decon_dosage
    - snps_to_test
    - fix_decon_dosage
    - run_decon_qtl
"""

TYPES = ["end", "ex", "in", "ip", "mic", "opc", "per", "pg", "rg"]


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/{file}",
            file=[
                "all_assoc.txt.gz",
                "significant_assoc.txt",
                "significant_feature_count.txt",
            ],
            cell_type=TYPES,
            num_hcp=np.arange(10, 101, 10),
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/all_assoc.txt.gz",
            zip, 
            cell_type=TYPES, 
            num_hcp=[100, 90, 90, 80, 80, 80, 70, 80, 100],
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/sig_pheno.txt",
            zip, 
            cell_type=TYPES, 
            num_hcp=[100, 90, 90, 80, 80, 80, 70, 80, 100],
        ),
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/decon/deconvolutionResults.csv",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/decon/predictedExpressionLevels.txt",


################################### Cell type specific ###################################
rule ct_cov:
    input:
        ct_expr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/{cell_type}.bed.gz",
        picard=(
            "/u/project/gandalm/cindywen/isoform_twas/picard/picard_QC_compiled.tsv"
        ),
        geno_pc="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/pca/data.ref.eigenvec",
        meta="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/metadata_inferSex.tsv",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/{cell_type}_{num_hcp}HCP_cov.txt",
    resources:
        mem_gb=8,
        time_min=60,
    params:
        num_gpc=5,
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/",
        prefix="{cell_type}",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/3.6.0
        Rscript scripts/decon_cov.R \
            --num_hcp {wildcards.num_hcp} \
            --expr {input.ct_expr} \
            --picard {input.picard} \
            --geno_pc {input.geno_pc} \
            --meta {input.meta} \
            --num_gpc {params.num_gpc} \
            --outdir {params.outdir} \
            --prefix {params.prefix}
        """


rule ct_fastqtl_nominal:
    input:
        expr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/{cell_type}.bed.gz",
        geno="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.vcf.gz",
        cov="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/{cell_type}_{num_hcp}HCP_cov.txt",
        rel_file="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/chunk{chunk}.txt.gz",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/chunk{chunk}.txt.gz.done",
    resources:
        mem_gb=4,
        num_cores=4,
        time_min=120,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp",
    shell:
        """
        mkdir -p {params.outdir}
        {config[FASTQTL]} --vcf {input.geno} \
                          --bed {input.expr} \
                          --cov {input.cov} \
                          --out {output[0]} \
                          --chunk {wildcards.chunk} 100 \
                          --seed 123 \
                          --window 1e6 \
                          --normal \
                          --exclude-samples {input.rel_file}
        touch {output[0]}.done
        """


rule ct_merge_nominal:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{{cell_type}}_nominal_{{num_hcp}}hcp/chunk{chunk}.txt.gz",
            chunk=np.arange(1, 101, 1),
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{{cell_type}}_nominal_{{num_hcp}}hcp/chunk{chunk}.txt.gz.done",
            chunk=np.arange(1, 101, 1),
        ),
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/all.chunks.txt.gz",
    resources:
        mem_gb=4,
        num_cores=6,
        time_min=240,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/",
    shell:
        """
        zcat {params.outdir}chunk*.txt.gz | gzip -c > {params.outdir}all.chunks.txt.gz
        """


rule ct_call_nominal:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/all.chunks.txt.gz",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{{cell_type}}_nominal_{{num_hcp}}hcp/{file}",
            file=[
                "all_assoc.txt.gz",
                "significant_assoc.txt",
                "significant_feature_count.txt",
            ],
        ),
    resources:
        mem_gb=4,
        time_min=160,
        num_cores=6,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_nominal_{num_hcp}hcp/",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/3.6.0

        Rscript scripts/call_nominal.R \
            --outdir {params.outdir}
        gzip {params.outdir}all_assoc.txt
        """

# Celine ran this in bash script
rule ct_fastqtl_perm:
    input:
        expr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/{cell_type}.bed.gz",
        geno="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.vcf.gz",
        cov="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/{cell_type}_{num_hcp}HCP_cov.txt",
        rel_file="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/chunk{chunk}.txt.gz",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/chunk{chunk}.txt.gz.done",
    resources:
        mem_gb=6,
        num_cores=4,
        time_min=480,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp",
    shell:
        """
        mkdir -p {params.outdir}
        {config[FASTQTL]} --vcf {input.geno} \
                          --bed {input.expr} \
                          --cov {input.cov} \
                          --permute 1000 10000 \
                          --out {output[0]} \
                          --chunk {wildcards.chunk} 300 \
                          --seed 123 \
                          --window 1e6 \
                          --normal \
                          --exclude-samples {input.rel_file}
        touch {output[0]}.done
        """

rule ct_merge_perm:
    input:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{{cell_type}}_perm_{{num_hcp}}hcp/chunk{chunk}.txt.gz",
            chunk=np.arange(1, 301, 1),
        ),
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{{cell_type}}_perm_{{num_hcp}}hcp/chunk{chunk}.txt.gz.done",
            chunk=np.arange(1, 301, 1),
        ),
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/all.chunks.txt.gz",
    resources:
        mem_gb=6,
        time_min=120,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/",
    shell:
        """
        zcat {params.outdir}chunk*.txt.gz | gzip -c > {params.outdir}all.chunks.txt.gz
        """

rule ct_call_perm:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/all.chunks.txt.gz",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{{cell_type}}_perm_{{num_hcp}}hcp/{file}",
            file=["all_assoc.txt.gz", "sig_pheno.txt"],
        ),
    resources:
        mem_gb=6,
        time_min=120,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/{cell_type}_perm_{num_hcp}hcp/",
        script="scripts/call_perm.R",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript {params.script} \
            --outdir {params.outdir}
        gzip {params.outdir}all_assoc.txt
        """

################################### Cell group interaction ###################################
# Change SNP position to ID in susie dosage file
rule make_decon_dosage:
    input:
        bim="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.bim.header",
        susie_dosage="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.tsv.gz",
    output:
        expand(
            "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/{file}",
            file=[
                "filtered.hg19.sorted.removeGeneOutlier.dose.decon.tsv.temp",
                "filtered.hg19.sorted.removeGeneOutlier.dose.decon.tsv.gz",
            ],
        ),
    resources:
        mem_gb=4,
        num_cores=4,
        time_min=120,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load htslib

        # Paste files together. Note can do this because they are in the same order
        paste {input.bim} <(zcat {input.susie_dosage}) > {output[0]}

        # Select columns
        cut -f2,11-639 {output[0]} | bgzip > {output[1]}
        """

# not using this
# rule remove_id:
#     input:
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.tsv.gz",
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/sample_list.tsv",
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.noHeader.tsv.gz",
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/sample_list_row_forDecon.tsv.gz",
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.usethis.tsv.gz",
#     resources:
#         mem_gb=4,
#         num_cores=4,
#         time_min=120,
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load htslib/1.9

#         # Decon QTL doesn't want ID in column name; just sample IDs
#         zcat {input[0]} | sed '1d' | gzip > {output[0]}
#         {config[CSVTK]} transpose {input[1]} -T | gzip  > {output[1]}
#         zcat {output[1]} {output[0]} | bgzip > {output[2]}
#         """


# Filter for SNPs to test in dosage file (orignal dosage file is to big to handle in R);
# header is dropped;
# in next step, need to change sample name order to match expression and cell count files;
# SNP IDs should be row names;
rule snps_to_test:
    input:
        gene_snp="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/gene_snp_file.txt",
        dosage="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.tsv.gz",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.snpsToTest.tsv.gz",
    resources:
        mem_gb=4,
        time_min=120,
        num_cores=4,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load htslib
        # join -1 2 -2 1 -o auto --nocheck-order <(sort -k2 {input.gene_snp}) <(sort -k1 <(zcat {input.dosage})) | bgzip > {output[0]}
        # grep -f <(awk '{{print $2}}' {input.gene_snp} ) <(zcat {input.dosage}) | bgzip > {output[0]}
        awk '{{if(NR==FNR){{h[$2] = $1}}else{{if($1 in h){{print $0}}}}}}' {input.gene_snp} <(zcat {input.dosage}) | bgzip > {output[0]}
        """


# Sample names have to be in the same order in dosage, expression, and cell counts
# had problem running the software with dosage file in bgzipped tsv.gz
rule fix_decon_dosage:
    input:
        dosage="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.snpsToTest.tsv.gz",
        sample="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/sample_list.tsv",
        prop="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/cellcounts.tsv",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.snpsToTest.fixed.txt",
    resources:
        mem_gb=4,
        time_min=120,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        module load htslib
        Rscript scripts/fix_decon_dosage.R \
            --dosage {input.dosage} \
            --sample {input.sample} \
            --cellcount {input.prop} \
            --out {output[0]}
        """


# n=629, relatives removed
rule run_decon_qtl:
    input:
        prop="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/cellcounts.tsv",
        expr="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/bulk_expr.tsv",
        dosage="/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.dose.decon.snpsToTest.fixed.txt",
        gene_snp="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/deconv/gene_snp_file.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/decon/deconvolutionResults.csv",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/decon/predictedExpressionLevels.txt",
    resources:
        time_min=120,
        mem_gb=4,
        num_cores=4,
    params:
        outdir="/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/decon",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load java/jre-1.8.0_281
        mkdir -p {params.outdir}
        java -jar {config[DECON-QTL]} \
            --cellcount {input.prop} \
            --expression {input.expr} \
            --genotype {input.dosage} \
            --snpsToTest {input.gene_snp} \
            --outfolder {params.outdir} \
            --outputPredictedExpression
            # --test_run
        """
