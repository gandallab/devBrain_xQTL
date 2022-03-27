from os.path import join
import os
import numpy as np
import pandas as pd
import sys


configfile: "config.yaml"


"""
rules:
    - extract_gwas_loci: for each GWAS locus, extract sum stats of all tested varaints in this locus (1Mb window around index SNP)
    - filter_1kg_variants: not all GWAS variants are in 1KG, filter for common variants, need to make LD later
    - locus_egene: get list of nominal eGene with at least one of the locus variants is eQTL
    - extract_cis_assoc: for each eGene, extract all cis variants associations
    - zscore: for each eGene, extract the overlapping varaints with GWAS; write zscore file for both GWAS and eQTL; sort variants by order in BIM/coordinates (eQTL sum stats is actually already sorted by coordinates, match with BIM just in case...; write done file for this rule and the above, because output files have gene, which is hard to specify...)
    - ld: calculate LD files for GWAS and eQTL
    - run_ecaviar: max_causal=2
    - analyze_ecaviar: extarct CLPP>0.01 
"""

GWAS_DIC = {s["trait"]: s["file"] for s in config["GWAS_LIST"]}
LOCI_DIC = {s["trait"]: s["locus"] for s in config["LOCI_LIST"]}

GWAS_LOCI_TABLE = pd.DataFrame(
    [(t, l) for t in LOCI_DIC.keys() for l in LOCI_DIC[t]], columns=["trait", "locus"]
)


def get_locus_chr(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return chromosome


def get_locus_start(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    start = LOCI_TABLE.loc[int(wildcards.locus), "START"]
    return start


def get_locus_end(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    end = LOCI_TABLE.loc[int(wildcards.locus), "END"]
    return end


def get_1kg_eur_bim_chr(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
        "LOCUS", drop=True
    )
    chromosome = LOCI_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
        + ".bim"
    )


def get_1kg_eur_bfile_chr(wildcards):
    LOCI_TABLE = pd.read_table(str(wildcards.trait) + "_1Mb.txt").set_index(
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
            "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/analyze.done",
            zip,
            trait=GWAS_LOCI_TABLE.trait.values,
            locus=GWAS_LOCI_TABLE.locus.values,
        ),


rule extract_gwas_loci:
    input:
        "{trait}_1Mb.txt",
        gwas_file=lambda wildcards: GWAS_DIC[wildcards.trait],
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas.txt",
    resources:
        mem_gb=4,
        time_min=60,
    params:
        # chromosome=lambda wildcards: LOCI_TABLE.loc[int(wildcards.locus), "CHR"],
        # start=lambda wildcards: LOCI_TABLE.loc[int(wildcards.locus), "START"],
        # end=lambda wildcards: LOCI_TABLE.loc[int(wildcards.locus), "END"],
        chromosome=get_locus_chr,
        start=get_locus_start,
        end=get_locus_end,
    shell:
        """
        mkdir -p /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/locus_{wildcards.locus}/
        awk '{{if($7=="{params.chromosome}" && $8>={params.start} && $8<={params.end}) print}}' <(zcat {input.gwas_file}) > {output[0]}
        """


rule filter_1kg_variants:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas.txt",
        get_1kg_eur_bim_chr,
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        awk 'FNR==NR{{a[$2]=$2; next}}; $1 in a {{print}}' {input[1]} {input[0]} > {output[0]}
        """


rule locus_egene:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/gwas_in_bim.txt",
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/locus_egene.txt",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --trait {wildcards.trait}
        """


rule extract_cis_assoc:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/gtex.allpairs.txt.gz",
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/locus_egene.txt",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/extract_cis_assoc.done",
    resources:
        mem_gb=6,
        time_min=240,
    shell:
        """
        cat {input[1]} | while read gene
        do
            awk -v a="${{gene}}" '$1 == a {{print}}' <(zcat {input[0]}) > /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/locus_{wildcards.locus}/${{gene}}_all_pairs.txt
        done
        touch {output[0]}
        """


rule zscore:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/extract_cis_assoc.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/zscore.done",
    resources:
        mem_gb=6,
        time_min=240,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        cat /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/locus_{wildcards.locus}/locus_egene.txt | while read gene
        do
            Rscript scripts/zscore.R \
                --locus {wildcards.locus} \
                --gene ${{gene}} \
                --trait {wildcards.trait}
        done
        touch {output[0]}
        """


# matrix is r
rule ld:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/zscore.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/ld.done",
    resources:
        mem_gb=4,
        time_min=60,
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/locus_{wildcards.locus}

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
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/ld.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/ecaviar.done",
    resources:
        mem_gb=8,
        num_cores=2,
        time_min=720,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/locus_{wildcards.locus}

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


# here parallel with loci, not genes; loci with many genes take too long to finish; separately ran some loci with genes in parallel
# while read locus gene;do qsub run_ecaviar.sh $locus $gene ;done < run_ecaviar.txt


rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/ecaviar.done",
    output:
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/locus_{locus}/analyze.done",
    resources:
        mem_gb=4,
        time_min=60,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/locus_{wildcards.locus}

        cat locus_egene.txt | while read gene
        do
            Rscript /u/project/gandalm/cindywen/isoform_twas/colocal/code/scripts/analyze.R \
                --locus {wildcards.locus} \
                --gene ${{gene}} \
                --trait {wildcards.trait}
        done
        touch {output[0]}
        """


# do this manually
# rule concat:
#     input:
#         expand(
#             "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{{trait}}/locus_{locus}/analyze.done",
#             locus=LOCI_TABLE.index,
#         ),
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{trait}/CLPP_sig.txt",
#     resources:
#         mem_gb=4,
#         time_min=60,
#     shell:
#         """
#         cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/{wildcards.trait}/
#         awk 'FNR==1 && NR!=1{{next;}}{{print}}' */*sig.txt > {output[0]}
#         """


# see locuszoom_working.txt
# rule locuszoom:
#     input:
#         "/u/project/gandalm/cindywen/isoform_twas/colocal/results/scz_locus_{locus}/{target_gene}.all_assoc_METAL.txt",
#         "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeGeneOutlier.vcf.gz"
#     output:
#         "/u/project/gandalm/cindywen/isoform_twas/colocal/figures/locuszoom_scz_locus{locus}_{target_gene}_rs7508148.pdf",
#     resources:
#         mem_gb=4,
#         time_min=60,
#     params:
#         "/u/project/gandalm/cindywen/isoform_twas/colocal/figures/locuszoom_scz_locus{locus}_{target_gene}",
#     shell:
#         """
#         . /u/local/Modules/default/init/modules.sh
#         module load R/4.1.0-BIO
#         module load python/2.7.18
#         module load htslib/1.12
#         module load plink/1.90b624

#         # find the FDR closet to 0.05, then signif cutoff is -log10(npval)
#         signif=3.41365
#         /u/project/gandalm/shared/apps/locuszoom/bin/locuszoom \
#             theme=publication \
#             --cache None \
#             --no-date \
#             --plotonly \
#             --gene-table gencode \
#             --build hg19 \
#             --metal {input[0]} \
#             --refsnp rs7508148 \
#             --flank 50kb \
#             --ld-vcf {input[1]}
#             --prefix {params.out} \
#             signifLine=\"${{signif}}\" \
#             signifLineColor=\"gray\" \
#             signifLineWidth=\"2\" \
#             showRecomb=TRUE \
#             width=5 \
#             height=5 \
#             showGenes=FALSE
#         """
