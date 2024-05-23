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

def get_1kg_eur_bfile_chr(wildcards):
    gwas = wildcards.gwas
    TEST_TABLE = pd.read_table(f"{gwas}_loci.tsv").set_index("locus", drop=True)
    chromosome = TEST_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
    )

LOCUS_TABLE_MB = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_gene_list_MB.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_MB.append([gwas, row.gene, row.locus])

LOCUS_TABLE_MB = pd.DataFrame(
    LOCUS_TABLE_MB, 
    columns = ['gwas', 'gene', 'locus']
)

LOCUS_TABLE_THISTLE = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_probe_list.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_THISTLE.append([gwas, row.gene, row.locus])

LOCUS_TABLE_THISTLE = pd.DataFrame(
    LOCUS_TABLE_THISTLE, 
    columns = ['gwas', 'gene', 'locus']
)

LOCUS_TABLE_eqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_1MB_locus_fetal_eqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_eqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_eqtl = pd.DataFrame(
    LOCUS_TABLE_eqtl, 
    columns = ['gwas', 'gene', 'locus']
)

LOCUS_TABLE_isoqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_1MB_locus_fetal_isoqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_isoqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_isoqtl = pd.DataFrame(
    LOCUS_TABLE_isoqtl, 
    columns = ['gwas', 'gene', 'locus']
)

LOCUS_TABLE_sqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_1MB_locus_fetal_sqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_sqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_sqtl = pd.DataFrame(
    LOCUS_TABLE_sqtl, 
    columns = ['gwas', 'gene', 'locus']
)

LOCUS_TABLE_tri1_eqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_fetal_tri1_eqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_tri1_eqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_tri1_eqtl = pd.DataFrame(
    LOCUS_TABLE_tri1_eqtl, 
    columns = ['gwas', 'gene', 'locus']
)
LOCUS_TABLE_tri2_eqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_fetal_tri2_eqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_tri2_eqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_tri2_eqtl = pd.DataFrame(
    LOCUS_TABLE_tri2_eqtl, 
    columns = ['gwas', 'gene', 'locus']
)
LOCUS_TABLE_tri1_isoqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_fetal_tri1_isoqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_tri1_isoqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_tri1_isoqtl = pd.DataFrame(
    LOCUS_TABLE_tri1_isoqtl, 
    columns = ['gwas', 'gene', 'locus']
)
LOCUS_TABLE_tri2_isoqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_fetal_tri2_isoqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_tri2_isoqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_tri2_isoqtl = pd.DataFrame(
    LOCUS_TABLE_tri2_isoqtl, 
    columns = ['gwas', 'gene', 'locus']
)

LOCUS_TABLE_tri1_sqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_fetal_tri1_sqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_tri1_sqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_tri1_sqtl = pd.DataFrame(
    LOCUS_TABLE_tri1_sqtl, 
    columns = ['gwas', 'gene', 'locus']
)
LOCUS_TABLE_tri2_sqtl = []
for gwas in LOCI_DIC.keys():
    df = pd.read_table(f"{gwas}_locus_fetal_tri2_sqtl.txt")
    for i, row in df.iterrows():
        LOCUS_TABLE_tri2_sqtl.append([gwas, row.gene, row.locus])

LOCUS_TABLE_tri2_sqtl = pd.DataFrame(
    LOCUS_TABLE_tri2_sqtl, 
    columns = ['gwas', 'gene', 'locus']
)
rule all:
    input:
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_MB_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_MB.gwas.values,
        #     locus = LOCUS_TABLE_MB.locus.values,
        #     gene = LOCUS_TABLE_MB.gene.values,
        # ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_thistle_ecaviar_col",
        #     zip,
        #     gwas = LOCUS_TABLE_THISTLE.gwas.values,
        #     locus = LOCUS_TABLE_THISTLE.locus.values,
        #     gene = LOCUS_TABLE_THISTLE.gene.values,
        # ),
        expand(
            "../out_1MB/{gwas}/locus{locus}/{gene}_fetal_eqtl.ld",
            zip,
            gwas = LOCUS_TABLE_eqtl.gwas.values,
            locus = LOCUS_TABLE_eqtl.locus.values,
            gene = LOCUS_TABLE_eqtl.gene.values,
        ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_fetal_tri1_eqtl_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_tri1_eqtl.gwas.values,
        #     locus = LOCUS_TABLE_tri1_eqtl.locus.values,
        #     gene = LOCUS_TABLE_tri1_eqtl.gene.values,
        # ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_fetal_tri2_eqtl_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_tri2_eqtl.gwas.values,
        #     locus = LOCUS_TABLE_tri2_eqtl.locus.values,
        #     gene = LOCUS_TABLE_tri2_eqtl.gene.values,
        # ),
        expand(
            "../out_1MB/{gwas}/locus{locus}/{gene}_fetal_isoqtl.ld",
            zip,
            gwas = LOCUS_TABLE_isoqtl.gwas.values,
            locus = LOCUS_TABLE_isoqtl.locus.values,
            gene = LOCUS_TABLE_isoqtl.gene.values,
        ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_fetal_tri1_isoqtl_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_tri1_isoqtl.gwas.values,
        #     locus = LOCUS_TABLE_tri1_isoqtl.locus.values,
        #     gene = LOCUS_TABLE_tri1_isoqtl.gene.values,
        # ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_fetal_tri2_isoqtl_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_tri2_isoqtl.gwas.values,
        #     locus = LOCUS_TABLE_tri2_isoqtl.locus.values,
        #     gene = LOCUS_TABLE_tri2_isoqtl.gene.values,
        # ),
        expand(
            "../out_1MB/{gwas}/locus{locus}/{gene}_fetal_sqtl.ld",
            zip,
            gwas = LOCUS_TABLE_sqtl.gwas.values,
            locus = LOCUS_TABLE_sqtl.locus.values,
            gene = LOCUS_TABLE_sqtl.gene.values,
        ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_fetal_tri1_sqtl_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_tri1_sqtl.gwas.values,
        #     locus = LOCUS_TABLE_tri1_sqtl.locus.values,
        #     gene = LOCUS_TABLE_tri1_sqtl.gene.values,
        # ),
        # expand(
        #     "../out/{gwas}/locus{locus}/{gene}_fetal_tri2_sqtl_ecaviar.done",
        #     zip,
        #     gwas = LOCUS_TABLE_tri2_sqtl.gwas.values,
        #     locus = LOCUS_TABLE_tri2_sqtl.locus.values,
        #     gene = LOCUS_TABLE_tri2_sqtl.gene.values,
        # ),


allowed_annot_fetal = "(eqtl|isoqtl|sqtl)"
allowed_annot_tri = "(tri1_eqtl|tri2_eqtl|tri1_isoqtl|tri2_isoqtl|tri1_sqtl|tri2_sqtl)"


rule ld_MB:
    input:
        "{gwas}_locus_gene_list_MB.txt",
        "../out/{gwas}/locus{locus}/{gene}_MB_snps.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_MB.ld",
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile {params[0]} \
            --r \
            --matrix \
            --extract {input[1]} \
            --out ../out/{wildcards.gwas}/locus{wildcards.locus}/{wildcards.gene}_MB
        """

rule ld_probe:
    input:
        "{gwas}_locus_probe_list.txt",
        "../out/{gwas}/locus{locus}/{gene}_thistle_snps.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_thistle.ld",
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile {params[0]} \
            --r \
            --matrix \
            --extract {input[1]} \
            --out ../out/{wildcards.gwas}/locus{wildcards.locus}/{wildcards.gene}_thistle
        """

rule ld_fetal:
    input:
        "{gwas}_1MB_locus_fetal_{annot}.txt",
        "../out_1MB/{gwas}/locus{locus}/{gene}_fetal_{annot}_snps.txt",
    output:
        "../out_1MB/{gwas}/locus{locus}/{gene}_fetal_{annot}.ld",
        "../out_1MB/{gwas}/locus{locus}/{gene}_fetal_{annot}_gwas.ld",
    params:
        get_1kg_eur_bfile_chr,
    wildcard_constraints:
        annot=allowed_annot_fetal,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile {params[0]} \
            --r \
            --matrix \
            --extract {input[1]} \
            --out ../out_1MB/{wildcards.gwas}/locus{wildcards.locus}/{wildcards.gene}_fetal_{wildcards.annot}_gwas

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
                --r \
                --matrix \
                --extract {input[1]} \
                --out ../out_1MB/{wildcards.gwas}/locus{wildcards.locus}/{wildcards.gene}_fetal_{wildcards.annot}
        """


rule ld_fetal_tri:
    input:
        "{gwas}_locus_fetal_{annot}.txt",
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_snps.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}.ld",
    params:
        get_1kg_eur_bfile_chr,
    wildcard_constraints:
        annot=allowed_annot_tri,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile {params[0]} \
            --r \
            --matrix \
            --extract {input[1]} \
            --out ../out/{wildcards.gwas}/locus{wildcards.locus}/{wildcards.gene}_fetal_{wildcards.annot}
        """

rule ecaviar_MB:
    input:
        "../out/{gwas}/locus{locus}/{gene}_MB.ld",
        "../out/{gwas}/locus{locus}/{gene}_MB_zscore.txt",
        "../out/{gwas}/locus{locus}/{gene}_MB_gwas_zscore.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_MB_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd ../out/{wildcards.gwas}/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_MB_ecaviar \
                -l {wildcards.gene}_MB.ld \
                -z {wildcards.gene}_MB_gwas_zscore.txt \
                -l {wildcards.gene}_MB.ld \
                -z {wildcards.gene}_MB_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule ecaviar_probe:
    input:
        "../out/{gwas}/locus{locus}/{gene}_thistle.ld",
        "../out/{gwas}/locus{locus}/{gene}_thistle_zscore.txt",
        "../out/{gwas}/locus{locus}/{gene}_thistle_gwas_zscore.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_thistle_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd ../out/{wildcards.gwas}/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_thistle_ecaviar \
                -l {wildcards.gene}_thistle.ld \
                -z {wildcards.gene}_thistle_gwas_zscore.txt \
                -l {wildcards.gene}_thistle.ld \
                -z {wildcards.gene}_thistle_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule ecaviar_fetal:
    input:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}.ld",
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_gwas.ld",
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_zscore.txt",
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_gwas_zscore.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_ecaviar_col",
    wildcard_constraints:
        annot=allowed_annot_fetal,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd ../out/{wildcards.gwas}/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_fetal_{wildcards.annot}_ecaviar \
                -l {wildcards.gene}_fetal_{wildcards.annot}_gwas.ld \
                -z {wildcards.gene}_fetal_{wildcards.annot}_gwas_zscore.txt \
                -l {wildcards.gene}_fetal_{wildcards.annot}.ld \
                -z {wildcards.gene}_fetal_{wildcards.annot}_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule ecaviar_fetal_tri:
    input:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}.ld",
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_zscore.txt",
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_gwas_zscore.txt",
    output:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_ecaviar_col",
    wildcard_constraints:
        annot=allowed_annot_tri,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd ../out/{wildcards.gwas}/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_fetal_{wildcards.annot}_ecaviar \
                -l {wildcards.gene}_fetal_{wildcards.annot}.ld \
                -z {wildcards.gene}_fetal_{wildcards.annot}_gwas_zscore.txt \
                -l {wildcards.gene}_fetal_{wildcards.annot}.ld \
                -z {wildcards.gene}_fetal_{wildcards.annot}_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule ecaviar_sig_fetal:
    input:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_ecaviar_col",
    output:
        "../out/{gwas}/locus{locus}/{gene}_fetal_{annot}_ecaviar.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/analyze.R \
            --gwas {wildcards.gwas} \
            --locus {wildcards.locus} \
            --annot fetal_{wildcards.annot} \
            --gene {wildcards.gene}
        touch {output[0]}
        """