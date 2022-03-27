library(plinkQC)
package.dir <- find.package('plinkQC')
#indir <- "/u/project/gandalm/shared/GenomicDatasets/FetalBrain/Walker2019_eQTL/geno"
#indir <- "/u/project/gandalm/shared/GenomicDatasets/FetalBrain/Obrien2018_eQTL/geno"
#indir <- "/u/project/gandalm/shared/GenomicDatasets/FetalBrain/HDBR/geno"
indir <- "/u/project/gandalm/shared/GenomicDatasets/FetalBrain/Werling2020_eQTL/geno/"
#indir <- "/u/project/gandalm/shared/GenomicDatasets/FetalBrain/LIBD_1and2/geno/phase1and2/position.filtered"
#indir <- "/u/project/gandalm/shared/GenomicDatasets/FetalBrain/LIBD_1and2/geno/phase2only/position.filtered"
qcdir <- "/u/home/c/cindywen/project-gandalm/isoform_twas/plinkQC/werling"
name <- "PEC_HSB_Yale-UCSF_WGS"
#name <- "libd.phase2only.1M.fixed.final"
path2plink <- "/u/project/gandalm/shared/apps/plink"

# now skipping individual level QC

# Per-marker quality control
fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
path2plink=path2plink,
verbose=TRUE, interactive=FALSE,
showPlinkOutput=TRUE,
do.check_snp_missingness = TRUE, lmissTh = 0.05,
do.check_hwe = TRUE, hweTh = 1e-06, 
do.check_maf = TRUE, mafTh = 0.01)

# save summary graph
pdf("file.pdf")
print(fail_markers)
dev.off()

# overview of failing markers
# overview_marker <- overviewPerMarkerQC(fail_markers, interactive=FALSE)
# Error in `.rowNamesDF<-`(x, value = value) :
#   missing values in 'row.names' are not allowed
