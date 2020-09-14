#https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html

data=./data.hg19
#ref=/u/project/gandalm/shared/refGenomes/1000genomes/chrs/kgp
#ref=/u/home/c/cindywen/project-gandalm/isoform_twas/genotype/ref/kgp
ref=/u/project/gandalm/shared/refGenomes/g1000/Phase3_ALL/ALL.autosomes.phase3
outdir=./pca
mkdir -p ${outdir}
mkdir -p ${outdir}/plink_log

# 0. In tutorial, but not doing: filter reference and study data for non A-T or G-C SNPs

# 1. Prune study data
# Filter for variants in LD with r^2>.2 in a 50kb window
plink --bfile ${data} \
--indep-pairwise 50 5 0.2 \
--allow-extra-chr \
--out ${outdir}/data 
mv ${outdir}/data.log ${outdir}/plink_log/

plink --bfile ${data} \
--extract ${outdir}/data.prune.in \
--allow-extra-chr \
--make-bed \
--out ${outdir}/data.pruned
mv ${outdir}/data.pruned.log ${outdir}/plink_log/

# 2. Filter reference data for the same SNP set as in study
# Reduce reference to the same size
plink --bfile ${ref} \
--extract ${outdir}/data.prune.in \
--make-bed \
--out ${outdir}/ref.pruned
mv ${outdir}/ref.pruned.log ${outdir}/plink_log/

# 3. Check and correct chromosome mismatch
# Check that the variant IDs of the reference data have the same chromosome ID as the study data, 
# update the pruned reference dataset
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
($2 in a && a[$2] != $1) {print a[$2],$2}' \
${outdir}/data.pruned.bim ${outdir}/ref.pruned.bim | \
sed -n '/^[XY]/!p' > ${outdir}/ref.toUpdateChr

## remove duplicates in ref.pruned
#plink --bfile ${outdir}/ref.pruned --write-snplist --out ${outdir}/ref.pruned
#cat ${outdir}/ref.pruned.snplist | sort | uniq -d > ${outdir}/ref.pruned.duplicated.snplist
#plink --bfile ${outdir}/ref.pruned \
#    --exclude ${outdir}/ref.pruned.duplicated.snplist --make-bed \
#    --out ${outdir}/ref.pruned.duplicated_excluded
#

# Correct chromosome mismatch
plink --bfile ${outdir}/ref.pruned \
--update-chr ${outdir}/ref.toUpdateChr 1 2 \
--make-bed \
--out ${outdir}/ref.updateChr
mv ${outdir}/ref.updateChr.log ${outdir}/plink_log/

# 4. Check and correct Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
($2 in a && a[$2] != $4) {print a[$2],$2}' \
${outdir}/data.pruned.bim ${outdir}/ref.pruned.bim > \
${outdir}/ref.toUpdatePos

# Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
${outdir}/data.pruned.bim ${outdir}/ref.pruned.bim > \
${outdir}/ref.toFlip

# Update position and allele flips
plink --bfile ${outdir}/ref.updateChr \
--update-map ${outdir}/ref.toUpdatePos 1 2 \
--flip ${outdir}/ref.toFlip \
--make-bed \
--out ${outdir}/ref.flipped
mv ${outdir}/ref.flipped.log ${outdir}/plink_log/

# 5. Remove mismatches in ref
# Any alleles not matching after allele flipping are identified and removed from the reference
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
${outdir}/data.pruned.bim ${outdir}/ref.flipped.bim > \
${outdir}/ref.mismatch

plink --bfile ${outdir}/ref.flipped \
--exclude ${outdir}/ref.mismatch \
--make-bed \
--out ${outdir}/ref.clean
mv ${outdir}/ref.clean.log ${outdir}/plink_log/

# 6. Merge study genotypes and ref
plink --bfile ${outdir}/data.pruned \
--bmerge ${outdir}/ref.clean.bed \
${outdir}/ref.clean.bim \
${outdir}/ref.clean.fam \
--allow-extra-chr \
--make-bed \
--out ${outdir}/data.merge.ref

# Will encounter error if multiallelic not removed

# # remove multiallelic snps from missnp list
# # different alleles in data pruned and ref.clean.bed
# plink --bfile ${outdir}/ref.clean \
# --exclude ${outdir}/data.merge.ref-merge.missnp \
# --make-bed \
# --out ${outdir}/ref.clean.missnp.removed
# 
# plink --bfile ${outdir}/data.pruned \
# --exclude ${outdir}/data.merge.ref-merge.missnp \
# --make-bed \
# --allow-extra-chr \
# --out ${outdir}/data.pruned.missnp.removed
# 
# 
# plink --bfile ${outdir}/data.pruned.missnp.removed \
# --bmerge ${outdir}/ref.clean.missnp.removed.bed \
# ${outdir}/ref.clean.missnp.removed.bim \
# ${outdir}/ref.clean.missnp.removed.fam \
# --allow-extra-chr \
# --make-bed \
# --out ${outdir}/data.merge.ref.clean
# 
# #mv ${outdir}/data.merge.ref.clean.log ${outdir}/plink_log/
# 


# 6. PCA on the merged data
plink --bfile ${outdir}/data.merge.ref \
--pca \
--allow-extra-chr \
--out ${outdir}/data.ref
#mv ${outdir}/data.ref.log ${outdir}/plink_log
mv ${outdir}/*.log ${outdir}/plink_log
