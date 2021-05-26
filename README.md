# Fetal brain mega QTL (05-26-2021)
## Quantifications
### Gene (`working/gene/`)
- [x] Estimated counts (`"countsFromAbundance="lengthScaledTPM"`): `gene.noVersion.scaled.counts.tsv`
- [x] TPM: `gene.noVersion.TPM.tsv`
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: `gene.counts.scaled.normalized.bed.gz`
### Isoform (`working/transcript/`)
- [x] Estimated counts (`"countsFromAbundance="lengthScaledTPM"`): `tx.counts.scaled.tsv`
- [x] TPM: `tx.TPM.tsv`
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: `tx.counts.scaled.normalized.bed.gz`
### Splicing (`working/splicing/`)
- [x] Intron excision ratios: filtered, normalized, ComBat: `leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.bed.gz`

## Genotype (`working/geno/`)
- [x] Imputed, filtered for R2>.3, intersection across studies, exlcude duplicates, map to rsID (still in hg38): `merge.reheader.chr*_rsid` (For Chao Chen's EAS analysis; TODO: make filtered EAS file)
***Please note some TOPMED variants are not mapped to rsID. They are still in hg38 coordinates as ID (chr:pos). Their coordinates are hg19 in the final genotype file.***

## Covariates (`working/covariates/`)

## cis-eQTL (`output/cis_eQTL_nominal_permutation_conditional_finemapping/`)
* Combined (`shared_90HCP/`)
- [x] Nominal pass: eGene and eQTL: `significant_assoc.txt`
- [x] Permutation pass: eGene and top eQTL: `sig_pheno_info.txt`
- [ ] Permutation pass: full list of eQTL with p-value passing nominal threshold: 
- [x] Conditional pass: top eQTL per rank: `conditional_top_variants.txt`
- [x] SuSiE finemapping; all variants in non-low purity CS: `mixed_ciseqtl_90hcp_perm_purity_filtered.txt.gz`
* EUR (`EUR_50HCP/`)
- [x] Nominal pass: eGene and eQTL: `significant_assoc.txt`
- [x] Permutation pass: eGene and top eQTL: `sig_pheno.txt`
* AMR (`AMR_15HCP/`)
- [x] Nominal pass: eGene and eQTL: `significant_assoc.txt`
- [x] Permutation pass: eGene and top eQTL: `sig_pheno.txt`
* AFR (`AFR_25HCP/`)
- [x] Nominal pass: eGene and eQTL: `significant_assoc.txt`
- [x] Permutation pass: eGene and top eQTL: `sig_pheno.txt`

## cis-isoQTL (`output/cis_isoQTL_nominal_permutation_conditional_finemapping/`)
* Combined (`shared_70HCP/`)
- [x] Nominal pass: isoform and isoQTL: `significant_assoc.txt`
- [x] Permutation pass: isoform, top isoQTL, isoGene: `sig_pheno_gene_info.txt`
- [ ] Permutation pass: full list of isoQTL with p-value passing nominal threshold: 
- [x] Conditional pass: top isoQTL per rank: `conditional_top_variants.txt`
- [x] SuSiE finemapping; all variants in non-low purity CS: `mixed_cisisoqtl_70hcp_perm_purity_filtered.txt.gz`
* EUR (`EUR_60HCP/`)
- [x] Nominal pass: isoform and isoQTL: `significant_assoc.txt`
- [x] Permutation pass: isoform, top isoQTL, isoGene: `sig_pheno_gene.txt`
* AMR (`AMR_15HCP/`)
- [x] Nominal pass: isoform and isoQTL: `significant_assoc.txt`
- [x] Permutation pass: isoform, top isoQTL, isoGene: `sig_pheno_gene.txt`
* AFR (`AFR_20HCP/`)
- [x] Nominal pass: isoform and isoQTL: `significant_assoc.txt`
- [x] Permutation pass: isoform, top isoQTL, isoGene: `sig_pheno_gene.txt`

## cis-sQTL (`output/cis_sQTL_nominal_permutation_conditional_finemapping/`)
* Combined (`shared_35HCP/`)
- [x] Nominal pass: intron and sQTL: `significant_assoc.txt`
- [x] Permutation pass: intron, top sQTL, sGene: `sig_pheno_gene_info.txt`
- [ ] Permutation pass: full list of sQTL with p-value passing nominal threshold: 
- [ ] Conditional pass: top sQTL per rank: 
- [ ] SuSiE finemapping; all variants in non-low purity CS: 
* EUR (`EUR_25HCP/`)
- [x] Nominal pass: intron and sQTL: `significant_assoc.txt`
- [x] Permutation pass: intron, top sQTL, sGene:: `sig_pheno_gene.txt`
* AMR (`AMR_10HCP/`)
- [x] Nominal pass: intron and sQTL: `significant_assoc.txt`
- [x] Permutation pass: intron, top sQTL, sGene:: `sig_pheno_gene.txt`
* AFR (`AFR_12HCP/`)
- [x] Nominal pass: intron and sQTL: `significant_assoc.txt`
- [x] Permutation pass: intron, top sQTL, sGene:: `sig_pheno_gene.txt`

## trans-eQTL (`output/trans_eQTL/`)