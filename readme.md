# Fetal brain mega QTL
### Quantifications
#### Gene level 
- [x] Estimated counts ("countsFromAbundance="lengthScaledTPM"): gene.noVersion.scaled.counts.tsv
- [x] TPM: gene.noVersion.TPM.tsv
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: gene.counts.scaled.normalized.bed.gz
- [x] QTLtools BED file
#### Isoform level
- [x] Estimated counts ("countsFromAbundance="lengthScaledTPM"): tx.counts.scaled.tsv
- [x] TPM: tx.TPM.tsv
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: tx.counts.scaled.normalized.bed.gz
#### Splicing (re-running STAR and Leafcutter)
- [ ] Normalized intron excision ratios: 
### Genotype
- [x] Mixed, EUR, AMR, AFR: QC-ed and imputed genotype data
- [x] Mixed, EUR, AMR, AFR: genotype PCs
- [x] Imputed, filtered for R2>.3, and intersection across studies: merge.reheader.vcf.gz
### Covariates
- [x] Compiled PicardTools QC metrics: picard_QC_compiled.tsv
- [x] Metadata with inferred sex: metadata_inferSex.tsv. 
- [x] Mixed ancestry genotype PC: mixed.data.ref.eigenvec
- [x] EUR/AMR/AFR geontype PC: eur/amr/afr.data.ref.eigenvec
- [x] Covariates files for QTL mapping: shared/EUR/AMR/AFR

***For all QTL summary statistics, please note: some variant ID are in hg38 chr:pos:REF:ALT format. Please use the last three columns for hg19 coordinates.***
### cis-eQTL 
#### Mixed
- [x] List of eGene from permutations and top eQTL per eGene: eqtl_mixed_permutations_90hcp_egenes_with_coord.txt
- [x] Full list of eQTL from permutations: eqtl_mixed_permutations_90hcp_eqtl_chunk_concat_with_coord.txt
- [x] List of eGene from permutations with ENSG ID, gene name and gene type: eqtl_mixed_permutations_90hcp_egene_name_type.txt
- [ ] List of all gene-SNP pair from nominals pass: 
- [ ] List of all eGene and top eQTL per rank from conditional pass:
#### Ancestries
* EUR
- [x] List of eGene from permutations and top eQTL per eGene: eqtl_eur_permutations_60hcp_egenes_with_coord.txt
- [ ] Full list of eQTL from permutations: 
- [ ] List of eGene from permutations with ENSG ID, gene name and gene type: 
- [ ] List of all gene-SNP pair from nominals pass: 
* AMR
- [x] List of eGene from permutations and top eQTL per eGene: eqtl_amr_permutations_20hcp_egenes_with_coord.txt
- [ ] Full list of eQTL from permutations: 
- [ ] List of eGene from permutations with ENSG ID, gene name and gene type: 
- [ ] List of all gene-SNP pair from nominals pass: 
* AFR
- [x] List of eGene from permutations and top eQTL per eGene: eqtl_afr_permutations_25hcp_egenes_with_coord.txt
- [ ] Full list of eQTL from permutations: 
- [ ] List of eGene from permutations with ENSG ID, gene name and gene type: 
- [ ] List of all gene-SNP pair from nominals pass: 
### cis-isoQTL 
#### Mixed
- [x] List of isoTx from permutations and top isoQTL per isoTx: isoqtl_mixed_permutations_80hcp_isotx_with_coord.txt
- [ ] Full list of isoQTL from permutations: 
- [x] List of top isoQTL per isoGene: isoqtl_mixed_permutations_80hcp_top_isoqtl_per_isogene_with_coord.txt
- [x] List of isoGene from permutations with ENSG ID, gene name and gene type: isoqtl_mixed_permutations_80hcp_isogene_name_type.txt
- [ ] List of all Tx-SNP pair from nominals pass: 
#### Ancestries
* EUR
- [x] List of isoTx from permutations and top isoQTL per isoTx: isoqtl_eur_permutations_60hcp_isotx_with_coord.txt
- [ ] Full list of isoQTL from permutations: 
- [ ] List of top isoQTL per isoGene:  
- [ ] List of isoGene from permutations with ENSG ID, gene name and gene type: 
- [ ] List of all Tx-SNP pair from nominals pass: 
* AMR
- [x] List of isoTx from permutations and top isoQTL per isoTx: isoqtl_amr_permutations_10hcp_isotx_with_coord.txt
- [ ] Full list of isoQTL from permutations: 
- [ ] List of top isoQTL per isoGene:  
- [ ] List of isoGene from permutations with ENSG ID, gene name and gene type: 
- [ ] List of all Tx-SNP pair from nominals pass: 
* AFR
- [x] List of isoTx from permutations and top isoQTL per isoTx: isoqtl_afr_permutations_20hcp_isotx_with_coord.txt
- [ ] Full list of isoQTL from permutations: 
- [ ] List of top isoQTL per isoGene:  
- [ ] List of isoGene from permutations with ENSG ID, gene name and gene type: 
- [ ] List of all Tx-SNP pair from nominals pass: 
### cis-sQTL (Re-running STAR and Leafcutter)
#### Mixed
- [ ] List of sQTL/intron from permutations and top sSNP per sQTL/intron:
- [ ] Full list of sQTL from permutations:
- [ ] List of all intron-SNP pair from nominal pass: 
#### Ancestries
* EUR
* AMR
* AFR
