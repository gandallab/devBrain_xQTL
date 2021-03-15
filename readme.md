# Fetal brain mega QTL
### Quantifications
#### Gene level (working/gene/)
- [x] Estimated counts ("countsFromAbundance="lengthScaledTPM"): gene.noVersion.scaled.counts.tsv
- [x] TPM: gene.noVersion.TPM.tsv
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: gene.counts.scaled.normalized.bed.gz
#### Isoform level (working/transcript/)
- [x] Estimated counts ("countsFromAbundance="lengthScaledTPM"): tx.counts.scaled.tsv
- [x] TPM: tx.TPM.tsv
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: tx.counts.scaled.normalized.bed.gz
#### Splicing (working/splicing/)
- [ ] Normalized intron excision ratios: 

### Genotype (working/geno/)
- [x] Mixed, EUR, AMR, AFR: QC-ed and imputed genotype data
- [x] Imputed, filtered for R2>.3, intersection across studies, exlcude duplicates, map to rsID (still in hg38): merge.reheader.chr*_rsid
*** Please note some TOPMED variants are not mapped to rsID. They are still in hg38 coordinates as ID (chr:pos). Their coordinates are are hg19 in the final genotype file. ***

### Covariates (working/covariates/)
- [x] Compiled PicardTools QC metrics: picard_QC_compiled.tsv
- [x] Metadata with inferred sex: metadata_inferSex.tsv 
- [x] Genotype PCs: shared/EUR/AMR/AFR

### cis-eQTL (output/cis_eQTL_nominal_and_permutation)
* Mixed
- [x] List of eGene from permutations and top eQTL per eGene: 
- [ ] Full list of eQTL from permutations: 
- [ ] List of all eGene and top eQTL per rank from conditional pass:
* EUR
- [x] List of eGene from permutations and top eQTL per eGene: 
- [ ] Full list of eQTL from permutations: 
- [ ] List of all eGene and top eQTL per rank from conditional pass:
* AMR
- [x] List of eGene from permutations and top eQTL per eGene: 
- [ ] Full list of eQTL from permutations: 
* AFR
- [x] List of eGene from permutations and top eQTL per eGene: 
- [ ] Full list of eQTL from permutations: 

### cis-isoQTL (output/cis_isoQTL_nominal_and_permutation)
* Mixed
- [ ] List of isoTx from permutations and top isoQTL per isoTx: 
- [ ] Full list of isoQTL from permutations: 
* EUR
- [ ] List of isoTx from permutations and top isoQTL per isoTx: 
- [ ] Full list of isoQTL from permutations: 
* AMR
- [ ] List of isoTx from permutations and top isoQTL per isoTx: 
- [ ] Full list of isoQTL from permutations: 
* AFR
- [ ] List of isoTx from permutations and top isoQTL per isoTx: 
- [ ] Full list of isoQTL from permutations: 

### cis-sQTL (output/cis_sQTL_nominal_and_permutation)
* Mixed
- [ ] List of sQTL/intron from permutations and top sSNP per sQTL/intron:
- [ ] Full list of sQTL from permutations:
- [ ] List of all intron-SNP pair from nominal pass: 
* EUR
* AMR
* AFR
