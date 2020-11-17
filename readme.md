# Fetal brain mega QTL

[Google Drive](https://drive.google.com/drive/u/1/folders/1-O4BE5-_3xmMhwoX_DKbSwi2hZw7DgGw)

### Quantifications
#### Gene level 
  * Estimated counts ("countsFromAbundance="lengthScaledTPM"): working/gene.noVersion.scaled.counts.tsv
  * TPM: working/gene.noVersion.TPM.tsv
#### Isoform level
  * Estimated counts ("countsFromAbundance="lengthScaledTPM"): working/tx.counts.scaled.tsv
  * TPM: working/tx.TPM.tsv
#### Splicing
  * Normalized intron excision ratios: working/fetalALL_perind.counts.noChr.gz.qqnorm_allChr_subj.tsv

### Covariates
* Compiled PicardTools QC metrics: working/picard_QC_compiled.tsv
* Metadata with inferred sex: working/metadata_inferSex.tsv. 
* Mixed ancestry genotype PC: working/mixed.data.ref.eigenvec
* EUR geontype PC: working/eur.data.ref.eigenvec
* AMR geontype PC: working/amr.data.ref.eigenvec
* AFR geontype PC: working/afr.data.ref.eigenvec


**For all QTL summary statistics, please note: some variant ID are in hg38 chr:pos:REF:ALT format. Please use the last three columns for hg19 coordinates.**
### eQTL summary statistics 
* output/eQTL/
#### Full Dataset
  * List of eGene from permutations and top eQTL per eGene: eqtl_mixed_permutations_90hcp_egenes_with_coord.txt
  * Full list of eQTL from permutations: eqtl_mixed_permutations_90hcp_eqtl_chunk_concat_with_coord.txt
  * List of eGene from permutations with ENSG ID, gene name and gene type: eqtl_mixed_permutations_90hcp_egene_name_type.txt
  * List of all eGene-eQTL pair from nominals pass: eqtl_mixed_nominals_90hcp_egenes_eqtl_with_coord.txt
#### Ancestry specific
  1. EUR
  * List of eGene from permutations and top eQTL per eGene: eqtl_eur_permutations_60hcp_egenes_with_coord.txt
  * Full list of eQTL from permutations: **TODO**
  * List of eGene from permutations with ENSG ID, gene name and gene type: **TODO**
  * List of all eGene-eQTL pair from nominals pass: eqtl_eur_nominals_60hcp_egenes_eqtl_with_coord.txt
  2. AMR
  * List of eGene from permutations and top eQTL per eGene: eqtl_amr_permutations_20hcp_egenes_with_coord.txt
  * Full list of eQTL from permutations: **TODO**
  * List of eGene from permutations with ENSG ID, gene name and gene type: **TODO**
  * List of all eGene-eQTL pair from nominals pass: eqtl_amr_nominals_20hcp_egenes_eqtl_with_coord.txt
  3. AFR
  * List of eGene from permutations and top eQTL per eGene: eqtl_afr_permutations_25hcp_egenes_with_coord.txt
  * Full list of eQTL from permutations: **TODO**
  * List of eGene from permutations with ENSG ID, gene name and gene type: **TODO**
  * List of all eGene-eQTL pair from nominals pass: eqtl_afr_nominals_25hcp_egenes_eqtl_with_coord.txt


### isoQTL summary statistics
* output/isoQTL/
#### Full Dataset
  * List of isoTx from permutations and top isoQTL per isoTx: isoqtl_mixed_permutations_80hcp_isotx_with_coord.txt
  * Full list of isoQTL from permutations: **TODO**
  * List of top isoQTL per isoGene: isoqtl_mixed_permutations_80hcp_top_isoqtl_per_isogene_with_coord.txt
  * List of isoGene from permutations with ENSG ID, gene name and gene type: isoqtl_mixed_permutations_80hcp_isogene_name_type.txt
  * List of all isoTx-isoQTL pair from nominals pass: isoqtl_mixed_nominals_80hcp_isotx_isoqtl_with_coord.txt
#### Ancestry specific: 
  1. EUR
  * List of isoTx from permutations and top isoQTL per isoTx: isoqtl_eur_permutations_60hcp_isotx_with_coord.txt
  * Full list of isoQTL from permutations: **TODO**
  * List of top isoQTL per isoGene: **TODO**
  * List of isoGene from permutations with ENSG ID, gene name and gene type: **TODO**
  * List of all isoTx-isoQTL pair from nominals pass: isoqtl_eur_nominals_60hcp_isotx_isoqtl_with_coord.txt
  2. AMR
  * List of isoTx from permutations and top isoQTL per isoTx: isoqtl_amr_permutations_10hcp_isotx_with_coord.txt
  * Full list of isoQTL from permutations: **TODO**
  * List of top isoQTL per isoGene: **TODO**
  * List of isoGene from permutations with ENSG ID, gene name and gene type: **TODO**
  * List of all isoTx-isoQTL pair from nominals pass: isoqtl_amr_nominals_10hcp_isotx_isoqtl_with_coord.txt
  3. AFR
  * List of isoTx from permutations and top isoQTL per isoTx: isoqtl_afr_permutations_20hcp_isotx_with_coord.txt
  * Full list of isoQTL from permutations: **TODO**
  * List of top isoQTL per isoGene: **TODO**
  * List of isoGene from permutations with ENSG ID, gene name and gene type: **TODO**
  * List of all isoTx-isoQTL pair from nominals pass: isoqtl_afr_nominals_20hcp_isotx_isoqtl_with_coord.txt
 
 
### sQTL summary statistics
* output/sQTL/
#### Full Dataset:
  * List of sQTL/intron from permutations and top sSNP per sQTL/intron: sqtl_mixed_permutations_40hcp_sqtl_with_coord.txt
  * List of all sQTL-sSNP pair from nominal pass: sqtl_mixed_nominals_40hcp_sqtl_ssnp_with_coord.txt
#### Ancestry specific: 
  1. EUR
  * List of all sQTL-sSNP pair from nominal pass: sqtl_eur_nominals_25hcp_sqtl_ssnp.txt (TODO: add hg19 coordinates?)
  2. AFR
  * List of all sQTL-sSNP pair from nominal pass: sqtl_amr_nominals_6hcp_sqtl_ssnp.txt (TODO: add hg19 coordinates?)
  3. AMR
  * List of all sQTL-sSNP pair from nominal pass: sqtl_afr_nominals_14hcp_sqtl_ssnp.txt (TODO: add hg19 coordinates?)
  
  

