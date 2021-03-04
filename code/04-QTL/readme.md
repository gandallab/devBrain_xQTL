## cis-e/iso/sQTL mapping 
- basics.ipynb: plot data basics, age, sex, infer NA sex
- cis-eQTL and cis-isoQTL follow a similar pipeline. See Snakefile
```
rules:
    expr_prep: expression filter, VST normalization, outlier detection, ComBat, create bed file
    ancestry_expr: subset ancestries
    cov: generate covariates files
    ancestry_cov:
    fastqtl_nominal:
    ancestry_fastqtl_nominal:
    call_nominal: identify significant nominal associations
```
  