#! /usr/bin/env Rscript

library(tidyverse)
library(data.table)

setwd("~/project-gandalm/isoform_twas/sLDSC/")
# gwa <- fread("/u/project/gandalm/shared/GWAS/SCZ.Pardinas.PGC.2018/CLOZUK_rsFiltered_SNP_P.txt", data.table = F)
gwa <- fread("/u/project/gandalm/shared/GWAS/SCZ.PGC3.2021/wave3_v3/PGC3_SCZ_wave3.european.autosome.public.v3.tsv", data.table = F)
cat("Finished reading GWAS\n")

eqtl_set <- read.table("data/mixed_top_eqtl.txt", header = F, stringsAsFactors = F)
isoqtl_set <- read.table("data/mixed_top_isoqtl_grp_perm.txt", header = F, stringsAsFactors = F)
sqtl_set <- read.table("data/mixed_top_sqtl_grp_perm.txt", header = F, stringsAsFactors = F)

eqtl_set <- eqtl_set %>% inner_join(gwa, by = c("V1"="ID"))
isoqtl_set <- isoqtl_set %>% inner_join(gwa, by = c("V1"="ID"))
sqtl_set <- sqtl_set %>% inner_join(gwa, by = c("V1"="ID"))

eqtl_qq <- data.frame("lobs" = -(log10(sort(eqtl_set$PVAL))), 
                      "lexp" = -(log10( c(1:nrow(eqtl_set)) / (nrow(eqtl_set)+1) )),
                      "group" = "cis-eQTL")
isoqtl_qq <- data.frame("lobs" = -(log10(sort(isoqtl_set$PVAL))), 
                        "lexp" = -(log10( c(1:nrow(isoqtl_set)) / (nrow(isoqtl_set)+1) )),
                        "group" = "cis-isoQTL")
sqtl_qq <- data.frame("lobs" = -(log10(sort(sqtl_set$PVAL))), 
                      "lexp" = -(log10( c(1:nrow(sqtl_set)) / (nrow(sqtl_set)+1) )),
                      "group" = "cis-sQTL")
observed <- sort(gwa$PVAL)
lobs <- -(log10(observed))
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))
gwas_qq <- data.frame("lobs" = lobs, "lexp" = lexp, "group" = "SCZ GWAS")                                          
cat("Finished df prep\n")

qq_df <- rbind(eqtl_qq, isoqtl_qq, sqtl_qq, gwas_qq)
colors <- c("cis-eQTL" = "#648FFF", "cis-isoQTL" = "#DC267F", "cis-sQTL" = "#FFB000", "SCZ GWAS" = "#999999")

ggplot(qq_df, aes(color = group, x = lexp, y = lobs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 1, color = "grey") +
    geom_point(size = 2.5, alpha = 0.8, shape = 20) +
    scale_color_manual(name = "", values = colors) +
    # ggtitle("SCZ GWAS QQ-plot", subtitle = "") +
    xlab("-log10(expected P-value)") +
    ylab("-log10(observed P-value)") +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = c(0.15, 0.85),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14))
ggsave("figures/qqplot.pdf", width =6, height = 6)
