#! /usr/bin/env Rscript
suppressMessages(library(argparser))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(wesanderson))

p <- arg_parser("QTL effect size")
p <- add_argument(p, "--pop1", help="pop1 perm sig")
p <- add_argument(p, "--pop2", help="pop2 nominal all assoc")
p <- add_argument(p, "--out", help="out PNG")
p <- add_argument(p, "--pop1_group", help="")
p <- add_argument(p, "--pop2_group", help="")
args <- parse_args(p)

# Read all_assoc files
pop1 <- fread(args$pop1, data.table = F)
pop2 <- fread(args$pop2, data.table = F)

# Gene-snp pairs
pop1 <- pop1 %>% unite("gene_snp", pid, sid, sep = "-", remove = FALSE)
pop2 <- pop2 %>% unite("gene_snp", pid, sid, sep = "-", remove = FALSE)
count <- sum(unique(pop1$gene_snp) %in% unique(pop2$gene_snp))
cat(count, "shared pairs found!\n")

# Merge
shared <- pop1 %>% inner_join(pop2, by = "gene_snp")
dim(shared)
shared <- shared %>% mutate(group = ifelse(fdr <= 0.05, "sig", "non-sig"))

cols <- c("sig" = "#FF0000", "non-sig" = "#F98400")
p <- ggplot(shared, aes(x = slope.x, y = slope.y)) +
    geom_point(aes(color = group), alpha = 1, size = 1) +
    theme_classic() +
    labs(x = paste0(args$pop1_group, " eGene-eQTL effect size"), 
         y = paste0(args$pop2_group, " effect size")) +
#     scale_color_manual(values = wes_palette("Darjeeling1"), name = "") +
    scale_color_manual(values = cols, 
                       name = "", 
                       labels = c(paste0("Nominally significant in ", args$pop2_group), 
                                  paste0("Nominally non-significant in ", args$pop2_group))) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = c(0.35, 0.95)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 1.2) +
    geom_smooth(method = 'lm', aes(color = group), fullrange = TRUE)

ggsave(args$out, p, width = 6, height = 6)
