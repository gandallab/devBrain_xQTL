#! /usr/bin/env Rscript
suppressMessages(library(argparser))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
# suppressMessages(library(wesanderson))

p <- arg_parser("Nominal QTL effect size")
p <- add_argument(p, "--pop1", help="pop1 nominal all assoc")
p <- add_argument(p, "--pop2", help="pop2 nominal all assoc")
p <- add_argument(p, "--out", help="out PNG")
p <- add_argument(p, "--pop12_group", help="string")
p <- add_argument(p, "--pop1_group", help="")
p <- add_argument(p, "--pop2_group", help="")
args <- parse_args(p)

# Read all_assoc files
all_assoc_pop1_nominal <- fread(args$pop1, data.table = F)
all_assoc_pop2_nominal <- fread(args$pop2, data.table = F)

# # Comment these for full file; first test one chunk
# colnames(all_assoc_pop1_nominal) <- c("pid","sid","dist","npval","slope")
# colnames(all_assoc_pop2_nominal) <- c("pid","sid","dist","npval","slope")
# all_assoc_pop1_nominal$fdr <- p.adjust(all_assoc_pop1_nominal$npval, method = 'fdr')
# all_assoc_pop2_nominal$fdr <- p.adjust(all_assoc_pop2_nominal$npval, method = 'fdr')

# Gene-snp pairs
all_assoc_pop1_nominal <- all_assoc_pop1_nominal %>% unite("gene_snp", pid, sid, sep = "-", remove = FALSE)
all_assoc_pop2_nominal <- all_assoc_pop2_nominal %>% unite("gene_snp", pid, sid, sep = "-", remove = FALSE)

# Significant subset
pop1_sig_nominal <- all_assoc_pop1_nominal %>% filter(fdr <= 0.05)
pop2_sig_nominal <- all_assoc_pop2_nominal %>% filter(fdr <= 0.05)

# Significant gene-QTL assoc in both pop1 and pop2, with their slopes
pop1_pop2_beta <- pop1_sig_nominal %>% 
    inner_join(pop2_sig_nominal, by = "gene_snp") %>%
    select(gene_snp, slope.x, slope.y)

# Significant in pop1 only
pop1_only_beta <- pop1_sig_nominal %>% 
    filter(!gene_snp %in% pop2_sig_nominal$gene_snp) %>% 
    select(gene_snp, slope)
names(pop1_only_beta)[2] <- "slope.x"
pop1_only_beta <- pop1_only_beta %>% 
    left_join(all_assoc_pop2_nominal, by = "gene_snp") %>% 
    mutate(slope.y = slope) %>% 
    select(gene_snp, slope.x, slope.y)

# Significant in pop2 only
pop2_only_beta <- pop2_sig_nominal %>% 
    filter(!gene_snp %in% pop1_sig_nominal$gene_snp) %>% 
    select(gene_snp, slope)
names(pop2_only_beta)[2] <- "slope.y"
pop2_only_beta <- pop2_only_beta %>% 
    left_join(all_assoc_pop1_nominal, by = "gene_snp") %>% 
    mutate(slope.x = slope) %>% 
    select(gene_snp, slope.x, slope.y)

# 
df <- rbind(pop1_pop2_beta, pop1_only_beta, pop2_only_beta)
df$group <- c(rep(args$pop12_group, nrow(pop1_pop2_beta)),
              rep(args$pop1_group, nrow(pop1_only_beta)),
              rep(args$pop2_group, nrow(pop2_only_beta)))
df$group <- factor(df$group, levels = c(args$pop1_group, args$pop2_group, args$pop12_group))

cols <- c( "group1" = "#00A08A", 
           "group2" = "#F2AD00", 
           "group12" = "#FF0000" )
names(cols) <- c(args$pop1_group, args$pop2_group, args$pop12_group)

p <- ggplot(df, aes(x = slope.x, y = slope.y)) +
    geom_point(aes(color = group), alpha = .4, size = .5) +
    theme_classic() +
    labs(x = paste0(substring(args$pop1_group, 1, 3), " effect size"), 
         y = paste0(substring(args$pop2_group, 1, 3), " effect size")) +
    # scale_color_manual(values = wes_palette("Darjeeling1"), name = "Significant cis-eQTL") +
    scale_color_manual(values = cols, name = "Significant cis-eQTL") +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 1.2) +
    geom_smooth(method = 'lm', data = subset(df, group == args$pop12_group), aes(color = group), fullrange = TRUE)
ggsave(args$out, p, width = 8, height = 6)
