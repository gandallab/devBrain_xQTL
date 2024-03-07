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

# eQTL slope correlation
shared_sig <- shared %>% filter(group == "sig")
cor <- cor(shared_sig$slope.x, shared_sig$slope.y, method = "spearman")

# TODO: geom point alpha = .8, size = 1.2, shape = 20
# check qq plot code
cols <- c("sig" = "#046C9A", "non-sig" = "#ABDDDE")
p <- ggplot(shared, aes(x = slope.x, y = slope.y)) +
    geom_point(aes(color = group), alpha = .6, size = .8) +
    theme_classic() +
    labs(x = paste0(args$pop1_group, " effect size"), 
         y = paste0(args$pop2_group, " effect size"),
	 title = paste0(args$pop1_group, "-", args$pop2_group, " eQTL"),
	 subtitle = paste0("Spearman correlation = ", cor)) +
#     scale_color_manual(values = wes_palette("Darjeeling1"), name = "") +
    scale_color_manual(values = cols, 
                       name = "", 
                       labels = c(paste0("Nominally significant in ", args$pop2_group), 
                                  paste0("Nominally non-significant in ", args$pop2_group))) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
	      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 16, hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = c(0.35, 0.95)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey", size = 1.2) +
    geom_smooth(method = 'lm', aes(color = group), fullrange = TRUE) 
    # geom_rug(data = shared[shared$group == "sig",], alpha = 0.8, color = "#046C9A", aes(x = NULL)) +
    # geom_rug(data = shared[shared$group == "non-sig",], alpha = 0.8, color = "#ABDDDE", aes(x = NULL))

ggsave(args$out, p, width = 6, height = 6)
