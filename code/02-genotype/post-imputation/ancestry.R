# Plot genotype PCA and infer data ancestry

#library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("~/Desktop/final/data/")
#pcainput = read.table("data.ref.eigenvec")
pcainput <- read.table("ancestry/afr.ref.eigenvec")
onekgpop <- read.table("kgp.pop", header = TRUE)

pcainlist <- pcainput$V2 %in% onekgpop$IID
onekglist <- onekgpop$IID %in% pcainput$V2
# check ref samples have same IID in kgp.pop and merged data
sum(as.vector(pcainput[pcainlist,]$V2) != as.vector(onekgpop[onekglist,]$IID))

################## 1. Assign populations
pcainput$suppop <- "Data"
pcainput$Pop <- "Data"
onekgpop$suppop <- NA
for(i in 1: dim(onekgpop)[1]) {
  if (onekgpop[i,"Pop"]%in% c("GBR", "IBS", "CEU", "TSI", "FIN")) {
    onekgpop[i,"suppop"] <- "EUR"
  }
  if (onekgpop[i,"Pop"]%in% c("MXL", "PUR", "CLM", "PEL")) {
    onekgpop[i,"suppop"] <- "AMR"
  }
  if (onekgpop[i,"Pop"]%in% c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")) {
    onekgpop[i,"suppop"] <- "AFR"
  }
  if (onekgpop[i,"Pop"]%in% c("GIH", "PJL", "BEB", "STU", "ITU")) {
    onekgpop[i,"suppop"] <- "SEA"
  }
  if (onekgpop[i,"Pop"]%in% c("CHB", "JPT", "CHS", "CDX", "KHV")) {
    onekgpop[i,"suppop"] <- "EA"
  }
}
pcainput[pcainlist,]$Pop <- as.vector(onekgpop[onekglist,]$Pop)
pcainput[pcainlist,]$suppop <- as.vector(onekgpop[onekglist,]$suppop)

# color_scale2 <- c("#696969", "#8d43b0", "#74aff3", "#e9d69a", "#53c6ef", "#ecbe8f")
# 
# color_scale <- c("#696969", "#74aff3", "#e9d69a", "#53c6ef", "#ecbe8f", "#5ddff2",
#                 "#e8ab90", "#5ac7dc", "#ecabb7", "#82ebde", "#e0add9", "#a7eabe",
#                 "#b7b3ea", "#cae8ad", "#8dc0f2", "#bfbb81", "#d8c5f3", "#9ac48e",
#                 "#9db3d6", "#e4e9b7", "#a6ceee", "#c3bc96", "#97e0eb", "#86ceac",
#                 "#bcedd8",  "#6ecac6")

poplevels <- c("GBR", "IBS", "CEU", "TSI", "FIN",               # EUR
              "MXL", "PUR", "CLM", "PEL",                      # Admix American
              "YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB", # African
              "GIH", "PJL", "BEB", "STU", "ITU",               # SEA
              "CHB", "JPT", "CHS", "CDX", "KHV",               # EA
              "Data"
)

poplevels2 <- c("EUR", "AMR", "AFR", "SEA", "EA", "Data")

pcainput$Pop <- factor(pcainput$Pop, levels = poplevels)
pcainput$suppop <- factor(pcainput$suppop, levels = poplevels2)

################## 2. Plotting
# 2-1: general plot of PC1 vs PC2
# p <- ggscatter(pcainput, x = "V3", y = "V4",
#           color = "suppop",
#           ellipse = TRUE,
#           xlab = "PC1",
#           ylab = "PC2") 
# 
# ggexport(p, filename = "./../figures/ancestry/PC1vsPC2.2.pdf",  width = 15, height = 10)

p <- ggplot(pcainput, aes(x=V3, y=V4, color=suppop)) +
  geom_point(size=.5) +
  stat_ellipse(geom="polygon", alpha=.1, aes(fill=suppop)) +
  labs(title="Genotype PCA\n(African Data, N=152)",
       x="PC1", y="PC2") +
  theme_light() +
  theme(axis.text = element_text(size=10), 
        axis.title = element_text(size=12),
        plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
  scale_fill_discrete("Super Population") + 
  scale_colour_discrete("Super Population") 

p

ggsave("./../figures/ancestry/afr-PC1vsPC2.png", p, width=8, height=6)


# 2-2: PC1vsPC2
p2 <- ggplot(pcainput, aes(x = V3, y = V4, color = suppop, shape = Pop)) + 
  geom_point(size=.7) +
  scale_shape_manual(values = c(0:18,33,20:25,19)) +
  labs(x = "PC1", y = "PC2", title = "PCA PC1 VS PC2") +
  theme_light() +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        plot.title = element_text(size=16, face="bold", hjust = 0.5)) 
p2

ggsave("./../figures/ancestry/PCA1v2.png", p2, width = 9, height = 6)


# 2-3: plot PC1 only
p3 <- ggplot(pcainput, aes(x = Pop, y = V3, shape = Pop, color = suppop)) + 
  geom_point(size=.5) +
  scale_shape_manual(values = c(0:18,33,20:25,19)) +
  labs(y = "PC1", title = "PCA PC1") +
  theme_light() +
  theme(axis.text=element_text(size=10), 
        axis.text.x=element_text(angle = 45, size=8, margin = margin(0.4, unit = "cm")),
        axis.title=element_text(size=12), 
        plot.title = element_text(size=16, face="bold", hjust = 0.5))

p3

ggsave("./../figures/ancestry/PC1.png", p3, width = 9, height = 6)

# 2-4: plot PC2 only
p4 <- ggplot(pcainput, aes(x = Pop, y = V4, shape = Pop, color = suppop)) + 
  geom_point(size=.5) +
  scale_shape_manual(values = c(0:18,33,20:25,19)) +
  labs(y = "PC2", title = "PCA PC2") +
  theme_light() +
  theme(axis.text=element_text(size=10), 
        axis.text.x=element_text(angle = 45, size=8, margin = margin(0.4, unit = "cm")),
        axis.title=element_text(size=12), 
        plot.title = element_text(size=16, face="bold", hjust = 0.5))

ggsave("./../figures/ancestry/PC2.png", p4, width = 9, height = 6)


# 2-5: PC3, PC4, PC5
p5 <- ggplot(pcainput, aes(x = Pop, y = V5, shape = Pop, color = suppop)) + 
  geom_point(size=.5) +
  scale_shape_manual(values = c(0:18,33,20:25,19)) +
  labs(y = "PC3", title = "PCA PC3") +
  theme_light() +
  theme(axis.text=element_text(size=10), 
        axis.text.x=element_text(angle = 45, size=8, margin = margin(0.4, unit = "cm")),
        axis.title=element_text(size=12), 
        plot.title = element_text(size=16, face="bold", hjust = 0.5))

ggsave("./../figures/ancestry/PC3.png", p5, width = 9, height = 6)

p6 <- ggplot(pcainput, aes(x = Pop, y = V6, shape = Pop, color = suppop)) + 
  geom_point(size=.5) +
  scale_shape_manual(values = c(0:18,33,20:25,19)) +
  labs(y = "PC4", title = "PCA PC4") +
  theme_light() +
  theme(axis.text=element_text(size=10), 
        axis.text.x=element_text(angle = 45, size=8, margin = margin(0.4, unit = "cm")),
        axis.title=element_text(size=12), 
        plot.title = element_text(size=16, face="bold", hjust = 0.5))

ggsave("./../figures/ancestry/PC4.png", p6, width = 9, height = 6)

p7 <- ggplot(pcainput, aes(x = Pop, y = V7, shape = Pop, color = suppop)) + 
  geom_point(size=.5) +
  scale_shape_manual(values = c(0:18,33,20:25,19)) +
  labs(y = "PC5", title = "PCA PC5") +
  theme_light() +
  theme(axis.text=element_text(size=10), 
        axis.text.x=element_text(angle = 45, size=8, margin = margin(0.4, unit = "cm")),
        axis.title=element_text(size=12), 
        plot.title = element_text(size=16, face="bold", hjust = 0.5))

ggsave("./../figures/ancestry/PC5.png", p7, width = 9, height = 6)




################## 3: Ancestry calling: k nearest neighbors classification
library(class)
pca <- pcainput[,3:24]
rownames(pca) <- pcainput[,2]
pca_scale <- as.data.frame(scale(pca[,1:20], center = TRUE, scale = TRUE))

# 1kg as training set, data as testing set
n_ref <- 2504
n_data<-dim(pca_scale)[1]-n_ref

train_set <- pca_scale[1:n_ref,]
test_set <- pca_scale[(n_ref+1):dim(pca_scale)[1],]
train_label <- pca[1:n_ref,21] #suppop

# initialize K as the square root of number of sample in training set
K <- floor(sqrt(n_ref))
knn <- knn(train=train_set, test=test_set, cl=train_label, k=K)

table(knn)
# knn
# EUR  AMR  AFR  SEA   EA Data 
# 293  156  152   28   25    0 

df_count <- data.frame("Ancestry" = c("EUR","AMR","AFR","SEA","EA"),
                      "Count" = c(sum(knn=="EUR"),
                                  sum(knn=="AMR"),
                                  sum(knn=="AFR"),
                                  sum(knn=="SEA"),
                                  sum(knn=="EA")))
df_count$Ancestry <- factor(df_count$Ancestry, levels = c("EUR","AMR","AFR","SEA","EA"))

p_ancestry <- ggplot(df_count, aes(x=Ancestry, y=Count)) +
  geom_bar(stat = "identity", aes(fill=Ancestry)) +
  labs(x="Ancestry", y="Number of Subjects", title="Data Ancestry") +
  geom_text(aes(label=Count),vjust=-1, size=6) +
  theme_light() +
  theme(axis.text = element_text(size=10),
        plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
  ylim(0,350)

p_ancestry

ggsave("./../figures/ancestry/anecstry.png", p_ancestry, width = 6, height = 6)






################## 4: Write to list
eur<-amr<-afr<-sea<-ea<-c()
pcainput[] <- lapply(pcainput, as.character)
for(i in 1:n_data) {
  if (knn[i]=="EUR") {
    eur <- append(eur,pcainput[n_ref+i,2])
  }
  if (knn[i]=="AFR") {
    afr <- append(afr,pcainput[n_ref+i,2])
  }
  if (knn[i]=="AMR") {
    amr <- append(amr,pcainput[n_ref+i,2])
  }
  if (knn[i]=="SEA") {
    sea <- append(sea,pcainput[n_ref+i,2])
  }
  if (knn[i]=="EA") {
    ea <- append(ea,pcainput[n_ref+i,2])
  }
}
write.table(data.frame(eur), "./ancestry/list.eur.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(data.frame(afr), "./ancestry/list.afr.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(data.frame(amr), "./ancestry/list.amr.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(data.frame(sea), "./ancestry/list.sea.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(data.frame(ea), "./ancestry/list.ea.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')





##################    hierarchical clustering
# pca <- pcainput[,-1]
# rownames(pca) <- pcainput[,1]
# pca_data_list <- pca$Pop%in%c("data")
# pca_data <- pca[pca_data_list,]
# pca_data <- pca_data[,1:20]
# pca_data_scale <- as.data.frame(scale(pca_data, center = TRUE, scale = TRUE))
# 
# dist_mat <- dist(pca_data_scale, method = 'euclidean')
# hclust <- hclust(dist_mat, method = 'complete')
# # dend <- as.dendrogram(hclust)
# # par(cex=0.4)
# # plot(dend, main = "hclust complete")
# 
# groups<-cutree(hclust, k=5)
# x<-cbind(pca_data_scale,groups)
# x1<- subset(x, groups==1)
# x2<- subset(x, groups==2)
# x3<- subset(x, groups==3)
# x4<- subset(x, groups==5)
# x5<- subset(x, groups==5)
### results seem off!

##################    k-means clustering, using top ? PC
### problem: cannot select number of PC to use based on cluster mean PC1 and PC2!!!
# pca <- pcainput[,-1]
# rownames(pca) <- pcainput[,1]
# pca_data_list <- pca$Pop%in%c("data")
# pca_data <- pca[pca_data_list,]
# pca_data <- pca_data[,1:2]
# #not scaling since it's already PC...?
# #pca_data_scale <- as.data.frame(scale(pca_data, center = TRUE, scale = TRUE))
# set.seed(20)
# clusters <- kmeans(pca_data, centers = 5, nstart = 3)
# str(clusters)
# x <- cbind(pca_data,clusters$cluster)
# # $ size        : int [1:5] 56 125 110 25 195
# x1<- subset(x, clusters$cluster==1)
# x2<- subset(x, clusters$cluster==2)
# x3<- subset(x, clusters$cluster==3)
# x4<- subset(x, clusters$cluster==4)
# x5<- subset(x, clusters$cluster==5)
# 
# mean_1 <- c(mean(x1$V3), mean(x1$V4))
# mean_2 <- c(mean(x2$V3), mean(x2$V4))
# mean_3 <- c(mean(x3$V3), mean(x3$V4))
# mean_4 <- c(mean(x4$V3), mean(x4$V4))
# mean_5 <- c(mean(x5$V3), mean(x5$V4))
# means <- list(mean_1,mean_2,mean_3,mean_4,mean_5)
# df <- data.frame(matrix(ncol = 2, nrow = 5))
# colnames(df) <- c("PC1","PC2")
# for(i in 1:5) {
#   df[i,1] <- means[[i]][1]
#   df[i,2] <- means[[i]][2]
# }
# 
# p <- ggscatter(pcainput, x = "V3", y = "V4",
#               color = "suppop",
#               ellipse = TRUE,
#               xlab = "PC1",
#               ylab = "PC2") + geom_point(df, mapping = aes(x=PC1,y=PC2))
# ggexport(p, filename = "k-means_cluster_PC_mean_top2.pdf", width = 15, height = 10, units = "in")
# 
# df_count <- data.frame(matrix(ncol=2, nrow = 5))
# colnames(df_count) <- c("Ancestry","Count")
# df_count[,1] <- c("EUR","AFR","AMR","SEA","EA")
# df_count[,2] <- c(195,110,56,125,25)
# df_count$Ancestry <- factor(df_count$Ancestry, levels = c("EUR","AFR","AMR","SEA","EA"))
# p <- ggplot(df_count, aes(Ancestry, Count)) +
#   geom_bar(stat = "identity", aes(fill=Ancestry)) +
#   xlab("Ancestry") +
#   ylab("Number of samples")
# ggexport(p, filename = "ancestry_top2PC.pdf",width = 10, height = 10, units = "in")


##################    ancestry calling based on top PC threshold
# select EUR samples
# EUR_list <- pcainput$suppop %in% c("EUR") #503
# pc_threshold <- list()
# for(i in 1:10) {
#   pc_threshold[[i]] <- range(pcainput[EUR_list,i+1]) 
# }
# 
# eursamples = filter(pcainput, Pop == "data" & 
#                      V3 >= pc_threshold[[1]][1] & V3 <= pc_threshold[[1]][2] &
#                      V4 >= pc_threshold[[2]][1] & V4 <= pc_threshold[[2]][2] &
#                      V5 >= pc_threshold[[3]][1] & V5 <= pc_threshold[[3]][2] &
#                      V6 >= pc_threshold[[4]][1] & V6 <= pc_threshold[[4]][2] &
#                      V7 >= pc_threshold[[5]][1] & V7 <= pc_threshold[[5]][2] &
#                      V8 >= pc_threshold[[6]][1] & V8 <= pc_threshold[[6]][2] &
#                      V9 >= pc_threshold[[7]][1] & V9 <= pc_threshold[[7]][2] &
#                       V10 >= pc_threshold[[8]][1] & V10 <= pc_threshold[[8]][2] 
#                       #V11 >= pc_threshold[[9]][1] & V11 <= pc_threshold[[9]][2]
#                       #V12 >= pc_threshold[[10]][1] & V12 <= pc_threshold[[10]][2]
#                       )
# nrow(eursamples) #120
# write.table(eursamples[, 1], "data.eur", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
# 
# # select AMR samples
# AMR_list <- pcainput$suppop %in% c("AMR")
# pc1_threshold <- range(pcainput[AMR_list,]$V3) 
# pc2_threshold <- range(pcainput[AMR_list,]$V4) 
# pc3_threshold <- range(pcainput[AMR_list,]$V5)
# pc4_threshold <- range(pcainput[AMR_list,]$V6)
# pc5_threshold <- range(pcainput[AMR_list,]$V7)
# pc6_threshold <- range(pcainput[AMR_list,]$V8)
# pc7_threshold <- range(pcainput[AMR_list,]$V9)
# amrsamples = filter(pcainput, Pop == "data"  & 
#                       V3 >= pc1_threshold[1] & V3 <= pc1_threshold[2] &
#                       V4 >= pc2_threshold[1] & V4 <= pc2_threshold[2] &
#                       V5 >= pc3_threshold[1] & V5 <= pc3_threshold[2] &
#                       V6 >= pc4_threshold[1] & V6 <= pc4_threshold[2] &
#                       V7 >= pc5_threshold[1] & V7 <= pc5_threshold[2] &
#                       V8 >= pc6_threshold[1] & V8 <= pc6_threshold[2] &
#                       V9 >= pc7_threshold[1] & V9 <= pc7_threshold[2])
# nrow(amrsamples) #173
# write.table(amrsamples[, 1], "data.amr", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
# 
# # select AFR samples
# AFR_list <- pcainput$suppop %in% c("AFR")
# pc1_threshold <- range(pcainput[AFR_list,]$V3) 
# pc2_threshold <- range(pcainput[AFR_list,]$V4) 
# pc3_threshold <- range(pcainput[AFR_list,]$V5)
# pc4_threshold <- range(pcainput[AFR_list,]$V6)
# pc5_threshold <- range(pcainput[AFR_list,]$V7)
# pc6_threshold <- range(pcainput[AFR_list,]$V8)
# pc7_threshold <- range(pcainput[AFR_list,]$V9)
# afrsamples = filter(pcainput, Pop == "data"  & 
#                       V3 >= pc1_threshold[1] & V3 <= pc1_threshold[2] &
#                       V4 >= pc2_threshold[1] & V4 <= pc2_threshold[2] &
#                       V5 >= pc3_threshold[1] & V5 <= pc3_threshold[2] &
#                       V6 >= pc4_threshold[1] & V6 <= pc4_threshold[2] &
#                       V7 >= pc5_threshold[1] & V7 <= pc5_threshold[2] &
#                       V8 >= pc6_threshold[1] & V8 <= pc6_threshold[2] &
#                       V9 >= pc7_threshold[1] & V9 <= pc7_threshold[2])
# nrow(afrsamples) #240
# write.table(afrsamples[, 1], "data.afr", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
# 
# # select SEA samples
# SEA_list <- pcainput$suppop %in% c("SEA")
# pc1_threshold <- range(pcainput[SEA_list,]$V3) 
# pc2_threshold <- range(pcainput[SEA_list,]$V4) 
# pc3_threshold <- range(pcainput[SEA_list,]$V5)
# pc4_threshold <- range(pcainput[SEA_list,]$V6)
# pc5_threshold <- range(pcainput[SEA_list,]$V7)
# pc6_threshold <- range(pcainput[SEA_list,]$V8)
# pc7_threshold <- range(pcainput[SEA_list,]$V9)
# seasamples = filter(pcainput, Pop == "data"  & 
#                       V3 >= pc1_threshold[1] & V3 <= pc1_threshold[2] &
#                       V4 >= pc2_threshold[1] & V4 <= pc2_threshold[2] &
#                       V5 >= pc3_threshold[1] & V5 <= pc3_threshold[2] &
#                       V6 >= pc4_threshold[1] & V6 <= pc4_threshold[2] &
#                       V7 >= pc5_threshold[1] & V7 <= pc5_threshold[2] &
#                       V8 >= pc6_threshold[1] & V8 <= pc6_threshold[2] &
#                       V9 >= pc7_threshold[1] & V9 <= pc7_threshold[2])
# nrow(seasamples) #2
# write.table(seasamples[, 1], "data.sea", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
# 
# # select EA samples
# EA_list <- pcainput$suppop %in% c("EA")
# pc1_threshold <- range(pcainput[EA_list,]$V3) 
# pc2_threshold <- range(pcainput[EA_list,]$V4) 
# pc3_threshold <- range(pcainput[EA_list,]$V5)
# pc4_threshold <- range(pcainput[EA_list,]$V6)
# pc5_threshold <- range(pcainput[EA_list,]$V7)
# pc6_threshold <- range(pcainput[EA_list,]$V8)
# pc7_threshold <- range(pcainput[EA_list,]$V9)
# easamples = filter(pcainput, Pop == "data" & 
#                      V3 >= pc1_threshold[1] & V3 <= pc1_threshold[2] &
#                      V4 >= pc2_threshold[1] & V4 <= pc2_threshold[2] &
#                      V5 >= pc3_threshold[1] & V5 <= pc3_threshold[2] &
#                      V6 >= pc4_threshold[1] & V6 <= pc4_threshold[2] &
#                      V7 >= pc5_threshold[1] & V7 <= pc5_threshold[2] &
#                      V8 >= pc6_threshold[1] & V8 <= pc6_threshold[2] &
#                      V9 >= pc7_threshold[1] & V9 <= pc7_threshold[2])
# nrow(easamples) #1
# write.table(easamples[, 1], "data.ea", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
# 
# 
# df <- data.frame(
#   pop=c("EUR","AMR","AFR","SEA","EA"),
#   count=c(120,173,240,2,1)
# )
# p <- ggplot(data=df, aes(x=pop, y=count)) +
#   geom_bar(stat = "identity") +
#   scale_x_discrete(limits=c("EUR", "AMR", "AFR", "SEA", "EA"))
# p
# ggsave("ancestry.pdf", width = 8, height = 8, units = "in")

