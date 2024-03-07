


ancestry = read.table('~/Downloads/somalier-ancestry.somalier-ancestry.txt')


download.file('https://ars.els-cdn.com/content/image/1-s2.0-S2211124720303673-mmc2.xlsx', destfile = '~/Downloads/werling_meta.xlsx')
werling = readxl::read_xlsx('~/Downloads/werling_meta.xlsx',sheet = 2)

werling.adult= werling$Braincode[!werling$AgeUnits=="PCW"]

to_keep = to_keep[!to_keep %in% werling.adult]
to_keep = to_keep[!grepl("_R0", to_keep)]
ancestry.fetal = ancestry[ancestry$individualID %in% to_keep,]

df=melt(ancestry.fetal[,c(14,grep("_prob", colnames(ancestry.fetal)))])

df$individualID = factor(df$individualID, levels=ancestry.fetal$individualID[order(ancestry.fetal$EUR_prob,ancestry.fetal$AMR_prob, ancestry.fetal$AFR_prob)])
ggplot(df,aes(x=individualID, y=value,fill=variable)) + geom_bar(stat='identity') + theme_bw() + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),            panel.grid.major.x = element_blank())