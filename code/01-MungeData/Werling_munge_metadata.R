#Werling_munge_metadata.R


library(readxl)


download.file('https://ars.els-cdn.com/content/image/1-s2.0-S2211124720303673-mmc2.xlsx',
              destfile = "rawdata/werling_meta_1-s2.0-S2211124720303673-mmc2.xlsx")

meta= read_xlsx("./rawdata/werling_meta_1-s2.0-S2211124720303673-mmc2.xlsx",sheet=2)


meta = meta[meta$AgeUnits=="PCW",]



#Subject metadata
df = data.frame(Subject=meta$Braincode, Study = meta$SampleGroup, Age=-(40-meta$Age)*7/365,
                Sex = meta$SexForAnalysis, Race=meta$AncestryReported, SNPChip=meta$WGSID)
clip = pipe("pbcopy","w"); write.table(df, file=clip, sep = '\t', col.names=F, row.names = FALSE); close(clip)


#Ssample metadata