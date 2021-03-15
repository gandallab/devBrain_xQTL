#! /usr/bin/env Rscript

library(argparser)
library(SeqVarTools)
library(Biobase)
library(data.table)
library(dplyr)

argp <- arg_parser("Generate AnnotatedDataFrame format of phenotype")
argp <- add_argument(argp, "--ancestry", help="Inferred ancestry")
argp <- add_argument(argp, "--out_file", help="Phenotype annotated RData")
argv <- parse_args(argp)

# Load data
pop <- fread(argv$ancestry, data.table=F )

# Make data for AnnotatedDataFrame
colnames(pop)<-c("sample.id","ethnicity") # note this has to be "sample.id", used in TopmedPipeline built-in functions
metaData <- data.frame(labelDescription=c("sample.id", "ethnicity"))
annot=AnnotatedDataFrame(data=pop, varMetadata=metaData)
save(annot,file=argv$out_file)
