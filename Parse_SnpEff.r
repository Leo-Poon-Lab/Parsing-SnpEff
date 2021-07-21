#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.csv"
}
library(tidyverse)

# input_vcf <- "../data/delta_hk.snpeff.vcf"
input_vcf <- args[1]

df_vcf_ann <- read_tsv(input_vcf, comment = "#", col_names = F)
df_vcf_ann$mut <- sapply(df_vcf_ann$X8, function(x){
	strsplit(x, ";")[[1]][1]
}, USE.NAMES = F)
df_vcf_ann$effect <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2] # only consider the nearest match
	strsplit(tmp, "|", fixed = T)[[1]][2]
}, USE.NAMES = F)
df_vcf_ann$impact <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][3]
}, USE.NAMES = F)
df_vcf_ann$gene <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][4]
}, USE.NAMES = F)
df_vcf_ann$gene_name <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][5]
}, USE.NAMES = F)
df_vcf_ann$mut_cdna <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][10]
}, USE.NAMES = F)
df_vcf_ann$mut_aa <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][11]
}, USE.NAMES = F)
df_vcf_ann$pos_cdna <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][12]
}, USE.NAMES = F)
df_vcf_ann$pos_aa <- sapply(df_vcf_ann$X8, function(x){
	tmp <- strsplit(x, ";")[[1]][2]
	strsplit(tmp, "|", fixed = T)[[1]][14]
}, USE.NAMES = F)

df_vcf_ann$in_gene <- df_ann$effect %in% c("missense_variant", "synonymous_variant")

df_vcf_ann %>% select(-X3, -X6, -X7, -X8) %>% write_csv(args[2])
