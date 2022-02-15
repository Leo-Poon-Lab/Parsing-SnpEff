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

file.copy(input_vcf, "./tmp.vcf", overwrite = T)
# the following works for MacOS
system("sed -i ' ' 's/Ala/A/'g tmp.vcf
sed -i ' ' 's/Arg/R/'g tmp.vcf
sed -i ' ' 's/Asn/N/'g tmp.vcf
sed -i ' ' 's/Asp/D/'g tmp.vcf
sed -i ' ' 's/Cys/C/'g tmp.vcf
sed -i ' ' 's/Glu/E/'g tmp.vcf
sed -i ' ' 's/Gln/Q/'g tmp.vcf
sed -i ' ' 's/Gly/G/'g tmp.vcf
sed -i ' ' 's/His/H/'g tmp.vcf
sed -i ' ' 's/His/H/'g tmp.vcf
sed -i ' ' 's/Ile/I/'g tmp.vcf
sed -i ' ' 's/Leu/L/'g tmp.vcf
sed -i ' ' 's/Lys/K/'g tmp.vcf
sed -i ' ' 's/Met/M/'g tmp.vcf
sed -i ' ' 's/Phe/F/'g tmp.vcf
sed -i ' ' 's/Pro/P/'g tmp.vcf
sed -i ' ' 's/Ser/S/'g tmp.vcf
sed -i ' ' 's/Thr/T/'g tmp.vcf
sed -i ' ' 's/Trp/W/'g tmp.vcf
sed -i ' ' 's/Tyr/Y/'g tmp.vcf
sed -i ' ' 's/Val/V/'g tmp.vcf")

df_vcf_ann <- read_tsv("./tmp.vcf", comment = "#", col_names = F)
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

df_vcf_ann$in_gene <- df_vcf_ann$effect %in% c("missense_variant", "synonymous_variant")

df_vcf_ann %>% select(-X3, -X6, -X7, -X8) %>% write_csv(args[2])
file.remove("./tmp.vcf")
file.remove("./tmp.vcf ")
