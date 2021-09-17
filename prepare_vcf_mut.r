prepare_vcf <- function(mut_nt, outfile, ref_name="MN908947.3") { 
	df_vcf <- tibble(Gene = rep(ref_name, length(mut_nt)))
	df_vcf$Pos <- sapply(mut_nt, function(x){
		gsub("\\D", "", x)
	})
	df_vcf$QUAL <- "."
	df_vcf$ori <- sapply(mut_nt, function(x){
		strsplit(gsub("\\d", "", x), "")[[1]][1]
	})
	df_vcf$alle <- sapply(mut_nt, function(x){
		strsplit(gsub("\\d", "", x), "")[[1]][2]
	})
	df_vcf$dep <- 100
	df_vcf$pass <- "PASS"
	df_vcf$info <- paste0("MUT=", mut_nt)

	df_vcf <- df_vcf %>% arrange(as.numeric(Pos))
	df_vcf <- df_vcf %>% filter((!is.na(Pos)) & Pos != "")
	# df_vcf <- df_vcf %>% filter(alle != "-")
	df_vcf$alle[df_vcf$alle=="-"] <- "<DEL>"
	outfile <- outfile
	write_tsv(df_vcf, outfile, col_names = F)
	header <- c('##fileformat=VCFv4.0', '##reference=reference.fasta', '##INFO=<ID=MUT,Number=1,Type=String,Description="Mutation">', '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO')
	writeLines(c(header, readLines(outfile)), outfile)
}