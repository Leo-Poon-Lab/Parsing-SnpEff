# Parsing-SnpEff

We use [SnpEff](https://pcingola.github.io/SnpEff/se_inputoutput/#ann-field-vcf-output-files) to annotate vcf files, then parse the results to regular tables with R scripts.

## Install and download database
Download SnpEff via this [link](https://pcingola.github.io/SnpEff/download/).
Download SARS-CoV-2 database via: 
``` sh
java -jar ~/softwares/snpEff/snpEff.jar download -c ~/softwares/snpEff/snpEff.config -v MN908947.3`
```

## Annotation
Run SnpEff as:
``` sh
# run snpeff annotation on vcf
java -jar ~/softwares/snpEff/snpEff.jar MN908947.3 -fastaProt ../results/example.snpeff.vcf.faa -csvStats ../results/example.snpeff.vcf.stats ../results/example.vcf > ../results/example.snpeff.vcf
# run bcftools csq to link consecutive SNPs on the same codon (BCSQ field)
bcftools csq --force --phase a -f reference.fasta -g ../data/genes.gff -Ov variants.snpeff.vcf -o variants.snpeff.csq.vcf

```

## Parsing the results
Run the R scripts:
``` sh
Rscript ./Parse_SnpEff.r ../results/example.snpeff.vcf ../results/example.snpeff.csv
```
