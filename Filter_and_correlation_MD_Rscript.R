library(dplyr)

############FILTER GENES
#----------------FILTER THE WORST VARIANTS
#---------------PRINT OUT THE TABLE VARIANTS WITH CADD>25, LoF, or 5/6 damaging 
filter <- (data$Alt_allele_count < 590)
RareVariants <- data[filter,]

filter <- ((RareVariants$PHRED>25 | RareVariants$effect == "LoF" | RareVariants$Prediction == "5 of 6 Predicted as Damaging" |
              RareVariants$Prediction == "6 of 6 Predicted as Damaging"))
Bad_variants <- RareVariants[filter,]
Bad_variants <- Bad_variants[complete.cases(Bad_variants$pos),]

View(Bad_variants)

#Filter further to Pho
filter <- (Bad_variants$Pho > 0.75)
Bad_variants2 <- Bad_variants[filter,]

# Count number of families with a damaging variants
genes <- Bad_variants %>% 
  group_by(fam, gene) %>%
  summarise(HGVS_c = paste0(HGVS_c, collapse = ","), 
            HGVS_p = paste0(HGVS_p, collapse = ","), 
            effect = paste0(effect, collapse = ","),
            PHRED = paste0(PHRED, collapse = ", "),
            Prediction = paste0(Prediction, collapse = ", "),
            Pho = paste0(Pho, collapse = ", "))

resbad = genes %>% 
  group_by(gene) %>%
  summarise(families = paste0(fam, collapse = ","), 
            count = n(),
            HGVS_c = paste0(HGVS_c, collapse = ", "),
            HGVS_p = paste0(HGVS_p, collapse = ", "),
            effect = paste0(effect, collapse = ", "),
            PHRED = paste0(PHRED, collapse = ", "),
            Prediction = paste0(Prediction, collapse = ", "),
            Pho = paste0(Pho, collapse = ", ")) %>%
  arrange(desc(count))
View(resbad)


# Count number of families with damaging variants and Pho>0.75
genes <- Bad_variants2 %>% 
  group_by(fam, gene) %>%
  summarise(HGVS_c = paste0(HGVS_c, collapse = ","), 
            HGVS_p = paste0(HGVS_p, collapse = ","), 
            effect = paste0(effect, collapse = ","),
            PHRED = paste0(PHRED, collapse = ", "),
            Prediction = paste0(Prediction, collapse = ", "),
            Pho = paste0(Pho, collapse = ", "))

resbad2 = genes %>% 
  group_by(gene) %>%
  summarise(families = paste0(fam, collapse = ","), 
            count = n(),
            HGVS_c = paste0(HGVS_c, collapse = ", "),
            HGVS_p = paste0(HGVS_p, collapse = ", "),
            effect = paste0(effect, collapse = ", "),
            PHRED = paste0(PHRED, collapse = ", "),
            Prediction = paste0(Prediction, collapse = ", "),
            Pho = paste0(Pho, collapse = ", ")) %>%
  arrange(desc(count))
View(resbad2)

##########################################
######  CORRELATION ANALYSiS     #########
##########################################


#Overview 
#use GnomAD to get ENST for each gene. Then Use biomart and add length to the list of genes, matching the ENST. 
#This allows for the maximum of genes to be assigned a length).

#Upload files
biomartGenes <- read.delim("mart_export ENST.txt")
colnames(biomartGenes) <- c("length","gene","RefSeq","ENST", "ENSG")

metrics <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene.txt")

ranked_genes_corr <- res1
ranked_genes_corr$families <- NULL
ranked_genes_corr$HGVS_c <- NULL
ranked_genes_corr$HGVS_p <- NULL
ranked_genes_corr$effect <- NULL
ranked_genes_corr <- as.data.frame(ranked_genes_corr)

#Add ENST from the GnomAD data matching the gene name 
for (i in 1:nrow(ranked_genes_corr)){
  gene <- ranked_genes_corr[i,"gene"]
  x <- metrics[metrics$gene == gene,"transcript"]
  if(length(x) == 0) {
    x <- 0
  }
  ranked_genes_corr[i,"ENST"] <- x
}

#Use Biomart data and add the gene length (5' and 3' UTR included), matching the ENST.
for (i in 1:nrow(ranked_genes_corr)){
  ENST <- ranked_genes_corr[i,"ENST"]
  x <- biomartGenes[biomartGenes$ENST == ENST,"length"]
  if (length(x) == 0) {
    x <- 0
  }
  ranked_genes_corr[i,"length"] <- x
}

#Add the number of variants observed per gene
#The final table has a list of genes, family counts, gene length, ENST, and variant count. 
#This part adds the variant count to the final table (ranked_genes_corr)
#Number of variants for correlation

res3 = variants %>% 
  group_by(gene) %>%
  summarise(HGVS_c = paste0(HGVS_c, collapse = ", "), count = n()) %>%
  arrange(desc(count))


allgenes <- res3

for (i in 1:nrow(ranked_genes_corr)){
  gene <- ranked_genes_corr[i,"gene"]
  x <- allgenes[allgenes$gene == gene,"count"]
  if (length(x) == 0) {
    x <- 0
  }
  ranked_genes_corr[i,"Variant_count"] <- x
}

#Eliminate the unmatched
ranked_genes_corr <- ranked_genes_corr[which(ranked_genes_corr$length > 0),]

#Run the correlation
summary(ranked_genes_corr)

cor1 <- cor(x = ranked_genes_corr$count, y = ranked_genes_corr$length)
cor1

cor2 <- cor(x = ranked_genes_corr$length, y = ranked_genes_corr$Variant_count)
cor2
