library(readxl)
library(dplyr)
library(ggplot2)

files = list.files(pattern="*-simple.xlsx")
files

#####----FILE WITH THE ALL VARIANTS FROM 29 FAMILIES

read_files <- function(x){
  fam_name = strsplit(strsplit(x, "-")[[1]][1], " ")[[1]][1]
  dat <- read_excel(x, col_names = F, skip = 1)[, 1:20 ]
  colnames(dat) <- c("pos","ref/alt","gene","seq_ontology","effect","Transcript", "HGVS_c","HGVS_p", "AF", "Alt_allele_count", 
                     "Homozygote_count", "Hemizygote_count", "Prediction", "PHRED", "HGMD", 
                     "Disorder", "Inheritance", "Clinical_significance", "Aggregate", "Pho")
  dat$fam <- fam_name
  return(as.data.frame(dat))
}

data <- lapply(files, FUN = read_files)
data <- do.call(rbind, data)
data

# Count number of families with a gene
genes <- data %>% 
  group_by(fam, gene) %>%
  summarise(HGVS_c = paste0(HGVS_c, collapse = ","), 
            HGVS_p = paste0(HGVS_p, collapse = ","), 
            effect = paste0(effect, collapse = ","))

res1 = genes %>% 
  group_by(gene) %>%
  summarise(families = paste0(fam, collapse = ","), 
            count = n(),
            HGVS_c = paste0(HGVS_c, collapse = ", "),
            HGVS_p = paste0(HGVS_p, collapse = ", "),
            effect = paste0(effect, collapse = ", ")) %>%
  arrange(desc(count))
View(res1)

# Count variants
variants <- data %>% 
  group_by(fam, gene) %>%
  count(HGVS_c, HGVS_p, seq_ontology, effect)

res2 = variants %>% 
  group_by(gene, HGVS_c, HGVS_p, seq_ontology, effect) %>%
  summarise(families = paste0(fam, collapse = ", "), count = n()) %>%
  arrange(desc(count))
#arrange(gene,desc(count))
View(res2)

