library(readxl)
library(dplyr)

# List Excel files with all variants from the 29 families
#  - Output from Varseq with relevant columns selected
files = list.files(pattern="*-simple.xlsx")
files

# Function to read Excel files
read_files <- function(x){
  # Extract family name from file name
  fam_name = strsplit(strsplit(x, "-")[[1]][1], " ")[[1]][1]
  # Read Excel file, skipping the first line, and selecting the first 20 columns
  dat <- read_excel(x, col_names = F, skip = 1)[, 1:20 ]
  # naming the columns
  colnames(dat) <- c("pos","ref/alt","gene","seq_ontology","effect",
                     "Transcript", "HGVS_c","HGVS_p", "AF", "Alt_allele_count",
                     "Homozygote_count", "Hemizygote_count", "Prediction",
                     "PHRED", "HGMD", "Disorder", "Inheritance",
                     "Clinical_significance", "Aggregate", "Pho")
  # Adding the family name
  dat$fam <- fam_name
  # Return data frame for one family
  return(as.data.frame(dat))
}

# Call function to read all data frames into list
data <- lapply(files, FUN = read_files)
# Bind the data frames together to form one big data frame
data <- do.call(rbind, data)
# Show data
data

# Count number of families with a variant in a specific gene
# - First: summarise all variants within a family i.e. keep only one row for
#   each gene with variants for each family
genes <- data %>%
  group_by(fam, gene) %>%
  summarise(HGVS_c = paste0(HGVS_c, collapse = ","),
            HGVS_p = paste0(HGVS_p, collapse = ","),
            effect = paste0(effect, collapse = ","))

# - Second: Count for each gene, the number of families with at least one
#   variant in that gene, and sort the counts in descending order.
res1 = genes %>%
  group_by(gene) %>%
  summarise(families = paste0(fam, collapse = ","),
            count = n(),
            HGVS_c = paste0(HGVS_c, collapse = ", "),
            HGVS_p = paste0(HGVS_p, collapse = ", "),
            effect = paste0(effect, collapse = ", ")) %>%
  arrange(desc(count))

# View the resulting data frame
View(res1)

# Count the number of times a unique variant is seen in the families
# - First: summarise all unique variants within a family i.e. keep only one row
#   for each unique variant for each family
variants <- data %>%
  group_by(fam, gene) %>%
  count(HGVS_c, HGVS_p, seq_ontology, effect)

# - Second: Count for each unique variant, the number of families with at least
#   one variant in that gene, and sort the counts in descending order.
res2 = variants %>%
  group_by(gene, HGVS_c, HGVS_p, seq_ontology, effect) %>%
  summarise(families = paste0(fam, collapse = ", "), count = n()) %>%
  arrange(desc(count))

# View the resulting data frame
View(res2)
