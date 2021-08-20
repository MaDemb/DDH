DESCRIPTION OF SCRIPTS AT https://github.com/MaDemb/DDH
The two scripts contain R code used to analyze the variants found in an analysis of 29 families with hip dysplasia.

The file "Count_genes_LvBA_Rscripts.R" contains the code created by Lars van Brakel Andersen and it is used to process the list of variants obtained from VarSeq software for each family.
The code provide a list of genes and a list of variants occuring in all families.

The file "Filter_and_correlation_MD_Rscript.R" contains the code created by Maja Dembic and it is to filter the variants and to perform correlation analysis between gene legth and family count, and gene length and variant count per each gene. 
Maximum length (in bp) of each gene was obtained using BioMart - Ensembl to export data. Family count and variant count was added through our scripts.