library("ChemmineR") 
library(plotly)
setwd('/Users/bia/Desktop/bioinformatics2/drug-protein-interaction/') 

# load data
dti_original <- read.table("data_project/DTI_100molecules.txt", header = FALSE)
dti[is.na(dti_original)] <- 0

# pki > 7 = interaction
pki_interest = 7
interactions_count <- rowSums(dti > pki_interest, na.rm = TRUE)
hist(interactions_count, main="Interactions distribution")

#p <- plot_ly(z = volcano, type = "contour")

