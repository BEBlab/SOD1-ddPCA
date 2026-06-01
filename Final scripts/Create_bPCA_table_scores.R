library(openxlsx)

path_suppl <- "C:/Users/tquiroga/OneDrive - IBEC/Projects/From owncloud/SOD1 binding PCA/SOD1 manuscript/Supplementary datasets"

#load data from the 3 bPCA SOD1 libraries
l1_bpca <- read.csv("bPCA_L1_table.csv")
l2_bpca <- read.csv("bPCA_L2_table.csv")
l3_bpca <- read.csv("bPCA_L3_table.csv")

all_bpca <- rbind(l1_bpca, l2_bpca, l3_bpca)
all_bpca$low_sigma <- as.character(all_bpca$low_sigma)

write.xlsx(all_bpca, file = file.path(path_suppl, "Dataset S2.xlsx"), rowNames = FALSE)
