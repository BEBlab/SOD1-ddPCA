library(openxlsx)

path_suppl <- "C:/Users/tquiroga/OneDrive - IBEC/Projects/From owncloud/SOD1 binding PCA/SOD1 manuscript/Supplementary datasets"

#load data from the 3 bPCA SOD1 libraries
l1_apca <- read.csv("aPCA_L1_table.csv")
l2_apca <- read.csv("aPCA_L2_table.csv")
l3_apca <- read.csv("aPCA_L3_table.csv")

all_apca <- rbind(l1_apca, l2_apca, l3_apca)
all_apca$low_sigma <- as.character(all_apca$low_sigma)

write.xlsx(all_apca, file = file.path(path_suppl, "Dataset S2_apca.xlsx"), rowNames = FALSE)
