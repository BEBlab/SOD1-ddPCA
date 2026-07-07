#This script allows to quantify the minimal distance between residues in chain A and chain B. The SOD1 PDB strcuture I use is 2V0A. To measure distances I use the alpha carbon of each residue

library(bio3d)
library(dplyr)
library(tidyr)

pdb_id_or_path <- "2V0A"  
chain1 <- "A"
chain2 <- "F"
output_prefix <- "SOD1_CA_dist"
contact_cutoff <- 5 #I define that atoms are in contact if their distance is < 5 Angstroms         
sod1_res_range <- 1:153     

pdb <- read.pdb(pdb_id_or_path)
atoms <- pdb$atom

ca_chain1 <- atoms %>% filter(chain == chain1 & resno %in% sod1_res_range)
ca_chain2 <- atoms %>% filter(chain == chain2 & resno %in% sod1_res_range)

if(nrow(ca_chain1) == 0 | nrow(ca_chain2) == 0){
  stop("invalid CA")
}

#Calculate distances between alpha carbons
res_pairs <- expand.grid(i = 1:nrow(ca_chain1), j = 1:nrow(ca_chain2))
res_pair_dist <- res_pairs %>%
  rowwise() %>%
  mutate(
    resno_A = ca_chain1$resno[i],
    resname_A = ca_chain1$resname[i],
    resno_B = ca_chain2$resno[j],
    resname_B = ca_chain2$resname[j],
    dist_CA = sqrt(sum((as.numeric(ca_chain1[i, c("x","y","z")]) -
                          as.numeric(ca_chain2[j, c("x","y","z")]))^2)),
    contact = ifelse(dist_CA <= contact_cutoff, "Yes", "No")
  ) %>%
  ungroup()



min_dist_A <- res_pair_dist %>%
  group_by(resno_A) %>%
  summarize(min_dist_to_B = min(dist_CA, na.rm = TRUE), .groups = "drop")

#Get minimal distance of each residue in chain B with respect residues in chain A
min_dist_B <- res_pair_dist %>%
  group_by(resno_B) %>%
  summarize(min_dist_to_A = min(dist_CA, na.rm = TRUE), .groups = "drop")

write.csv(min_dist_A, file = "min_atom_distance_SOD1.csv", row.names = FALSE)

write.csv(res_pair_dist, paste0(output_prefix, "_CA_dist_filtered.csv"), row.names = FALSE)

contacts_df <- res_pair_dist %>% filter(contact == "Yes")
contacts_df



##### Contacts within one chain ####
ca_chain <- atoms %>%
  filter(chain == chain & resno %in% sod1_res_range)

res_pairs <- expand.grid(i = 1:nrow(ca_chain), j = 1:nrow(ca_chain)) %>%
  filter(i < j)  

res_pair_dist <- res_pairs %>%
  rowwise() %>%
  mutate(
    resno_1 = ca_chain$resno[i],
    resname_1 = ca_chain$resname[i],
    resno_2 = ca_chain$resno[j],
    resname_2 = ca_chain$resname[j],
    dist_CA = sqrt(sum((as.numeric(ca_chain[i, c("x","y","z")]) -
                          as.numeric(ca_chain[j, c("x","y","z")]))^2)),
    contact = ifelse(dist_CA <= contact_cutoff, "Yes", "No")
  ) %>%
  ungroup() %>%
  select(resno_1, resname_1, resno_2, resname_2, dist_CA, contact)

head(res_pair_dist)
