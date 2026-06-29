#### Add side_chain and location data #####
#first I load the rsa data, predicted by pop music
#add amino acid orientation category
sasa <- read.csv("SOD1_rSASA.csv", sep = ";")
sasa$rSASA.monomer <- str_replace_all(sasa$rSASA.monomer, "[^0-9.,-]", "")  
sasa$rSASA.dimer <- str_replace_all(sasa$rSASA.dimer, "[^0-9.,-]", "")  
sasa$Pos <- str_replace_all(sasa$Pos, "[^0-9.,-]", "")  

sasa$rSASA.monomer <- as.numeric(sasa$rSASA.monomer)
sasa$rSASA.dimer <- as.numeric(sasa$rSASA.dimer)
sasa$Pos <- as.numeric(sasa$Pos)


sasa$drSASA <- sasa$rSASA.monomer - sasa$rSASA.dimer

# predictions <- read.csv("predictions_SOD1.csv", sep = ";")
##2. Extract positions of each mutation (for all subs, ins and dels)
sasa <- sasa %>%
  mutate(side_chain_monomer = case_when(
    Pos %in% c("5", "7", "50", "51", "52", "53", "54", "113", "114", "148", "150", "151", "152", "153") ~ "dimer interface",
    rSASA.monomer <=25 ~ "core",
    rSASA.monomer >25 ~ "surface")) %>% 
  mutate(side_chain_dimer = case_when(
    Pos %in% c("5", "7", "50", "51", "52", "53", "54", "113", "114", "148", "150", "151", "152", "153") ~ "dimer interface",
    rSASA.dimer <=25 ~ "core",
    rSASA.dimer >25 ~ "surface")) %>% 
  mutate(side_chain = case_when(
    Pos %in% c("5", "7", "50", "51", "52", "53", "54", "113", "114", "148", "150", "151", "152", "153") ~ "dimer interface",
    rSASA.dimer <=25 ~ "core",
    rSASA.dimer >25 ~ "surface"))

sasa <- sasa %>%
  mutate(side_chain = case_when(
    Pos %in% c(63, 71, 80, 83) ~ "Zn binding",
    TRUE ~ side_chain)) %>%
  mutate(side_chain = recode(side_chain,
                             "binding interface" = "dimer interface")) %>%
  mutate(location = case_when(
    Pos %in% c(3:9, 15:22, 29:36, 41:48, 86:89, 95:101, 116:120, 143:148) ~ "β-sheet",
    Pos %in% c(1,2, 10:14, 23:28, 37:40, 49:85, 90:94, 102:115, 121:142, 149:153) ~ "loop")) %>%
  mutate(s_location = case_when(
    Pos %in% c(63, 71, 80, 83) ~ "Zn binding residues",
    Pos %in% c(1,2) ~ "N term",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(3:9) ~ "β1",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(15:22) ~ "β2",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(29:36) ~ "β3",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(41:48) ~ "β4",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(83:89) ~ "β5",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(95:101) ~ "β6",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(116:120) ~ "β7",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(143:148) ~ "β8",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(10:14) ~ "loop 1",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(23:28) ~ "loop 2",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(37:40) ~ "loop 3",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(49:82) ~ "Zn binding loop",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(90:94) ~ "loop 5",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(102:115) ~ "loop 6",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(121:142) ~ "Electrostatic loop",
    !(Pos %in% c(71, 80, 83)) & Pos %in% c(149:153) ~ "loop 8"
  ))


sasa <- sasa %>%
  distinct(Pos, side_chain, location, s_location, drSASA, rSASA.monomer, rSASA.dimer, .keep_all = FALSE)

side_chain_vec <- setNames(sasa$side_chain, sasa$Pos)
location_vec <- setNames(sasa$location, sasa$Pos)
s_location_vec <- setNames(sasa$s_location, sasa$Pos)
drSASA_vec <- setNames(sasa$drSASA, sasa$Pos)
rSASA.monomer_vec <- setNames(sasa$rSASA.monomer, sasa$Pos)
rSASA.dimer_vec <- setNames(sasa$rSASA.dimer, sasa$Pos)



#load all SOD1 normalized data
SOD1 <- read.csv("SOD1_abundance_binding_scores.csv")

SOD1_subsitutions <- SOD1 %>%
  filter(mutation_type == "subs") %>%
  mutate(WT_AA = str_extract(ID, "^[A-za-z]"),
         Pos = str_extract(ID, "\\d+"),
         Mut_AA = str_extract(ID, "[A-za-z]$")) 

WT_AA_vec <- setNames(SOD1_subsitutions$WT_AA, SOD1_subsitutions$Pos)


SOD1_insertions <- SOD1 %>%
  filter(mutation_type == "ins") %>%
  mutate(Pos = str_extract(ID, "\\d+"),
         WT_AA = WT_AA_vec[as.character(Pos)],
         Mut_AA = str_extract(ID, "[A-za-z]$"))



SOD1_deletions <- SOD1 %>%
  filter(mutation_type == "del") %>%
  mutate(Pos = str_extract(ID, "\\d+"),
         WT_AA = WT_AA_vec[as.character(Pos)],
         Mut_AA = str_extract(ID, "^[A-za-z]"))


SOD1_subsitutions$side_chain <- side_chain_vec[as.character(SOD1_subsitutions$Pos)]
SOD1_subsitutions$location <- location_vec[as.character(SOD1_subsitutions$Pos)]
SOD1_subsitutions$s_location <- s_location_vec[as.character(SOD1_subsitutions$Pos)]
SOD1_subsitutions$drSASA <- drSASA_vec[as.character(SOD1_subsitutions$Pos)]
SOD1_subsitutions$rSASA.monomer <- rSASA.monomer_vec[as.character(SOD1_subsitutions$Pos)]
SOD1_subsitutions$rSASA.dimer <- rSASA.dimer_vec[as.character(SOD1_subsitutions$Pos)]

SOD1_subsitutions <- SOD1_subsitutions %>%
  mutate(interface = case_when(
    side_chain == "dimer interface" ~ TRUE,
    TRUE ~ FALSE
  ))



SOD1_insertions$side_chain <- side_chain_vec[as.character(SOD1_insertions$Pos)]
SOD1_insertions$location <- location_vec[as.character(SOD1_insertions$Pos)]
SOD1_insertions$s_location <- s_location_vec[as.character(SOD1_insertions$Pos)]
SOD1_insertions$drSASA <- drSASA_vec[as.character(SOD1_insertions$Pos)]
SOD1_insertions$rSASA.monomer <- rSASA.monomer_vec[as.character(SOD1_insertions$Pos)]
SOD1_insertions$rSASA.dimer <- rSASA.dimer_vec[as.character(SOD1_insertions$Pos)]

SOD1_insertions <- SOD1_insertions %>%
  mutate(interface = case_when(
    side_chain == "dimer interface" ~ TRUE,
    TRUE ~ FALSE
  ))



SOD1_deletions$side_chain <- side_chain_vec[as.character(SOD1_deletions$Pos)]
SOD1_deletions$location <- location_vec[as.character(SOD1_deletions$Pos)]
SOD1_deletions$s_location <- s_location_vec[as.character(SOD1_deletions$Pos)]
SOD1_deletions$drSASA <- drSASA_vec[as.character(SOD1_deletions$Pos)]
SOD1_deletions$rSASA.monomer <- rSASA.monomer_vec[as.character(SOD1_deletions$Pos)]
SOD1_deletions$rSASA.dimer <- rSASA.dimer_vec[as.character(SOD1_deletions$Pos)]

SOD1_deletions <- SOD1_deletions %>%
  mutate(interface = case_when(
    side_chain == "dimer interface" ~ TRUE,
    TRUE ~ FALSE
  ))

SOD1_allmut <- bind_rows(SOD1_subsitutions, SOD1_insertions, SOD1_deletions)
write.csv(SOD1_allmut, file = "SOD1_allmut_no_modelled.csv")