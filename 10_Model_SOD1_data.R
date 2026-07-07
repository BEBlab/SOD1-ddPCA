#load required packages
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)


#Add structural information#####
sasa <- read.csv("SOD1_rSASA.csv", sep = ";")
sasa$rSASA.monomer <- str_replace_all(sasa$rSASA.monomer, "[^0-9.,-]", "")  
sasa$rSASA.dimer <- str_replace_all(sasa$rSASA.dimer, "[^0-9.,-]", "")  
sasa$Pos <- str_replace_all(sasa$Pos, "[^0-9.,-]", "")  

sasa$rSASA.monomer <- as.numeric(sasa$rSASA.monomer)
sasa$rSASA.dimer <- as.numeric(sasa$rSASA.dimer)
sasa$Pos <- as.numeric(sasa$Pos)


sasa$drSASA <- sasa$rSASA.monomer - sasa$rSASA.dimer

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
    Pos %in% c("1","3","4","6","9","16","18","19","48","49","55","56","57","59","60","61","106","112","114","115","116","146","147","149") ~ "second_shell",
    rSASA.dimer <=25 ~ "core",
    rSASA.dimer >25 ~ "surface"))

sasa <- sasa %>%
  mutate(location = case_when(
    Pos %in% c(3:9, 15:22, 29:36, 41:48, 86:89, 95:101, 116:120, 143:148) ~ "β-sheet",
    Pos %in% c(1,2, 10:14, 23:28, 37:40, 49:85, 90:94, 102:115, 121:142, 149:153) ~ "loop")) %>%
  mutate(s_location = case_when( #specific_location
    Pos %in% c(71, 80, 83) ~ "Zn binding residues",
    Pos %in% c(46, 48, 120) ~ "Cu binding",
    Pos == 63 ~ "Zn and Cu binding",
    Pos %in% c(1,2) ~ "N term",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(3:9) ~ "β1",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(15:22) ~ "β2",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(29:36) ~ "β3",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(41:48) ~ "β4",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(83:89) ~ "β5",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(95:101) ~ "β6",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(116:120) ~ "β7",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(143:148) ~ "β8",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(10:14) ~ "loop 1",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(23:28) ~ "loop 2",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(37:40) ~ "loop 3",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(49:82) ~ "Zn binding loop",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(90:94) ~ "loop 5",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(102:115) ~ "loop 6",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(121:142) ~ "Electrostatic loop",
    !(Pos %in% c(46, 48, 63, 71, 80, 83, 120)) & Pos %in% c(149:153) ~ "loop 8"
  ))


sasa <- sasa %>%
  distinct(Pos, side_chain, location, s_location, drSASA, rSASA.monomer, rSASA.dimer, .keep_all = FALSE)

side_chain_vec <- setNames(sasa$side_chain, sasa$Pos)
location_vec <- setNames(sasa$location, sasa$Pos)
s_location_vec <- setNames(sasa$s_location, sasa$Pos)
drSASA_vec <- setNames(sasa$drSASA, sasa$Pos)
rSASA.monomer_vec <- setNames(sasa$rSASA.monomer, sasa$Pos)
rSASA.dimer_vec <- setNames(sasa$rSASA.dimer, sasa$Pos)



#load SOD1 normalized data
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

#Model data####
#1. check correlation between abundance and binding
r_corr <- SOD1_allmut %>%
  summarise(
    R = cor(abundance_score, binding_score, method = "pearson"),
    p_value = cor.test(abundance_score, binding_score, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p < 0.001 ")
  )

p_corr <- ggplot(data = SOD1_allmut, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 3, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Normalized abundance score", y = "Normalized heterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_corr, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE) 
p_corr

ggsave(p_corr, file = "normal_correlation_apca_boca.tiff", width = 10, height = 8)

#label by stop-like mutations
SOD1_allmut <- SOD1_allmut %>% 
  mutate(stop_like = case_when(
    category_fdr_apca == "stop-like abundance" & category_fdr_bpca == "stop-like heterodimerization" ~ "TRUE",
    TRUE ~ "FALSE"
  ))

p_corr_stop <- ggplot(data = SOD1_allmut, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 3, alpha = 0.2, aes(color = stop_like)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Normalized abundance score", y = "Normalized heterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_corr, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE) 
p_corr_stop

ggsave(p_corr_stop, file = "correlation_apca_boca_label_by_stoplike.tiff", width = 10, height = 8)

#model the data
focal_points <- data.frame(
  abundance_score = c(-1, 0),
  binding_score   = c(-1, 0),
  weights = c(1e04, 1e04)
)

SOD1_allmut <- SOD1_allmut %>%
  mutate(weights = 1) %>%
  bind_rows(focal_points)

model <- loess(binding_score ~ abundance_score,
               data = SOD1_allmut,
               span = 0.5,
               weights = SOD1_allmut$weights)

# loess predictions
pred <- data.frame(
  abundance_score = seq(min(SOD1_allmut$abundance_score),
                        max(SOD1_allmut$abundance_score),
                        length.out = 200)
)
pred$binding_pred <- predict(model, newdata = pred)

SOD1_allmut$fitted <- fitted(model)
SOD1_allmut$residuals <- SOD1_allmut$binding_score - SOD1_allmut$fitted

mad_resid <- mad(SOD1_allmut$residuals, na.rm = TRUE)
SOD1_allmut$outlier <- abs(SOD1_allmut$residuals) > 2 * mad_resid


r_corr <- SOD1_allmut %>%
  summarise(
    R = cor(abundance_score, binding_score, method = "pearson"),
    p_value = cor.test(abundance_score, binding_score, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p < 0.001 ")
  )

p_loess <- ggplot(SOD1_allmut, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = residuals), size = 2.5, shape = 16, alpha = 1) +
  scale_color_gradient2(
    name = "Heterodimerization residual",                     
    low = "#B99C33", mid = "grey95", high = "#AC7299",
    midpoint = 0,
    breaks=c(-0.4, -0.2, 0, 0.2, 0.4,  0.6), labels = c(-0.4, -0.2, 0 , 0.2, 0.4, 0.6))+
  geom_line(data = pred, aes(y = binding_pred), color = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  coord_cartesian(xlim = range(SOD1_allmut$abundance_score),
                  ylim = range(SOD1_allmut$binding_score)) +
  labs(x = "Normalized\nabundance score", y = "Normalized\nheterodimerization score") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = "none") +
  geom_text(data = r_corr, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 2, size = 10, inherit.aes = FALSE) 

p_loess

ggsave(p_loess, file = "correlation_apca_boca_label_by_residual_unfiltered.tiff", width = 10, height = 8)

#label by interface outlier
SOD1_allmut <- SOD1_allmut %>% 
  mutate(outlier_interface = case_when(
    interface == T & outlier == T ~ "TRUE",
    T ~ "FALSE"
  ))

SOD1_allmut <- SOD1_allmut %>% 
  mutate(outlier_type = case_when(
    interface == T & outlier == T ~ "interface outlier",
    interface == F & outlier == T ~ "non-interface outlier",
    T ~ "FALSE"
  ))

SOD1_allmut <- SOD1_allmut %>%
  dplyr::arrange(outlier_type == "interface outlier")

p_corr <- ggplot(data = SOD1_allmut, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 2.5, alpha = 0.8, aes(color = outlier_type)) +
  scale_color_manual(values = c("interface outlier" = "#2A5783", "non-interface outlier" = "#B9DDF1", "FALSE" = "lightgrey")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Normalized\nabundance score", y = "Normalized\nheterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.position = "top") +
  geom_text(data = r_corr, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 2, size = 10, inherit.aes = FALSE) 
p_corr
ggsave(p_corr, file = "correlation_apca_boca_label_by_outlier.tiff", width = 10, height = 8)


#Define residual categories####
#I define significant residuals (i.e outliers) as mutations whose residual module is higher than mad_resid
SOD1_allmut <- SOD1_allmut %>%
  mutate(
    residual_intensity = case_when(
      abs(residuals) > 2 * mad_resid ~ "high",
      TRUE ~ "low"
    ),
    residual_category = case_when( #add direction of residual
      residual_intensity == "high" & residuals > 0 ~ "high_positive",
      residual_intensity == "high" & residuals < 0 ~ "high_negative",
      residual_intensity == "low"                  ~ "low"
    ))

#Save modeled data####
SOD1_subsitutions_modeled <- SOD1_allmut %>%
  filter(mutation_type %in% c("subs", "stop"))

SOD1_insertions_modeled <- SOD1_allmut %>%
  filter(mutation_type == "ins")

SOD1_deletions_modeled <- SOD1_allmut %>% 
  filter(mutation_type == "del")


#save modeled data
write.csv(SOD1_subsitutions_modeled,file="SOD1_subs_modelled.csv", row.names = FALSE)
write.csv(SOD1_insertions_modeled,file ="SOD1_ins_modelled.csv", row.names = FALSE)
write.csv(SOD1_deletions_modeled,file = "SOD1_del_modelled.csv", row.names = FALSE)
write.csv(SOD1_allmut,file="SOD1_allmut_modelled.csv", row.names = FALSE)

#Save to color 3D structure####
#save data for the python script
SOD_final_subs_3D <- SOD1_subsitutions_modelled %>%
  select(ID, abundance_score, binding_score, category_fdr_apca, category_fdr_bpca, WT_AA, Pos, Mut_AA, residuals)

SOD_final_ins_3D <- SOD1_insertions_modelled %>%
  select(ID, abundance_score, binding_score, category_fdr_apca, category_fdr_bpca, WT_AA, Pos, Mut_AA, residuals)

SOD_final_del_3D <- SOD1_deletions_modelled %>%
  select(ID, abundance_score, binding_score, category_fdr_apca, category_fdr_bpca, WT_AA, Pos, Mut_AA, residuals)

write.csv(SOD_final_subs_3D,file="SOD1_subs_color3D.csv", row.names = FALSE)
write.csv(SOD_final_ins_3D,file ="SOD1_ins_color3D.csv", row.names = FALSE)
write.csv(SOD_final_del_3D,file = "SOD1_del_color3D.csv", row.names = FALSE)




