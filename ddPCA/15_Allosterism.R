library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(ggtext)
library(ggnewscale)
library(minpack.lm)

#load modeled data
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")
SOD1_insertions <-read.csv("SOD1_ins_modelled.csv")
SOD1_deletions <- read.csv("SOD1_del_modelled.csv")

SOD1_allmut <- bind_rows(SOD1_subsitutions, SOD1_insertions, SOD1_deletions)

#adjust data
SOD1_allmut <- SOD1_allmut %>% 
  mutate(s_location= factor(s_location, levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                                   "β4", "Zn binding loop", "β5", "loop 5", "β6", "loop 6",
                                                   "β7", "Electrostatic loop", "β8", "loop 8", "Zn binding residues", "Cu binding residues"))) %>% 
  mutate(outlier_interface = case_when(
    outlier == "TRUE" & interface =="TRUE" ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(outlier_nointerface = case_when(
    outlier == "TRUE" & interface =="FALSE" ~ TRUE,
    TRUE ~ FALSE
  )) 


#Identification of allosteric sites####
#load atom distance data (run the "Calculate atom distances" to generate this dataset)
min_distance <- read.csv("min_atom_distance_SOD1.csv")
min_distance <- rename(min_distance, Pos = resno_A)

SOD1_subsitutions$Pos <- as.numeric(SOD1_subsitutions$Pos)

SOD1_subsitutions_dist <- left_join(SOD1_subsitutions, min_distance, by = "Pos")
SOD1_subsitutions <- left_join(SOD1_subsitutions, min_distance, by = "Pos")

SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>% 
  mutate(WT_Pos = paste0(WT_AA, Pos, ""))

#per position
#calculate de median residual in interface residues
SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>%
  group_by(Pos) %>%
  mutate(median_residual_abs = median(abs(residuals), na.rm = TRUE),
         median_residual = median(residuals),
         residual_type = case_when(
           median_residual > 0 ~ "positive",
           median_residual < 0 ~ "negative")) %>%
  select(median_residual_abs, median_residual, Pos, WT_Pos, interface, min_dist_to_B, location, residual_type, s_location)



SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>% 
  mutate(residue_type = case_when(
    Pos %in% c(5,7,50,51,52,53,54,113,114,148,150,151,152,153) ~"interface",
    Pos %in% c(1,3,4,6,9,16,18,19,48,49,55,56,57,59,60,61,106,112,114,115,116,146,147,149) ~ "second_shell",
    TRUE ~ "rest"
    
  ))

SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>%
  distinct(Pos, .keep_all = TRUE)
SOD1_subsitutions_dist <- na.omit(SOD1_subsitutions_dist)

interface_subs <- SOD1_subsitutions %>% 
  filter(interface == TRUE)
interface <- SOD1_allmut %>% 
  filter(interface == T)

median_res_interface <- median(abs(interface$residuals))



p_min_dist_pos_subs <- ggplot(SOD1_subsitutions_dist, aes(x = min_dist_to_B, y = median_residual_abs)) +
  geom_point(aes(color = interface, shape = residual_type), size = 6, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5783", "FALSE" = "#81B1D6"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(
    data = subset(
      SOD1_subsitutions_dist,
      median_residual_abs > median_res_interface | min_dist_to_B < 3.881280),
    aes(label = Pos, color = interface),
    size = 6
  )  +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median heterodimerization residual|") +
  theme(axis.text = element_text(size = 33),
        axis.title = element_text(size = 36),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 3.881280, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_res_interface, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_pos_subs

#fit exponential model (is there an exponential allosteric decay?)
fit_exp <- nlsLM(
  median_residual_abs ~ a * exp(-b * min_dist_to_B) + c,
  data = SOD1_subsitutions_dist,
  start = list(
    a = max(SOD1_subsitutions_dist$median_residual_abs),
    b = 0.5,
    c = min(SOD1_subsitutions_dist$median_residual_abs)
  )
)

pred_df <- data.frame(
  min_dist_to_B = seq(
    min(SOD1_subsitutions_dist$min_dist_to_B),
    max(SOD1_subsitutions_dist$min_dist_to_B),
    length.out = 200
  )
)

pred_df$pred <- predict(fit_exp, newdata = pred_df)



p_min_dist_pos_subs <- ggplot(SOD1_subsitutions_dist, aes(x = min_dist_to_B, y = median_residual_abs)) +
  geom_point(aes(color = residual_type, shape = interface), size = 7, alpha = 1) +
  scale_color_manual(
    values = c("positive" = "#7C4D79", "negative" = "#AD9024"),
    labels = c("positive" = "HD residual > 0", "negative" = "HD residual < 0")) +
  scale_shape_manual(
    values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  geom_text_repel(
    data = subset(
      SOD1_subsitutions_dist,
      median_residual_abs > median_res_interface | min_dist_to_B < 3.881280),
    aes(label = Pos, color = residual_type),
    size = 6
  )  +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median heterodimerization residual|") +
  theme(axis.text = element_text(size = 33),
        axis.title = element_text(size = 36),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 3.881280, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_res_interface, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  ) +
  geom_line(
    data = pred_df,
    aes(x = min_dist_to_B, y = pred),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 1.5
  )

p_min_dist_pos_subs

ggsave(p_min_dist_pos_subs, file = "allosteric_sites_subs.tiff", width = 10, height = 10, dpi = 300)

#by mutations
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(residual_type = case_when(
    residuals < 0 ~ "negative",
    residuals > 0 ~ "positive"
  )) %>% 
  mutate(abs_residual = abs(residuals))


p_min_dist_pos_subs_mut <- ggplot(SOD1_subsitutions, aes(x = min_dist_to_B, y = abs_residual)) +
  geom_point(aes(color = residual_type, shape = interface), size = 6, alpha = 0.5) +
  scale_color_manual(
    values = c("positive" = "#7C4D79", "negative" = "#AD9024"),
    labels = c("positive" = "HD residual > 0", "negative" = "HD residual < 0")) +
 scale_shape_manual(
    values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Heterodimerization\nresidual|") +
  theme(axis.text = element_text(size = 33),
        axis.title = element_text(size = 36),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 3.881280, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_res_interface, linetype = "dashed", color = "black", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  ) 

p_min_dist_pos_subs_mut
ggsave(p_min_dist_pos_subs_mut, file = "allosteric_sites_subs_mut.tiff", width = 12, height = 10, dpi = 300)


fit_pos <- nlsLM(
  abs_residual ~ a * exp(-b * min_dist_to_B) + c,
  data = SOD1_subsitutions %>% filter(residual_type == "positive"),
  start = list(
    a = max(SOD1_subsitutions$abs_residual),
    b = 0.5,
    c = min(SOD1_subsitutions$abs_residual)
  )
)

fit_neg <- nlsLM(
  abs_residual ~ a * exp(-b * min_dist_to_B) + c,
  data = SOD1_subsitutions %>% filter(residual_type == "negative"),
  start = list(
    a = max(SOD1_subsitutions$abs_residual),
    b = 0.5,
    c = min(SOD1_subsitutions$abs_residual)
  )
)


pred_pos <- data.frame(
  min_dist_to_B = seq(min(SOD1_subsitutions$min_dist_to_B),
                      max(SOD1_subsitutions$min_dist_to_B),
                      length.out = 200),
  residual_type = "positive"
)
pred_pos$fit <- predict(fit_pos, newdata = pred_pos)

pred_neg <- data.frame(
  min_dist_to_B = seq(min(SOD1_subsitutions$min_dist_to_B),
                      max(SOD1_subsitutions$min_dist_to_B),
                      length.out = 200),
  residual_type = "negative"
)
pred_neg$fit <- predict(fit_neg, newdata = pred_neg)

pred_all <- rbind(pred_pos, pred_neg)



p_mut_type_allos <- ggplot(SOD1_subsitutions, aes(x = min_dist_to_B, y = abs_residual)) +
  geom_point(aes(color = residual_type, shape = interface), size = 6, alpha = 0.5) +
  scale_color_manual(
    values = c("positive" = "#7C4D79", "negative" = "#AD9024"),
    labels = c("positive" = "HD residual > 0", "negative" = "HD residual < 0")) +
  scale_shape_manual(
    values = c("TRUE" = 17, "FALSE" = 16),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  facet_wrap(~residual_type) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Heterodimerization\nresidual|") +
  theme(axis.text = element_text(size = 33),
        axis.title = element_text(size = 36),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 3.881280, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = median_res_interface, linetype = "dashed", color = "black", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  ) +
  geom_line(
    data = pred_all,
    aes(x = min_dist_to_B, y = fit),
    color = "red",
    linewidth = 1.5
  )

p_mut_type_allos


#Analysis of allosteric positions####
allosteric_positions <- SOD1_subsitutions_dist %>% 
  filter(interface == "FALSE" & median_residual_abs > median_res_interface)
allosteric_positions$allosteric_pos <- T
allosteric_positions <- left_join(SOD1_subsitutions, allosteric_positions, by = "Pos")
allosteric_positions <- allosteric_positions %>% 
  filter(allosteric_pos == T)

allosteric_positions <- allosteric_positions %>% 
  mutate(s_location.x = recode(s_location.x,
                               "Zn binding residues" = "Zn binding loop"))
allosteric_positions$s_location.x <- factor(allosteric_positions$s_location.x,
                                            levels = c("β2", "β5", "β6", "β7", "β8", "Zn binding loop", "loop 5", "loop 6"))

allosteric_positions <- allosteric_positions %>%
  mutate(Pos = factor(Pos, levels = sort(unique(Pos))))
allosteric_positions <- allosteric_positions %>% 
  mutate(residual_in_site = case_when(
    residuals > median_res_interface ~ "positive",
    residuals < -median_res_interface ~ "negative",
    TRUE ~ "non-allosteric"
  ))

#plot the distriubution of residuals at each allosteric position
p_loc_allos_sites <- ggplot(allosteric_positions, aes(x = Pos, y = residuals)) +
  geom_violin(position = position_dodge(), scale = "width", size = 1) +
  geom_jitter(position = position_jitter(), aes(fill = residual_in_site), size = 4, shape = 21, stroke = 0.5, color = "black")+
  scale_color_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024", "non-allosteric" = "#EDEDED")) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024", "non-allosteric" = "#EDEDED")) +
  
  #scale_x_discrete(labels = c("Zn binding loop" = "Zn binding\nloop")) +
  geom_hline(yintercept = median_res_interface, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = -median_res_interface, linetype = "dashed", size = 0.8) +
  theme_classic() +
  labs(x = "Allosteric position", y = "Heterodimerization residual") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
p_loc_allos_sites

ggsave(p_loc_allos_sites, file = "violin_allos_sites_distribution.tiff", width = 14, height = 7)

#highlight pathogenic mutations in allosteric sites (only showing negative allosteric mutations)
allosteric_positions_clinical <- allosteric_positions %>% 
  mutate(classification = case_when(
    ID %in% c("G16S", "G16C", "H71T", "H80T") ~ "pathogenic",
    TRUE ~"rest"
  )) %>% 
  filter(Pos %in% c(16,61,64,71,79,80,82,83))

allosteric_positions_clinical$WT_Pos <- factor(allosteric_positions_clinical$WT_Pos,
                                               levels = c("G16", "G61", "F64", "H71", "R79", "H80", "G82", "D83"))

p_loc_allos_sites_patho <- ggplot(allosteric_positions_clinical, aes(x = WT_Pos, y = residuals)) +
  geom_violin(position = position_dodge(), scale = "width", size = 0.5) +
  geom_jitter(position = position_jitter(), aes(fill = residual_in_site, shape = classification, size = classification, alpha = classification, stroke = classification),   color = "black")+
  scale_color_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024", "non-allosteric" = "#EDEDED")) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024", "non-allosteric" = "#EDEDED")) +
  scale_shape_manual(values = c("pathogenic" = 23, "rest" = 21)) +
  scale_size_manual(values = c("pathogenic" = 6, "rest" = 3)) +
  scale_alpha_manual(values = c("pathogenic" = 1, "rest" = 0.5)) +
  scale_discrete_manual(
    aesthetics = "stroke",
    values = c("pathogenic" = 2, "rest" = 0.5)
  ) +
  geom_hline(yintercept = median_res_interface, linetype = "dashed", size = 0.8) +
  geom_hline(yintercept = -median_res_interface, linetype = "dashed", size = 0.8) +
  theme_classic() +
  labs(x = "Negative allosteric site", y = "Heterodimerization\nresidual") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        legend.position = "none")
p_loc_allos_sites_patho

ggsave(p_loc_allos_sites_patho, file = "violin_allos_sites_distribution_patho.tiff", width = 13, height = 5)

##Abundance and binding scores of allosteric positions####
allos_ab <- allosteric_positions %>% 
  rename(score = abundance_score) %>% 
  select(WT_Pos, Pos, score)
allos_ab$assay <- "apca"

allos_bind <- allosteric_positions %>% 
  rename(score = binding_score) %>% 
  select(WT_Pos, Pos, score)
allos_bind$assay <- "bpca"

allos_res <- allosteric_positions %>% 
  rename(score = residuals) %>% 
  select(WT_Pos,Pos, score)
allos_res$assay <- "residual"

allosteric_pos_scores <- rbind(allos_ab, allos_bind, allos_res)

allosteric_pos_scores <- allosteric_pos_scores %>% 
  group_by(WT_Pos, assay, Pos) %>% 
  summarise(median = median(score))

#analyze significant differences between abundance and binding in allosteric sites
wilcox_allosteric <- allosteric_positions %>%
  group_by(WT_Pos) %>%
  summarise(
    p_value = wilcox.test(
      abundance_score,
      binding_score,
      paired = TRUE
    )$p.value
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      p_adj < 0.05 ~ TRUE,
      TRUE ~ FALSE))

signif_allos <- wilcox_allosteric %>% filter(p_adj < 0.05)
signif_allos <- signif_allos$WT_Pos

allosteric_pos_scores <- allosteric_pos_scores %>% 
  mutate(signif_0.05 = WT_Pos %in% signif_allos)



allosteric_pos_scores$WT_Pos <- factor(allosteric_pos_scores$WT_Pos, levels = c("G16", "G61", "F64", "H71", "R79", "H80",
                                                                                "G82", "D83", "L84", "G85", "V87", "A89", "G93",
                                                                                "V97", "I104", "S105", "L106", "G108", "C111", "I112",
                                                                                "L117", "V119", "A145", "G147", "I149") )


p_allos_scores <- ggplot() +
  geom_tile(data = subset(allosteric_pos_scores, assay %in% c("apca", "bpca")),aes(x = WT_Pos, y = assay, fill = median),color = "black") +
  scale_fill_gradient2(low = "darkorange",mid = "grey95",high = "darkgreen",midpoint = 0) +
  new_scale_fill() +
  geom_tile(data = subset(allosteric_pos_scores, assay == "residual"),aes(x = WT_Pos, y = assay, fill = median),color = "black") +
  scale_fill_gradient2(low = "#AD9024",mid = "grey95",high = "#7C4D79",midpoint = 0) +
  scale_y_discrete(labels = c("apca" = "Abundance", "bpca" = "Heterodimerization", "residual" = "HD residual")) +
  theme_classic() +
  labs(x = "Allosteric site", y = "") +
  theme(axis.text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 28),
        legend.position = "none")
p_allos_scores



p_allos <- ggarrange(p_loc_allos_sites,
                     p_allos_scores,
                             nrow = 2, align = "v", heights = c(2,1))
p_allos

ggsave(p_allos, file = "allosteric_sites_vioplots.tiff", width = 16, height = 8)



#check relative FDRs by position
allosteric_pos_median_scores <- allosteric_positions %>% 
  group_by(Pos) %>% 
  summarise(median_apca = median(abundance_score),
            median_bpca = median(binding_score))

#asign relative FDR by using the synonymous and stops thresholds
syn_thres_apca <- -0.19
syn_thres_bpca <- -0.05

stop_thres_apca <- -0.87
stop_thres_bpca <- -0.86

allosteric_pos_median_scores <- allosteric_pos_median_scores %>% 
  mutate(fdr_pos_apca = case_when(
    median_apca > stop_thres_apca & median_apca < syn_thres_apca ~ "low",
    median_apca > syn_thres_apca ~ "wt-like",
    TRUE ~ "stop-like"
  )) %>% 
  mutate(fdr_pos_bpca = case_when(
    median_bpca > stop_thres_bpca & median_bpca < syn_thres_bpca ~ "low",
    median_bpca > syn_thres_bpca ~ "wt-like",
    TRUE ~ "stop-like"
  ))


#abundance and binding scores of allosteric sites
allosteric_pos_scores <- bind_rows(allos_ab, allos_bind, allos_res)
allosteric_pos_scores <- allosteric_pos_scores %>% 
  filter(assay != "residual")

allosteric_pos_scores <- allosteric_pos_scores %>% 
  group_by(assay, Pos) %>% 
  summarise(score = score)


p_allosites <- ggplot(allosteric_pos_scores, aes(x = Pos, y = score, fill = assay, color = score)) +
  geom_violin(position = position_dodge(), scale = "width", size = 1) +
  geom_jitter(position = position_jitterdodge(), aes(color = score), size = 3)+
  scale_fill_manual(values = c("apca" = "white", "bpca" = "white")) +
  scale_color_gradient2(low = "darkorange", mid = "lightgrey", high = "darkgreen", midpoint = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
  theme_classic() +
  labs(x = "Allosteric position", y = "Score") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p_allosites

apca_fdr <- allosteric_pos_median_scores %>% 
  mutate(assay = "apca",
         fdr = fdr_pos_apca) %>% 
  select(Pos, assay, fdr)

bpca_fdr <- allosteric_pos_median_scores %>% 
  mutate(assay = "bpca",
         fdr = fdr_pos_bpca) %>% 
  select(Pos, assay, fdr)


allosteric_pos_median_scores <- rbind(apca_fdr, bpca_fdr)

p_allos_fdr <- ggplot() +
  geom_tile(data = subset(allosteric_pos_median_scores, assay == "apca"),aes(x = Pos, y = assay, fill = fdr), color = "black") +
  scale_fill_manual(values = c("wt-like" = "grey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  new_scale_fill() +
  geom_tile(data = subset(allosteric_pos_median_scores, assay == "bpca"),aes(x = Pos, y = assay, fill = fdr), color = "black") +
  scale_y_discrete(labels = c( "apca" = "Abundance", "bpca" = "Binding")) +
  scale_fill_manual(values = c("wt-like" = "grey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  theme_classic() +
  labs(x = "Allosteric site", y = "") +
  theme(axis.text = element_text(size = 25),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 28),
        legend.position = "none")

p_allos_fdr

p_allos_fdr_vio <- ggarrange(p_allosites,
                             p_allos_fdr,
                             nrow = 2, align = "v", heights = c(2,1))
p_allos_fdr_vio


ggsave(p_allos_fdr_vio, file = "allosteric_sites_fdr.tiff", width = 17, height = 9)

#Analazye enrichment of allosteric mutations####
##Substitutions####
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(residual_type = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  )) %>% 
  mutate(abs_residual = abs(residuals)) %>% 
  mutate(allosteric_mutation = case_when(
    abs_residual > median_res_interface & interface == FALSE ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(allosteric_mutation_type = case_when(
    allosteric_mutation == T & residual_type == "positive" ~ "positive",
    allosteric_mutation == T & residual_type == "negative" ~ "negative",
    T ~ "non_allosteric"
  ))

SOD1_subsitutions$allosteric_mutation_type <- factor(SOD1_subsitutions$allosteric_mutation_type,
                                                     levels = c("non_allosteric", "negative", "positive"))
allosteric_mutations <- SOD1_subsitutions %>% 
  filter(allosteric_mutation == T)


#not normalized to number of positions
p_side_allosterism <- ggplot(allosteric_mutations, aes(x=side_chain)) +
  geom_bar(aes(fill = location), width = 0.7, linewidth=0.5, color = "black") +
  scale_fill_manual(values = c("loop" = "#49525E", "β-sheet" = "lightgrey")) +
  labs(x = "Side-chain", y = "Number of allosteric mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22)) +
  coord_flip()
p_side_allosterism

ggsave(p_side_allosterism, file = "p_n_allosteric_mut_side_chain_subs.tiff", width = 10, height = 4)


allosteric_mutations_subs <- SOD1_subsitutions %>%
  dplyr::count(allosteric_mutation_type, name = "n")

p_pie_allosterism_subs <- ggplot(allosteric_mutations_subs, aes(x = "", y = n, fill = allosteric_mutation_type)) +
  geom_col(width = 1, color = "black") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024","non_allosteric" = "#EDEDED"),
                    labels = c("positive" = "Positive allostery", "negative" = "Negative allostery", "non_allosteric" = "Non-allosteric")) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Type of mutation") +
  theme(legend.position = "none")
p_pie_allosterism_subs

global_allosterism <- mean(SOD1_subsitutions$allosteric_mutation)

N_total <- nrow(SOD1_subsitutions)
K_total <- sum(SOD1_subsitutions$allosteric_mutation)

SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(is_allosteric_site = case_when(
    Pos %in% c(16,61,64,71,79,80,82,83,84,85,87,89,93,97,104,105,106,108,111,112,117,119,145,147,149) ~ T,
    T~F
  )) %>% 
  mutate(WT_Pos = paste0(WT_AA, Pos, ""))

per_pos_counts <- SOD1_subsitutions %>%
  group_by(WT_Pos, is_allosteric_site) %>%
  summarise(
    k_pos = sum(allosteric_mutation),   
    n_pos = n(),                 
    .groups = "drop"
  )

#Fisher exact test
fisher_pos <- function(k_pos, n_pos, K_total, N_total) {
  
  mat <- matrix(
    c(
      k_pos,
      n_pos - k_pos,
      K_total - k_pos,
      (N_total - n_pos) - (K_total - k_pos)
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  ft <- fisher.test(mat, alternative = "greater")
  
  c(
    OR = unname(ft$estimate),
    p_value = ft$p.value
  )
}


fisher_results <- per_pos_counts %>%
  rowwise() %>%
  mutate(
    test = list(fisher_pos(k_pos, n_pos, K_total, N_total)),
    OR = test["OR"],
    p_value = test["p_value"]
  ) %>%
  ungroup() %>%
  select(-test)

#FDR correction
fisher_results <- fisher_results %>%
  filter(is_allosteric_site == T)%>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")) %>% 
  mutate(
    sig = case_when(#p-value correction
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""))

fisher_results$WT_Pos <- factor(fisher_results$WT_Pos, levels = c("G16", "G61", "F64", "H71", "R79", "H80",
                                                                                "G82", "D83", "L84", "G85", "V87", "A89", "G93",
                                                                                "V97", "I104", "S105", "L106", "G108", "C111", "I112",
                                                                                "L117", "V119", "A145", "G147", "I149") )
#odds ratio of allosteric sites
p_odds <- ggplot(fisher_results, aes(x = WT_Pos, y = OR)) +
  geom_col(fill = "grey", color = "black") +
  geom_text(aes(label = sig), size = 8) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 1) +
  labs(x = "Allosteric site", y = "Odds ratio") +
  theme_classic() +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
p_odds

ggsave(p_odds, file = "p_odds_allosteric_sites.tiff", width = 14, height = 8)


#now per type of allosteric site
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(allosteric_site_type = case_when(
    Pos %in% c(16,61,64,71,79,80,82,83) ~ "negative",
    Pos %in% c(84,85,87,89,93,97,104,105,106,108,111,112,117,119,145,147,149) ~ "positive",
    T ~ "no"
  ))

per_type_counts <- SOD1_subsitutions %>%
  group_by(allosteric_site_type) %>%
  summarise(
    k_pos = sum(allosteric_mutation),   
    n_pos = n(),                  
    .groups = "drop"
  )

fisher_allos_type <- per_type_counts %>%
  rowwise() %>%
  mutate(
    test = list(fisher_pos(k_pos, n_pos, K_total, N_total)),
    OR = test["OR"],
    p_value = test["p_value"]
  ) %>%
  ungroup() %>%
  select(-test)

fisher_allos_type <- fisher_allos_type %>%
  #filter(is_allosteric_site == T)%>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")) %>% 
  mutate(
    sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""))

fisher_allos_type$allosteric_site_type <- factor(fisher_allos_type$allosteric_site_type,
                                                 levels = c("positive", "negative", "no"))

p_odds_allos_type <- ggplot(fisher_allos_type, aes(x = allosteric_site_type, y = OR)) +
  geom_col(fill = "grey", color = "black", width = 0.5) +
  geom_text(aes(label = sig), size = 8) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 1) +
  scale_x_discrete(labels = c("positive" = "Positive\nallosteric sites", "negative" = "Negative\nallosteric sites", "no" = "Non-allosteric\nsites")) +
  labs(x = "", y = "Odds ratio") +
  theme_classic() +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25))

p_odds_allos_type

ggsave(p_odds_allos_type, file = "p_odds_allosteric_types.tiff", width = 10, height = 8)


##Insertions####
interface_ins <- SOD1_insertions %>% 
  filter(interface == TRUE) 

SOD1_insertions <- SOD1_insertions %>% 
  mutate(residual_type = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  )) %>% 
  mutate(abs_residual = abs(residuals)) %>% 
  mutate(allosteric_mutation = case_when(
    abs_residual > median_res_interface & interface == FALSE ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(allosteric_mutation_type = case_when(
    allosteric_mutation == T & residual_type == "positive" ~ "positive",
    allosteric_mutation == T & residual_type == "negative" ~ "negative",
    T ~ "non_allosteric"
  ))

SOD1_insertions$allosteric_mutation_type <- factor(SOD1_insertions$allosteric_mutation_type,
                                                   levels = c("non_allosteric", "negative", "positive"))
allosteric_mutations_ins <- SOD1_insertions %>% 
  filter(allosteric_mutation == TRUE)


p_side_allosterism_ins <- ggplot(allosteric_mutations_ins, aes(x=side_chain)) +
  geom_bar(aes(fill = location), width = 0.5, linewidth=0.5, color = "black") +
  scale_fill_manual(values = c("loop" = "#49525E", "β-sheet" = "lightgrey")) +
  labs(x = "Side-chain", y = "Number of allosteric mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22)) +
  coord_flip()
p_side_allosterism_ins

allosteric_mutations_ins <- SOD1_insertions %>%
  dplyr::count(allosteric_mutation_type, name = "n")

p_pie_allosterism_ins <- ggplot(allosteric_mutations_ins, aes(x = "", y = n, fill = allosteric_mutation_type)) +
  geom_col(width = 1, color = "black") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024","non_allosteric" = "#EDEDED"),
                    labels = c("positive" = "Positive allostery", "negative" = "Negative allostery", "non_allosteric" = "Non-allosteric")) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Type of mutation") +
  theme(legend.position = "none")
p_pie_allosterism_ins

##Deletions####
interface_del <- SOD1_deletions %>% 
  filter(interface == TRUE) 

SOD1_deletions <- SOD1_deletions %>% 
  mutate(residual_type = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  )) %>% 
  mutate(abs_residual = abs(residuals)) %>% 
  mutate(allosteric_mutation = case_when(
    abs_residual > median_res_interface & interface == FALSE ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(allosteric_mutation_type = case_when(
    allosteric_mutation == T & residual_type == "positive" ~ "positive",
    allosteric_mutation == T & residual_type == "negative" ~ "negative",
    T ~ "non_allosteric"
  ))

SOD1_deletions$allosteric_mutation_type <- factor(SOD1_deletions$allosteric_mutation_type,
                                                  levels = c("non_allosteric", "negative", "positive"))
allosteric_mutations_del <- SOD1_deletions %>% 
  filter(allosteric_mutation == T )

allosteric_mutations_del <- SOD1_deletions %>%
  dplyr::count(allosteric_mutation_type, name = "n")

p_pie_allosterism_del <- ggplot(allosteric_mutations_del, aes(x = "", y = n, fill = allosteric_mutation_type)) +
  geom_col(width = 1, color = "black") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024","non_allosteric" = "#EDEDED"),
                    labels = c("positive" = "Positive allostery", "negative" = "Negative allostery", "non_allosteric" = "Non-allosteric")) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Type of mutation") +
  theme(legend.position = "none")
p_pie_allosterism_del

p_pie_all <- ggarrange(p_pie_allosterism_subs,
                       p_pie_allosterism_ins,
                       p_pie_allosterism_del,
                       ncol = 3, common.legend = F)

p_pie_all

ggsave(p_pie_all, file = "pie_n_allosteric_mutations.tiff", width = 10, height = 8)


#Second shell allosteric mutations####
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(second_shell = case_when(
    Pos %in% c(1,3,4,6,9,16,18,19,48,49,55,56,57,59,60,61,106,112,114,115,116,146,147,149) ~ TRUE,
    TRUE ~ FALSE
  ))

n_second <- length(c(1,3,4,6,9,16,18,19,48,49,55,56,57,59,60,61,106,112,114,115,116,146,147,149))
n_total <- length(unique(SOD1_subsitutions$Pos))
n_interface <- length(unique(interface$Pos))
n_total_clean <- n_total - n_interface

n_rest <- n_total- n_second

norm_df <- data.frame(
  second_shell = c(TRUE, FALSE),
  n_pos = c(n_second, n_rest)
)


allosteric_mutations_2nd_shell <- SOD1_subsitutions %>%
  filter(allosteric_mutation == TRUE) %>%
  dplyr::count(second_shell, name = "n") %>%
  left_join(norm_df, by = "second_shell") %>%
  mutate(n_norm = n / n_pos,
         prop = n_norm/sum(n_norm))

# allosteric_mutations_2nd_shell <- SOD1_subsitutions %>% 
#   group_by(second_shell) %>% 
#   count(allosteric_mutation) %>% 
#   mutate(prop = n/sum(n))


allosteric_mutations_2nd_shell$second_shell <- factor(allosteric_mutations_2nd_shell$second_shell, levels = 
                                                        c("TRUE", "FALSE"))

tab <- table(SOD1_subsitutions$second_shell,
             SOD1_subsitutions$allosteric_mutation)
chi_test <- chisq.test(tab)
p_val <- chi_test$p.value
p_label <- paste0("p = ", signif(p_val, 3))

p_2nd_shell <- ggplot(allosteric_mutations_2nd_shell, aes(x = second_shell, y = prop)) +
  geom_col(fill = "grey", color = "black", width = 0.4) +
  scale_x_discrete(labels = c( "FALSE" = "Rest", "TRUE" = "2nd shell"),) +
  labs(x = "", y = "Proportion of\n allosteric mutations") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 18),
        legend.position = "none") +
  annotate("text",
           x = 1.5, 
           y = max(allosteric_mutations_2nd_shell$prop) * 1.1,
           label = p_label,
           size = 9)
p_2nd_shell
#1.3 fold

ggsave(p_2nd_shell, file = "allosteric_mut_2ndshell_rest.tiff", width = 8, height = 6)


#Allosterism by mutation type####
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(Pos_mut = paste0(Pos, Mut_AA))

SOD1_insertions <- SOD1_insertions %>% 
  mutate(Pos_mut = paste0(Pos, Mut_AA))

subs_ins <- inner_join(SOD1_subsitutions, SOD1_insertions, by = "Pos_mut")
subs_ins <- rename(subs_ins, residual_subs = residuals.x, residual_ins = residuals.y)

#scatter plot
r_labels_sub_ins <- subs_ins %>%
  summarise(
    R = cor(residual_subs, residual_ins, method = "pearson"),
    p_value = cor.test(residual_subs, residual_ins, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_cor_subs_ins <- ggplot(subs_ins, aes(x = residual_subs, y = residual_ins)) +
  geom_point(size = 4, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions\nheterodimerization residual", y = "Insertions\nheterodimerization residual", color = "") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        legend.position = "top",
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(
    override.aes = list(size = 10))) +
  geom_text(data = r_labels_sub_ins, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_cor_subs_ins


#subs vs dels residuals
median_subs_residual <- SOD1_subsitutions %>% 
  group_by(Pos, is_allosteric_site, allosteric_site_type) %>% 
  summarise(median_residual = median(residuals))

subs_del <- inner_join(median_subs_residual, SOD1_deletions, by = "Pos")
subs_del <- rename(subs_del, residual_subs = median_residual, residual_del = residuals)

#scatter plot
r_labels_sub_del <- subs_del %>%
  ungroup() %>%
  summarise(
    R = cor(residual_subs, residual_del, use = "complete.obs"),
    p_value = cor.test(residual_subs, residual_del)$p.value
  ) %>% 
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )

p_cor_subs_del <- ggplot(subs_del, aes(x = residual_subs, y = residual_del)) +
  geom_point(size = 4, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nheterodimerization residual", y = "Deletions\nheterodimerization residual", color = "") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        legend.position = "top",
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(
    override.aes = list(size = 10))) +
  geom_text(data = r_labels_sub_del, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_cor_subs_del

#ins vs dels residuals
median_ins_residual <- SOD1_insertions %>% 
  group_by(Pos) %>% 
  summarise(median_residual = median(residuals)) %>% 
  mutate(allosteric_position = case_when(
    abs(median_residual) > median_res_interface ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(allosteric_site_type = case_when(
    median_residual > 0 ~ "positive",
    median_residual < 0 ~ "negative",
    TRUE ~ "0"
  ))

ins_del <- inner_join(median_ins_residual, SOD1_deletions, by = "Pos")
ins_del <- rename(ins_del, residual_ins = median_residual, residual_del = residuals)
ins_del <- na.omit(ins_del)

#scatter plot
r_labels_ins_del <- ins_del %>%
  summarise(
    R = cor(residual_ins, residual_del, method = "pearson"),
    p_value = cor.test(residual_ins, residual_del, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_cor_ins_del <- ggplot(ins_del, aes(x = residual_ins, y = residual_del)) +
  geom_point(size = 4, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Insertions median\nheterodimerization residual", y = "Deletions\nheterodimerization residual") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        legend.position = "top",
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(
    override.aes = list(size = 10))) +
  geom_text(data = r_labels_ins_del, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_cor_ins_del


p_all_comparisons <- ggarrange(p_cor_subs_ins,
                               p_cor_subs_del,
                               p_cor_ins_del,
                               ncol = 3, align = "hv")
p_all_comparisons

ggsave(p_all_comparisons, file = "p_correlation_residual_mut_type.tiff", width = 20, height = 6)

subs_ins <- rename(subs_ins, allosteric_mutation_subs = allosteric_mutation.x, allosteric_mutation_ins = allosteric_mutation.y,
                   allosteric_mutation_type_subs = allosteric_mutation_type.x, allosteric_mutation_type_ins = allosteric_mutation_type.y)

subs_ins <- subs_ins %>% 
  mutate(residual_comparison = case_when(
    allosteric_mutation_type_subs == "non_allosteric" & allosteric_mutation_type_ins == "non_allosteric" ~ "non_allosteric",
    allosteric_mutation_type_subs == "positive" & allosteric_mutation_type_ins == "positive" ~ "same_allostery",
    allosteric_mutation_type_subs == "negative" & allosteric_mutation_type_ins == "negative" ~ "same_allostery",
    allosteric_mutation_type_subs == "positive" & allosteric_mutation_type_ins == "negative" ~ "opposite_allostery",
    allosteric_mutation_type_subs == "negative" & allosteric_mutation_type_ins == "positive" ~ "opposite_allostery",
    TRUE ~"rest"
  ))

subs_ins_filtered <- subs_ins %>%
  filter(residual_comparison %in% c("same_allostery", "opposite_allostery")) %>%
  dplyr::count(residual_comparison) %>%
  mutate(perc = n / sum(n) * 100)

p_subs_ins_allosteric <- ggplot(subs_ins_filtered, aes(x = residual_comparison, y = perc)) +
  geom_col(width = 0.5) +
  scale_x_discrete(labels = c("opposite_allostery" = "Opposite allosteric\neffect", "same_allostery" = "Same allosteric\neffect")) +
  theme_classic() +
  labs(x = "", y = "Percentage of mutations (%)", title = "Substitutions vs insertions") +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
p_subs_ins_allosteric





subs_del <- subs_del %>% 
  mutate(residual_comparison = case_when(
    allosteric_site_type == "no" & allosteric_mutation_type == "non_allosteric" ~ "non_allosteric",
    allosteric_site_type == "positive" & allosteric_mutation_type == "positive" ~ "same_allostery",
    allosteric_site_type == "negative" & allosteric_mutation_type == "negative" ~ "same_allostery",
    allosteric_site_type == "positive" & allosteric_mutation_type == "negative" ~ "opposite_allostery",
    allosteric_site_type == "negative" & allosteric_mutation_type == "positive" ~ "opposite_allostery",
    TRUE ~"rest"
  ))

subs_del_filtered <- subs_del %>% 
  filter(residual_comparison %in% c("same_allostery", "opposite_allostery")) %>% 
  dplyr::count(residual_comparison) %>%
  mutate(perc = n/sum(n) * 100)

p_subs_del_allosteric <- ggplot(subs_del_filtered, aes(x = residual_comparison, y = n)) +
  geom_col(width = 0.5) +
  scale_x_discrete(labels = c("opposite_allostery" = "Opposite allosteric\neffect", "same_allostery" = "Same allosteric\neffect")) +
  theme_classic() +
  labs(x = "", y = "Nuumber of positions", title = "Substitutions vs deletions\n(by position)") +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
p_subs_del_allosteric



ins_del <- ins_del %>% 
  mutate(residual_comparison = case_when(
    allosteric_site_type == FALSE & allosteric_mutation_type == "non_allosteric" ~ "non_allosteric",
    allosteric_site_type == "positive" & allosteric_mutation_type == "positive" ~ "same_allostery",
    allosteric_site_type == "negative" & allosteric_mutation_type == "negative" ~ "same_allostery",
    allosteric_site_type == "positive" & allosteric_mutation_type == "negative" ~ "opposite_allostery",
    allosteric_site_type == "negative" & allosteric_mutation_type == "positive" ~ "opposite_allostery",
    TRUE ~"rest"
  ))

ins_del_filtered <- ins_del %>% 
  filter(residual_comparison %in% c("same_allostery", "opposite_allostery")) %>% 
  dplyr::count(residual_comparison) %>%
  mutate(perc = n/sum(n) * 100)

p_ins_del_allosteric <- ggplot(ins_del_filtered, aes(x = residual_comparison, y = n)) +
  geom_col(width = 0.5) +
  scale_x_discrete(labels = c("opposite_allostery" = "Opposite allosteric\neffect", "same_allostery" = "Same allosteric\neffect")) +
  theme_classic() +
  labs(x = "", y = "Nuumber of positions", title = "Insertions vs deletions\n(by position)") +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5))

p_ins_del_allosteric


ggsave(p_subs_ins_allosteric, file = "allosteruc_effects_subs_ins.tiff", width = 8, height = 6)
ggsave(p_subs_del_allosteric, file = "allosteruc_effects_subs_del.tiff", width = 8, height = 6)
ggsave(p_ins_del_allosteric, file = "allosteruc_effects_del_ins.tiff", width = 8, height = 6)

#Allosteric pockets####
#load pockets prediction by CASTp
pockets_pred <- read.csv("SOD1_pockets.csv", sep = ";")

#join allosteric sites with allosteric pockets data frames
pockets_assignment <- inner_join(pockets_pred, SOD1_subsitutions, by = "Pos")
pockets_assignment <- pockets_assignment %>% 
  mutate(Mut_aa = str_extract(ID, "[A-za-z]$"))


heatmap_pockets_subs <- pockets_assignment %>%
  group_by(Pocket) %>%
  mutate(Pos_fac = factor(Pos, levels = sort(unique(Pos)))) %>%
  ungroup()
heatmap_pockets_subs$Pos <- as.character(heatmap_pockets_subs$Pos)

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

min_res<-min(heatmap_pockets_subs$residuals, na.rm = TRUE)
max_res<-max(heatmap_pockets_subs$residuals, na.rm = TRUE)
cols_res <- c(colorRampPalette(c( "#B99C33", "grey95"))((-min_res/(-min_res+max_res)*100)-0.5), colorRampPalette("grey95")(1),
              colorRampPalette(c("grey95",  "#AC7299"), bias=1)((max_res/(-min_res+max_res)*100)-0.5))

heatmap_pockets_subs <- heatmap_pockets_subs %>%
  mutate(
    Pos_label = ifelse(
      is_allosteric_site,
      paste0("<span style='color:red;'>", Pos_fac, "</span>"),
      as.character(Pos_fac)
    )
  )

p_heatmap_pockets_subs<-ggplot(heatmap_pockets_subs) +
  geom_tile(aes(x = WT_Pos,
                y = factor(Mut_aa, levels = rev(vectorAA)), fill = residuals),
            color = 'white', size = 0.1) +
  facet_wrap(~Pocket, ncol = 5, scales = "free_x", labeller = as_labeller(function(x) paste("Pocket", x))) +
  theme_minimal()+
  theme()+
  labs(x="SOD1 position", y="Mutant amino acid", fill="HD residual")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.8), 
        strip.text = element_text(size = 22),
        axis.text.x = ggtext::element_markdown(size = 14, angle = 90,
                                               hjust = 0.5, vjust = 0.5)) +
  scale_fill_gradientn(colours=cols_res, limits=c(min_res,max_res), na.value = "grey") +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5, direction = "vertical", title.position = "top", ticks = TRUE, ticks.colour = "black", frame.colour = "white"))
p_heatmap_pockets_subs

ggsave(p_heatmap_pockets_subs, file = "pockets_heatmaps_residuals.tiff", width = 15, height = 18, dpi = 300)


#define pockets enriched with allosteric mutations. They should be enriched with the same allostery type
pockets_assignment_counts <- pockets_assignment %>% 
  group_by(Pocket, residual_type) %>% 
  dplyr::count(allosteric_mutation) %>%
  mutate(prop = n/sum(n))


pockets_assignment_counts <- pockets_assignment_counts %>% 
  filter(allosteric_mutation == T)

pockets_assignment_counts$Pocket <- factor(pockets_assignment_counts$Pocket)
pockets_assignment_counts <- pockets_assignment_counts %>% 
  mutate(allosteric_sites = case_when(
    Pocket %in% c("5", "6", "10", "11", "15", "20", "28") ~ T,
    T ~ F
  ))

pockets_with_allosteric_sites <- pockets_assignment_counts %>% 
  filter(allosteric_sites == T)

p_pockets_allosteric_sites <- ggplot(pockets_with_allosteric_sites, aes(x = Pocket, y = prop, fill = residual_type)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024")) +
  labs(x = "Predicted SOD1 pocket\nwith allosteric sites", y = "Proportion of\n allosteric mutations", fill = "Type of allostery") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "top")

p_pockets_allosteric_sites
ggsave(p_pockets_allosteric_sites, file = "pockets_allosteric_sites.tiff", width = 7, height = 5)

pockets_without_allosteric_sites <- pockets_assignment_counts %>% 
  filter(allosteric_sites == F)

p_pockets_no_allosteric_sites <- ggplot(pockets_without_allosteric_sites, aes(x = Pocket, y = prop, fill = residual_type)) +
  geom_col(color = "black", width = 0.6) +
  scale_fill_manual(values = c("positive" = "#7C4D79","negative" = "#AD9024")) +
  labs(x = "Predicted SOD1 pocket", y = "Proportion of\n allosteric mutations", fill = "Type of allostery") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "top") 

p_pockets_no_allosteric_sites

ggsave(p_pockets_no_allosteric_sites, file = "pockets_allosteric_mutations.tiff", width = 9, height = 5)
