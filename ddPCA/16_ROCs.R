library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(isotree)
library(pROC)
library(PRROC)

#load SOD1 dms data, predictors and clinvar data
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")

#load predictors
alphamiss <- read.csv("am_scores.csv", sep = ";")
esm <- read.csv("esm_scores.csv", sep = ",")
eve <- read.csv("eve_scores.csv", sep = ",")
revel <- read.csv("revel_benign_set.csv")

#load clinical information
clinvar <- read.delim("clinvar_result (6).txt")
gnomad <- read.csv("gnomad_SOD1.csv", sep = ";")

#adjust predictors/clinical data to work in the same format
alphamiss$ID <- sapply(alphamiss$ID, decrement_id_pos, USE.NAMES = FALSE)
alphamiss <- rename(alphamiss, am_score = pathogenicity_score)
alphamiss <- select(alphamiss, ID, am_score)

esm <- rename(esm, ID = variant, esm_score = score)
esm$ID <- sapply(esm$ID, decrement_id_pos, USE.NAMES = FALSE)
esm <- select(esm, ID, esm_score)

eve$ID <- gsub("-", "", eve$ID)
eve <- select(eve, ID, eve_score)

revel$ID <- sapply(revel$ID, decrement_id_pos, USE.NAMES = FALSE)
revel <- select(revel, ID, REVEL)
revel <- dplyr::distinct(revel)


SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(wt_aa = str_extract(ID,"^[A-za-z]"),
         mut_aa = str_extract(ID, "[A-za-z]$")) %>% 
  mutate(syn = case_when(
    wt_aa == mut_aa ~ T,
    T~F)) %>% 
  filter( syn == F)
SOD1_subsitutions <- select(SOD1_subsitutions, ID, abundance_score, binding_score, residuals, category_fdr_apca, category_fdr_bpca)

clinvar <- clinvar %>%
  mutate(ID = str_replace(Protein.change, 
                          "(?<=^[A-Z])\\d+(?=[A-Z]$)", 
                          function(x) as.character(as.numeric(x) - 1)))
clinvar <- clinvar %>%
  rename(Classification = Germline.classification) %>%
  select(ID, Classification)

clinvar <- clinvar %>%
  mutate(Classification = recode(Classification, #for simplicity, I define Pathogenic/Likely pathogenic as Pathogenic
                                 "Pathogenic/Likely pathogenic" = "Pathogenic",
                                 "Likely pathogenic" = "Pathogenic")) %>% 
  mutate(multiple_submitters = case_when(
    ID %in% c("K3E", "A4S", "A4T", "A4V", "G12R", "V14L", "V14M", "G16A", "F20L",
              "Q22L", "G37R", "L38V", "G41S", "G41D", "H43R", "H46R", "V47A", "H48R",
              "N65S", "G72S", "D76Y", "L84F", "G85S", "G85R", "N86S", "V87M", "V87A",
              "A89V", "G93R", "G93D", "A95T", "E100K", "E100G", "S105L", "L106V", "I113T",
              "R115G", "V119L", "D124G", "D124V", "T137I", "N139K", "L144S", "L144F", "G147D",
              "V148I", "V148G", "I149T") ~ TRUE,
    TRUE ~ FALSE
  ))

gnomad <- rename(gnomad, ID = Protein.Consequence)
gnomad$ID <- sapply(gnomad$ID, convert_hgvsp_to_id)
gnomad$ID <- sapply(gnomad$ID, decrement_id_pos)
gnomad <- select(gnomad, ID, Allele.Count, Allele.Frequency, ClinVar.Germline.Classification, VEP.Annotation)

#see allele count distribution
p_allele_count <- ggplot(gnomad, aes(x = Allele.Count)) +
  geom_histogram(bins = 100)
p_allele_count


#ROCs####
#with VAMP-seq####
#here I use the reference set what was used in Axakova et al. AJHG, 2025 (i.e. benign controls = variants in gnomad with no clinical annotations)
invitae_set <- read.csv("20250212_refSetInvitae_toSave.csv")
invitae_set <- invitae_set %>% 
  mutate(label = case_when(
    referenceSet == "Positive" ~ 1,
    referenceSet == "Negative" ~ 0
  ))

#work in the same ID format
invitae_set$ID <- sapply(invitae_set$hgvsp, convert_hgvsp_to_id, USE.NAMES = FALSE)
invitae_set$ID <- sapply(invitae_set$ID, decrement_id_pos, USE.NAMES = FALSE)

invitae_set <- full_join(gnomad, invitae_set, by = "ID")
#invitae_set <- na.omit(invitae_set)

#I want to check if the benign set can be refined by excluding extremely weird variants
#calculate quantile 0.85 of allele count in the negative reference
negative_reference <- invitae_set %>% 
  filter(referenceSet == "Negative")
negative_reference <- select(negative_reference, -label)

write.csv(negative_reference, file = "unfiltered_benign_set_SOD1.csv")


q90_benign_allelecount <- quantile(negative_reference$Allele.Count, probs = c(0.85))
q90_benign_allelecount

#discard negative reference variants with allele count lower than 5
negative_reference_filt <- negative_reference %>% 
  filter(Allele.Count >=q90_benign_allelecount)

negative_reference_disc <- negative_reference %>% 
  filter(Allele.Count < q90_benign_allelecount)

write.csv(negative_reference_filt, file = "filtered_benign_set_SOD1.csv")

invitae_set_filtered <- invitae_set %>%
  filter(!ID %in% negative_reference_disc$ID)
invitae_set <- inner_join(invitae_set, SOD1_subsitutions, by = "ID")


#load VAMP-seq data
fritz_ab <- read.csv("2025_merged_map_unfloored_abundance.csv")
fritz_ab$hgvs_pro <- sapply(fritz_ab$hgvs_pro, convert_hgvsp_to_id, USE.NAMES = FALSE)
fritz_ab$hgvs_pro <- sapply(fritz_ab$hgvs_pro, decrement_id_pos, USE.NAMES = FALSE)
fritz_ab$ID <- fritz_ab$hgvs_pro

fritz_ab <- fritz_ab %>% 
  mutate(stop = grepl("\\*", ID)) %>% 
  filter(stop == F)
fritz_ab <- fritz_ab %>% 
  mutate(wt_aa = str_extract(ID,"^[A-za-z]"),
         mut_aa = str_extract(ID, "[A-za-z]$")) %>% 
  mutate(syn = case_when(
    wt_aa == mut_aa ~ T,
    T~F))

#join predictors with SOD1 aPCA data
SOD1_predictions <- SOD1_subsitutions %>%
  left_join(alphamiss, by = "ID") %>%
  left_join(esm, by = "ID") %>%
  left_join(eve, by = "ID")

SOD1_predictions <- select(SOD1_predictions, ID, abundance_score, binding_score, am_score, esm_score, eve_score, residuals, category_fdr_apca, category_fdr_bpca)

#and with Fritz reference and VAMP-seq datasets
SOD1_predictions <- left_join(SOD1_predictions, invitae_set, by = "ID")
SOD1_predictions <- left_join(SOD1_predictions, fritz_ab, by = "ID")

SOD1_predictions <- rename(SOD1_predictions, abundance_score = abundance_score.x, binding_score = binding_score.x, vampseq = score, residuals = residuals.x)

# SOD1_predictions_filter_benign <- SOD1_subsitutions %>%
#   left_join(alphamiss, by = "ID") %>%
#   left_join(esm, by = "ID") %>%
#   left_join(eve, by = "ID")
# 
# SOD1_predictions_filter_benign <- select(SOD1_predictions_filter_benign, ID, abundance_score, binding_score, am_score, esm_score, eve_score)
# SOD1_predictions_filter_benign <- left_join(SOD1_predictions_filter_benign, invitae_set_filtered, by = "ID")
# SOD1_predictions_filter_benign <- left_join(SOD1_predictions_filter_benign, fritz_ab, by = "ID")

# SOD1_predictions_clean <- na.omit(SOD1_predictions_filter_benign)
# SOD1_predictions_clean <- SOD1_predictions_clean %>% 
#   mutate(abundance_score_prc = -abundance_score,
#          binding_score_prc = -binding_score,
#          vampseq_prc = -score)

#generate the ROC
mylist_roc <- list()
auc_labels <- c()

#ab_pca
apca_roc <- SOD1_predictions[!is.na(SOD1_predictions$abundance_score) & !is.na(SOD1_predictions$label), ]
glm1 <- glm(label ~abundance_score, data = apca_roc, family = binomial)
roc1 <- roc(apca_roc$label, glm1$fitted.values)
mylist_roc[["AbundancePCA"]] <- roc1
auc_labels <- c(auc_labels, paste0("Abundance PCA (AUC = ", round(roc1$auc, 2), ")"))


#bpca
bpca_roc <- SOD1_predictions[!is.na(SOD1_predictions$binding_score) & !is.na(SOD1_predictions$label), ]
glm2 <- glm(label ~binding_score, data = bpca_roc, family = binomial)
roc2 <- roc(bpca_roc$label, glm2$fitted.values)
mylist_roc[["Heterodimerization"]] <- roc2
auc_labels <- c(auc_labels, paste0("Heterodimerization (AUC = ", round(roc2$auc, 2), ")"))


#alphamissense
am_roc <- SOD1_predictions[!is.na(SOD1_predictions$am_score) & !is.na(SOD1_predictions$label), ]
glm3 <- glm(label ~am_score, data = am_roc, family = binomial)
roc3 <- roc(am_roc$label, glm3$fitted.values)
mylist_roc[["Alphamissense"]] <- roc3
auc_labels <- c(auc_labels, paste0("Alphamissense (AUC = ", round(roc3$auc, 2), ")"))

#ESM
esm_roc <- SOD1_predictions[!is.na(SOD1_predictions$esm_score) & !is.na(SOD1_predictions$label), ]
glm4 <- glm(label ~esm_score, data = esm_roc, family = binomial)
roc4 <- roc(esm_roc$label, glm4$fitted.values)
mylist_roc[["ESM-b1"]] <- roc4
auc_labels <- c(auc_labels, paste0("ESM-b1 (AUC = ", round(roc4$auc, 2), ")"))

#Eve
eve_roc <- SOD1_predictions[!is.na(SOD1_predictions$eve_score) & !is.na(SOD1_predictions$label), ]
glm5 <- glm(label ~eve_score, data = eve_roc, family = binomial)
roc5 <- roc(eve_roc$label, glm5$fitted.values)
mylist_roc[["EVE"]] <- roc5
auc_labels <- c(auc_labels, paste0("EVE (AUC = ", round(roc5$auc, 2), ")"))


#vampseq
vampseq_roc <- SOD1_predictions[!is.na(SOD1_predictions$vampseq) & !is.na(SOD1_predictions$label), ]
glm7 <- glm(label ~vampseq, data = vampseq_roc, family = binomial)
roc7 <- roc(vampseq_roc$label, glm7$fitted.values)
mylist_roc[["VAMP-seq"]] <- roc7
auc_labels <- c(auc_labels, paste0("VAMP-seq (AUC = ", round(roc7$auc, 2), ")"))

auc_values <- sapply(mylist_roc, function(x) x$auc)
scores_ordered <- names(sort(auc_values, decreasing = TRUE))
auc_labels_ordered <- paste0(scores_ordered, " (AUC = ", round(sort(auc_values, decreasing = TRUE), 2), ")")

colors_ordered <- c("#E9AF1D", "#0A1829", "#81B357", "#1F733D", "purple", "#ACC7A0") 


g_roc <- ggroc(mylist_roc[scores_ordered], legacy.axes = TRUE, size = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 28),
        legend.position = c(0.7, 0.3),  
        legend.key.height = unit(1.2, "cm"),
        legend.background = element_rect(
          fill = alpha("white", 0.8),  
          colour = "black",            
          linewidth = 0.5              
        ),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        title = element_text(size = 20)) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  geom_abline(linetype = "dashed") +
  scale_color_manual(
    values = colors_ordered,
    labels = auc_labels_ordered,
    guide = guide_legend(
      override.aes = list(color = colors_ordered, size = 3)  
    )
  )

g_roc

ggsave(g_roc, file = "roc_predictors.tiff", width = 12, height = 10)


#with clinical controls####
controls_ddpca <- read.csv("SOD1_DMS_clinical_ctrl_final.csv", sep = ";")
controls_ddpca$Class <- factor(controls_ddpca$Class,
                               levels = c("Synonymous", "Negative", "Negative_filt", "Likely benign", "Pathogenic", "Likely pathogenic"))

##scores distribution of controls####
controls_vampseq <- left_join(fritz_ab, controls_ddpca, by = "ID")
controls_am <- left_join(alphamiss, controls_ddpca, by = "ID")
controls_esm <- left_join(esm, controls_ddpca, by = "ID")
controls_eve <- left_join(eve, controls_ddpca, by = "ID")

controls_ddpca_long  <- controls_ddpca %>%
  pivot_longer(
    cols = c(abundance_score, binding_score),
    names_to = "assay",
    values_to = "score"
  ) %>% 
  filter(Class != "VOUS")

controls_ddpca_long <- na.omit(controls_ddpca_long)

p_controls_dist_ddPCA <- ggplot(controls_ddpca_long, aes(x = Class, y = score, color = assay)) +
  geom_violin(scale = "width", size = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 1), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("abundance_score" = "#0A1829",
                                "binding_score" = "#E9AF1D")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "", y = "Normalized score") +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 22),
    legend.position = "none"
    
  )
p_controls_dist_ddPCA

controls_vampseq  <- controls_vampseq %>%
  pivot_longer(
    cols = score,
    names_to = "assay",
    values_to = "score"
  ) %>% 
  filter(Class != "VOUS")

controls_vampseq <- na.omit(controls_vampseq)

p_controls_dist_vampseq <- ggplot(controls_vampseq, aes(x = Class, y = score, color = assay)) +
  geom_violin(scale = "width", size = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 1), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("score" = "darkgrey")) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "", y = "VAMP-seq score") +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 22),
    legend.position = "none"
  )
p_controls_dist_vampseq

controls_am_long  <- controls_am %>%
  pivot_longer(
    cols = am_score,
    names_to = "assay",
    values_to = "score"
  ) %>% 
  filter(Class != "VOUS")


controls_am_long <- na.omit(controls_am_long)

p_controls_dist_am <- ggplot(controls_am_long, aes(x = Class, y = score, color = assay)) +
  geom_violin(scale = "width", size = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 1), size = 3, alpha = 0.6) +
  scale_color_manual(values = c("am_score" = "darkgrey")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "", y = "Alphamissense score") +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 22),
    legend.position = "none"
    
  )
p_controls_dist_am

p_controls_dist <- ggarrange(p_controls_dist_ddPCA,
                             p_controls_dist_vampseq,
                             p_controls_dist_am,
                             ncol=3)
p_controls_dist

ggsave(p_controls_dist, file = "clinical_ctrls_scores_distribution.tiff", width = 22, height = 6)


#use our clinical control strategy
control_4 <- controls_ddpca %>% 
  filter(Class %in% c("Pathogenic", "Likely pathogenic", "Negative_filt", "Synonymous", "Likely benign")) %>% 
  mutate(
    label = case_when(
      Class %in% c("Pathogenic", "Likely pathogenic") ~ 1,
      Class %in% c("Negative_filt", "Synonymous", "Likely benign") ~ 0))

control_4_dms <- inner_join(controls_vampseq, control_4, by = "ID")

#crate the ROC
mylist_roc <- list()
auc_labels <- c()

#ab_pca
apca_roc <- control_4[!is.na(control_4$abundance_score) & !is.na(control_4$label), ]
glm1 <- glm(label ~abundance_score, data = apca_roc, family = binomial)
roc1 <- roc(apca_roc$label, glm1$fitted.values)
mylist_roc[["abundance"]] <- roc1
auc_labels <- c(auc_labels, paste0("abundance (AUC = ", round(roc1$auc, 2), ")"))

#bpca
bpca_roc <- control_4[!is.na(control_4$binding_score) & !is.na(control_4$label), ]
glm2 <- glm(label ~binding_score, data = bpca_roc, family = binomial)
roc2 <- roc(bpca_roc$label, glm2$fitted.values)
mylist_roc[["heterodimerization"]] <- roc2
auc_labels <- c(auc_labels, paste0("heterodimerization (AUC = ", round(roc2$auc, 2), ")"))

auc_values <- sapply(mylist_roc, function(x) x$auc)
scores_ordered <- names(sort(auc_values, decreasing = TRUE))
auc_labels_ordered <- paste0(scores_ordered, " (AUC = ", round(sort(auc_values, decreasing = TRUE), 2), ")")

roc.test(roc1, roc2, method = "delong")

g_roc_ctrl_4 <- ggroc(mylist_roc[scores_ordered], legacy.axes = TRUE, size = 4) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 30),
        legend.position = c(0.63, 0.2),  
        legend.key.height = unit(2, "cm"),
        legend.background = element_rect(
          fill = alpha("white", 0.8),  
          colour = "black",            
          linewidth = 0.5              
        ),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        title = element_text(size = 20)) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  geom_abline(linetype = "dashed", size = 1) +
  scale_color_manual(
    values = c("heterodimerization" = "#E9AF1D", "abundance" = "#0A1829"),
    labels = auc_labels_ordered)

g_roc_ctrl_4
ggsave(g_roc_ctrl_4, file = "roc_ctrl_4.tiff", width = 12, height = 10)


#caclulate AUCs for combined scores
weights <- seq(0, 1, by = 0.05)

results <- map_dfr(weights, function(w) {
  
  score <- w * control_4$binding_score +
    (1 - w) * control_4$abundance_score
  
  roc_obj <- roc(control_4$label, score, quiet = TRUE)
  
  tibble(
    weight_binding = w,
    AUC = as.numeric(auc(roc_obj))
  )
})

p_binding_weights <- ggplot(results, aes(x = weight_binding, y = AUC)) +
  geom_line(size = 2, color = "red") +
  theme_classic() +
  labs(x = "Heterodimerization weight", title = "Combined score = w x HD score + (1-w) x abundance score") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))

p_binding_weights
ggsave(p_binding_weights, file = "binding_weights_auc.tiff", width = 11, height = 9)


#BS3_Supporting/PS3_Moderate evidence
controls_ddpca <- controls_ddpca %>% 
  mutate(evidence_apca = case_when(
    category_fdr_apca %in% c("low", "Stop-like") ~ "PS3_Moderate",
    category_fdr_apca == "WT-like" ~ "BS3_Supporting",
    category_fdr_apca == "high" ~ "Indeterminate")) %>% 
  mutate(evidence_bpca = case_when(
    category_fdr_bpca %in% c("low", "Stop-like") ~ "PS3_Moderate",
    category_fdr_bpca == "WT-like" ~ "BS3_Supporting",
    category_fdr_bpca == "high" ~ "Indeterminate"
  ))

controls_ddpca_evidence <- controls_ddpca %>%
  dplyr::count(evidence_apca, evidence_bpca) %>%
  mutate(perc = n / sum(n) * 100)


p_bar_evidence_apca <- ggplot(controls_ddpca_evidence, aes(x = evidence_apca, y = n)) +
  geom_col() +
  coord_flip()

p_bar_evidence_apca

p_bar_evidence_bpca <- ggplot(controls_ddpca_evidence, aes(x = evidence_bpca, y = n)) +
  geom_col() +
  coord_flip()

p_bar_evidence_bpca


p_lluvial_evi <- ggplot(controls_ddpca_evidence,
                              aes(axis1 = evidence_apca,
                                  axis2 = evidence_bpca,
                                  y = n)) +
  geom_alluvium(aes(fill = evidence_apca), width = 1/12, curve_type = "sigmoid") +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black", size = 1) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
  scale_x_discrete(limits = c("Abundance", "Hetero-\ndimerization"), expand = c(.05, .05)) +
  scale_fill_manual(
    values = c(
      "BS3_Supporting" = "#469A4D",
      "Indeterminate"   = "darkgrey",
      "PS3_Moderate"   = "darkorange")) +
  labs(y = "Number of mutations") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.position = "none")

p_lluvial_evi

ggsave(p_lluvial_evi, file = "alluvial_evidence_allmut.tiff", width = 12, height = 8)


clinvar_evidence <- inner_join(clinvar, controls_ddpca, by = "ID") %>%
  filter(Classification == "Uncertain significance") %>%
  mutate(Classification = "VUS")

clinvar_evidence <- clinvar_evidence %>%
  dplyr::count(Classification, evidence_apca, evidence_bpca, name = "n")

clinvar_evidence <- clinvar_evidence %>%
  mutate(
    Classification = factor(Classification, levels = "VUS"),
    evidence_apca = factor(evidence_apca),
    evidence_bpca = factor(evidence_bpca)
  )

clinvar_evidence$evidence_apca <- dplyr::recode(
  clinvar_evidence$evidence_apca,
  "PS3_Moderate" = "PS3\nModerate",
  "BS3_Supporting" = "BS3\nSupporting"
)

clinvar_evidence$evidence_bpca <- dplyr::recode(
  clinvar_evidence$evidence_bpca,
  "PS3_Moderate" = "PS3\nModerate",
  "BS3_Supporting" = "BS3\nSupporting"
)


p_lluvial_vus <- ggplot(clinvar_evidence,
                        aes(axis1 = Classification,
                            axis2 = evidence_apca,
                            axis3 = evidence_bpca,
                            y = n)) +
  geom_alluvium(aes(fill = evidence_apca), width = 1/12, curve_type = "sigmoid", alpha = 0.6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black", size = 1) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
  scale_x_discrete(limits = c("VUS", "Abundance", "Hetero-\ndimerization"),
                   labels = c("VUS" = "", "Abundance" = "Abundance", "Hetero-\ndimerization" = "Hetero-\ndimerization"),expand = c(.05, .05)) +
  scale_fill_manual(
    values = c("BS3\nSupporting" = "#469A4D","Indeterminate"   = "darkgrey","PS3\nModerate"   = "darkorange","VUS" = "#95BDED")) +
  labs(y = "Number of mutations") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 25),
        axis.title.y = element_text(size = 26),
        legend.position = "none")

p_lluvial_vus
ggsave(p_lluvial_vus, file = "alluvial_evidence_VUS.tiff", width = 13, height = 5)

clinvar_evidence_conf <- inner_join(clinvar, controls_ddpca, by = "ID") %>%
  filter(Classification == "Conflicting classifications of pathogenicity") %>%
  mutate(Classification = "Conf. classif.")

clinvar_evidence_conf <- clinvar_evidence_conf %>%
  dplyr::count(Classification, evidence_apca, evidence_bpca, name = "n")


clinvar_evidence_conf <- clinvar_evidence_conf %>%
  mutate(
    Classification = factor(Classification, levels = "Conf. classif."),
    evidence_apca = factor(evidence_apca),
    evidence_bpca = factor(evidence_bpca)
  )



p_lluvial_vus_conf <- ggplot(clinvar_evidence_conf,
                        aes(axis1 = Classification,
                            axis2 = evidence_apca,
                            axis3 = evidence_bpca,
                            y = n)) +
  geom_alluvium(aes(fill = evidence_apca), width = 1/12, curve_type = "sigmoid", alpha = 0.6) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black", size = 1) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
  scale_x_discrete(limits = c("Conf. classif.", "Abundance", "Hetero-\ndimerization"),
                   labels = c("Conf. classif." = "", "Abundance" = "Abundance", "Hetero-\ndimerization" = "Hetero-\ndimerization"),expand = c(.05, .05)) +
  scale_fill_manual(
    values = c(
      "BS3_Supporting" = "#469A4D",
      "Indeterminate"   = "darkgrey",
      "PS3_Moderate"   = "darkorange",
      "Conf. classif." = "darkblue")) +
  labs(y = "Number of mutations") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.position = "none")

p_lluvial_vus_conf
ggsave(p_lluvial_vus_conf,file = "alluvial_evidence_conf.tiff", width = 17, height = 5)