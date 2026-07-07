library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)


flim <- read.csv("FLIM.csv")
SOD_final <- read.csv("SOD1_abundance_binding_scores.csv")

flim <- rename(flim, ID = treatment, flim_score = mean)

#normalize flim
# wt_flim = 1.78
# flim <- flim %>%
#   mutate(flim_score = flim_score - wt_flim) %>% 
#   mutate(flim_score = -flim_score)

SOD1_flim_dms <- inner_join(flim, SOD_final, by = "ID")

SOD1_flim_dms <- SOD1_flim_dms %>%
  mutate(category_fdr_apca = recode(category_fdr_apca,
                                    "high abundance" = "high",
                                    "WT-like abundance" = "WT-like",
                                    "low abundance" = "low",
                                    "stop-like abundance" = "stop-like"))


cor_test <-cor.test(SOD1_flim_dms$abundance_score, SOD1_flim_dms$flim_score, method = "pearson")

cor_coef <- round(cor_test$estimate, 2) 
p_value <- format(cor_test$p.value, digits = 3, scientific = TRUE)




p_corr_apca_flim <-ggplot(SOD1_flim_dms, aes(x = abundance_score, y = flim_score, label = ID, fill = category_fdr_apca)) +
  geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
  scale_fill_manual(values = c("high" = "darkgreen", "WT-like" = "darkgrey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 1.78, linetype = "dashed", size = 1) + #wt flim score
  theme_classic() +
    labs(x = "Abundance score\n(yeast)",
    y = "AcGFP lifetime (ns)\n(HEK293T)")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 2.4, label = paste("R =", cor_coef), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 2.35, label = paste("p =", p_value), size = 10, hjust = 0)


p_corr_apca_flim

ggsave(p_corr_apca_flim, file = "flim_apca_corr.tiff", width = 9, height = 7, dpi = 300)

# cor_test_bpca <-cor.test(SOD1_flim_dms$binding_score, SOD1_flim_dms$flim_score, method = "pearson")
# 
# cor_coef_bpca <- round(cor_test_bpca$estimate, 2) 
# p_value_bpca <- format(cor_test_bpca$p.value, digits = 3, scientific = TRUE)
# 
# p_corr_bpca_flim <-ggplot(SOD1_flim_dms, aes(x = binding_score, y = flim_score, label = ID, fill = category_fdr_bpca)) +
#   geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
#   scale_fill_manual(values = c("high heterodimerization" = "darkgreen", "WT-like heterodimerization" = "darkgrey", "low heterodimerization" = "darkorange")) +
#   geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
#   geom_errorbar(aes(xmin = binding_score - binding_sigma, xmax = binding_score + binding_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
#   geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
#   geom_hline(yintercept = 1.78, linetype = "dashed", size = 1) +
#   theme_classic() +
#   labs(x = "Heterodimerization score\n(yeast)",
#        y = "AcGFP lifetime (ns)\n(HEK293T)")+
#   theme(axis.text = element_text(size = 24),
#         axis.title = element_text(size = 28),
#         legend.key = element_rect(fill = NA, color = "black"),
#         legend.title = element_text(size = 26),
#         legend.title.align = 0.5,
#         legend.text = element_text(size = 20),
#         legend.key.width = unit(2, "cm"),
#         legend.position = "none") +
#   geom_text_repel(size = 6) +
#   annotate("text", x = -1, y = 2.4, label = paste("R =", cor_coef_bpca), size = 10, hjust = 0) +
#   annotate("text", x = -1, y = 2.35, label = paste("p =", p_value_bpca), size = 10, hjust = 0)
# 
# 
# p_corr_bpca_flim

#In vitro validation####
invitro <- read.csv("invitro_data.csv", sep = ";")
invitro$mean_score <- rowMeans(invitro[, c("score_r1", "score_r2", "score_r3")], na.rm = TRUE)

survival <- invitro %>% 
  filter(assay == "survival") %>% 
  rename(survival_score = mean_score)

survival$ID <- factor(survival$ID, levels = unique(survival$ID))
survival <- survival %>%
  rowwise() %>%
  mutate(
    sd_surv = sd(c(score_r1, score_r2, score_r3), na.rm = TRUE),
    sem_surv = sd_surv / sqrt(3)
  )

survival_long <- survival %>%
  pivot_longer(
    cols = c(score_r1, score_r2, score_r3),
    names_to = "replicate",
    values_to = "score"
  )

p_survival <- ggplot(survival_long, aes(x = ID, y = score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +
  geom_point(position = position_jitter(width = 0.1),
             size = 3, color = "black") +
  theme_classic() +
  labs(x = "", y = "Survival\nscore") +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))


p_survival


inclusions <- invitro %>% 
  filter(assay == "inclusions") %>% 
  rename(inclusions_score = mean_score)

inclusions$ID <- factor(inclusions$ID, levels = unique(inclusions$ID))
inclusions <- inclusions %>%
  rowwise() %>%
  mutate(
    sd_inc = sd(c(score_r1, score_r2, score_r3), na.rm = TRUE),
    sem_inc = sd_inc / sqrt(3)
  )

inclusions_long <- inclusions %>%
  pivot_longer(
    cols = c(score_r1, score_r2, score_r3),
    names_to = "replicate",
    values_to = "score"
  )

p_inclusions <- ggplot(inclusions_long, aes(x = ID, y = score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +
  geom_point(position = position_jitter(width = 0.1),
             size = 3, color = "black") +
  theme_classic() +
  labs(x = "", y = "Aggregation\nscore") +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))


p_inclusions

intensity <- invitro %>% 
  filter(assay == "intensity") %>% 
  rename(intensity_score = mean_score)

intensity$ID <- factor(intensity$ID, levels = unique(intensity$ID))
intensity <- intensity %>%
  rowwise() %>%
  mutate(
    sd_int = sd(c(score_r1, score_r2, score_r3), na.rm = TRUE),
    sem_int = sd_int / sqrt(3)
  )

intensity_long <- intensity %>%
  pivot_longer(
    cols = c(score_r1, score_r2, score_r3),
    names_to = "replicate",
    values_to = "score"
  )

p_intensity <- ggplot(intensity_long, aes(x = ID, y = score)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +
  geom_point(position = position_jitter(width = 0.1),
             size = 3, color = "black") +
  theme_classic() +
  labs(x = "", y = "Intensity\nscore") +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

p_intensity

p_all_phenotypes <- ggarrange(p_inclusions,
                              p_survival,
                              p_intensity,
                              nrow = 3, align = "hv")

p_all_phenotypes
ggsave(p_all_phenotypes, file = "val_phen_individual_measurementes.tiff", width = 20, height = 12)


survival <- survival %>%
  rowwise() %>%
  mutate(
    sd_surv = sd(c(score_r1, score_r2, score_r3), na.rm = TRUE),
    sem_surv = sd_surv / sqrt(3)
  ) %>%
  ungroup()

inclusions <- inclusions %>%
  rowwise() %>%
  mutate(
    sd_inc = sd(c(score_r1, score_r2, score_r3), na.rm = TRUE),
    sem_inc = sd_inc / sqrt(3)
  ) %>%
  ungroup()

intensity <- intensity %>%
  rowwise() %>%
  mutate(
    sd_int = sd(c(score_r1, score_r2, score_r3), na.rm = TRUE),
    sem_int = sd_int / sqrt(3)
  ) %>%
  ungroup()


invitro <- inner_join(survival, inclusions, by = "ID")
invitro <- inner_join(invitro, intensity, by = "ID")
invitro <- left_join(invitro, flim, by = "ID")

invitro <- select(invitro, ID, survival_score, sem_surv, inclusions_score, sem_inc, intensity_score, sem_int, flim_score)


dms_inv <- inner_join(invitro, SOD_final, by = "ID")
dms_inv <- select(dms_inv, ID, abundance_score, binding_score, survival_score, sem_surv, inclusions_score, sem_inc, intensity_score, sem_int, flim_score, category_fdr_apca, category_fdr_bpca)



dms_inv_norm <- dms_inv %>% 
  mutate(survival_score = survival_score -1,
         inclusions_score = inclusions_score -1,
         intensity_score = intensity_score -1,
         flim_score = flim_score - 1.78)

dms_inv_norm <- dms_inv_norm %>%
  mutate(
    abundance_norm = abundance_score/max(abs(c(min(abundance_score, na.rm = TRUE),
                                             max(abundance_score, na.rm = TRUE))),
                                       na.rm = TRUE),
    survival_norm = survival_score/max(abs(c(min(survival_score, na.rm = TRUE),
                                             max(survival_score, na.rm = TRUE))),
                                       na.rm = TRUE),
    inclusions_norm = inclusions_score/max(abs(c(min(inclusions_score, na.rm = TRUE),
                                                 max(inclusions_score, na.rm = TRUE))),
                                           na.rm = TRUE),
    intensity_norm = intensity_score/max(abs(c(min(intensity_score, na.rm = TRUE),
                                               max(intensity_score, na.rm = TRUE))),
                                         na.rm = TRUE),
    flim_norm = flim_score/max(abs(c(min(flim_score, na.rm = TRUE),
                                     max(flim_score, na.rm = TRUE))),
                               na.rm = TRUE)
  ) %>% 
  select(ID,abundance_norm, intensity_norm, survival_norm, inclusions_norm, flim_norm, category_fdr_apca)


dms_inv_long <- dms_inv_norm %>%
  pivot_longer(
    cols = -c(ID, category_fdr_apca),
    names_to = "assay",
    values_to = "score"
  ) %>%
  mutate(Pos = as.numeric(str_extract(ID, "\\d+")))

dms_inv_long <- dms_inv_long %>%
  arrange(Pos)
dms_inv_long$ID <- factor(dms_inv_long$ID, levels = unique(dms_inv_long$ID))
dms_inv_long$assay <- factor(dms_inv_long$assay, levels = c("abundance_norm", "flim_norm", "inclusions_norm", "survival_norm", "intensity_norm"))


p_heatmap <- ggplot(dms_inv_long, aes(x = ID, y = assay, fill = score)) +
  geom_tile(color = "black", size = 0.8) +
  scale_y_discrete(labels = c("abundance_norm" = "Abundance", "flim_norm" = "AcGFP lifetime", "inclusions_norm" = "Aggregation", 
                              "survival_norm" = "Survival", "intensity_norm" = "Intensity")) +
  labs(x = "", y = "", fill = "Score") +
  scale_fill_gradient2(low = "#AD9024", mid = "#F4F8FB", high = "#2A5783", breaks = c(-1, 0, 1), na.value = "darkgrey") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 25, hjust = 1, vjust = 0.5, angle = 90),
        axis.text.y = element_text(size = 22),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 35, hjust = 0.5),
        legend.position = "bottom") +
  guides(fill = guide_colorbar(
    barwidth = 14, 
    barheight = 1.7,
    title.position = "top",
    title.hjust = 0.5
  )) + 
  coord_flip()
p_heatmap


ggsave(p_heatmap, file ="invitro_heatmap.tiff", width = 4, height = 19)

agg_ab_surv <- dms_inv_long %>% 
  filter(assay %in% c("abundance_norm", "inclusions_norm", "survival_norm"))

#correlation plots####
#abundance vs aggregation with % cells (not normalized)
cor_ab_aggre <-cor.test(dms_inv$abundance_score, dms_inv$inclusions_score, method = "pearson")

cor_coef_1 <- round(cor_ab_aggre$estimate, 2) 
p_value_1 <- format(cor_ab_aggre$p.value, digits = 3, scientific = TRUE)
conf_int_1 <- cor_ab_aggre$conf.int

p_corr_apca_aggre <-ggplot(dms_inv_norm, aes(x = abundance_score, y = inclusions_score, label = ID, fill = category_fdr_apca)) +
  geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
  scale_fill_manual(values = c("high abundance" = "darkgreen", "WT-like abundance" = "lightgrey", "low abundance" = "darkorange", "stop-like abundance" = "#DD4636")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score\n(yeast)",
       y = "% of cells with aggregates\n(HEK293T)")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6, max.overlaps = 6) +
  annotate("text", x = -1.3, y = 20, label = paste("R =", cor_coef_1), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 18, label = paste("p =", p_value_1), size = 10, hjust = 0)
p_corr_apca_aggre

ggsave(p_corr_apca_aggre, file = "corr_aggregation_apca.tiff", width = 9, height = 7)
#abundance vs aggregation
cor_ab_aggre <-cor.test(dms_inv_norm$abundance_norm, dms_inv_norm$inclusions_norm, method = "pearson")

cor_coef_1 <- round(cor_ab_aggre$estimate, 2) 
p_value_1 <- format(cor_ab_aggre$p.value, digits = 3, scientific = TRUE)
conf_int_1 <- cor_ab_aggre$conf.int

p_corr_apca_aggre <-ggplot(dms_inv_norm, aes(x = abundance_norm, y = inclusions_norm, label = ID, fill = category_fdr_apca)) +
  geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
  scale_fill_manual(values = c("high abundance" = "darkgreen", "WT-like abundance" = "lightgrey", "low abundance" = "darkorange", "stop-like abundance" = "#DD4636")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score\n(yeast)",
       y = "Aggregation score\n(HEK293T)")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6, max.overlaps = 6) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_1), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.92, label = paste("p =", p_value_1), size = 10, hjust = 0)
p_corr_apca_aggre

ggsave(p_corr_apca_aggre, file = "corr_aggregation_apca.tiff", width = 9, height = 7)


#abundance vs survival
cor_ab_surv <-cor.test(dms_inv_norm$abundance_norm, dms_inv_norm$survival_norm, method = "pearson")

cor_coef_2 <- round(cor_ab_surv$estimate, 2) 
p_value_2 <- format(cor_ab_surv$p.value, digits = 3, scientific = TRUE)

p_corr_apca_surv <-ggplot(dms_inv_norm, aes(x = abundance_norm, y = survival_norm, label = ID, fill = category_fdr_apca)) +
  geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
  scale_fill_manual(values = c("high abundance" = "darkgreen", "WT-like abundance" = "darkgrey", "low abundance" = "darkorange", "stop-like abundance" = "#DD4636")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score\n(yeast)",
       y = "Survival score\n(NSC-34)")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6, max.overlaps = 7) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_2), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_2), size = 10, hjust = 0)
p_corr_apca_surv

ggsave(p_corr_apca_surv, file = "corr_survival_apca.tiff", width = 9, height = 7)


#abundance vs intensity
cor_ab_int <-cor.test(dms_inv_norm$abundance_norm, dms_inv_norm$intensity_norm, method = "pearson")

cor_coef_3 <- round(cor_ab_int$estimate, 2) 
p_value_3 <- format(cor_ab_int$p.value, digits = 3, scientific = TRUE)

p_corr_apca_int<-ggplot(dms_inv_norm, aes(x = abundance_norm, y = intensity_norm, label = ID, fill = category_fdr_apca)) +
  geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
  scale_fill_manual(values = c("high abundance" = "darkgreen", "WT-like abundance" = "darkgrey", "low abundance" = "darkorange", "stop-like abundance" = "#DD4636")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score\n(yeast)",
       y = "Intensity score\n(NSC-34)")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6, max.overlaps = 5) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_3), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_3), size = 10, hjust = 0)
p_corr_apca_int

ggsave(p_corr_apca_int, file = "corr_intensity_apca.tiff", width = 9, height = 7)


#Corelation between mammalian cells assays####
##FLIM vs inclusions####
cor_flim_aggre <-cor.test(dms_inv_norm$flim_norm, dms_inv_norm$inclusions_norm, method = "pearson")

cor_coef_4 <- round(cor_flim_aggre$estimate, 2) 
p_value_4 <- format(cor_flim_aggre$p.value, digits = 3, scientific = TRUE)

p_corr_flim_aggre <-ggplot(dms_inv_norm, aes(x = flim_norm, y = inclusions_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "FLIM score",
       y = "Aggregation score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_4), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_4), size = 10, hjust = 0)
p_corr_flim_aggre

ggsave(p_corr_flim_aggre, file = "corr_flim_aggre.tiff", width = 8, height = 6)

##FLIM vs survival####
cor_flim_surv <-cor.test(dms_inv_norm$flim_norm, dms_inv_norm$survival_norm, method = "pearson")

cor_coef_5 <- round(cor_flim_surv$estimate, 2) 
p_value_5 <- format(cor_flim_surv$p.value, digits = 3, scientific = TRUE)

p_corr_flim_surv <-ggplot(dms_inv_norm, aes(x = flim_norm, y = survival_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "FLIM score",
       y = "Survival score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_5), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_5), size = 10, hjust = 0)
p_corr_flim_surv

ggsave(p_corr_flim_surv, file = "corr_flim_surv.tiff", width = 8, height = 6)


##FLIM vs intensity####
cor_flim_int <-cor.test(dms_inv_norm$flim_norm, dms_inv_norm$intensity_norm, method = "pearson")

cor_coef_6 <- round(cor_flim_int$estimate, 2) 
p_value_6 <- format(cor_flim_int$p.value, digits = 3, scientific = TRUE)

p_corr_flim_int <-ggplot(dms_inv_norm, aes(x = flim_norm, y = intensity_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "FLIM score",
       y = "Intensity score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_6), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_6), size = 10, hjust = 0)
p_corr_flim_int

ggsave(p_corr_flim_int, file = "corr_flim_int.tiff", width = 8, height = 6)



##Inclusions vs survival####
cor_aggre_surv <-cor.test(dms_inv_norm$inclusions_norm, dms_inv_norm$survival_norm, method = "pearson")

cor_coef_7 <- round(cor_aggre_surv$estimate, 2) 
p_value_7 <- format(cor_aggre_surv$p.value, digits = 3, scientific = TRUE)
conf_int_7 <- cor_aggre_surv$conf.int

p_corr_aggre_surv <-ggplot(dms_inv_norm, aes(x = inclusions_norm, y = survival_norm, label = ID, fill = category_fdr_apca)) +
  geom_point(size = 5, alpha = 1, color = "black", stroke = 1, shape = 21) +
  scale_fill_manual(values = c("high abundance" = "darkgreen", "WT-like abundance" = "darkgrey", "low abundance" = "darkorange", "Stop-like abundance" = "#DD4636")) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Aggregation score",
       y = "Survival score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 4) +
  annotate("text", x = 0.6, y = 0.15, label = paste("R =", cor_coef_7), size = 10, hjust = 0) +
  annotate("text", x = 0.6, y = 0.05, label = paste("p =", p_value_7), size = 10, hjust = 0)
p_corr_aggre_surv

ggsave(p_corr_aggre_surv, file = "corr_aggre_surv.tiff", width = 8, height = 6)

##Intensity vs survival####
cor_int_surv <-cor.test(dms_inv_norm$intensity_norm, dms_inv_norm$survival_norm, method = "pearson")

cor_coef_8 <- round(cor_int_surv$estimate, 2) 
p_value_8 <- format(cor_int_surv$p.value, digits = 3, scientific = TRUE)

p_corr_int_surv <-ggplot(dms_inv_norm, aes(x = intensity_norm, y = survival_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Intensity score",
       y = "Survival score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_8), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_8), size = 10, hjust = 0)
p_corr_int_surv

ggsave(p_corr_int_surv, file = "corr_int_surv.tiff", width = 8, height = 6)

##Inclusions vs intensity####
cor_aggre_int <-cor.test(dms_inv_norm$inclusions_norm, dms_inv_norm$intensity_norm, method = "pearson")

cor_coef_9 <- round(cor_aggre_int$estimate, 2) 
p_value_9 <- format(cor_aggre_int$p.value, digits = 3, scientific = TRUE)

p_corr_aggre_int <-ggplot(dms_inv_norm, aes(x = inclusions_norm, y = intensity_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Aggregation score",
       y = "Intensity score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_9), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_9), size = 10, hjust = 0)
p_corr_aggre_int

ggsave(p_corr_aggre_int, file = "corr_aggre_int.tiff", width = 8, height = 6)

#Correlation matrix (mammalian cells data)####
mammalian_correlations <- data.frame(
  assay_x = c("FLIM", "FLIM", "FLIM", "FLIM", 
              "Aggregation", "Aggregation", "Aggregation",
              "Survival", "Survival",
              "Intensity"),
  assay_y = c("Abundance", "Aggregation", "Survival", "Intensity", 
              "Abundance", "Survival", "Intensity",
              "Abundance", "Intensity",
              "Abundance"),
  r = c(cor_coef, cor_coef_4, cor_coef_5, cor_coef_6,
        cor_coef_1, cor_coef_7, cor_coef_9,
        cor_coef_2, cor_coef_8,
        cor_coef_3)
)

p_pearson_heatmap <- ggplot(mammalian_correlations, aes(x = assay_x, y = assay_y, fill = r)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(r, 2)), size = 6) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  theme_minimal() +
  labs(x = "", y = "")
p_pearson_heatmap


df_full <- mammalian_correlations |>
  bind_rows(mammalian_correlations |> rename(assay_x = assay_y, assay_y = assay_x))

p_pearson_heatmap <- ggplot(df_full, aes(x = assay_x, y = assay_y, fill = r)) +
  geom_tile(color = "black", size = 0.5) +
  geom_text(aes(label = round(r, 2)), size = 8) +
  scale_fill_gradient2(
    low = "#A70C20", mid = "#F5F9FC", high = "#26456E",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson's R"
  ) +
  labs(x ="", y = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        legend.position = "top")

p_pearson_heatmap

ggsave(p_pearson_heatmap, file = "pearson_matrix_validarion.tiff", width = 10, height = 8)

#Visualize pearson's correlations####
vars <- c("flim_norm", "inclusions_norm", "survival_norm", "intensity_norm")

res <- lapply(vars, function(v) {
  test <- cor.test(
    dms_inv_norm$abundance_norm,
    dms_inv_norm[[v]],
    method = "pearson"
  )
  
  data.frame(
    assay = v,
    r = unname(test$estimate),
    lower = test$conf.int[1],
    upper = test$conf.int[2],
    p = test$p.value
  )
})

pearson_coeff <- do.call(rbind, res)

pearson_coeff$assay <- factor(
  pearson_coeff$assay,
  levels = c("flim_norm", "inclusions_norm", "survival_norm", "intensity_norm")
)

p_coef <- ggplot(pearson_coeff, aes(x = assay, y = r)) +
  geom_col(fill = "lightgrey", color = "black", width = 0.3, linewidth = 0.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Pearson's R") +
  scale_x_discrete(labels = c("flim_norm" = "FLIM", "inclusions_norm" = "Aggregation", "survival_norm" = "Survival", "intensity_norm" = "Intensity")) +
  scale_y_continuous(breaks = c(-0.8, -0.4, 0, 0.4, 0.8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22))
p_coef

ggsave(p_coef, file = "pearson_coeff_invitro.tiff", width = 8, height = 6)
