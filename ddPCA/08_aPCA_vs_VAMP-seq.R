#This script compares the SOD1 aPCA data with the SOD1 VAMP-seq data (Axakova et. al, AMJHG 2025)
library(ggplot2)
library(dplyr)

#load VAMP-seq and aPCA data
vampseq <- read.csv("merged_map_unfloored.csv")
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")

#adjust IDs to work in the same format
vampseq$ID <- sapply(vampseq$hgvs_pro, convert_hgvsp_to_id)
vampseq$ID <- sapply(vampseq$ID, decrement_id_pos)
vampseq$vampseq_score <- vampseq$score

both_scores <- inner_join(SOD1_subsitutions, vampseq, by = "ID")

#see scatter plot
cor_test <-cor.test(both_scores$abundance_score, both_scores$vampseq_score, method = "spearman")

cor_coef <- round(cor_test$estimate, 2) 
p_value <- format(cor_test$p.value, digits = 3, scientific = TRUE)

p_cor_both <- ggplot(both_scores, aes(x = abundance_score, y = vampseq_score)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  theme_classic() +
  labs(x = "aPCA score", y = "VAMP-seq score\n(Axakova et al.)") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22)) +
  annotate("text", x = -1.2, y = 2, label = paste("Spearman’s ρ =", cor_coef), size = 10, hjust = 0)
p_cor_both

ggsave(p_cor_both, file = "correlation_vampseq_apca.tiff", width = 10, height = 8)

#matrix scatter plot by secondary structure
both_scores$s_location <- factor(both_scores$s_location,
                                 levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                            "β4", "Zn binding loop", "β5", "loop 5", "β6", "loop 6",
                                            "β7", "Electrostatic loop", "β8", "loop 8", "Zn binding residues"))

cor_sec_struc <- both_scores %>%
  group_by(s_location) %>%
  summarise(
    r = cor(abundance_score, vampseq_score, method = "spearman"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Spearman’s ρ = ", round(r, 2))
  )


p_cor_sec_struc <- ggplot(both_scores, aes(x = abundance_score, y = vampseq_score)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  theme_classic() +
  facet_wrap(~s_location) +
  labs(x = "aPCA score", y = "VAMP-seq score\n(Axakova et al.)") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        strip.text = element_text(size = 15)) +
  geom_text(
    data = cor_sec_struc,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.1,
    size = 4
  )
p_cor_sec_struc

ggsave(p_cor_sec_struc, file = "correlation_vampseq_apca_sec_struc.tiff", width = 14, height = 10)


p_spearman_sec <- ggplot(cor_sec_struc, aes(x = s_location, y = r)) +
  geom_col(fill = "grey", color = "black") +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  theme_classic() +
  labs(x = "", y = "Spearman’s ρ\n(aPCA vs VAMP-seq)") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_spearman_sec

#scatter plot by beta strands and loops
cor_loc <- both_scores %>%
  group_by(location) %>%
  summarise(
    r = cor(abundance_score, vampseq_score, method = "spearman"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "Spearman’s ρ =", round(r, 2)
    )
  )

p_cor_loc <- ggplot(both_scores, aes(x = abundance_score, y = vampseq_score)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  theme_classic() +
  facet_wrap(~location) +
  labs(x = "aPCA score", y = "VAMP-seq score\n(Axakova et al.)") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22)) +
  geom_text(
    data = cor_loc,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.1,
    size = 6
  )
p_cor_loc

p_spearman_loc <- ggplot(cor_loc, aes(x = location, y = r)) +
  geom_col(fill = "grey", color = "black", width = 0.5) +
  theme_classic() +
  labs(x = "", y = "Spearman’s ρ\n(aPCA vs VAMP-seq)") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22))
p_spearman_loc

p_spearman_comb <- ggarrange(p_spearman_sec,
                             p_spearman_loc,
                             ncol = 2, align = "h")
p_spearman_comb

ggsave(p_spearman_comb, file = "spearman_corr_vampseq_apca.tiff", width = 18, height = 8)


assays <- c("apca", "vampseq")
mutations <- c(5980, 2906)

comparison_mut <- data.frame(
  assay = assays,
  n_mut = mutations
)

p_comp_mut <- ggplot(comparison_mut, aes(x = assay, y = n_mut)) +
  geom_col(width = 0.5, color = "black") +
  theme_classic() +
  scale_x_discrete(labels = c("apca" = "Abundance PCA", "vampseq" = "VAMP-seq\n(Axakova et al.)")) +
  labs(x="", y = "Number of mutations") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22))
p_comp_mut

ggsave(p_comp_mut, file = "nmut_vampseq_apca.tiff", width = 8, height = 6)


#VAMP-SEQ vs in vitro data####
flim <- read.csv("FLIM.csv")

flim <- rename(flim, ID = treatment, flim_score = mean)

#normalize flim
wt_flim = 1.78
flim <- flim %>%
  mutate(flim_score = flim_score - wt_flim) %>% 
  mutate(flim_score = -flim_score)

SOD1_flim_vampseq <- inner_join(flim, vampseq, by = "ID")

cor_flim_vampseq <-cor.test(SOD1_flim_vampseq$vampseq_score, SOD1_flim_vampseq$flim_score, method = "pearson")

coef_flim_vampseq <- round(cor_flim_vampseq$estimate, 2) 
p_value <- format(cor_flim_vampseq$p.value, digits = 3, scientific = TRUE)

p_corr_vampseq_flim <-ggplot(SOD1_flim_vampseq, aes(x = vampseq_score, y = flim_score, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_errorbar(aes(xmin = vampseq_score - se, xmax = vampseq_score + se), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "VAMP-seq score\n(Axakova et al.)",
       y = "FLIM score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = -0.03, label = paste("R =", coef_flim_vampseq), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = -0.06, label = paste("p =", p_value), size = 10, hjust = 0)


p_corr_vampseq_flim

ggsave(p_corr_vampseq_flim, file = "flim_vampseq_corr.tiff", width = 8, height = 6, dpi = 300)


invitro <- read.csv("invitro_data.csv", sep = ";")
invitro$mean_score <- rowMeans(invitro[, c("score_r1", "score_r2", "score_r3")], na.rm = TRUE)

survival <- invitro %>% 
  filter(assay == "survival") %>% 
  rename(survival_score = mean_score)

inclusions <- invitro %>% 
  filter(assay == "inclusions") %>% 
  rename(inclusions_score = mean_score)

intensity <- invitro %>% 
  filter(assay == "intensity") %>% 
  rename(intensity_score = mean_score)



invitro <- inner_join(survival, inclusions, by = "ID")
invitro <- inner_join(invitro, intensity, by = "ID")
invitro <- left_join(invitro, flim, by = "ID")

invitro <- select(invitro, ID, survival_score, inclusions_score, intensity_score, flim_score)
vampseq_inv <- inner_join(invitro, vampseq, by = "ID")
vampseq_inv <- select(vampseq_inv, ID, vampseq_score, survival_score, inclusions_score, intensity_score, flim_score)

vampseq_inv_norm <- vampseq_inv %>% 
  mutate(survival_score = survival_score -1,
         inclusions_score = inclusions_score -1,
         intensity_score = intensity_score -1)

vampseq_inv_norm <- vampseq_inv_norm %>%
  mutate(
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
  select(ID,vampseq_score, intensity_norm, survival_norm, inclusions_norm, flim_norm)



#abundance vs aggregation
cor_ab_aggre <-cor.test(vampseq_inv_norm$vampseq_score, vampseq_inv_norm$inclusions_norm, method = "pearson")

cor_coef_1 <- round(cor_ab_aggre$estimate, 2) 
p_value_1 <- format(cor_ab_aggre$p.value, digits = 3, scientific = TRUE)

p_corr_vampseq_aggre <-ggplot(vampseq_inv_norm, aes(x = vampseq_score, y = inclusions_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "VAMP-seq score\n(Axakova et al.)",
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
  annotate("text", x = 0, y = 1, label = paste("R =", cor_coef_1), size = 10, hjust = 0) +
  annotate("text", x = 0, y = 0.92, label = paste("p =", p_value_1), size = 10, hjust = 0)
p_corr_vampseq_aggre

ggsave(p_corr_vampseq_aggre, file = "corr_aggregation_vampseq.tiff", width = 8, height = 6)

#abundance vs survival
cor_ab_surv <-cor.test(vampseq_inv_norm$vampseq_score, vampseq_inv_norm$survival_norm, method = "pearson")

cor_coef_2 <- round(cor_ab_surv$estimate, 2) 
p_value_2 <- format(cor_ab_surv$p.value, digits = 3, scientific = TRUE)

p_corr_vampseq_surv <-ggplot(vampseq_inv_norm, aes(x = vampseq_score, y = survival_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "VAMP-seq score\n(Axakova et al.)",
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
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_2), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_2), size = 10, hjust = 0)
p_corr_vampseq_surv

ggsave(p_corr_vampseq_surv, file = "corr_survival_vampseq.tiff", width = 8, height = 6)


#abundance vs intensity
cor_ab_int <-cor.test(vampseq_inv_norm$vampseq_score, vampseq_inv_norm$intensity_norm, method = "pearson")

cor_coef_3 <- round(cor_ab_int$estimate, 2) 
p_value_3 <- format(cor_ab_int$p.value, digits = 3, scientific = TRUE)

p_corr_vampseq_int<-ggplot(vampseq_inv_norm, aes(x = vampseq_score, y = intensity_norm, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "VAMP-seq score\n(Axakova et al.)",
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
  annotate("text", x = -1.3, y = 1, label = paste("R =", cor_coef_3), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.82, label = paste("p =", p_value_3), size = 10, hjust = 0)
p_corr_vampseq_int

ggsave(p_corr_vampseq_int, file = "corr_intensity_vampseq.tiff", width = 8, height = 6)

assays <- c( "flim","aggregation", "survival", "intensity")
coeficients <- c(0.79, -0.75, 0.57, 0.82)

pearson_coeff <- data.frame(
  assay = assays,
  coefficient = coeficients
)

pearson_coeff$assay <- factor(pearson_coeff$assay, levels = c("flim", "aggregation", "survival", "intensity"))

p_coef <- ggplot(pearson_coeff, aes(x = assay, y = coefficient)) +
  geom_col(fill = "lightgrey", color = "black", width = 0.3, linewidth = 0.8) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Pearson's R") +
  scale_x_discrete(labels = c("flim" = "FLIM", "aggregation" = "Aggregation", "survival" = "Cell\nsurvival", "intensity" = "Intensity")) +
  scale_y_continuous(breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 22))
p_coef

ggsave(p_coef, file = "pearson_coeff_invitro.tiff", width = 8, height = 6)


#Correlation of validation variants (aPCA vs VAMP-seq)####
vampseq_apca_invitro <- inner_join(vampseq_inv_norm, SOD1_subsitutions, by = "ID")

cor_ab_vampseq <-cor.test(vampseq_apca_invitro$abundance_score, vampseq_apca_invitro$vampseq_score, method = "pearson")

cor_coef_4 <- round(cor_ab_vampseq$estimate, 2) 
p_value_4 <- format(cor_ab_vampseq$p.value, digits = 3, scientific = TRUE)

p_corr_vampseq_apca<-ggplot(vampseq_apca_invitro, aes(x = abundance_score, y = vampseq_score, label = ID)) +
  geom_point(size = 5, shape=16, alpha = 1, aes(color = category_fdr_apca)) +
  #geom_errorbar(aes(ymin = flim_score - sem, ymax = flim_score + sem), width = 0.02, alpha = 0.7, color = "darkgrey") +
  #geom_errorbar(aes(xmin = abundance_score - abundance_sigma, xmax = abundance_score + abundance_sigma), width = 0.02, alpha = 0.7, color = "darkgrey") +
  scale_color_manual(values = c("high abundance" = "darkgreen", "WT-lke abundance" = "lightgrey", "low abundance" = "darkorange", "stop-like abundance" = "#DD4636")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(y = "VAMP-seq score\n(Axakova et al.)",
       x = "aPCA score")+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.key = element_rect(fill = NA, color = "black"),
        legend.title = element_text(size = 26),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.key.width = unit(2, "cm"),
        legend.position = "none") +
  geom_text_repel(size = 6) +
  annotate("text", x = -1.3, y = 0.9, label = paste("R =", cor_coef_4), size = 10, hjust = 0) +
  annotate("text", x = -1.3, y = 0.75, label = paste("p =", p_value_4), size = 10, hjust = 0)
p_corr_vampseq_apca


ggsave(p_corr_vampseq_apca, file = "corr_apca_vampseq_validation_mutants.tiff", width = 8, height = 6)

