library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(ggalluvial)
library(ggnewscale)

#load modelled data
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")
SOD1_insertions <-read.csv("SOD1_ins_modelled.csv")
SOD1_deletions <- read.csv("SOD1_del_modelled.csv")

SOD1_allmut <- bind_rows(SOD1_subsitutions, SOD1_insertions, SOD1_deletions)
SOD1_allmut$abs_residual <- abs(SOD1_allmut$residuals)

#### Substiutions ####
#calculate the median residual by position of all interface residues
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
  )) %>% 
  mutate(outlier_interface = case_when(
    outlier == "TRUE" & interface =="TRUE" ~ TRUE,
    TRUE ~ FALSE
  ))



SOD1_allmut_filt <- na.omit(SOD1_allmut)

r_corr <- SOD1_allmut_filt %>%
  summarise(
    R = cor(abundance_score, binding_score, method = "pearson"),
    p_value = cor.test(abundance_score, binding_score, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p < 2.2e-16")
  )

SOD1_allmut <- SOD1_allmut %>% 
  mutate(outlier_type = case_when(
    outlier_interface == T ~ "outlier interface",
    outlier_nointerface == T ~ "outlier no interface",
    TRUE ~ "no outlier"
  ))
SOD1_allmut <- SOD1_allmut %>%
  dplyr::arrange(outlier_type == "outlier interface")

p_corr <- ggplot(data = SOD1_allmut, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 2.5, alpha = 0.8, aes(color = outlier_type)) +
  scale_color_manual(values = c("outlier interface" = "#2A5783", "outlier no interface" = "#B9DDF1", "no outlier" = "lightgrey")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Normalized\nabundance score", y = "Normalized\nheterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.position = "none") +
  geom_text(data = r_corr, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE) 
p_corr
ggsave(p_corr, path = path_fig4, file = "correlation_apca_boca_label_by_outlier.tiff", width = 10, height = 8)

SOD1_subsitutions <- SOD1_allmut %>% 
  filter(mutation_type == "subs")

SOD1_insertions <- SOD1_allmut %>% 
  filter(mutation_type == "ins")

SOD1_deletions<- SOD1_allmut %>% 
  filter(mutation_type == "del")

interface <- SOD1_allmut %>% 
  filter(interface == T)


SOD1_allmut$mutation_type <- factor(SOD1_allmut$mutation_type, 
                                    levels = c("subs", "ins", "del"))
prop_fimpact <- SOD1_allmut %>% 
  group_by(interface, mutation_type, residual_category) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(interface, mutation_type) %>%         
  mutate(prop = n / sum(n),
         perc = n / sum(n) * 100) %>% 
  ungroup()

prop_fimpact$interface <- factor(prop_fimpact$interface, 
                                 levels = c("TRUE", "FALSE"))

prop_fimpact$residual_category <- factor(prop_fimpact$residual_category, 
                                         levels = c("high_positive", "high_negative", "low"))
prop_fimpact <- na.omit(prop_fimpact)

p_prop_impact_specific <- ggplot(prop_fimpact, aes(x = mutation_type, y = perc, fill = residual_category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 1, width = 0.5) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  scale_x_discrete(labels = c("subs" = "Subtitutions", "ins" = "Insertions", "del" = "Deletions")) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            color = "black", size = 7) +
  facet_wrap(~interface, ncol = 1, labeller = labeller(interface = c("TRUE" = "Dimer interface", "FALSE" = "Non-dimer interface"))) +
  labs(y = "Percentage of mutations (%)", x = "") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 28),
    axis.title = element_text(size = 33),
    axis.text.x = element_text(size = 30, hjust = 0.5),
    axis.text.y = element_text(size = 28),
    strip.text = element_text(size = 30)) +
  guides(fill = guide_legend(override.aes = list(size = 14)))

p_prop_impact_specific

ggsave(p_prop_impact_specific, file = "perc_signif_residuals_allmut.tiff", width = 10, height = 8)


SOD1_subsitutions$interface <- factor(SOD1_subsitutions$interface,
                                        levels = c("TRUE", "FALSE"))

test_allsubs <- wilcox.test(abs_residual ~ interface, data = SOD1_subsitutions)

  
p_residual_dist <- ggplot(SOD1_subsitutions, aes(x = interface, y = abs_residual, fill = interface)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1")) +
  scale_x_discrete(labels = c("FALSE" = "Non-dimer interface", "TRUE" = "Dimer interface")) +
  theme_classic() +
  labs(x = "", y = "|Heterodimerization residual|", title = "Substitutions") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "none") +
  annotate(
    "text",
    x = 1.5,
    y = max(SOD1_subsitutions$abs_residual) * 1.1,
    label = paste0("p = ", signif(test_allsubs$p.value, 3)),
    size = 6
  )
p_residual_dist

test_allsubs <- wilcox.test(abs_residual ~ interface, data = SOD1_subsitutions)


p_residual_density_subs <- ggplot(SOD1_subsitutions, aes(x = abs_residual, fill = interface)) +
  geom_density(alpha = 0.9) +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1"),
                    labels = c("TRUE" = "Dimer interface", "FALSE" = "Non-dimer interface")) +
  labs(x = "|Heterodimerization residual|", y = "Density", title = "Substitutions", fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
       axis.title = element_text(size = 24),
       plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
       legend.text = element_text(size = 18),
       legend.position = c(0.8, 0.8)) +
  annotate(
    "text",
    x = 0.2,
    y =15,
    label = paste0("p = ", signif(test_allsubs$p.value, 3)),
    size = 6
  )
  
p_residual_density_subs

p_residual_subs_mut <- ggarrange(p_residual_dist,
                                 p_residual_density_subs,
                                 nrow =  2, align = "v")

p_residual_subs_mut

ggsave(p_residual_subs_mut, file = "residuals_subs_by_mut.tiff", width = 10, height = 10)

SOD1_median_res <- SOD1_subsitutions %>% 
  group_by(interface, Pos) %>% 
  summarise(median_res = median(abs_residual))

SOD1_median_res$interface <- factor(SOD1_median_res$interface,
                                        levels = c("TRUE", "FALSE"))


test_median_subs <- wilcox.test(median_res ~ interface, data = SOD1_median_res)

p_median_dist <- ggplot(SOD1_median_res, aes(x = interface, y = median_res, color = interface)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_jitter(width = 0.1, alpha = 1, size = 4) +
  scale_color_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1")) +
  scale_x_discrete(labels = c( "FALSE" = "Non-dimer interface", "TRUE" = "Dimer interface")) +
  theme_classic() +
  labs(x = "", y = "|Median of binding\nresidual by position|", title = "Substitutions") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "none") +
  annotate(
    "text",
    x = 1.5,
    y = max(SOD1_median_res$median_res) * 1.1,
    label = paste0("p = ", signif(test_median_subs$p.value, 3)),
    size = 6
  )
p_median_dist


p_residual_density_subs_median <- ggplot(SOD1_median_res, aes(x = median_res, fill = interface)) +
  geom_density(alpha = 0.9) +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1"),
                    labels = c("TRUE" = "Dimer interface", "FALSE" = "Non-dimer interface")) +
  labs(x = "|Median of heterodimerization residual|", y = "Density", title = "Substitutions", fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 22),
        legend.position = c(0.8, 0.5)) +
  annotate(
    "text",
    x = 0.1,
    y = 20,
    label = paste0("p = ", signif(test_median_subs$p.value, 3)),
    size = 10
  )

p_residual_density_subs_median

ggsave(p_residual_density_subs_median, file = "density_residuals_int_vs_non_int_subs.tiff", width = 8, height = 6)

p_residual_subs_pos <- ggarrange(p_median_dist,
                                 p_residual_density_subs_median,
                                 nrow = 2, align = "v")

p_residual_subs_pos
ggsave(p_residual_subs_pos, file = "residuals_subs_by_pos.tiff", width = 10, height = 10)


SOD1_subs_outlier_counts <- SOD1_subsitutions %>%
  group_by(interface) %>%
  dplyr::count(outlier) %>%
  mutate(
    prop = n / sum(n),
    perc = prop * 100
  ) 

SOD1_subs_outlier_counts <- SOD1_subs_outlier_counts %>% 
  filter(outlier == T)


p_outliers <- ggplot(SOD1_subs_outlier_counts, aes(x = interface, y = perc, fill = interface)) +
  geom_col(color = "black", width = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1")) +
  scale_x_discrete(labels = c( "FALSE" = "Non-dimer interface", "TRUE" = "Dimer interface")) +
  labs(x = "", y = "% of significant residuals", title = "Substitutions") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.position = "none") 

p_outliers

ggsave(p_outliers, file = "subs_outlier_prop.tiff", width = 8, height = 6)

##### ALL MUTATIONS####
SOD1_allmut$interface <- factor(SOD1_allmut$interface,
                                          levels = c("TRUE", "FALSE"))

test_allmut <- wilcox.test(abs_residual ~ interface, data = SOD1_allmut)


p_residual_dist_all <- ggplot(SOD1_allmut, aes(x = interface, y = abs_residual, fill = interface)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1")) +
  scale_x_discrete(labels = c("FALSE" = "Non-dimer interface", "TRUE" = "Dimer interface")) +
  theme_classic() +
  labs(x = "", y = "|Binding residual|", title = "All") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "none") +
  annotate(
    "text",
    x = 1.5,
    y = 0.55,
    label = paste0("p = ", signif(test_allmut$p.value, 3)),
    size = 6
  )
p_residual_dist_all


p_residual_density_all <- ggplot(SOD1_allmut, aes(x = abs_residual, fill = interface)) +
  geom_density(alpha = 0.9) +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1"),
                    labels = c("TRUE" = "Dimer interface", "FALSE" = "Non-dimer interface")) +
  labs(x = "|Binding residual", y = "Density", title = "All", fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.position = c(0.8, 0.8))

p_residual_density_all

p_residual_all_mut <- ggarrange(p_residual_dist_all,
                                p_residual_density_all,
                                 nrow =  2, align = "v")

p_residual_all_mut

ggsave(p_residual_all_mut, file = "residuals_allmut_by_mut.tiff", width = 10, height = 10)

SOD1_median_res <- SOD1_allmut %>% 
  group_by(interface, Pos) %>% 
  summarise(median_res = median(abs_residual))

SOD1_median_res$interface <- factor(SOD1_median_res$interface,
                                        levels = c("TRUE", "FALSE"))


test_median_all <- wilcox.test(median_res ~ interface, data = SOD1_median_res)

p_median_dist_all_pos <- ggplot(SOD1_median_res, aes(x = interface, y = median_res, color = interface)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_jitter(width = 0.1, alpha = 1, size = 4) +
  scale_color_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1")) +
  scale_x_discrete(labels = c( "FALSE" = "Non-dimer interface", "TRUE" = "Dimer interface")) +
  theme_classic() +
  labs(x = "", y = "|Median of binding\nresidual by position|", title = "All") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "none") +
  annotate(
    "text",
    x = 1.5,
    y = 0.25,
    label = paste0("p = ", signif(test_median_all$p.value, 3)),
    size = 6
  )
p_median_dist_all_pos


p_residual_density_all_median <- ggplot(SOD1_median_res, aes(x = median_res, fill = interface)) +
  geom_density(alpha = 0.9) +
  scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1"),
                    labels = c("TRUE" = "Dimer interface", "FALSE" = "Non-dimer interface")) +
  labs(x = "|Median of binding residual|", y = "Density", title = "All mut", fill = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 24),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 18),
        legend.position = c(0.8, 0.8)) +
  annotate(
    "text",
    x = 0.15,
    y = 18,
    label = paste0("p = ", signif(test_median_all$p.value, 3)),
    size = 10
  )

p_residual_density_all_median
ggsave(p_residual_density_all_median, file = "density_residuals_int_vs_non_int_all.tiff", width = 8, height = 6)

p_residual_all_pos <- ggarrange(p_median_dist_all_pos,
                                p_residual_density_all_median,
                                 nrow = 2, align = "v")

p_residual_all_pos
ggsave(p_residual_all_pos, file = "residuals_allmut_by_pos.tiff", width = 10, height = 10)

# 
# SOD1_all_outlier_counts <- SOD1_allmut %>% 
#   group_by(interface) %>% 
#   count(outlier) %>% 
#   mutate(prop = n/sum(n))
# 
# SOD1_all_outlier_counts <- SOD1_all_outlier_counts %>% 
#   filter(outlier == T)
# 
# p_outliers <- ggplot(SOD1_all_outlier_counts, aes(x = interface_new, y = prop, fill = interface_new)) +
#   geom_col(color = "black", width = 0.6) +
#   scale_fill_manual(values = c("TRUE" = "#2A5783", "FALSE" = "#B9DDF1")) +
#   scale_x_discrete(labels = c( "FALSE" = "Non-dimer interface", "TRUE" = "Dimer interface")) +
#   labs(x = "", y = "Proportion of\n significant residuals", title = "All") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 22),
#         axis.title = element_text(size = 24),
#         plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
#         legend.text = element_text(size = 18),
#         legend.position = "none")
# 
# p_outliers
# 
# ggsave(p_outliers, file = "allmut_outlier_prop.tiff", width = 8, height = 6)


#Mutational impact of interface mutations in aPCA and bPCA####
interface <- interface %>% 
  mutate(
    category_fdr_apca = factor(category_fdr_apca,
                               levels = c("high", "WT-like", "low", "stop-like")),
    category_fdr_bpca = factor(category_fdr_bpca,
                               levels = c("high", "WT-like", "low", "stop-like"))
  )

median_residual_interface <- median(interface$abs_residual)


interface_high_res <- interface %>% 
  filter(residual_intensity == "high")

labels_cat <- c(
  "high" = "high",
  "WT-like" = "WT-like",
  "low" = "low",
  "stop-like" = "stop-like"
)

p_lluvial_fdr_int_highres <- ggplot(interface_high_res,
                                     aes(axis1 = category_fdr_apca,
                                         axis2 = category_fdr_bpca,
                                         y = 1)) +
  geom_alluvium(aes(fill = category_fdr_apca), width = 1/12) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black") +
  geom_label(stat = "stratum", aes(label = labels_cat[after_stat(stratum)]), size = 6) +
  scale_x_discrete(limits = c("Abundance\nFDR = 0.1", "Hetero-\ndimerization\nFDR = 0.1"), expand = c(.05, .05)) +
  scale_fill_manual(
    values = c(
      "high" = "darkgreen",
      "WT-like"   = "lightgrey",
      "low"   = "darkorange",
      "stop-like"    = "#DD4636")) +
  labs(y = "Number of mutations", title = "Dimer interface mutations\n(significant residual)") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        plot.title  = element_text(size = 26, face = "bold", hjust = 0.5),
        panel.grid = element_blank(),
        legend.position = "none")

p_lluvial_fdr_int_highres

ggsave(p_lluvial_fdr_int_highres, file = "alluvial_plot_interface_highres.tiff", width = 11, height = 8)


# heatmap_res_median <- interface %>%
#   group_by(mutation_type,Pos) %>%
#   summarise(residual_median = median(residuals, na.rm = TRUE)) %>%
#   ungroup() %>%
#   mutate(row = "row")
# 
# 
# heatmap_res_median$Pos <- factor(
#   heatmap_res_median$Pos,
#   levels = as.character(sort(as.numeric(unique(heatmap_res_median$Pos))))
# )
# heatmap_res_median$mutation_type <- factor(heatmap_res_median$mutation_type,
#                                            levels = c("subs", "ins", "del"))
# 
# p_heatmap_residual_interface<- ggplot(data = heatmap_res_median, aes(x = Pos, y = row, fill=residual_median)) +
#   geom_tile(size = 0.2, color = "black") +
#   scale_fill_gradient2(
#     low = "#B99C33",
#     mid = "grey95",
#     high = "#AC7299",
#     midpoint = 0,
#     limits = c(min(heatmap_res_median$residual_median), 
#                max(heatmap_res_median$residual_median))) +
#   facet_wrap(~mutation_type, nrow = 3, labeller = as_labeller(c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions"))) +
#   theme_classic() +
#   labs(x = "Dimer interface position", y = "", fill = "Median residual") +
#   theme(axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.text = element_text(size = 20, face = "bold"),
#         axis.text.x = element_text(size = 18),
#         axis.title.x = element_text(size = 24),
#         legend.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         legend.text = element_text(size=18),
#         legend.position = "right"
#   )
# p_heatmap_residual_interface

interface_subs <- interface %>% 
  filter(mutation_type == "subs")
interface_subs$Pos <- factor(interface_subs$Pos)

interface_subs <- interface_subs %>% 
  mutate(residual_type_position = case_when(
    Pos %in% c(5,50,51,52,53,54,150) ~ "negative",
    Pos %in% c(7,113,114,148,151,152,153) ~ "positive"
  )) %>% 
  mutate(WT_pos = paste0(WT_AA, Pos, ""))


interface_subs_wilcox <- interface_subs %>%
  group_by(Pos) %>%
  summarise(
    p_value = wilcox.test(residuals, mu = 0)$p.value,
    median = median(residuals),
    n = n()
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significant_pos = case_when(
           p_adj < 0.05 ~ TRUE,
           TRUE ~FALSE
         ))

interface_subs <- interface_subs %>%
  left_join(interface_subs_wilcox %>% select(Pos, significant_pos),
            by = "Pos")

interface_subs$WT_pos <- factor(interface_subs$WT_pos,
                                    levels = c("V5","V7","F50","G51","D52","N53","T54", "I113","G114", "R115","V148","G150","I151","A152","Q153"))

p_int_subs <- ggplot(interface_subs, aes(x = WT_pos, y = residuals)) +
  geom_violin(aes(fill = residual_type_position)) +
  geom_jitter(position = position_dodge(), color = "black", alpha = 0.5) +
  scale_fill_manual(values = c(positive = "#7C4D79", "negative" = "#AD9024")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Dimer interface position", y = "Heterodimerization residual", title = "Subtitutions") +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.position = "none") +
  geom_text(
    data = subset(interface, significant_pos == TRUE),
    aes(x = WT_pos, y = max(residuals) * 1.1, label = "*"),
    color = "black",
    size = 8
  ) +
  coord_flip()

p_int_subs

ggsave(p_int_subs, file = "interface_subs_residuals_vioplot.tiff", width = 6, height = 10)

interface_ins <- interface %>% 
  filter(mutation_type == "ins")


interface_ins <- interface_ins %>% 
  mutate(residual_type_position = case_when(
    Pos %in% c(5,50,51,52,53,54,150) ~ "negative",
    Pos %in% c(7,113,114,148,151,152,153) ~ "positive"
  )) %>% 
  mutate(WT_pos = paste0(WT_AA, Pos, ""))

interface_ins_wilcox <- interface_ins %>%
  group_by(Pos) %>%
  summarise(
    p_value = wilcox.test(residuals, mu = 0)$p.value,
    median = median(residuals),
    n = n()
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significant_pos = case_when(
           p_adj < 0.05 ~ TRUE,
           TRUE ~FALSE
         ))

interface_ins <- interface_ins %>%
  left_join(interface_ins_wilcox %>% select(Pos, significant_pos),
            by = "Pos")


interface_ins$WT_pos <- factor(interface_ins$WT_pos,
                                    levels = c("V5","V7","F50","G51","D52","N53","T54", "I113","G114", "R115","V148","G150","I151","A152","Q153"))

p_int_ins <- ggplot(interface_ins, aes(x = WT_pos, y = residuals)) +
  geom_violin(aes(fill = residual_type_position)) +
  geom_jitter(position = position_dodge(), color = "black", alpha = 0.5) +
  scale_fill_manual(values = c(positive = "#7C4D79", "negative" = "#AD9024")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Dimer interface position", y = "Heterodimerization residual", title = "Insertions") +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.position = "none") +
  geom_text(
    data = subset(interface_ins, significant_pos == TRUE),
    aes(x = WT_pos, y = max(residuals) * 1.1, label = "*"),
    color = "black",
    size = 8
  )

p_int_ins

ggsave(p_int_ins, file = "interface_ins_residuals_vioplot.tiff", width = 10, height = 6)


interface_del <- interface %>% 
  filter(mutation_type == "del")

interface_del <- interface_del %>% 
  mutate(residual_type_position = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  )) %>% 
  mutate(WT_pos = paste0(WT_AA, Pos, ""))



interface_del$WT_pos <- factor(interface_del$WT_pos,
                                   levels = c("V5","V7","F50","G51","D52","N53","T54", "I113","G114", "R115","V148","G150","I151","A152","Q153"))

p_int_del <- ggplot(interface_del, aes(x = WT_pos, y = residuals)) +
  geom_col(aes(fill = residual_type_position), size = 5, color = "black") +
  scale_fill_manual(values = c(positive = "#7C4D79", "negative" = "#AD9024")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Dimer interface position", y = "Heterodimerization residual", title = "Deletions") +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        legend.position = "none") +
  geom_text(
    data = subset(interface_del, outlier == TRUE),
    aes(x = WT_pos, y = max(residuals) * 1.1, label = "*"),
    color = "black",
    size = 8
  )


p_int_del

ggsave(p_int_del, file = "interface_del_residuals_vioplot.tiff", width = 10, height = 6)
