library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(ggbeeswarm)
library(ggforce)

#load modelled data
SOD1_allmut <- read.csv("SOD1_allmut_modelled.csv")
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")
SOD1_insertions <-read.csv("SOD1_ins_modelled.csv")
SOD1_deletions <- read.csv("SOD1_del_modelled.csv")



#adjust data
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(s_location= factor(s_location, levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                                     "β4", "Zn binding loop", "β5", "loop 5", "β6", "loop 6",
                                                     "β7", "Electrostatic loop", "β8", "loop 8", "Zn binding residues", "Cu binding residues")))

SOD1_insertions <- SOD1_insertions %>% 
  mutate(s_location= factor(s_location, levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                                   "β4", "Zn binding loop", "β5", "loop 5", "β6", "loop 6",
                                                   "β7", "Electrostatic loop", "β8", "loop 8", "Zn binding residues", "Cu binding residues")))


SOD1_insertions <- na.omit(SOD1_insertions)

SOD1_deletions <- SOD1_deletions %>% 
  mutate(s_location= factor(s_location, levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                                   "β4", "Zn binding loop", "β5", "loop 5", "β6", "loop 6",
                                                   "β7", "Electrostatic loop", "β8", "loop 8", "Zn binding residues", "Cu binding residues")))


SOD1_allmut <- SOD1_allmut %>% 
  mutate(s_location= factor(s_location, levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                                   "β4", "Zn binding loop", "β5", "loop 5", "β6", "loop 6",
                                                   "β7", "Electrostatic loop", "β8", "loop 8", "Zn binding residues", "Cu binding residues")))


SOD1_allmut$mutation_type <- factor(SOD1_allmut$mutation_type, levels = c("subs", "ins", "del"))
SOD1_allmut$category_fdr_apca <- factor(SOD1_allmut$category_fdr_apca, levels = c("high", "WT-like", "low", "stop-like"))
SOD1_allmut$category_fdr_bpca <- factor(SOD1_allmut$category_fdr_bpca, levels = c("high", "WT-like", "low", "stop-like"))

SOD1_subsitutions$category_fdr_apca <- factor(SOD1_subsitutions$category_fdr_apca, levels = c("high", "WT-like", "low", "stop-like"))
SOD1_subsitutions$category_fdr_bpca <- factor(SOD1_subsitutions$category_fdr_bpca, levels = c("high", "WT-like", "low", "stop-like"))

SOD1_insertions$category_fdr_apca <- factor(SOD1_insertions$category_fdr_apca, levels = c("high", "WT-like", "low", "stop-like"))
SOD1_insertions$category_fdr_bpca <- factor(SOD1_insertions$category_fdr_bpca, levels = c("high", "WT-like", "low", "stop-like"))

SOD1_deletions$category_fdr_apca <- factor(SOD1_deletions$category_fdr_apca, levels = c("high", "WT-like", "low", "stop-like"))
SOD1_deletions$category_fdr_bpca <- factor(SOD1_deletions$category_fdr_bpca, levels = c("high", "WT-like", "low", "stop-like"))

#FDR categories proportions####
counts_fdr_apca <- SOD1_allmut %>% 
  group_by(mutation_type, category_fdr_apca) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(mutation_type) %>% 
  mutate(perc = n / sum(n) * 100) %>% 
  ungroup()

counts_fdr_apca <- na.omit(counts_fdr_apca)

p_bar_fdr_apca <- ggplot(counts_fdr_apca, 
                         aes(x = mutation_type, y = perc, fill = category_fdr_apca)) +
  geom_col(color = "black", width = 0.9) +
  scale_fill_manual(values = c("high" = "darkgreen",
                               "WT-like" = "lightgrey",
                               "low" = "darkorange",
                               "stop-like" = "#DD4636")) +
  scale_x_discrete(labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  theme_classic() +
  labs(x = "", y = "Percentage of mutations (%)", fill = "Abundance FDR = 0.1") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 28),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 22),
        legend.position = "none")
p_bar_fdr_apca
ggsave(p_bar_fdr_apca, file = "fdr_prop_abundance.tiff", width = 5, height = 8)


counts_fdr_bpca <- SOD1_allmut %>% 
  group_by(mutation_type, category_fdr_bpca) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(mutation_type) %>% 
  mutate(perc = n / sum(n) * 100) %>% 
  ungroup()

counts_fdr_bpca <- na.omit(counts_fdr_bpca)

p_bar_fdr_bpca <- ggplot(counts_fdr_bpca, 
                         aes(x = mutation_type, y = perc, fill = category_fdr_bpca)) +
  geom_col(color = "black", width = 0.9) +
  scale_fill_manual(values = c("high" = "darkgreen",
                               "WT-like" = "lightgrey",
                               "low" = "darkorange",
                               "stop-like" = "#DD4636")) +
  scale_x_discrete(labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  theme_classic() +
  labs(x = "", y = "Percentage of mutations (%)", fill = "Heterodimerization FDR = 0.1") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 28),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 22),
        legend.position = "none")
p_bar_fdr_bpca
ggsave(p_bar_fdr_bpca, file = "fdr_prop_binding.tiff", width = 5, height = 8)

#Mutational impact in strands and loops
comparisons <- list(
  c("subs", "ins"),
  c("subs", "del"),
  c("ins", "del"))

SOD1_allmut <- SOD1_allmut %>% 
  mutate(location = recode(location, "β-sheet" = "β-strand"))
SOD1_allmut <- na.omit(SOD1_allmut)

p_loc_apca <- ggplot(SOD1_allmut, aes(x = mutation_type, y = abundance_score, color = mutation_type)) +
  geom_violin(width = 1.0, size = 0.8, color = "black", alpha = 1, trim = TRUE, adjust = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(alpha = mutation_type), size = 2) +
  scale_alpha_manual(values = c("subs" = 0.2, "ins" = 0.2, "del" = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = c("subs" = "#0B2130", "ins" = "#1E6085", "del" = "#A1C8DC"),
                     labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  facet_wrap(~location, ncol = 2) +
  theme_classic() +
  labs(x = "", y = "Abundance score", color = "Mutation") +
  theme(
    axis.text = element_text(size = 30),
    axis.title.y = element_text(size = 33),
    strip.text = element_text(size = 31, margin = margin(b = 0.5, t = 0.5)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 38)
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 8)
  )) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE, color = "black", test = "wilcox.test", stat = "signif", y_position = c(0.25, 0.5, 0.75), size = 0.5, textsize = 6)

p_loc_apca

ggsave(p_loc_apca, file = "violin_loc_alltype_apca.tiff", width = 14, height = 7, dpi = 300)

#Mutational impact in side-chain
SOD1_allmut_noint <- SOD1_allmut %>% 
  filter(side_chain != "dimer interface")

SOD1_allmut_noint$side_chain <- factor(SOD1_allmut_noint$side_chain,
                                       levels = c("core", "surface", "Zn binding", "Cu binding", "Zn and Cu binding"))

p_side_chain_apca <- ggplot(SOD1_allmut_noint, aes(x = mutation_type, y = abundance_score, color = mutation_type)) +
  geom_violin(width = 1.0, size = 1.2, color = "black", alpha = 1, trim = TRUE, adjust = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2), aes(alpha = mutation_type), size = 2) +
  scale_alpha_manual(values = c("subs" = 0.2, "ins" = 0.2, "del" = 1)) +  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = c("subs" = "#0B2130", "ins" = "#1E6085", "del" = "#A1C8DC"),
                     labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  facet_wrap(~side_chain, nrow = 1) +
  theme_classic() +
  labs(x = "", y = "Abundance score", color = "Mutation") +
  theme(
    axis.text = element_text(size = 30),
    axis.title.y = element_text(size = 33),
    strip.text = element_text(size = 31, margin = margin(b = 0.5, t = 0.5)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 38)
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 8)
  )) +
  geom_signif(data = subset(SOD1_allmut_noint, !side_chain %in% "Zn and Cu binding"), 
    comparisons = list(c("subs", "ins"), c("subs", "del"), c("ins", "del")),
    map_signif_level = TRUE, test = "wilcox.test", y_position = c(0.25, 0.5, 0.75), size = 0.5, color = "black", textsize = 6) +
  geom_signif(data = subset(SOD1_allmut_noint, side_chain %in% "Zn and Cu binding"), #omit test with single point
    comparisons = list(c("subs", "ins")),
    map_signif_level = TRUE, test = "wilcox.test", y_position = 0.25, size = 0.5, color = "black", textsize = 6)
 

p_side_chain_apca

ggsave(p_side_chain_apca, file = "violin_side_chain_alltype_apca.tiff", width = 20, height = 7, dpi = 300)

p_loc_bpca <- ggplot(SOD1_allmut, aes(x = mutation_type, y = binding_score, color = mutation_type)) +
  geom_violin(width = 1.0, size = 0.8, color = "black", alpha = 1, trim = TRUE, adjust = 0.5) +  
  geom_jitter(position = position_jitter(width = 0.2), aes(alpha = mutation_type), size = 2) +
  scale_alpha_manual(values = c("subs" = 0.2, "ins" = 0.2, "del" = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = c("subs" = "#0B2130", "ins" = "#1E6085", "del" = "#A1C8DC"),
                     labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  facet_wrap(~location, ncol = 2) +
  theme_classic() +
  labs(x = "", y = "Heterodimerization score", color = "Mutation") +
  theme(
    axis.text = element_text(size = 30),
    axis.title.y = element_text(size = 33),
    strip.text = element_text(size = 31, margin = margin(b = 0.5, t = 0.5)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 38)
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 8)
  )) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE, color = "black",
    test = "wilcox.test", stat = "signif", y_position = c(0.25, 0.5, 0.75), size = 0.5, textsize = 6)

p_loc_bpca
ggsave(p_loc_bpca, file = "violin_loc_alltype_bpca.tiff", width = 14, height = 7, dpi = 300)


#Mutational impact in side-chain
SOD1_allmut$side_chain <- factor(SOD1_allmut$side_chain,
                                       levels = c("dimer interface", "second_shell", "core", "surface",  "Zn binding", "Cu binding", "Zn and Cu binding"))

p_side_chain_bpca <- ggplot(SOD1_allmut, aes(x = mutation_type, y = binding_score, color = mutation_type)) +
  geom_violin(width = 1.0, size = 1.2, color = "black", alpha = 1, trim = TRUE, adjust = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2), aes(alpha = mutation_type), size = 2) +
  scale_alpha_manual(values = c("subs" = 0.2, "ins" = 0.2, "del" = 1)) +  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = c("subs" = "#0B2130", "ins" = "#1E6085", "del" = "#A1C8DC"),
                     labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  facet_wrap(~side_chain, nrow = 1) +
  theme_classic() +
  labs(x = "", y = "Heterodimerization score", color = "Mutation") +
  theme(
    axis.text = element_text(size = 30),
    axis.title.y = element_text(size = 33),
    strip.text = element_text(size = 31, margin = margin(b = 0.5, t = 0.5)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 38)
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 8)
  )) +
  geom_signif(data = subset(SOD1_allmut, !side_chain %in% "Zn and Cu binding"),
    comparisons = list(c("subs", "ins"), c("subs", "del"), c("ins", "del")),
    map_signif_level = TRUE, test = "wilcox.test", y_position = c(0.25, 0.5, 0.75), size = 0.5, color = "black", textsize = 6) +
  geom_signif(data = subset(SOD1_allmut, side_chain %in% "Zn and Cu binding"),
    comparisons = list(c("subs", "ins")), map_signif_level = TRUE, test = "wilcox.test", y_position = 0.25, size = 0.5, color = "black", textsize = 6)

p_side_chain_bpca
ggsave(p_side_chain_bpca, file = "violin_side_chain_alltype_bpca.tiff", width = 24, height = 7, dpi = 300)


# Mutant aminoacids impact####
SOD1_allmut <- SOD1_allmut %>%
  mutate(Property_Mut = case_when(
    Mut_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    Mut_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    Mut_AA %in% c("K", "R", "H") ~ "Positive charge",
    Mut_AA %in% c("D", "E") ~ "Negative charge",
    Mut_AA %in% c("F", "Y", "W") ~ "Aromatic",
    Mut_AA %in% "G" ~ "Glycine",
    Mut_AA %in% "P" ~ "Proline"
  )) %>%
  mutate(Property_wt =  case_when(
    WT_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    WT_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    WT_AA %in% c("K", "R", "H") ~ "Positive charge",
    WT_AA %in% c("D", "E") ~ "Negative charge",
    WT_AA %in% c("F", "Y", "W") ~ "Aromatic",
    WT_AA %in% "G" ~ "Glycine",
    WT_AA %in% "P" ~ "Proline"))


p_vio_mutaa <- ggplot(SOD1_allmut, aes(x = Property_Mut, y = abundance_score)) +
  geom_violin(width = 0.7,size = 1,color = "black",alpha = 0.3,trim = TRUE,adjust = 0.5 ) +
  geom_sina(maxwidth = 0.4,   alpha = 0.7,size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~location) +
  theme_classic()
p_vio_mutaa


p_vio_wtaa <- ggplot(SOD1_allmut, aes(x = Property_wt, y = abundance_score)) +
  geom_violin(width = 0.7,size = 1,color = "black",alpha = 0.3,trim = TRUE,adjust = 0.5 ) +
  geom_sina(maxwidth = 0.4,   alpha = 0.7,size = 1.5) +
  facet_wrap(~location) +
  theme_classic()
p_vio_wtaa


#Correlation abundance - rSASA ####
#Substitutions
r_labels_sasa_subs <- SOD1_subsitutions %>%
  summarise(
    R = cor(abundance_score, rSASA.monomer, method = "pearson"),
    p_value = cor.test(abundance_score, rSASA.monomer, method = "pearson", exact = FALSE)$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_corr_sasa_ab <- ggplot(SOD1_subsitutions, aes(x = rSASA.monomer, y = abundance_score)) +
  geom_point(size = 3, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 25, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "rSASA SOD1 monomer", y = "Substitutions\nabundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_labels_sasa_subs, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_corr_sasa_ab

#Insertions
r_labels_sasa_ins <- SOD1_insertions %>%
  summarise(
    R = cor(abundance_score, rSASA.monomer, method = "pearson"),
    p_value = cor.test(abundance_score, rSASA.monomer, method = "pearson", exact = FALSE)$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_corr_sasa_ab_ins <- ggplot(SOD1_insertions, aes(x = rSASA.monomer, y = abundance_score)) +
  geom_point(size = 3, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 25, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "rSASA SOD1 monomer", y = "Insertions\nabundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_labels_sasa_ins, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_corr_sasa_ab_ins


#Deletions
r_labels_sasa_del <- SOD1_deletions %>%
  summarise(
    R = cor(abundance_score, rSASA.monomer, method = "pearson"),
    p_value = cor.test(abundance_score, rSASA.monomer, method = "pearson", exact = FALSE)$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )

p_corr_sasa_ab_del <- ggplot(SOD1_deletions, aes(x = rSASA.monomer, y = abundance_score)) +
  geom_point(size = 4, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 25, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "rSASA SOD1 monomer", y = "Deletions\nabundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_labels_sasa_del, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE) 
p_corr_sasa_ab_del

#Comparing median abundance impact in  subs and indels ####
#Substitutions
median_by_location <- SOD1_subsitutions %>% 
  group_by(s_location) %>% 
  summarise(median_loc = median(abundance_score))

p_median_loc_subs <- ggplot(median_by_location, aes(x = median_loc, y = s_location)) +
  geom_col(aes(fill = median_loc), color = "black", linewidth = 0.7, width = 0.8) +
  scale_fill_gradient2(
    low = "darkorange",
    mid = "grey95",
    high = "darkgreen",
    midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "solid", size = 1) +
  theme_classic() +
  labs(x = "Median abundance score", y = "", fill = "Median abundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        legend.position = "top")
p_median_loc_subs

median_by_location <- SOD1_subsitutions %>% 
  group_by(location) %>% 
  summarise(median_loc = median(abundance_score))

p_median_sloc_subs <- ggplot(median_by_location, aes(x = median_loc, y = location)) +
  geom_col(aes(fill = median_loc), color = "black", linewidth = 0.7, width = 0.8) +
  scale_fill_gradient2(
    low = "darkorange",
    mid = "grey95",
    high = "darkgreen",
    midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "solid", size = 1) +
  theme_classic() +
  labs(x = "Median abundance score", y = "", fill = "Median abundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        legend.position = "top")
p_median_sloc_subs

#Insertions
median_by_location <- SOD1_insertions %>% 
  group_by(s_location) %>% 
  summarise(median_loc = median(abundance_score))

p_median_loc_ins <- ggplot(median_by_location, aes(x = median_loc, y = s_location)) +
  geom_col(aes(fill = median_loc), color = "black", linewidth = 0.7, width = 0.8) +
  scale_fill_gradient2(
    low = "darkorange",
    mid = "grey95",
    high = "darkgreen",
    midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "solid", size = 1) +
  theme_classic() +
  labs(x = "Median abundance score", y = "", fill = "Median abundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        legend.position = "top")
p_median_loc_ins

median_by_s_location <- SOD1_insertions %>% 
  group_by(location) %>% 
  summarise(median_loc = median(abundance_score))

p_median_sloc_ins <- ggplot(median_by_s_location, aes(x = median_loc, y = location)) +
  geom_col(aes(fill = median_loc), color = "black", linewidth = 0.7, width = 0.8) +
  scale_fill_gradient2(
    low = "darkorange",
    mid = "grey95",
    high = "darkgreen",
    midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "solid", size = 1) +
  theme_classic() +
  labs(x = "Median abundance score", y = "", fill = "Median abundance score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        legend.position = "top")
p_median_sloc_ins


p_median_loc <- ggarrange(p_median_loc_subs,
                          p_median_loc_ins + theme(axis.text.y = element_blank()),
                          ncol = 2, common.legend = TRUE)
p_median_loc

ggsave(p_median_loc, file = "median_ab_subs_ins.tiff", width = 15, height = 10, dpi = 300)


# Correlation between type of mutations####
#adjust the dataset to be able to compare between mutation types
#create a column to join substitutions and insertions
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(Pos_mut = paste0(Pos, Mut_AA))

SOD1_insertions <- SOD1_insertions %>% 
  mutate(Pos_mut = paste0(Pos, Mut_AA))

subs_ins <- inner_join(SOD1_subsitutions, SOD1_insertions, by = "Pos_mut")
subs_ins <- rename(subs_ins, abundance_score_subs = abundance_score.x, abundance_score_ins = abundance_score.y)

#plot correlation
r_labels_sub_ins <- subs_ins %>%
  summarise(
    R = cor(abundance_score_subs, abundance_score_ins, method = "pearson"),
    p_value = cor.test(abundance_score_subs, abundance_score_ins, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_cor_subs_ins <- ggplot(subs_ins, aes(x = abundance_score_subs, y = abundance_score_ins)) +
  geom_point(aes(color = side_chain.x), size = 4, alpha = 0.5) + 
  scale_color_manual(values = c("darkred", "darkblue", "goldenrod", "darkgreen"))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions\nabundance score", y = "Insertions\nabundance score", color = "") +
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


r_labels_sub_ins_loc <- subs_ins %>%
  group_by(s_location.x) %>% 
  summarise(
    R = cor(abundance_score_subs, abundance_score_ins, method = "pearson"),
    p_value = cor.test(abundance_score_subs, abundance_score_ins, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )

p_cor_subs_ins_sloc <- ggplot(subs_ins, aes(x = abundance_score_subs, y = abundance_score_ins)) +
  geom_point(aes(color = side_chain.x), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c("darkred", "darkblue", "goldenrod", "darkgreen"))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  facet_wrap(~s_location.x) +
  theme_classic() +
  labs(x = "Substitutions abundance score", y = "Insertions abundance score", color = "") +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 23),
        legend.position = "top",
        legend.text = element_text(size = 20),
        strip.text = element_text(face = "bold")) +
  guides(color = guide_legend(
    override.aes = list(size = 10))) +
  geom_text(data = r_labels_sub_ins_loc, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 3, inherit.aes = FALSE)  
p_cor_subs_ins_sloc

ggsave(p_cor_subs_ins_sloc, file = "matrix_corr_subs_ins_byloc.tiff", width = 12, height = 8, dpi = 300)

#subs vs dels
#to compare the data, I will use the median abundanc score per position of substitutions
median_subs_pos <- SOD1_subsitutions %>% 
  group_by(Pos) %>% 
  summarise(median_abundance = median(abundance_score))

subs_del <- inner_join(median_subs_pos, SOD1_deletions, by = "Pos")
subs_del <- rename(subs_del, abundance_score_del = abundance_score)

#plot correlation
r_labels_sub_del <- subs_del %>%
  summarise(
    R = cor(median_abundance, abundance_score_del, method = "pearson"),
    p_value = cor.test(median_abundance, abundance_score_del, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_cor_subs_del <- ggplot(subs_del, aes(x = median_abundance, y = abundance_score_del)) +
  geom_point(aes(color = side_chain), size = 4, alpha = 0.8) + 
  scale_color_manual(values = c("darkred", "darkblue", "goldenrod", "darkgreen"))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score by position", y = "Deletions\nabundance score", color = "") +
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



#ins vs dels
#to compare the data, I will use the median abundanc score per position of insertions
median_ins_pos <- SOD1_insertions %>% 
  group_by(Pos) %>% 
  summarise(median_abundance = median(abundance_score))

ins_del <- inner_join(median_ins_pos, SOD1_deletions, by = "Pos")
ins_del <- rename(ins_del, abundance_score_del = abundance_score)

#plot correlation
r_labels_ins_del <- ins_del %>%
  summarise(
    R = cor(median_abundance, abundance_score_del, method = "pearson"),
    p_value = cor.test(median_abundance, abundance_score_del, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\n", "p = ", format(p_value, digits = 2, scientific = TRUE))
  )


p_cor_ins_del <- ggplot(ins_del, aes(x = median_abundance, y = abundance_score_del)) +
  geom_point(aes(color = side_chain), size = 4, alpha = 0.8) + 
  scale_color_manual(values = c("darkred", "darkblue", "goldenrod", "darkgreen"))+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Insertions median\nabundance score by position", y = "Deletions\nabundance score", color = "") +
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


x_min <- min(subs_ins$abundance_score_subs, subs_ins$abundance_score_del, ins_del$median_abundance, na.rm = TRUE)
x_max <- max(subs_ins$abundance_score_subs, subs_ins$abundance_score_del, ins_del$median_abundance, na.rm = TRUE)

y_min <- min(subs_ins$abundance_score_ins, subs_ins$abundance_score_del, ins_del$abundance_score_del, na.rm = TRUE)
y_max <- max(subs_ins$abundance_score_ins, subs_ins$abundance_score_del, ins_del$abundance_score_del, na.rm = TRUE)

p_cor_subs_ins <- p_cor_subs_ins + coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
p_cor_subs_del <- p_cor_subs_del + coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
p_cor_ins_del <- p_cor_ins_del + coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))

p_corr_all <- ggarrange(p_cor_subs_ins,
                        p_cor_subs_del,
                        p_cor_ins_del,
                        ncol = 3, common.legend = TRUE, align = "hv")

p_corr_all


ggsave(p_corr_all, file = "corr_subs_ins_del.tiff", width = 25, height = 7, dpi = 300)


# FDR categories by side-chain orientation####
df_prop <- SOD1_allmut %>%
  group_by(mutation_type, category_fdr_apca, side_chain) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(side_chain) %>%
  mutate(prop = n / sum(n))



df_prop$category_fdr_apca <- factor(
  df_prop$category_fdr_apca,
  levels = c("high", "WT-like", "low", "stop-like") 
)

df_prop$mutation_type <- factor(
  df_prop$mutation_type,
  levels = c("subs", "ins", "del") 
)

df_prop$side_chain <- factor(
  df_prop$side_chain,
  levels = c("core", "surface", "dimer interface", "Zn binding") 
)

df_prop <- df_prop %>% 
  filter(prop > 0.0009)

p_prop_fdr <- ggplot(data = df_prop, aes(x = side_chain, y = prop, fill = category_fdr_apca)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.8, linewidth = 1) +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = c("high" = "darkgreen", "WT-like" = "lightgrey", "low" = "darkorange", "Stop-like" = "darkred")) +
  scale_x_discrete(labels = c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions")) +
  facet_wrap(~mutation_type, labeller = as_labeller(c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions"))) +
  geom_text(aes(label = sprintf("%.0f", prop * 100)),
            position = position_stack(vjust = 0.5),
            color = "black", size = 10) +
  labs(x = "", y = "Percetage of mutations (%)", fill= "Abundance\nFDR = 0.1") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 33),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 30),
        legend.title = element_text(size = 36, face = "bold"),
        legend.text = element_text(size = 33),
        legend.position = "top")
p_prop_fdr

ggsave(p_prop_fdr, path = path_fig2, file = "fdr_by_side_chain.tiff", width = 15, height = 10, dpi = 300)


#Comparison of WT_aa and Mut_aa between subs and ins####
#adjust the data
SOD1_subsitutions <- SOD1_subsitutions %>%
  mutate(Property_Mut = case_when(
    Mut_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    Mut_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    Mut_AA %in% c("K", "R", "H") ~ "Positive charge",
    Mut_AA %in% c("D", "E") ~ "Negative charge",
    Mut_AA %in% c("F", "Y", "W") ~ "Aromatic",
    Mut_AA %in% "G" ~ "Glycine",
    Mut_AA %in% "P" ~ "Proline"
  )) %>%
  mutate(Property_wt =  case_when(
    WT_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    WT_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    WT_AA %in% c("K", "R", "H") ~ "Positive charge",
    WT_AA %in% c("D", "E") ~ "Negative charge",
    WT_AA %in% c("F", "Y", "W") ~ "Aromatic",
    WT_AA %in% "G" ~ "Glycine",
    WT_AA %in% "P" ~ "Proline"))

SOD1_insertions <- SOD1_insertions %>%
  mutate(Property_Mut = case_when(
    Mut_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    Mut_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    Mut_AA %in% c("K", "R", "H") ~ "Positive charge",
    Mut_AA %in% c("D", "E") ~ "Negative charge",
    Mut_AA %in% c("F", "Y", "W") ~ "Aromatic",
    Mut_AA %in% "G" ~ "Glycine",
    Mut_AA %in% "P" ~ "Proline"
  )) %>%
  mutate(Property_wt =  case_when(
    WT_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    WT_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    WT_AA %in% c("K", "R", "H") ~ "Positive charge",
    WT_AA %in% c("D", "E") ~ "Negative charge",
    WT_AA %in% c("F", "Y", "W") ~ "Aromatic",
    WT_AA %in% "G" ~ "Glycine",
    WT_AA %in% "P" ~ "Proline"))

SOD1_deletions <- SOD1_deletions %>%
  mutate(Property_Mut = case_when(
    Mut_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    Mut_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    Mut_AA %in% c("K", "R", "H") ~ "Positive charge",
    Mut_AA %in% c("D", "E") ~ "Negative charge",
    Mut_AA %in% c("F", "Y", "W") ~ "Aromatic",
    Mut_AA %in% "G" ~ "Glycine",
    Mut_AA %in% "P" ~ "Proline"
  )) %>%
  mutate(Property_wt =  case_when(
    WT_AA %in% c("A", "V", "L", "I", "M") ~ "Hydrophobic",  
    WT_AA %in% c("S", "T", "C", "N", "Q") ~ "Polar",
    WT_AA %in% c("K", "R", "H") ~ "Positive charge",
    WT_AA %in% c("D", "E") ~ "Negative charge",
    WT_AA %in% c("F", "Y", "W") ~ "Aromatic",
    WT_AA %in% "G" ~ "Glycine",
    WT_AA %in% "P" ~ "Proline"))


#to compare the impact of each mutant, I will compare the median score of each
median_mut_aa_subs <- SOD1_subsitutions %>% 
  group_by(Mut_AA, Property_Mut) %>% 
  summarise(median_mut_aa_subs = median(abundance_score))

median_mut_aa_ins <- SOD1_insertions %>% 
  group_by(Mut_AA, Property_Mut) %>% 
  summarise(median_mut_aa_ins = median(abundance_score))

median_mut_both <- inner_join(median_mut_aa_subs, median_mut_aa_ins, by = "Mut_AA")

cor_test <- cor.test(median_mut_both$median_mut_aa_subs,
                     median_mut_both$median_mut_aa_ins,
                     method = "pearson",
                     exact = FALSE)

r_value <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
label <- paste0(
  "R = ", r_value, "\n",
  "p = ", format(p_value, scientific = TRUE, digits = 2)
)



p_cor_median_mutaa <- ggplot(median_mut_both, aes(x = median_mut_aa_subs, y = median_mut_aa_ins)) +
  geom_point(aes(color = Property_Mut.x), size = 5) +
  scale_color_manual(values = c("Hydrophobic" = "blue", "Polar" = "orange", "Positive charge" = "red", "Negative charge" = "darkgreen",
                     "Aromatic" = "#B254A5", "Glycine" = "#52A8BC", "Proline" = "#515A66")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score", y = "Insertions median\nabundance score", color = "Amino acid type") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 26),
        legend.position = "top",
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 23)) +
  geom_text_repel(aes(label = Mut_AA, color = Property_Mut.x),
                  size = 5,
                  min.segment.length = 0,   
                  segment.size = 0.5,       
                  segment.alpha = 0.8,
                  max.overlaps = Inf) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = -0.5, y = -0.15, label = label, size = 10, hjust = 0)
p_cor_median_mutaa

#to compare the impact on the WT type, I will compare the median score of each
median_wt_aa_subs <- SOD1_subsitutions %>% 
  group_by(WT_AA, Property_wt) %>% 
  summarise(median_wt_aa_subs = median(abundance_score))

median_wt_aa_ins <- SOD1_insertions %>% 
  group_by(WT_AA, Property_wt) %>% 
  summarise(median_wt_aa_ins = median(abundance_score))

median_wt_both <- inner_join(median_wt_aa_subs, median_wt_aa_ins, by = "WT_AA")

cor_test_wt <- cor.test(median_wt_both$median_wt_aa_subs,
                        median_wt_both$median_wt_aa_ins,
                     method = "pearson",
                     exact = FALSE)

r_value <- round(cor_test_wt$estimate, 2)
p_value <- signif(cor_test_wt$p.value, 2)
label <- paste0(
  "R = ", r_value, "\n",
  "p = ", format(p_value, scientific = TRUE, digits = 2)
)



p_cor_median_wtaa <- ggplot(median_wt_both, aes(x = median_wt_aa_subs, y = median_wt_aa_ins)) +
  geom_point(aes(color = Property_wt.x), size = 5) +
  scale_color_manual(values = c("Hydrophobic" = "blue", "Polar" = "orange", "Positive charge" = "red", "Negative charge" = "darkgreen",
                                "Aromatic" = "#B254A5", "Glycine" = "#52A8BC", "Proline" = "#515A66")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score", y = "Insertions median\nabundance score", color = "Amino acid type") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 26)) +
  geom_text_repel(aes(label = WT_AA, color = Property_wt.x),
                  size = 5,
                  min.segment.length = 0,   
                  segment.size = 0.5,       
                  segment.alpha = 0.8,
                  max.overlaps = Inf) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = -0.75, y = -0.15, label = label, size = 10, hjust = 0)

p_cor_median_wtaa



p_cor_combined <- ggarrange(p_cor_median_mutaa, 
                            p_cor_median_wtaa,
                            ncol = 2, common.legend = TRUE)
p_cor_combined

ggsave(p_cor_combined, file = "cor_mut_wt_aa_subs_ins.tiff", width = 14, height = 6, dpi = 300)


#subs vs dels
median_mut_aa_del <- SOD1_deletions %>% 
  group_by(Mut_AA, Property_Mut) %>% 
  summarise(median_mut_aa_del = median(abundance_score))

median_mut_subs_del <- inner_join(median_mut_aa_subs, median_mut_aa_del, by = "Mut_AA")

cor_test_subs_del <- cor.test(median_mut_subs_del$median_mut_aa_subs,
                              median_mut_subs_del$median_mut_aa_del,
                     method = "pearson",
                     exact = FALSE)

r_value <- round(cor_test_subs_del$estimate, 2)
p_value <- signif(cor_test_subs_del$p.value, 2)
label <- paste0(
  "R = ", r_value, "\n",
  "p = ", format(p_value, scientific = TRUE, digits = 2)
)



p_cor_median_subs_del_mutaa <- ggplot(median_mut_subs_del, aes(x = median_mut_aa_subs, y = median_mut_aa_del)) +
  geom_point(aes(color = Property_Mut.x), size = 5) +
  scale_color_manual(values = c("Hydrophobic" = "blue", "Polar" = "orange", "Positive charge" = "red", "Negative charge" = "darkgreen",
                                "Aromatic" = "#B254A5", "Glycine" = "#52A8BC", "Proline" = "#515A66")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score", y = "Deletions median\nabundance score", color = "Amino acid type") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 26)) +
  geom_text_repel(aes(label = Mut_AA, color = Property_Mut.x),
                  size = 5,
                  min.segment.length = 0,   
                  segment.size = 0.5,       
                  segment.alpha = 0.8,
                  max.overlaps = Inf) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = -0.5, y = -0.5, label = label, size = 10, hjust = 0)
p_cor_median_subs_del_mutaa

#to compare the impact on the WT type, I will compare the median score of each
median_wt_aa_del <- SOD1_deletions %>% 
  group_by(WT_AA, Property_wt) %>% 
  summarise(median_wt_aa_del = median(abundance_score))

median_wt_subs_del <- inner_join(median_wt_aa_subs, median_wt_aa_del, by = "WT_AA")

cor_test_wt_subs_del<- cor.test(median_wt_subs_del$median_wt_aa_subs,
                        median_wt_subs_del$median_wt_aa_del,
                        method = "pearson",
                        exact = FALSE)

r_value <- round(cor_test_wt_subs_del$estimate, 2)
p_value <- signif(cor_test_wt_subs_del$p.value, 2)
label <- paste0(
  "R = ", r_value, "\n",
  "p = ", format(p_value, scientific = TRUE, digits = 2)
)



p_cor_median_wtaa_subsdel <- ggplot(median_wt_subs_del, aes(x = median_wt_aa_subs, y = median_wt_aa_del)) +
  geom_point(aes(color = Property_wt.x), size = 5) +
  scale_color_manual(values = c("Hydrophobic" = "blue", "Polar" = "orange", "Positive charge" = "red", "Negative charge" = "darkgreen",
                                "Aromatic" = "#B254A5", "Glycine" = "#52A8BC", "Proline" = "#515A66")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score", y = "Deletions median\nabundance score", color = "Amino acid type") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 26)) +
  geom_text_repel(aes(label = WT_AA, color = Property_wt.x),
                  size = 5,
                  min.segment.length = 0,   
                  segment.size = 0.5,       
                  segment.alpha = 0.8,
                  max.overlaps = Inf) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = -0.75, y = -0.15, label = label, size = 10, hjust = 0)

p_cor_median_wtaa_subsdel


p_cor_subs_del <- ggarrange(p_cor_median_subs_del_mutaa, 
                            p_cor_median_wtaa_subsdel,
                            ncol = 2, common.legend = TRUE)
p_cor_subs_del

ggsave(p_cor_subs_del, file = "cor_mut_wt_aa_subs_del.tiff", width = 14, height = 6, dpi = 300)


#ins vs dels
median_mut_ins_del <- inner_join(median_mut_aa_ins, median_mut_aa_del, by = "Mut_AA")

cor_test_ins_del <- cor.test(median_mut_ins_del$median_mut_aa_ins,
                             median_mut_ins_del$median_mut_aa_del,
                              method = "pearson",
                              exact = FALSE)

r_value <- round(cor_test_ins_del$estimate, 2)
p_value <- signif(cor_test_ins_del$p.value, 2)
label <- paste0(
  "R = ", r_value, "\n",
  "p = ", format(p_value, scientific = TRUE, digits = 2)
)



p_cor_median_ins_del_mutaa <- ggplot(median_mut_ins_del, aes(x = median_mut_aa_ins, y = median_mut_aa_del)) +
  geom_point(aes(color = Property_Mut.x), size = 5) +
  scale_color_manual(values = c("Hydrophobic" = "blue", "Polar" = "orange", "Positive charge" = "red", "Negative charge" = "darkgreen",
                                "Aromatic" = "#B254A5", "Glycine" = "#52A8BC", "Proline" = "#515A66")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Insertions median\nabundance score", y = "Deletions median\nabundance score", color = "Amino acid type") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 26)) +
  geom_text_repel(aes(label = Mut_AA, color = Property_Mut.x),
                  size = 5,
                  min.segment.length = 0,   
                  segment.size = 0.5,       
                  segment.alpha = 0.8,
                  max.overlaps = Inf) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = -0.7, y = -0.5, label = label, size = 10, hjust = 0)
p_cor_median_ins_del_mutaa

#to compare the impact on the WT type, I will compare the median score of each
median_wt_ins_del <- inner_join(median_wt_aa_ins, median_wt_aa_del, by = "WT_AA")

cor_test_wt_ins_del<- cor.test(median_wt_ins_del$median_wt_aa_ins,
                               median_wt_ins_del$median_wt_aa_del,
                                method = "pearson",
                                exact = FALSE)

r_value <- round(cor_test_wt_ins_del$estimate, 2)
p_value <- signif(cor_test_wt_ins_del$p.value, 2)
label <- paste0(
  "R = ", r_value, "\n",
  "p = ", format(p_value, scientific = TRUE, digits = 2)
)



p_cor_median_wtaa_insdel <- ggplot(median_wt_ins_del, aes(x = median_wt_aa_ins, y = median_wt_aa_del)) +
  geom_point(aes(color = Property_wt.x), size = 5) +
  scale_color_manual(values = c("Hydrophobic" = "blue", "Polar" = "orange", "Positive charge" = "red", "Negative charge" = "darkgreen",
                                "Aromatic" = "#B254A5", "Glycine" = "#52A8BC", "Proline" = "#515A66")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Insertions median\nabundance score", y = "Deletions median\nabundance score", color = "Amino acid type") +
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 26)) +
  geom_text_repel(aes(label = WT_AA, color = Property_wt.x),
                  size = 5,
                  min.segment.length = 0,   
                  segment.size = 0.5,       
                  segment.alpha = 0.8,
                  max.overlaps = Inf) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = -0.75, y = -0.5, label = label, size = 10, hjust = 0)

p_cor_median_wtaa_insdel



p_cor_combined_subsins <- ggarrange(p_cor_median_mutaa, 
                           p_cor_median_wtaa,
                            ncol = 2, common.legend = TRUE)
p_cor_combined_subsins

p_cor_combined_subsdel <- ggarrange(p_cor_median_subs_del_mutaa, 
                                    p_cor_median_wtaa_subsdel,
                                    ncol = 2, common.legend = TRUE)
p_cor_combined_subsdel

p_cor_combined_insdel <- ggarrange(p_cor_median_ins_del_mutaa, 
                                   p_cor_median_wtaa_insdel,
                                    ncol = 2, common.legend = TRUE)
p_cor_combined_insdel


ggsave(p_cor_combined_subsins, file = "cor_mut_wt_aa_subs_ins.tiff", width = 14, height = 6, dpi = 300)
ggsave(p_cor_combined_subsdel,  file = "cor_mut_wt_aa_subs_del.tiff", width = 14, height = 6, dpi = 300)
ggsave(p_cor_combined_insdel, file = "cor_mut_wt_aa_del_ins.tiff", width = 14, height = 6, dpi = 300)


#Determining  mutants with higher impact and position more affected####
#Substitutions
median_mut_aa_subs <- SOD1_subsitutions %>% 
  group_by(Property_Mut, mutation_type) %>% 
  summarise(median_abundance = median(abundance_score))

median_mut_aa_subs <- median_mut_aa_subs %>%
  mutate(Property_Mut = factor(Property_Mut, levels = median_mut_aa_subs %>% arrange(desc(median_abundance)) %>% pull(Property_Mut)))

p_median_mut_subs <- ggplot(median_mut_aa_subs, aes(x = median_abundance, y = Property_Mut)) +
  geom_col(fill = "#49525E", color = "black", linewidth = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score", y = "") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 26))
p_median_mut_subs


#Insertions
median_mut_aa_ins <- SOD1_insertions %>% 
  group_by(Property_Mut, mutation_type) %>% 
  summarise(median_abundance = median(abundance_score))

median_mut_aa_ins <- median_mut_aa_ins %>%
  mutate(Property_Mut = factor(Property_Mut, levels = median_mut_aa_ins %>% arrange(desc(median_abundance)) %>% pull(Property_Mut)))

p_median_mut_ins <- ggplot(median_mut_aa_ins, aes(x = median_abundance, y = Property_Mut)) +
  geom_col(fill = "#49525E", color = "black", linewidth = 1) +
  theme_classic() +
  labs(x = "Insertions median\nabundance score", y = "") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 26))
p_median_mut_ins


#Deletions
median_mut_aa_del <- SOD1_deletions %>% 
  group_by(Property_Mut, mutation_type) %>% 
  summarise(median_abundance = median(abundance_score))

median_mut_aa_del <- median_mut_aa_del %>%
  mutate(Property_Mut = factor(Property_Mut, levels = median_mut_aa_del %>% arrange(desc(median_abundance)) %>% pull(Property_Mut)))

p_median_mut_del <- ggplot(median_mut_aa_del, aes(x = median_abundance, y = Property_Mut)) +
  geom_col(fill = "lightgrey", color = "black", linewidth = 1) +
  theme_classic() +
  labs(x = "Deletions median\nabundance score", y = "") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 26))
p_median_mut_del



#Substitutions
median_wt_aa_subs <- SOD1_subsitutions %>% 
  group_by(Property_wt, mutation_type) %>% 
  summarise(median_abundance = median(abundance_score))

median_wt_aa_subs <- median_wt_aa_subs %>%
  mutate(Property_wt = factor(Property_wt, levels = median_wt_aa_subs %>% arrange(desc(median_abundance)) %>% pull(Property_wt)))

p_median_WT_subs <- ggplot(median_wt_aa_subs, aes(x = median_abundance, y = Property_wt)) +
  geom_col(fill = "lightgrey", color = "black", linewidth = 1) +
  theme_classic() +
  labs(x = "Substitutions median\nabundance score", y = "") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 26))
p_median_WT_subs



#Insertions
median_wt_aa_ins <- SOD1_insertions %>% 
  group_by(Property_wt, mutation_type) %>% 
  summarise(median_abundance = median(abundance_score))

median_wt_aa_ins <- median_wt_aa_ins %>%
  mutate(Property_wt = factor(Property_wt, levels = median_wt_aa_ins %>% arrange(desc(median_abundance)) %>% pull(Property_wt)))

p_median_WT_ins <- ggplot(median_wt_aa_ins, aes(x = median_abundance, y = Property_wt)) +
  geom_col(fill = "lightgrey", color = "black", linewidth = 1) +
  theme_classic() +
  labs(x = "Insertions median\nabundance score", y = "") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 26))
p_median_WT_ins


#Deletions
#it is the same as in mutant

p_mut_all <- ggarrange(p_median_mut_subs,
                       p_median_mut_ins,
                       ncol = 2)
p_mut_all

p_wt_all <- ggarrange(p_median_WT_subs,
                      p_median_WT_ins,
                      p_median_mut_del,
                      ncol = 3)
p_wt_all


ggsave(p_mut_all, file = "median_mut_impact_subs_ins.tiff", width = 16, height = 8, dpi = 300)
ggsave(p_wt_all, file = "median_wt_impact_subs_ins_del.tiff", width = 22, height = 8, dpi = 300)
