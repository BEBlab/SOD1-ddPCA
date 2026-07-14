#required packages
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(ggbeeswarm)
library(ggalluvial)


#load modeled data
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")
SOD1_insertions <-read.csv("SOD1_ins_modelled.csv")
SOD1_deletions <- read.csv("SOD1_del_modelled.csv")

SOD1_allmut <- bind_rows(SOD1_subsitutions, SOD1_insertions, SOD1_deletions)
SOD1_allmut <- na.omit(SOD1_allmut)


SOD1_subsitutions <- SOD1_allmut %>% 
  filter(mutation_type == "subs")

 SOD1_subsitutions <- SOD1_subsitutions %>% 
   mutate(functional_impact = case_when(
     category_fdr_apca %in% c("high abundance", "WT-like abundance") & category_fdr_bpca %in% c("low heterodimerization", "stop-like heterodimerization") ~ "only binding",
     category_fdr_apca %in% c("low abundance", "stop-like abundance") & category_fdr_bpca %in% c("high heterodimerization", "WT-like heterodimerization") ~ "only abundance",
     category_fdr_apca %in% c("high abundance", "WT-like abundance") & category_fdr_bpca %in% c("high heterodimerization", "WT-like heterodimerization") ~ "wt in both",
     TRUE ~ "decrease both"
   ))
 SOD1_subsitutions$functional_impact <- factor(SOD1_subsitutions$functional_impact,
                                        levels = c("wt in both", "only abundance", "only binding", "decrease both"))

 SOD1_allmut <- SOD1_allmut %>% 
   mutate(functional_impact = case_when(
     category_fdr_apca %in% c("high abundance", "WT-like abundance") & category_fdr_bpca %in% c("low heterodimerization", "stop-like heterodimerization") ~ "only binding",
     category_fdr_apca %in% c("low abundance", "stop-like abundance") & category_fdr_bpca %in% c("high heterodimerization", "WT-like heterodimerization") ~ "only abundance",
     category_fdr_apca %in% c("high abundance", "WT-like abundance") & category_fdr_bpca %in% c("high heterodimerization", "WT-like heterodimerization") ~ "wt in both",
     TRUE ~ "decrease both"
   ))
 SOD1_allmut$functional_impact <- factor(SOD1_allmut$functional_impact,
                                               levels = c("wt in both", "only abundance", "only binding", "decrease both"))

#load clinvar and gnomad data
clinvar <- read.delim("clinvar_result (6).txt")
clinvar <- clinvar %>%
  mutate(ID = str_replace(Protein.change, 
                          "(?<=^[A-Z])\\d+(?=[A-Z]$)", 
                          function(x) as.character(as.numeric(x) - 1)))
clinvar <- clinvar %>%
  rename(Classification = Germline.classification) %>%
  select(ID, Classification)

clinvar$ID[64] <- "F64S"
clinvar_as_it_is <- clinvar

clinvar <- clinvar %>%
  mutate(Classification = recode(Classification,
                                 "Pathogenic/Likely pathogenic" = "Pathogenic",
                                 "Likely pathogenic" = "Pathogenic")) %>% 
  filter(ID != "G82R")

gnomad <- read.csv("gnomad_SOD1.csv", sep = ";")

aa_dict <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
  Glu="E", Gln="Q", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
  Ter="*", X="*")

aa_short <- c(
  A="A", R="R", N="N", D="D", C="C", E="E", Q="Q", G="G",
  H="H", I="I", L="L", K="K", M="M", F="F", P="P",
  S="S", T="T", W="W", Y="Y", V="V", "*"="*")

gnomad$Protein.Consequence <- sapply(gnomad$Protein.Consequence, convert_hgvsp_to_id, USE.NAMES = FALSE)
gnomad$Protein.Consequence <- sapply(gnomad$Protein.Consequence, decrement_id_pos, USE.NAMES = FALSE)
gnomad <- select(gnomad, Protein.Consequence, Allele.Frequency, ClinVar.Germline.Classification)
gnomad <- rename(gnomad, ID = Protein.Consequence, af = Allele.Frequency, Classification = ClinVar.Germline.Classification)
gnomad <- gnomad %>%
  mutate(Classification = case_when(
    ID == "E133d" ~ "Pathogenic/Likely pathogenic",
    TRUE ~ Classification
  ))

gnomad_patho <- gnomad %>% 
  filter(Classification == "Pathogenic/Likely pathogenic")

gnomad_clinvar <- full_join(gnomad, clinvar, by = "ID")
gnomad_clinvar_as_it_is <- gnomad_clinvar

gnomad_clinvar <- rename(gnomad_clinvar, Classification_gnomad = Classification.x, Classification = Classification.y)
gnomad_clinvar <- gnomad_clinvar %>%
  mutate(Classification = case_when(
    ID == "E133d" ~ "Pathogenic",
    ID == "I18d" ~ "Uncertain significance",
    TRUE ~ Classification
  ))

gnomad_clinvar_as_it_is <- rename(gnomad_clinvar_as_it_is, Classification_gnomad = Classification.x, Classification_clinvar = Classification.y)
gnomad_clinvar_as_it_is <- gnomad_clinvar_as_it_is %>%
  mutate(Classification_clinvar = case_when(
    ID == "E133d" ~ "Pathogenic",
    ID == "I18d" ~ "Uncertain significance",
    TRUE ~ Classification_clinvar
  ))

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

SOD1_allmut$category_fdr_apca <- factor(SOD1_allmut$category_fdr_apca,
                                levels = c("high", "WT-like", "low", "stop-like"))
SOD1_allmut$category_fdr_bpca <- factor(SOD1_allmut$category_fdr_bpca,
                                levels = c("high", "WT-like", "low", "stop-like"))

clinvar_table <- inner_join(SOD1_allmut, clinvar_as_it_is, by = "ID")
clinvar_table <- clinvar_table %>% 
  select(ID, Classification, category_fdr_apca, category_fdr_bpca) %>% 
  rename(clinical_classification = Classification, abundance_FDR = category_fdr_apca, heterodimerization_FDR = category_fdr_bpca)

write.csv(clinvar_table, file = "Dataset_S5.csv")

#join DMS data with clinvar data
SOD1_clinvar <- inner_join(SOD1_allmut, gnomad_clinvar, by = "ID")

pathogenic  <- SOD1_clinvar %>% 
  filter(Classification == "Pathogenic") %>% 
  distinct(ID, .keep_all = TRUE)

pathogenic_low_res <- pathogenic %>% 
  filter(residual_intensity == "low")

p_lluvial_fdr_patho <- ggplot(pathogenic,
                                 aes(axis1 = category_fdr_apca,
                                     axis2 = category_fdr_bpca,
                                     y = 1)) +
  geom_alluvium(aes(fill = category_fdr_apca), width = 1/12) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
  scale_x_discrete(limits = c("Abundance", "Hetero-\ndimerization"), expand = c(.05, .05)) +
  scale_fill_manual(
    values = c(
      "high" = "darkgreen",
      "WT-like"   = "lightgrey",
      "low"   = "darkorange",
      "stop-like"    = "#DD4636")) +
  labs(y = "Number of mutations", title = "P/LP") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title  = element_text(size = 26, face = "bold", hjust = 0.5),
        legend.position = "none")

p_lluvial_fdr_patho
ggsave(p_lluvial_fdr_patho, file = "alluvial_pathogenic.tiff", width = 8, height = 5)

vus  <- SOD1_clinvar %>% 
  filter(Classification == "Uncertain significance") %>% 
  distinct(ID, .keep_all = TRUE)


p_lluvial_fdr_vus <- ggplot(vus,
                               aes(axis1 = category_fdr_apca,
                                   axis2 = category_fdr_bpca,
                                   y = 1)) +
  geom_alluvium(aes(fill = category_fdr_apca), width = 1/12) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
  scale_x_discrete(limits = c("Abundance", "Hetero-\ndimerization"), expand = c(.05, .05)) +
  scale_fill_manual(
    values = c(
      "high" = "darkgreen",
      "WT-like"   = "lightgrey",
      "low"   = "darkorange",
      "stop-like"    = "#DD4636")) +
  labs(y = "Number of mutations", title = "VUS") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title  = element_text(size = 26, face = "bold", hjust = 0.5),
        legend.position = "none")

p_lluvial_fdr_vus

ggsave(p_lluvial_fdr_vus, file = "alluvial_vus.tiff", width = 8, height = 5)

#Violin plot of clinvar variants####
SOD1_long <- SOD1_allmut %>%
  pivot_longer(
    cols = c(`abundance_score`, `binding_score`),
    names_to = "assay",
    values_to = "score"
  )
apca <- SOD1_long %>% 
  filter(assay == "abundance_score")
apca$Classification <- "All aPCA substitutions"

bpca <- SOD1_long %>% 
  filter(assay == "binding_score")
bpca$Classification <- "All bPCA substitutions"

del_apca <- apca %>%
  filter(ID %in% c("E133d", "I18d")) %>%
  mutate(Classification = case_when(
    ID == "E133d" ~ "Pathogenic deletion aPCA",
    ID == "I18d" ~ "Uncertain significance deletion aPCA"
  ))

del_bpca <- bpca %>%
  filter(ID %in% c("E133d", "I18d")) %>%
  mutate(Classification = case_when(
    ID == "E133d" ~ "Pathogenic deletion bPCA",
    ID == "I18d" ~ "Uncertain significance deletion bPCA"
  ))

apca_clinvar <- SOD1_clinvar_long %>% 
  filter(assay == "abundance_score")

apca_clinvar <- apca_clinvar %>%
  mutate(Classification = recode(Classification,
                                   "Pathogenic" = "Pathogenic/Likely pathogenic aPCA",
                                   "Conflicting classifications of pathogenicity" = "Conflict. classif. aPCA",
                                   "Uncertain significance" = "Uncertain significance aPCA"))


bpca_clinvar <- SOD1_clinvar_long %>% 
  filter(assay == "binding_score")

bpca_clinvar <- bpca_clinvar %>%
  mutate(Classification = recode(Classification,
                                 "Pathogenic" = "Pathogenic/Likely pathogenic bPCA",
                                 "Conflicting classifications of pathogenicity" = "Conflict. classif. bPCA",
                                 "Uncertain significance" = "Uncertain significance bPCA"))
apca_clinvar <- apca_clinvar %>%
  filter(!ID %in% c("E133d", "I18d"))

bpca_clinvar <- bpca_clinvar %>%
  filter(!ID %in% c("E133d", "I18d"))


SOD1_violin <- bind_rows(
  apca, bpca,
  apca_clinvar, bpca_clinvar,
  del_apca, del_bpca
)
SOD1_violin$Classification <- factor(
  SOD1_violin$Classification,
  levels = c(
    "All aPCA substitutions",
    "All bPCA substitutions",
    "Pathogenic/Likely pathogenic aPCA",
    "Pathogenic/Likely pathogenic bPCA",
    "Uncertain significance aPCA",
    "Uncertain significance bPCA",
    "Conflict. classif. aPCA",
    "Conflict. classif. bPCA",
    "Pathogenic deletion aPCA",
    "Pathogenic deletion bPCA",
    "Uncertain significance deletion aPCA",
    "Uncertain significance deletion bPCA")
)

SOD1_violin_patho <- SOD1_violin %>% 
  filter(Classification %in% c("All aPCA substitutions", "All bPCA substitutions",
                               "Pathogenic/Likely pathogenic aPCA", "Pathogenic/Likely pathogenic bPCA",
                               "Pathogenic deletion aPCA", "Pathogenic deletion bPCA"))


p_violin_patho <- ggplot(SOD1_violin_patho, aes(x = Classification, y = score, alpha = Classification)) +
  geom_violin(scale = "width", trim = TRUE, aes(color = Classification), size = 1) +
  geom_jitter(
    data = SOD1_violin_patho %>% filter(!grepl("deletion", Classification)),
    aes(fill = Classification, color = Classification),
    shape = 21, width = 0.1, size = 4) +
  geom_jitter(
    data = SOD1_violin_patho %>% filter(grepl("deletion", Classification)),
    aes(fill = Classification, color = Classification),
    shape = 21, width = 0.1, size = 7) +
  scale_alpha_manual(values = c("All aPCA substitutions" = 0.05, "All bPCA substitutions" = 0.05, "Pathogenic/Likely pathogenic aPCA"  = 0.7, "Pathogenic/Likely pathogenic bPCA" = 0.7)) +
  scale_fill_manual(values = c(
    "All aPCA substitutions" = "#0A1829",
    "All bPCA substitutions" = "#E9AF1D",
    "Pathogenic/Likely pathogenic aPCA" = "#0A1829",
    "Pathogenic/Likely pathogenic bPCA" = "#E9AF1D",
    "Pathogenic deletion aPCA" = "#0A1829",
    "Pathogenic deletion bPCA" = "#E9AF1D"
  ))+
  scale_color_manual(values = c(
    "All aPCA substitutions" = "#0A1829",
    "All bPCA substitutions" = "#E9AF1D",
    "Pathogenic/Likely pathogenic aPCA" = "#0A1829",
    "Pathogenic/Likely pathogenic bPCA" = "#E9AF1D",
    "Pathogenic deletion aPCA" = "#0A1829",
    "Pathogenic deletion bPCA" = "#E9AF1D"
  )) +
  geom_line(data = SOD1_violin_patho %>% 
              filter(Classification %in% c("Pathogenic/Likely pathogenic aPCA",
                                           "Pathogenic/Likely pathogenic bPCA",
                                           "Pathogenic deletion aPCA",
                                           "Pathogenic deletion bPCA")),
            aes(group = ID),
            color = "lightgrey", alpha = 0.8) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  theme_classic(base_size = 14, base_line_size = 0.2) +
  labs(y = "Score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 30),
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
  )
p_violin_patho

ggsave(p_violin_patho, file="p_violin_pathogenic.tiff", width=10, height = 6, dpi = 300)

SOD1_violin_vus <- SOD1_violin %>% 
  filter(Classification %in% c("All aPCA substitutions", "All bPCA substitutions",
                               "Uncertain significance aPCA", "Uncertain significance bPCA",
                               "Uncertain significance deletion aPCA", "Uncertain significance deletion bPCA"))

p_violin_vus <- ggplot(SOD1_violin_vus, aes(x = Classification, y = score, alpha = Classification)) +
  geom_violin(scale = "width", trim = TRUE, aes(color = Classification), size = 1) +
  geom_jitter(
    data = SOD1_violin_vus %>% filter(!grepl("deletion", Classification)),
    aes(fill = Classification, color = Classification),
    shape = 21, width = 0.1, size = 4) +
  geom_jitter(
    data = SOD1_violin_vus %>% filter(grepl("deletion", Classification)),
    aes(fill = Classification, color = Classification),
    shape = 21, width = 0.1, size = 7) + 
  scale_alpha_manual(values = c("All aPCA substitutions" = 0.05, "All bPCA substitutions" = 0.05, "Uncertain significance aPCA"  = 0.7, "Uncertain significance bPCA" = 0.7)) +
  scale_fill_manual(values = c(
    "All aPCA substitutions" = "#0A1829",
    "All bPCA substitutions" = "#E9AF1D",
    "Uncertain significance aPCA" = "#0A1829",
    "Uncertain significance bPCA" = "#E9AF1D",
    "Uncertain significance deletion aPCA" = "#0A1829",
    "Uncertain significance deletion bPCA" = "#E9AF1D"
  ))+
  scale_color_manual(values = c(
    "All aPCA substitutions" = "#0A1829",
    "All bPCA substitutions" = "#E9AF1D",
    "Uncertain significance aPCA" = "#0A1829",
    "Uncertain significance bPCA" = "#E9AF1D",
    "Uncertain significance deletion aPCA" = "#0A1829",
    "Uncertain significance deletion bPCA" = "#E9AF1D"
  )) +
  geom_line(data = SOD1_violin_vus %>% 
              filter(Classification %in% c("Uncertain significance aPCA",
                                           "Uncertain significance bPCA",
                                           "Uncertain significance deletion aPCA",
                                           "Uncertain significance deletion bPCA")),
            aes(group = ID),
            color = "lightgrey", alpha = 0.8) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  theme_classic(base_size = 14, base_line_size = 0.2) +
  labs(y = "Score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 30),
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
  )


p_violin_vus
ggsave(p_violin_vus, file="p_violin_vus.tiff", width=10, height = 6, dpi = 300)


SOD1_violin_conf <- SOD1_violin %>% 
  filter(Classification %in% c("All aPCA substitutions", "All bPCA substitutions",
                               "Conflict. classif. aPCA", "Conflict. classif. bPCA"))

p_violin_conf <- ggplot(SOD1_violin_conf, aes(x = Classification, y = score)) +
  geom_violin(scale = "width", trim = TRUE, aes(color = Classification), size = 1) +
  geom_jitter(aes(fill = Classification, color = Classification), shape = 21, width = 0.1, size = 4, alpha = 0.6) +
  scale_fill_manual(values = c(
    "All aPCA substitutions" = "#0A1829",
    "All bPCA substitutions" = "#E9AF1D",
    "Conflict. classif. aPCA" = "#0A1829",
    "Conflict. classif. bPCA" = "#E9AF1D"
  ))+
  scale_color_manual(values = c(
    "All aPCA substitutions" = "#0A1829",
    "All bPCA substitutions" = "#E9AF1D",
    "Conflict. classif. aPCA" = "#0A1829",
    "Conflict. classif. bPCA" = "#E9AF1D"
  )) +
  geom_line(data = SOD1_violin_conf %>% 
              filter(Classification %in% c("Conflict. classif. aPCA",
                                           "Conflict. classif. bPCA")),
            aes(group = ID),
            color = "lightgrey", alpha = 0.8) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  theme_classic(base_size = 14, base_line_size = 0.2) +
  labs(y = "Normalized score") +
  theme(axis.title = element_text(size = 33),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 30),
        axis.line = element_line(size = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
  )


p_violin_conf

ggsave(p_violin_conf, file="p_violin_conf.tiff", width=10, height = 6, dpi = 300)


#Correlation between assays of clinvar variants####
r_labels_pathogenic <- pathogenic %>%
  summarise(
    R = cor(abundance_score, binding_score, method = "spearman"),
    p_value = cor.test(abundance_score, binding_score, method = "spearman")$p.value
  ) %>%
  mutate(
    label = paste0("rho = ", round(R, 2), "\n", "p < 2.2e-16")
  )

p_cor_clinvar <- ggplot(pathogenic, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 6, alpha = 1, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score", y = "Heterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_labels_pathogenic, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_cor_clinvar

r_labels_vus <- vus %>%
  summarise(
    R = cor(abundance_score, binding_score, method = "spearman"),
    p_value = cor.test(abundance_score, binding_score, method = "spearman")$p.value
  ) %>%
  mutate(
    label = paste0("rho = ", round(R, 2), "\n", "p < 2.2e-16")
  )

p_cor_clinvar_vus <- ggplot(vus, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 6, alpha = 1, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score", y = "Heterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_labels_vus, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_cor_clinvar_vus

conf_variants_scores <- SOD1_clinvar %>% 
  filter(Classification == "Conflicting classifications of pathogenicity")

r_labels_conf <- conf_variants_scores %>%
  summarise(
    R = cor(abundance_score, binding_score, method = "spearman"),
    p_value = cor.test(abundance_score, binding_score, method = "spearman")$p.value
  ) %>%
  mutate(
    label = paste0("rho = ", round(R, 2), "\n", "p < 2.2e-16 ")
  )

p_cor_clinvar_conf <- ggplot(conf_variants_scores, aes(x = abundance_score, y = binding_score)) +
  geom_point(size = 6, alpha = 1, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  labs(x = "Abundance score", y = "Heterodimerization score") +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30)) +
  geom_text(data = r_labels_conf, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1, size = 10, inherit.aes = FALSE)  
p_cor_clinvar_conf

p_cor_all <- ggarrange(p_cor_clinvar,
                       p_cor_clinvar_vus,
                       p_cor_clinvar_conf,
                       ncol = 3)
p_cor_all
ggsave(p_cor_all, file = "p_corr_clinvar_variants_scores.tiff", width = 23, height = 7)

#FDR and residual of clinvar variants####
pathogenic_counts <- pathogenic %>% 
  group_by(category_fdr_apca, category_fdr_bpca, residual_category) %>% 
  dplyr::count()

pathogenic_counts <- pathogenic_counts %>%
  pivot_longer(
    cols = c(category_fdr_apca, category_fdr_bpca, residual_category),
    names_to = "category",
    values_to = "value"
  ) %>%
  group_by(category, value) %>%
  summarise(n = sum(n), .groups = "drop")


pathogenic_counts_fdr <- pathogenic_counts %>%
  filter(category %in% c("category_fdr_apca", "category_fdr_bpca"))

pathogenic_counts_fdr$value <- factor(pathogenic_counts_fdr$value,
                                         levels = c("high", "WT-like", "low", "stop-like"))

pathogenic_counts_residual <- pathogenic_counts %>%
  filter(category == "residual_category")

order_levels <- c("residual_category", "category_fdr_bpca", "category_fdr_apca")

pathogenic_counts_fdr$category <- factor(
  pathogenic_counts_fdr$category,
  levels = order_levels
)

pathogenic_counts_residual$category <- factor(
  pathogenic_counts_residual$category,
  levels = order_levels
)

p_pathogenic_categories <- ggplot() +
  geom_bar(data = pathogenic_counts_fdr, aes(x = category, y = n, fill = value),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(name = "Abundance / Binding",
                    values = c("high" = "darkgreen",
                               "WT-like" = "lightgrey",
                               "low" = "darkorange",
                               "stop-like" = "#DD4636")) +
  ggnewscale::new_scale_fill() +
  geom_bar(data = pathogenic_counts_residual, aes(x = category, y = n, fill = value),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  scale_x_discrete(
    limits = order_levels,
    labels = c("category_fdr_apca" = "Abundance",
               "category_fdr_bpca" = "Heterodimerization",
               "residual_category" = "Heterodimerization\nresidual")
  ) +
  labs(x = "", y = "Number of mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(hjust = 0.5, vjust = 1),
        legend.position = "none") +
  coord_flip()
p_pathogenic_categories
ggsave(p_pathogenic_categories, file = "pathogenic_barplot.tiff", width = 10, height = 4)

pathogenic_counts_loc <- pathogenic %>%
  mutate(
    side_chain = case_when(
      side_chain %in% c("core", "surface") ~ "rest",
      TRUE ~ side_chain
    )
  ) %>%
  dplyr::count(residual_category, side_chain)

pathogenic_counts_loc$residual_category <- factor(pathogenic_counts_loc$residual_category,
                                      levels = c("high_negative", "high_positive", "low"))


order_levels <- c("rest", "second_shell", "dimer interface")

pathogenic_counts_loc$side_chain <- factor(
  pathogenic_counts_loc$side_chain,
  levels = order_levels
)

p_pathogenic_loc <- ggplot() +
  geom_bar(data = pathogenic_counts_loc, aes(x = side_chain, y = n, fill = residual_category),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  scale_x_discrete(
    limits = order_levels,
    labels = c("rest" = "rest", "second_shell" = "2nd shell", "dimer interface" = "dimer\ninterface")) +
  labs(x = "", y = "Number of mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(hjust = 0.5, vjust = 1),
        legend.position = "none") +
  coord_flip()
p_pathogenic_loc
ggsave(p_pathogenic_loc, file = "pathogenic_barplot_loc_byresidual.tiff", width = 10, height = 4)

p_pathogenic_bar <- ggarrange(p_pathogenic_categories,
                              p_pathogenic_loc,
                              nrow = 2, align = "hv")
p_pathogenic_bar
ggsave(p_pathogenic_bar, file = "pathogenic_barplot.tiff", width = 10, height = 8)

vus_counts <- vus %>% 
  group_by(category_fdr_apca, category_fdr_bpca, residual_category) %>% 
  dplyr::count()

vus_counts <- vus_counts %>%
  pivot_longer(
    cols = c(category_fdr_apca, category_fdr_bpca, residual_category),
    names_to = "category",
    values_to = "value"
  ) %>%
  group_by(category, value) %>%
  summarise(n = sum(n), .groups = "drop")

vus_counts_fdr <- vus_counts %>%
  filter(category %in% c("category_fdr_apca", "category_fdr_bpca"))
vus_counts_fdr$value <- factor(vus_counts_fdr$value,
                                      levels = c("high", "WT-like", "low", "stop-like"))

vus_counts_residual <- vus_counts %>%
  filter(category == "residual_category")

order_levels <- c("residual_category", "category_fdr_bpca", "category_fdr_apca")


vus_counts_fdr$category <- factor(
  vus_counts_fdr$category,
  levels = order_levels
)

vus_counts_residual$category <- factor(
  vus_counts_residual$category,
  levels = order_levels
)

p_vus_categories <- ggplot() +
  geom_bar(data = vus_counts_fdr, aes(x = category, y = n, fill = value),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(name = "Abundance / Binding",
                    values = c("high" = "darkgreen",
                               "WT-like" = "lightgrey",
                               "low" = "darkorange",
                               "stop-like" = "#DD4636")) +
  ggnewscale::new_scale_fill() +
  geom_bar(data = vus_counts_residual, aes(x = category, y = n, fill = value),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  scale_x_discrete(
    limits = order_levels,
    labels = c("category_fdr_apca" = "Abundance",
               "category_fdr_bpca" = "Heterodimerization",
               "residual_category" = "Heterodimerization\nresidual")
  ) +
  labs(x = "", y = "Number of mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(hjust = 0.5, vjust = 1),
        legend.position = "none") +
  coord_flip()
p_vus_categories
ggsave(p_vus_categories, file = "vus_barplot.tiff", width = 10, height = 4)


vus_counts_loc <- vus %>%
  mutate(
    side_chain = case_when(
      side_chain %in% c("core", "surface") ~ "rest",
      TRUE ~ side_chain
    )
  ) %>%
  dplyr::count(residual_category, side_chain)


vus_counts_loc$residual_category <- factor(vus_counts_loc$residual_category,
                                                  levels = c("high_negative", "high_positive", "low"))


order_levels <- c("rest", "second_shell", "dimer interface")

vus_counts_loc$side_chain <- factor(
  vus_counts_loc$side_chain,
  levels = order_levels
)

p_vus_loc <- ggplot() +
  geom_bar(data = vus_counts_loc, aes(x = side_chain, y = n, fill = residual_category),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  scale_x_discrete(
    limits = order_levels,
    labels = c("rest" = "rest", "second_shell" = "2nd shell", "dimer interface" = "dimer\ninterface")) +
  labs(x = "", y = "Number of mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(hjust = 0.5, vjust = 1),
        legend.position = "none") +
  coord_flip()
p_vus_loc

ggsave(p_vus_loc, file = "vus_barplot_loc_byresidual.tiff", width = 10, height = 4)

p_vus_bar <- ggarrange(p_vus_categories,
                              p_vus_loc,
                              nrow = 2, align = "hv")
p_vus_bar
ggsave(p_vus_bar, file = "vus_barplot.tiff", width = 10, height = 8)

conf  <- SOD1_clinvar %>% 
  filter(Classification == "Conflicting classifications of pathogenicity") %>% 
  distinct(ID, .keep_all = TRUE)

conf_counts <- conf %>% 
  group_by(category_fdr_apca, category_fdr_bpca, residual_category) %>% 
  dplyr::count()

conf_counts <- conf_counts %>%
  pivot_longer(
    cols = c(category_fdr_apca, category_fdr_bpca, residual_category),
    names_to = "category",
    values_to = "value"
  ) %>%
  group_by(category, value) %>%
  summarise(n = sum(n), .groups = "drop")

conf_counts_fdr <- conf_counts %>%
  filter(category %in% c("category_fdr_apca", "category_fdr_bpca"))

conf_counts_fdr$value <- factor(conf_counts_fdr$value,
                                 levels = c("high", "WT-like", "low", "stop-like"))

conf_counts_residual <- conf_counts %>%
  filter(category == "residual_category")

order_levels <- c("residual_category", "category_fdr_bpca", "category_fdr_apca")


conf_counts_fdr$category <- factor(
  conf_counts_fdr$category,
  levels = order_levels
)

conf_counts_residual$category <- factor(
  conf_counts_residual$category,
  levels = order_levels
)


p_conf_categories <- ggplot() +
  geom_bar(data = conf_counts_fdr, aes(x = category, y = n, fill = value),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(name = "Abundance / Binding",
                    values = c("high" = "darkgreen",
                               "WT-like" = "lightgrey",
                               "low" = "darkorange",
                               "stop-like" = "#DD4636")) +
  ggnewscale::new_scale_fill() +
  geom_bar(data = conf_counts_residual, aes(x = category, y = n, fill = value),
           stat = "identity", position = "stack", color = "black", width = 0.6) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  scale_x_discrete(
    limits = order_levels,
    labels = c("category_fdr_apca" = "Abundance",
               "category_fdr_bpca" = "Heterodimerization",
               "residual_category" = "Heterodimerization\nresidual")
  ) +
  labs(x = "", y = "Number of mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(hjust = 0.5, vjust = 1),
        legend.position = "none") +
  coord_flip()
p_conf_categories
ggsave(p_conf_categories, file = "conf_barplot.tiff", width = 10, height = 4)


#Pathogenic WT-like mutations####
wt_patho <- SOD1_clinvar_long %>% 
  filter(category_fdr_apca == "WT-like" & category_fdr_bpca == "WT-like" & Classification == "Pathogenic")

wt_patho <- wt_patho %>% 
  filter(category_fdr_apca != "high" & category_fdr_bpca != "high")

SOD1_subs_long <- SOD1_subsitutions %>%
  pivot_longer(
    cols = c(`abundance_score`, `binding_score`),
    names_to = "assay",
    values_to = "score"
  )

wt_patho_median_score <- SOD1_subs_long %>% 
  group_by(Pos, assay) %>% 
  summarise(
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  )
wt_patho_median_score <- inner_join(wt_patho, wt_patho_median_score, by = "Pos")

wt_patho_median_score <- wt_patho_median_score %>% 
  mutate(
    ID_num = as.numeric(str_extract(ID, "\\d+")),
    ID = fct_reorder(ID, ID_num)
  )


p_median_score <- ggplot(wt_patho_median_score, aes(x = ID, y = median_score, color = assay.y, group = ID)) +
  geom_line(color = "lightgrey", linewidth = 0.5) +
  geom_point(size = 7) +
  scale_color_manual(values = c("abundance_score" = "#0A1829", "binding_score" = "#E9AF1D"),
                     labels = c("abundance_score" = "Abundance", "binding_score" = "Binding")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  geom_hline(yintercept = -0.19, linetype = "dashed", color = "#0A1829", size = 1) +
  geom_hline(yintercept = -0.05, linetype = "dashed", color = "#E9AF1D", size = 1) +
  theme_classic() +
  labs(x = "P/LP WT-like mutation", y = "Median score at position\nof WT-like mutation", color = "") +
  theme(axis.title = element_text(size = 28),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 25),
        legend.position = "top",
        legend.text = element_text(size = 22)) +
  coord_flip()

p_median_score

#load predictors
alphamiss <- read.csv("am_scores.csv", sep = ";")
revel <- read.csv("revel_fullset.csv", sep = ",")
eve <- read.csv("eve_all.csv", sep = ";")
patho_info <- read.csv("Pathogenic_WT-like_variants.csv", sep = ";")
clinical_phenotypes <- read.csv("SOD1_clinical_phenotypes_dms.csv")

alphamiss <- alphamiss %>% 
  mutate(am_class = case_when(
    am_class == "Amb" ~ "VUS",
    am_class == "LPath" ~ "Pathogenic",
    am_class == "LBen" ~ "Benign")) %>% 
  rename(score = pathogenicity_score, classification = am_class) 
  

revel <- revel %>% 
  rename(score = REVEL) %>% 
  mutate(classification  = case_when(
    score >= 0.7 ~ "Pathogenic",
    score < 0.4 ~ "Benign",
    T ~ "VUS"
  )) %>% 
  select(ID, score, classification)

eve$ID <- paste0(eve$wt_aa, eve$position, eve$mt_aa)
eve <- eve %>% 
  rename(classification = EVE_classes_75_pct_retained_ASM, score = EVE_scores_ASM) %>% 
  select(ID, classification,score) %>% 
  mutate(classification = recode(classification,
                                 "Uncertain" = "VUS"))

alphamiss$predictor <- "alphamissense"
revel$predictor <- "revel"
eve$predictor <- "eve"

predictors <- rbind(alphamiss, revel, eve)
predictors$ID <- sapply(predictors$ID, decrement_id_pos)

penetrance <- clinical_phenotypes %>% 
  select(ID, penetrance) %>%
  mutate(penetrance = recode(penetrance, "high/incomplete" = "complete + high")) %>% 
  rename(category = penetrance) %>% 
  mutate(phenotype = "penetrance")


survival_time <- clinical_phenotypes %>% 
  select(ID, survival_time) %>% 
  rename(category = survival_time) %>% 
  mutate(phenotype = "survival_time")

progression <- clinical_phenotypes %>% 
  select(ID, progresion) %>% 
  rename(category = progresion) %>% 
  mutate(phenotype = "progression")

clinical_phenotypes <- rbind(penetrance, survival_time, progression)


patho_info <- patho_info %>% 
  mutate(Clinvar.classification = recode(Clinvar.classification,
                                         "Pathogenic" = "P",
                                         "Pathogenic/likely pathogenic" = "P/LP",
                                         "Likely pathogenic" = "LP"
                                         ))

wt_patho_median_score <- left_join(wt_patho_median_score, patho_info, by = "ID")
wt_patho_median_score <- left_join(wt_patho_median_score, gnomad, by = "ID")
wt_patho_median_score <- left_join(wt_patho_median_score, predictors, by = "ID")
wt_patho_median_score <- left_join(wt_patho_median_score, clinical_phenotypes, by = "ID")

wt_patho_median_score <- wt_patho_median_score %>% 
  mutate(Type_ALS = recode(Type_ALS,
                                     "familial" = "F",
                                     "familial/sporadic" = "F/S",
                                    "sporadic" = "S")) %>% 
  mutate(s_location = recode(s_location,
                           "Electrostatic loop" = "Elect.loop",
                           "Zn binding loop" = "Zn loop"))
  
wt_patho_median_score <- wt_patho_median_score %>% 
  mutate(ID_num = as.numeric(str_extract(ID, "\\d+")),
         ID = fct_reorder(ID, ID_num))



p_predictors <- ggplot(wt_patho_median_score, aes(x = predictor, y = ID, fill = classification)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Pathogenic" = "#D55E00", "VUS" = "darkblue", "Benign"= "#009E73"), na.value = "darkgrey") +
  theme_void() +
  theme(legend.position = "none")
p_predictors



p_clinical_phenotypes <- ggplot(wt_patho_median_score, aes(x = phenotype, y = ID, fill = category)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("slow"= "#2F6790", "incomplete" = "#2F6790", "long" = "#2F6790", "complete + high" = "#76669B", "short" = "#76669B"), na.value = "white") +
  theme_void() +
  scale_x_discrete(na.translate = FALSE) +
  theme(legend.position = "none")
p_clinical_phenotypes


p_als_type <- ggplot(wt_patho_median_score %>% distinct(ID, Type_ALS),
aes(x = 1, y = ID, label = Type_ALS)) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)
p_als_type

p_clin_text <- ggplot(
  wt_patho_median_score %>% distinct(ID, Clinvar.classification),
  aes(x = 1, y = ID, label = Clinvar.classification)
) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)
p_clin_text

wt_patho_median_score <- wt_patho_median_score %>%
  mutate(
    Clinvar.stars = as.integer(Clinvar.stars),        
    Clinvar.stars = pmax(Clinvar.stars, 0),          
    stars_label = strrep("★", Clinvar.stars)         
  )

p_stars_text <- ggplot(
  wt_patho_median_score %>% distinct(ID, Clinvar.stars),
  aes(
    x = 1,
    y = ID,
    label =
      strrep("★", Clinvar.stars))) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)

p_stars_text

p_af_text <- ggplot(
  wt_patho_median_score %>% distinct(ID, af.y),
  aes(
    x = 1,
    y = ID,
    label = signif(af.y, 2)
  )
) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)
p_af_text


p_loc <- ggplot(
  wt_patho_median_score %>% distinct(ID, s_location),
  aes(x = 1, y = ID, label = s_location)
) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)
p_loc

p_stars_text <- p_stars_text +
  theme(plot.margin = margin(5, -5, 5, 0))

p_af_text <- p_af_text +
  theme(plot.margin = margin(5, -5, 5, 0))


final_plot <- ggarrange(p_median_score, p_clin_text, p_stars_text, p_af_text, p_predictors, p_clinical_phenotypes, p_als_type, p_loc,
                        align = "h", nrow = 1,   widths = c(8, 1, 1, 1.4, 1, 1, 0.9, 2))
final_plot


ggsave(final_plot, file = "patho_wt_like_median_pos_effec.tiff", width = 20, height = 10)
ggsave(final_plot, file = "patho_wt_like_median_pos_effec_final.tiff", width = 20, height = 10)


count_fdr <- pathogenic %>% 
  mutate(fdr_both = case_when(
    category_fdr_apca == "WT-like" & 
      category_fdr_bpca == "WT-like" ~ "WT-like",
    TRUE ~ "non-WT-like"
  )) %>% 
  distinct(ID, .keep_all = TRUE)

counts_patho <- dplyr::count(count_fdr, fdr_both) %>% 
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%")
  )

pie_patho <- ggplot(counts_patho, aes(x = "", y = n, fill = fdr_both)) +
  geom_col(width = 1, color = "black", linewidth = 0.9) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "WT-like" = "lightgrey",
    "non-WT-like" = "darkorange"
    
  )) +
  geom_text(
    aes(label = paste0(label, "\n(n=", n, ")")),
    position = position_stack(vjust = 0.5),
    size = 11
  ) +
  theme_void() +
  theme(legend.position = "none")

pie_patho
ggsave(pie_patho, file = "pie_pathogenic_wtlike.tiff", width = 5, height = 5, dpi = 300)

#rSASA pathogenic####
sasa_patho <- count_fdr %>% 
  group_by(fdr_both) %>% 
  summarise(median_sasa = median(rSASA.dimer))

p_sasa_patho <- ggplot(sasa_patho, aes(x = fdr_both, y = median_sasa, fill = fdr_both)) +
  geom_col(linewidth = 0.9, color = "black", width = 0.5) +
  scale_fill_manual(values = c("WT-like" = "darkgrey", "non-WT-like" = "darkorange")) +
  theme_classic() +
  labs(x = "FDR = 0.1", y = "Median dimer rSASA") +
  theme(legend.position = "none",
        axis.title = element_text(size = 28),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 22))
p_sasa_patho

ggsave(p_sasa_patho, file = "barplot_sasa_pathogenic.tiff", width = 6,height = 6)

#VUS wt-like variants####
vus_wt <- SOD1_clinvar_long %>% 
  filter(category_fdr_apca == "WT-like" & category_fdr_bpca == "WT-like" & Classification == "Uncertain significance")

vus_median_score <- SOD1_subs_long %>% 
  group_by(Pos, assay) %>% 
  summarise(
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  )
vus_median_score <- inner_join(vus_wt, vus_median_score, by = "Pos")


vus_median_score$ID_sum <- sapply(vus_median_score$ID, sum_id_pos)
vus_median_score$hgvs_pro <- sapply(vus_median_score$ID_sum, id_to_hgvs_pro)


# vus_median_score <- vus_median_score %>% 
#   mutate(
#     Pos_num = as.numeric(str_extract(hgvs_pro, "\\d+")),
#     hgvs_pro = fct_reorder(hgvs_pro, Pos_num)
#   )
vus_median_score <- vus_median_score %>% 
  mutate(
    ID_num = as.numeric(str_extract(ID, "\\d+")),
    ID = fct_reorder(ID, ID_num)
  )

p_median_score_vus <- ggplot(vus_median_score, aes(x = ID, y = median_score, color = assay.y, group = ID)) +
  geom_line(color = "lightgrey", linewidth = 0.5) +
  geom_point(size = 7) +
  scale_color_manual(values = c("abundance_score" = "#0A1829", "binding_score" = "#E9AF1D"),
                     labels = c("abundance_score" = "Abundance", "binding_score" = "Binding")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  geom_hline(yintercept = -0.19, linetype = "dashed", color = "#0A1829", size = 1) +
  geom_hline(yintercept = -0.05, linetype = "dashed", color = "#E9AF1D", size = 1) +
  theme_classic() +
  labs(x = "VUS WT-like mutation", y = "Median score at position\nof WT-like mutation", color = "") +
  theme(axis.title = element_text(size = 28),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 25),
        legend.position = "top",
        legend.text = element_text(size = 22)) +
  coord_flip()

p_median_score_vus


#load predictors
vus_info <- read.csv("VUS_WT-like_variants.csv", sep = ";")

vus_median_score_plot <- left_join(vus_median_score, vus_info, by = "ID")
vus_median_score_plot <- left_join(vus_median_score_plot, gnomad, by = "ID")
vus_median_score_plot <- left_join(vus_median_score_plot, predictors, by = "ID")

vus_median_score_plot <- vus_median_score_plot %>% 
  mutate(
    ID_num = as.numeric(str_extract(ID, "\\d+")),
    ID = fct_reorder(ID, ID_num)
  )

p_predictors_vus <- ggplot(vus_median_score_plot, aes(x = predictor, y = ID, fill = classification)) +
  geom_tile(color = "black", linewidth = 1) +
  scale_fill_manual(values = c("Pathogenic" = "#D55E00", "VUS" = "darkblue", "Benign"= "#009E73"), na.value = "darkgrey") +
  theme_void() +
  theme(legend.position = "none")
p_predictors_vus

p_stars_text <- ggplot(
  vus_median_score_plot %>% distinct(ID, Clinvar.stars),
  aes(
    x = 1,
    y = ID,
    label =
      strrep("★", Clinvar.stars))) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)

p_stars_text

p_af_text <- ggplot(
  vus_median_score_plot %>% distinct(ID, af.x),
  aes(
    x = 1,
    y = ID,
    label = signif(af.x, 2)
  )
) +
  geom_text(size = 6, hjust = 0) +
  theme_void() +
  xlim(1, 2)
p_af_text

p_stars_text <- p_stars_text +
  theme(plot.margin = margin(5, -5, 5, 0))

p_af_text <- p_af_text +
  theme(plot.margin = margin(5, -5, 5, 0))


final_plot_vus <- ggarrange(p_median_score_vus, p_predictors_vus, p_stars_text, p_af_text,
                        align = "h", nrow = 1,   widths = c(6, 1, 0.8, 1))
final_plot_vus


ggsave(final_plot_vus, file = "vus_wt_like_median_pos_effec.tiff", width = 14, height = 10)
ggsave(final_plot_vus, file = "vus_wt_like_median_pos_effec_final.tiff", width = 14, height = 10)

count_fdr_vus <- vus %>% 
  mutate(fdr_both = case_when(
    category_fdr_apca == "WT-like" & 
    category_fdr_bpca == "WT-like" ~ "wt-like",
    TRUE ~ "non-wt-like"
  )) %>% 
  distinct(ID, .keep_all = TRUE)

counts_vus <- dplyr::count(count_fdr_vus, fdr_both) %>% 
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%")
  )

pie_vus <- ggplot(counts_vus, aes(x = "", y = n, fill = fdr_both)) +
  geom_col(width = 1, color = "black", linewidth = 0.9) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "wt-like" = "lightgrey",
    "non-wt-like" = "darkorange"
    
  )) +
  geom_text(
    aes(label = paste0(label, "\n(n=", n, ")")),
    position = position_stack(vjust = 0.5),
    size = 11
  ) +
  theme_void() +
  theme(legend.position = "none")

pie_vus
ggsave(pie_vus, file = "pie_vus_wtlike.tiff", width = 5, height = 5, dpi = 300)


  
#Concflicting wt-like####
count_fdr_conf <- SOD1_clinvar %>% 
  mutate(fdr_both = case_when(
    category_fdr_apca == "WT-like" & 
      category_fdr_bpca == "WT-like" ~ "wt-like",
    TRUE ~ "non-wt-like"
  )) %>% 
  filter(Classification == "Conflicting classifications of pathogenicity") %>% 
  distinct(ID, .keep_all = TRUE)

counts_conf <- dplyr::count(count_fdr_conf, fdr_both) %>% 
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%")
  )

pie_conf <- ggplot(counts_conf, aes(x = "", y = n, fill = fdr_both)) +
  geom_col(width = 1, color = "black", linewidth = 0.9) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "wt-like" = "lightgrey",
    "non-wt-like" = "darkorange"
  )) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 11
  ) +
  theme_void() +
  theme(legend.position = "none")

pie_conf
ggsave(pie_conf,  file = "pie_conf_wtlike.tiff", width = 5, height = 5, dpi = 300)
