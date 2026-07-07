library(ggbeeswarm)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

#load collected data from the literature
review <- read.csv("SOD1_variants_review.csv", sep = ";")
final_caco <- read.csv("Variantes_SOD1_final.csv", sep = ";")
als_type <- read.csv("type_ALS.csv", sep = ";")
penetrance <- read.csv("penetrance_variants_info.csv", sep = ";")

#load modeled data
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")

#organize data to work with only one data frame
als_type <- als_type %>% 
  filter(Gene == "SOD1")
als_type$Amino.acid.change.1 <- sapply(als_type$Amino.acid.change.1, convert_hgvsp_to_id)
als_type$Amino.acid.change.1 <- sapply(als_type$Amino.acid.change.1, decrement_id_pos)
als_type <- na.omit(als_type)
als_type <- als_type %>% 
  rename(ID = Amino.acid.change.1, als_type = FALS.SALS) %>% 
  select(ID, als_type)

als_type <- als_type %>%
  mutate(als_type = recode(als_type,
                           "SALS pronand" = "sALS",
                           "FALS and SALS probands (homozygous juvenile ALS J Med Genet. 1998 Feb;35(2):174.)" = "fALS/sALS",
                           "ALS proband" = "ALS",
                           "ALS probands" = "ALS",
                           "FALS" = "fALS",
                           "FALS and SALS" = "fALS/sALS",
                           "FALS and SALS (including recessive)" = "fALS/sALS",
                           "FALS and SALS (juvenile ALS)" = "fALS/sALS\n(juvenile)",
                           "FALS and SALS probands" = "fALS/sALS",
                           "FALS and SALS probands (homozygous juvenile)" = "fALS/sALS\n(homozygous juv.)",
                           "FALS proband" = "fALS",
                           "FALS probands" = "fALS",
                           "FALS(compound heterozygote)" = "fALS\n(compound hetzyg.)",
                           "SALS proband" = "sALS",
                           "SALS probands" = "sALS",
                           "SALS probands (including juvenile ALS)" = "sALS probands\n(juvenile)",
                           "FALS and SALS proband (including juvenile ALS)" = "fALS/sALS"))

final_caco <- final_caco %>% 
  mutate(Penetrancia = recode(Penetrancia,
                              "Alta" = "high",
                              "Baja" = "low", 
                              "Variable" = "variable")) %>% 
  mutate(Progresión = recode(Progresión,
                             "Rápida" = "fast",
                             "Lenta" = "slow")) %>% 
  mutate(Supervivencia.estimada = recode(Supervivencia.estimada,
                                         ">10 años" = "long",
                                         "<3 años" = "short",
                                         ">11 años" = "long")) %>% 
  rename(ID = Variante, penetrance = Penetrancia, progresion = Progresión, survival_time = Supervivencia.estimada) %>% 
  select(ID, penetrance, progresion, survival_time)


review$ID <- sapply(review$ID, convert_hgvsp_to_id)
review$ID <- sapply(review$ID, decrement_id_pos)

review <- review %>% 
  mutate(Penetrance = recode(Penetrance,
                             "Complete" = "complete + high",
                             "Greatly reduced" = "incomplete")) %>%
  mutate(Survival.time = recode(Survival.time,
                                "Long" = "long",
                                "Short" = "short",
                                "Longest (44-45 y)" = "long",
                                "Shortest (2-3 m)" = "short")) %>% 
  rename(survival_time = Survival.time, penetrance = Penetrance) %>% 
  select(ID, survival_time, penetrance)



penetrance <- penetrance %>% 
  mutate(ALS.type = recode(ALS.type,
                           "fALS and sALS" = "fALS/sALS",
                           "sALS and fALS" = "fALS/sALS")) %>% 
  mutate(Penetrance = recode(Penetrance,
                             "High" = "complete + high",
                             "Incomplete" = "incomplete",
                             "Incomplete/High" = "variable")) %>% 
  rename(als_type = ALS.type, penetrance = Penetrance) %>% 
  select(ID, als_type, penetrance)

sod1_phenotypes <- left_join(review, final_caco, by = "ID", suffix = c("_df1", "_df2"))
sod1_phenotypes <- sod1_phenotypes %>%
  mutate(
    survival_time = coalesce(survival_time_df1, survival_time_df2)
  ) %>% 
  select(ID, progresion, penetrance_df1, penetrance_df2, survival_time)

sod1_phenotypes <- left_join(sod1_phenotypes, penetrance, by = "ID")
sod1_phenotypes <- sod1_phenotypes %>%
  mutate(
    penetrance = case_when(
      ID == "D76V" ~ "complete + high",
      ID == "I149T" ~ "high/incomplete",
      ID == "N65S" ~ "high/incomplete",
      ID == "D90A" ~ "variable",
      ID == "V148G" ~ "complete + high",
      ID == "H46R" ~ "complete + high",
      ID == "H43R" ~ "complete + high",
      ID == "G93A" ~ "complete + high",
      ID == "G41S" ~ "complete + high",
      ID == "G16C" ~ "complete + high",
      ID == "G114A" ~ "complete + high",
      ID == "C6S" ~ "complete + high",
      ID == "C6F" ~ "complete + high",
      ID == "C6G" ~ "complete + high",
      ID == "D11Y" ~ "incomplete",
      ID == "D101H" ~ "complete + high",
      ID == "N86K" ~ "complete + high",
      ID == "P66S" ~ "incomplete",
      ID == "D76Y" ~ "incomplete",
      ID == "T137R" ~ "incomplete",
      ID == "D101G" ~ "complete + high",
      ID == "R115G" ~ "complete + high",
      ID == "A4V" ~ "complete + high",
      ID == "A4T" ~ "complete + high",
      ID == "G85R" ~ "complete + high",
      TRUE                          ~ penetrance
    )) %>% 
  select(ID, progresion, survival_time, als_type, penetrance)
      
      
sod1_phenotypes <- full_join(sod1_phenotypes, als_type, by = "ID", suffix = c("_df1", "_df2"))
sod1_phenotypes <- sod1_phenotypes %>%
  mutate(
    als_type = coalesce(als_type_df1, als_type_df2)
  ) %>% 
  select(ID, progresion, als_type, survival_time, penetrance)

sod1_phenotypes_dms <- inner_join(sod1_phenotypes, SOD1_subsitutions, by = "ID")
sod1_phenotypes_dms <- select(sod1_phenotypes_dms, ID, abundance_score, binding_score, residuals, residual_category, category_fdr_apca, category_fdr_bpca,
                              penetrance, als_type, survival_time, progresion)

write.csv(sod1_phenotypes, file = "SOD1_clinical_phenotypes_dms.csv")

sod1_phenotypes_dms <- sod1_phenotypes_dms %>% 
  mutate(residual_intensity = case_when(
    residual_category %in% c("high_positive", "high_negative") ~ "high",
    TRUE ~ "low"
  ))

#Analyze clinical phenotypes scores distribution####
##Penetrance####
sod1_penetrance <- select(sod1_phenotypes_dms, -als_type, -survival_time, -progresion)
sod1_penetrance <- na.omit(sod1_penetrance)
sod1_penetrance <- sod1_penetrance %>% 
  filter(penetrance %in% c("complete + high", "incomplete"))
apca <- sod1_penetrance %>%
  mutate(assay = "apca",
         fitness = abundance_score,
         fdr = category_fdr_apca)

bpca <- sod1_penetrance %>%
  mutate(assay = "bpca",
         fitness = binding_score,
         fdr = category_fdr_bpca)


sod1_penetrance <- rbind(apca, bpca)
sod1_penetrance$penetrance <- factor(sod1_penetrance$penetrance, levels = c("incomplete", "complete + high"))


# Abundance
tabla_apca <- table(
  sod1_penetrance$category_fdr_apca == "WT-like",
  sod1_penetrance$penetrance
)
ft_apca <- fisher.test(tabla_apca)

# Heterodimerization
tabla_bpca <- table(
  sod1_penetrance$category_fdr_bpca == "WT-like",
  sod1_penetrance$penetrance
)
ft_bpca <- fisher.test(tabla_bpca)

df_label <- data.frame(
  assay = c("apca", "bpca"),
  x = 1.5,
  y = max(sod1_penetrance$fitness, na.rm = TRUE),
  label = c(
    paste0("p = ", signif(ft_apca$p.value, 2),
           "\nOR = ", round(ft_apca$estimate, 2)),
    paste0("p = ", signif(ft_bpca$p.value, 2),
           "\nOR = ", round(ft_bpca$estimate, 2))
  )
)


labels_x <- c(
  "incomplete" = "incomplete",
  "complete + high" = "**complete +<br>high**"
)


p_penetrance<- ggplot() +
  geom_violin(data = sod1_penetrance,aes(x = penetrance, y = fitness, group = penetrance),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_penetrance,aes(x = penetrance, y = fitness, fill = fdr),width = 0.1,alpha = 0.9,size = 6, shape = 21, color = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = -0.9651394, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("high" = "darkgreen", "WT-like" = "darkgrey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  scale_x_discrete(labels = labels_x) +
  facet_wrap(~assay, ncol = 4, labeller = as_labeller(c("apca" = "Folded abundance", "bpca" = "Heterodimerization"))) +
  labs(x = "", y = "Score") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    axis.text.x = element_markdown(angle = 90, hjust = 0.5, vjust = 0.5, size = 25))
p_penetrance


ggsave(p_penetrance, file = "sod1_penetrance_dms.tiff", width = 10, height = 8, dpi = 300)

#penetrance by residual
sod1_penetrance_res <- select(sod1_phenotypes_dms, -als_type, -survival_time, -progresion)
sod1_penetrance_res$assay <- "Binding residual"
sod1_penetrance_res <- sod1_penetrance_res %>% 
  filter(penetrance %in% c("complete + high", "incomplete"))
sod1_penetrance_res$penetrance <- factor(sod1_penetrance_res$penetrance, levels = c("incomplete", "complete + high"))

p_penetrance_res<- ggplot() +
  geom_violin(data = sod1_penetrance_res,aes(x = assay, y = residuals),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_penetrance_res,aes(x = assay, y = residuals, fill = residual_category),width = 0.1,alpha = 0.9,size = 6, shape = 21, color = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  facet_wrap(~penetrance, ncol = 5) +
  labs(x = "", y = "Heterodimerization residual", color = "") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    axis.text.x = element_blank())
p_penetrance_res


##Survival####
sod1_survival<- select(sod1_phenotypes_dms, -als_type, -penetrance, -progresion)
sod1_survival <- na.omit(sod1_survival)

apca <- sod1_survival %>%
  mutate(assay = "apca",
         fitness = abundance_score,
         fdr = category_fdr_apca)

bpca <- sod1_survival %>%
  mutate(assay = "bpca",
         fitness = binding_score,
         fdr = category_fdr_bpca)


sod1_survival <- rbind(apca, bpca)

# Abundance
tabla_apca <- table(
  sod1_survival$category_fdr_apca == "WT-like",
  sod1_survival$survival_time
)
ft_apca <- fisher.test(tabla_apca)

# Heterodimerization
tabla_bpca <- table(
  sod1_survival$category_fdr_bpca == "WT-like",
  sod1_survival$survival_time
)
ft_bpca <- fisher.test(tabla_bpca)

df_label <- data.frame(
  assay = c("apca", "bpca"),
  x = 1.5,
  y = max(sod1_survival$fitness, na.rm = TRUE),
  label = c(
    paste0("p = ", signif(ft_apca$p.value, 2),
           "\nOR = ", round(ft_apca$estimate, 2)),
    paste0("p = ", signif(ft_bpca$p.value, 2),
           "\nOR = ", round(ft_bpca$estimate, 2))
  )
)

labels_x <- c(
  "long" = "long<br>(> 10 years)",
  "short" = "**short<br>(< 3 years)**"
)

p_survival <- ggplot() +
  geom_violin(data = sod1_survival,aes(x = survival_time, y = fitness, group = survival_time),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_survival,aes(x = survival_time, y = fitness, fill = fdr),width = 0.1,alpha = 0.9,size = 6, shape = 21, color = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = -0.9651394, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("high" = "darkgreen", "WT-like" = "darkgrey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  scale_x_discrete(labels = labels_x) +
  facet_wrap(~assay, ncol = 4, labeller = as_labeller(c("apca" = "Folded abundance", "bpca" = "Heterodimerization"))) +
  labs(x = "", y = "Score") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 25))

p_survival

ggsave(p_survival, file = "sod1_survival_dms.tiff", width = 10, height = 8, dpi = 300)

#survival by residual
sod1_survival_res <- select(sod1_phenotypes_dms, -als_type, -penetrance, -progresion)
sod1_survival_res$assay <- "Binding residual"
sod1_survival_res <- na.omit(sod1_survival_res)

p_survival_res<- ggplot() +
  geom_violin(data = sod1_survival_res,aes(x = assay, y = residuals),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_survival_res,aes(x = assay, y = residuals, fill = residual_category),width = 0.1,alpha = 0.9,size = 6, shape = 21, color = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  facet_wrap(~survival_time, ncol = 5) +
  labs(x = "", y = "Heterodimerization residual", color = "") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    axis.text.x = element_blank())
p_survival_res


##ALS type####
sod1_alstype<- select(sod1_phenotypes_dms, -survival_time, -penetrance, -progresion)
sod1_alstype <- na.omit(sod1_alstype)

apca <- sod1_alstype %>%
  mutate(assay = "apca",
         fitness = abundance_score,
         fdr = category_fdr_apca)

bpca <- sod1_alstype %>%
  mutate(assay = "bpca",
         fitness = binding_score,
         fdr = category_fdr_bpca)


sod1_alstype <- rbind(apca, bpca)

p_alstype<- ggplot() +
  geom_violin(data = sod1_alstype,aes(x = als_type, y = fitness, group = als_type),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_alstype,aes(x = als_type, y = fitness, color = fdr),width = 0.1,alpha = 0.9,size = 6) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = -0.9651394, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("high" = "darkgreen", "WT-like" = "darkgrey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  facet_wrap(~assay) +
  labs(x = "", y = "Normalized score", color = "FDR = 0.1") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 28),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p_alstype


##Progression####
sod1_progresion<- select(sod1_phenotypes_dms, -survival_time, -penetrance, -als_type)
sod1_progresion <- na.omit(sod1_progresion)

apca <- sod1_progresion %>%
  mutate(assay = "apca",
         fitness = abundance_score,
         fdr = category_fdr_apca)

bpca <- sod1_progresion %>%
  mutate(assay = "bpca",
         fitness = binding_score,
         fdr = category_fdr_bpca)


sod1_progresion <- rbind(apca, bpca)
sod1_progresion$progresion <- factor(sod1_progresion$progresion, levels = c("slow", "fast"))

# Abundance
tabla_apca <- table(
  sod1_progresion$category_fdr_apca == "WT-like",
  sod1_progresion$progresion
)
ft_apca <- fisher.test(tabla_apca)

# Heterodimerization
tabla_bpca <- table(
  sod1_progresion$category_fdr_bpca == "WT-like",
  sod1_progresion$progresion
)
ft_bpca <- fisher.test(tabla_bpca)

df_label <- data.frame(
  assay = c("apca", "bpca"),
  x = 1.5,
  y = max(sod1_progresion$fitness, na.rm = TRUE),
  label = c(
    paste0("p = ", signif(ft_apca$p.value, 2),
           "\nOR = ", round(ft_apca$estimate, 2)),
    paste0("p = ", signif(ft_bpca$p.value, 2),
           "\nOR = ", round(ft_bpca$estimate, 2))
  )
)

labels_x <- c(
  "slow" = "slow",
  "fast" = "**fast**"
)


p_progression<- ggplot() +
  geom_violin(data = sod1_progresion,aes(x = progresion, y = fitness, group = progresion),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_progresion,aes(x = progresion, y = fitness, fill = fdr),width = 0.1,alpha = 0.9,size = 6, shape = 21, color = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_hline(yintercept = -0.9651394, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("high" = "darkgreen", "WT-like" = "darkgrey", "low" = "darkorange", "stop-like" = "#DD4636")) +
  scale_x_discrete(labels = labels_x) +
  facet_wrap(~assay, ncol = 4, labeller = as_labeller(c("apca" = "Folded abundance", "bpca" = "Heterodimerization"))) +
  labs(x = "", y = "Score") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = 25))
p_progression


ggsave(p_progression, file = "sod1_progresion_dms.tiff", width = 10, height = 8, dpi = 300)

p_phenotypes <- ggarrange(p_penetrance,
                          p_survival,
                          p_progression,
                          ncol = 3, align = "h")

p_phenotypes
ggsave(p_phenotypes, file = "sod1_phenotypes_comb_dms.tiff", width = 28, height = 8, dpi = 300)


#progression by residual
sod1_progresion_res <- select(sod1_phenotypes_dms, -als_type, -penetrance, -survival_time)
sod1_progresion_res$assay <- "Binding residual"
sod1_progresion_res$progresion <- factor(sod1_progresion_res$progresion, levels = c("slow", "fast"))
sod1_progresion_res <- na.omit(sod1_progresion_res)


p_progression_res<- ggplot() +
  geom_violin(data = sod1_progresion_res,aes(x = assay, y = residuals),fill = NA,scale = "width",trim = TRUE,size = 0.8) +
  geom_jitter(data = sod1_progresion_res,aes(x = assay, y = residuals, fill = residual_category),width = 0.1,alpha = 0.9,size = 6, shape = 21, color = "black") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low" = "#EDEDED")) +
  facet_wrap(~progresion, ncol = 5) +
  labs(x = "", y = "Heterodimerization residual", color = "") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 25),
    axis.text.x = element_blank())
p_progression_res


p_progression_comb_res <- ggarrange(p_progression_res + theme(axis.text.x = element_blank()),
                                 p_bar_progression_res + theme(axis.text = element_blank(), axis.title = element_blank()),
                                 nrow = 2, align = "v", heights = c(2,1))
p_progression_comb_res

ggsave(p_progression_comb_res, file = "sod1_progression_dms_residual.tiff", width = 8, height = 5, dpi = 300)

p_phenotypes_res <- ggarrange(p_penetrance_res,
                          p_survival_res,
                          p_progression_res,
                          ncol = 3, align = "h")

p_phenotypes_res
ggsave(p_phenotypes_res, file = "sod1_phenotypes_res_comb_dms.tiff", width = 18, height = 6, dpi = 300)