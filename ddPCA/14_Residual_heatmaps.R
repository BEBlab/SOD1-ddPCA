library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(ggrepel)

#load modelled data
SOD1_subsitutions <- read.csv("SOD1_subs_modelled.csv")
SOD1_insertions <- read.csv("SOD1_ins_modelled.csv")
SOD1_deletions <- read.csv("SOD1_del_modelled.csv")


#Substitutions residual heatmap####
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(WT_AA = str_extract(ID, "^[A-za-z]"),
         Pos = str_extract(ID, "\\d+"),
         Mut_AA = str_extract(ID, "[A-za-z]$")) %>% 
  mutate(Pos = as.numeric(str_extract(ID, "\\d+"))) %>% 
  filter(mutation_type != "stop")


peptide_seq<-'ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")


#add syn

positions <- 1:153
syn.df <- data.frame(
  "WT_AA" = peptide_seq,
  "Mut_AA" = peptide_seq,
  "Pos" = positions,
  "residuals" = 0,
  "ID" = "syn"
)



heatmap_residual_subs<-rbind(SOD1_subsitutions[,c("WT_AA", "Mut_AA", "Pos", "residuals", "ID")], syn.df)


heatmap_residual_subs$label<-""
heatmap_residual_subs[heatmap_residual_subs$ID=="syn",]$label<-"WT"

heatmap_residual_subs <- heatmap_residual_subs %>%
  complete(Pos = 1:153, Mut_AA = vectorAA, fill = list(residuals = NA, WT_AA = NA, ID = NA, label = ""))


min<-min(heatmap_residual_subs$residuals, na.rm = TRUE)
max<-max(heatmap_residual_subs$residuals, na.rm = TRUE)
cols <- c(colorRampPalette(c( "#B99C33", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "#AC7299"), bias=1)((max/(-min+max)*100)-0.5))


p_heatmap_res_subs<-ggplot(heatmap_residual_subs)+
  geom_tile(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), fill=residuals), color='white', size=0.2)+
  scale_x_continuous(breaks=seq(1,153), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  theme()+
  labs(x="SOD1 WT amino acid position", y="Mutant amino acid", fill="Binding\nresidual")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.position = "none",
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-0.39, -0.2, 0, 0.2, 0.4,  0.6), labels = c(-0.4, -0.2, 0 , 0.2, 0.4, 0.7)) +
  guides(fill = guide_colorbar(
    barwidth = 1,   
    barheight = 5,
    direction = "vertical",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"))
p_heatmap_res_subs

#add side_chain information
sasa <- sasa %>%
  mutate(side_chain = case_when(
    Pos %in% c(63, 71, 80, 83) ~ "Zn binding"
    TRUE ~ side_chain))

sasa$row <- "row"
p_side_chain <- ggplot(sasa, aes(x = factor(Pos), y = row, fill = side_chain)) +
  geom_tile(color = "white", height = 0.1, size = 0.2) +
  scale_fill_manual(
    values = c(
      "core" = "#33608C",
      "surface" = "#C4D8F3",
      "dimer interface" = "#41B7C4",
      "Zn binding" = "#BB173A")
  ) +
  scale_y_discrete(expand = c(0,0)) +  # <--- elimina el espacio arriba y abajo
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0,0,0.5,0), 'cm')
  ) +
  labs(fill = "Side chain type")

p_side_chain

p_heatmap_side_chain <- ggarrange(p_side_chain + theme(axis.title = element_blank(), axis.text = element_blank()),
                                  p_heatmap_res_subs,
                                  heights = c(0.05, 1), nrow = 2, align = "v")
p_heatmap_side_chain



#2. Insertions
SOD1_insertions <- SOD1_allmut_filt %>%
  filter(mutation_type %in% c("ins", "stop"))

#add position 1. We have not mutated it but we show it as NA
aa_vector <- sort(unique(SOD1_insertions$Mut_AA))  

missing_pos1 <- data.frame(
  ID = paste0("ins1", aa_vector),
  Pos = 1,
  Mut_AA = aa_vector,
  residuals = NA  # valor NA para el heatmap
)

SOD1_insertions$Pos <- as.numeric(SOD1_insertions$Pos)
SOD1_insertions <- bind_rows(SOD1_insertions, missing_pos1)

positions <- 1:153  
SOD1_insertions <- SOD1_insertions %>%
  mutate(Pos = as.numeric(Pos)) %>%
  mutate(WT_AA = peptide_seq[Pos],
         Mut_AA = case_when(
           mutation_type == "stop" ~ "stop",
           TRUE ~ Mut_AA
         )) 


heatmap_df_apca_ins<-SOD1_insertions

heatmap_df_apca_ins <- heatmap_df_apca_ins %>%
  complete(Pos = 1:153, Mut_AA = vectorAA, fill = list(residuals = NA, WT_AA = NA, ID = NA, label = ""))

min_ins<-min(heatmap_df_apca_ins$residuals, na.rm = TRUE)
max_ins<-max(heatmap_df_apca_ins$residuals, na.rm = TRUE)
cols <- c(colorRampPalette(c( "#B99C33", "grey95"))((-min_ins/(-min_ins+max_ins)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "#AC7299"), bias=1)((max_ins/(-min_ins+max_ins)*100)-0.5))

p_heatmap_ins_residual <-ggplot(heatmap_df_apca_ins)+
  geom_tile(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), fill=residuals), color='white', size=0.2)+
  scale_x_continuous(breaks=positions, labels = positions, expand=c(0,0))+
  theme_minimal()+
  theme()+
  labs(x="Position of inserted position", y="Inserted amino acid", fill="Binding\nresidual")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.position = "none",
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-0.39, -0.2, 0, 0.2, 0.4,  0.6), labels = c(-0.4, -0.2, 0 , 0.2, 0.4, 0.7)) +
  guides(fill = guide_colorbar(
    barwidth = 1,   
    barheight = 5,
    direction = "vertical",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"))

p_heatmap_ins_residual

p_heatmap_side_chain_ins <- ggarrange(p_side_chain + theme(axis.title = element_blank(), axis.text = element_blank()),
                                      p_heatmap_ins_residual,
                                      heights = c(0.05, 1), nrow = 2, align = "v")
p_heatmap_side_chain_ins

#3. Deletions
SOD1_deletions <- SOD1_allmut_filt %>%
  filter(mutation_type=="del")


SOD1_deletions <- SOD1_deletions %>%
  separate_rows(ID, sep = ";") %>%
  mutate(
    del_pos = str_extract(ID, "\\d+"),
    del_aa = str_extract(ID, "^[A-Z]")
  )


min_del<-min(SOD1_deletions$residuals)
max_del<-max(SOD1_deletions$residuals)
cols <- c(colorRampPalette(c( "#B99C33", "grey95"))((-min_del/(-min_del+max_del)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "#AC7299"), bias=1)((max_del/(-min_del+max_del)*100)-0.5))

SOD1_deletions <- SOD1_deletions %>%
  distinct(ID, .keep_all = TRUE)



p_residual_deletions<- ggplot(SOD1_deletions, aes(x = factor(del_pos, levels = 1:153), y = residuals)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_point(data = SOD1_deletions,
             aes(fill = residuals), size = 8, shape = 21, stroke = 1.2) +
  labs(y = "Heterodimerization residual", x = "Position of deleted amino acid", fill = "Binding\nresidual") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.2),
    axis.line.x = element_line(color = "black", linewidth = 0.2),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size=14),
    legend.position = "none",
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-1, 0, 1), labels = c(-1, 0, 1)) +
  guides(fill = guide_colorbar(
    barwidth = 5,
    barheight = 1,
    direction = "horizontal",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"
  ))

p_residual_deletions

p_residual_deletions <- ggarrange(p_side_chain + theme(axis.title = element_blank(), axis.text = element_blank()),
                                  p_residual_deletions,
                                  heights = c(0.05, 1), nrow = 2, align = "v")
p_residual_deletions

#ggsave(p_heatmap_side_chain, path = path_fig4, file = "residual_heatmap_subs.tiff", width = 30, height = 10, dpi = 300)
#ggsave(p_heatmap_side_chain_ins, path = path_fig4, file = "residual_heatmap_ins.tiff", width = 30, height = 10, dpi = 300)
ggsave(p_residual_deletions, path = path_fig4, file = "residual_heatmap_del.tiff", width = 25, height = 8, dpi = 300)

p_heatmaps_subs_ins <- ggarrange(p_heatmap_side_chain,
                                 p_heatmap_side_chain_ins,
                                 nrow = 2, common.legend = T)

p_heatmaps_subs_ins
ggsave(p_heatmaps_subs_ins, path = path_fig4, file = "residual_heatmap_subsins.tiff", width = 30, height = 20, dpi = 300)


#Median residual heatmap####
SOD1_allmut_filt <- SOD1_allmut_filt %>% 
  mutate(second_shell = case_when(
    Pos %in% c(1,2,3,6,8,9,15,16,18,19,20,21,32,33,34,47,48,49,55,56,59,60,62,64,106,107,108,109,110,112,116,117,146,147) ~ TRUE,
    TRUE ~ FALSE
  ))

heatmap_res_median <- SOD1_allmut_filt %>%
  group_by(Pos, mutation_type, interface, second_shell) %>%
  summarise(residual_median = median(residuals, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(row = "row")  

# SOD1_subsitutions <- SOD1_allmut_filt %>% 
#   filter(mutation_type == "subs")


heatmap_res_median$mutation_type <- factor(heatmap_res_median$mutation_type,
                                           levels = c("subs", "ins", "del"))

heatmap_res_median <- na.omit(heatmap_res_median)

heatmap_res_median$Pos <- as.numeric(heatmap_res_median$Pos)
heatmap_res_median$Pos_range <- cut(
  heatmap_res_median$Pos,
  breaks = seq(0, max(heatmap_res_median$Pos), by = 10),
  right = FALSE,
  include.lowest = TRUE
)

p_heatmap_median_res<- ggplot(data = heatmap_res_median, aes(x = Pos, y = row, fill=residual_median)) +
  geom_tile(size = 0.1, width = 1, height = 2, color = "darkgrey") +
  geom_tile(data = subset(heatmap_res_median, interface == TRUE),aes(x = Pos, y = row),
            fill = NA,color = "black",size = 0.9, width = 1, height = 2) +
  geom_text(
    data = subset(heatmap_res_median, second_shell == TRUE),
    aes(x = Pos, y = row),
    label = "*",
    color = "#555F6A",
    size = 5
  ) +
  scale_fill_gradient2(
    low = "#B99C33",
    mid = "grey95",
    high = "#AC7299",
    midpoint = 0,
    limits = c(min(heatmap_res_median$residual_median), 
               max(heatmap_res_median$residual_median))) +
  theme_classic() +
  labs(x = "SOD1 WT AA position", y = "", fill = "Median of binding residual") +
  facet_wrap(~mutation_type, nrow = 3, labeller = as_labeller(c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions"))) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        strip.background = element_rect(linewidth = 0.1),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        legend.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size=18),
        legend.position = "top") +
  guides(fill = guide_colorbar(
    barwidth = 12,   
    barheight = 2,
    direction = "horizontal",
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE, ticks.colour = "white", frame.colour = "white"))

p_heatmap_median_res

ggsave(p_heatmap_median_res, path = path_fig4, file = "median_residual_heatmaps.tiff", width = 18, height = 8)




#### Analyse dimer interface ####
#1. First, I want to determine which are the interface residues that only affect binding (they have negative binding and ~wt-like abundance)
#plot only interface residues
interface_subs <- SOD1_subsitutions %>% 
  filter(interface == TRUE)

median_residual_int_subs <- median(abs(interface_subs$residuals))


SOD1_side_chains_all_int <- SOD1_side_chains_filt %>%
  filter(interface == TRUE)
SOD1_side_chains_all_int$mutation_type <- factor(SOD1_side_chains_filt$mutation_type,
                                                 levels = c("subs", "ins", "del"))


SOD1_subsitutions <- SOD1_side_chains_all_int %>%
  filter(mutation_type == "subs")

SOD1_insertions <- SOD1_side_chains_all_int %>%
  filter(mutation_type == "ins")

SOD1_deletions <- SOD1_side_chains_all_int %>%
  filter(mutation_type == "del")

p_loess_interface_subs <- ggplot(SOD1_subsitutions, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = outlier), size = 5, shape = 16) +
  scale_color_manual(values = c("TRUE" = "#B254A5", "FALSE" = "#D5D5D5"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  coord_cartesian(xlim = range(SOD1_side_chains_all_int$abundance_score),
                  ylim = range(SOD1_side_chains_all_int$binding_score)) +
  labs(x = "Normalized abundance score", y = "Normalized binding score", color = "")+ 
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 33),
        legend.position = "top") +
  guides(color = guide_legend(
    override.aes = list(size = 13)))

p_loess_interface_subs

p_loess_interface_ins <- ggplot(SOD1_insertions, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = outlier), size = 5, shape = 16) +
  scale_color_manual(values = c("TRUE" = "#B254A5", "FALSE" = "#D5D5D5"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  coord_cartesian(xlim = range(SOD1_side_chains_all_int$abundance_score),
                  ylim = range(SOD1_side_chains_all_int$binding_score)) +
  labs(x = "Normalized abundance score", y = "Normalized binding score", color = "")+ 
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 33),
        legend.position = "top") +
  guides(color = guide_legend(
    override.aes = list(size = 13)))

p_loess_interface_ins

p_loess_interface_del <- ggplot(SOD1_deletions, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = outlier), size = 5, shape = 16) +
  scale_color_manual(values = c("TRUE" = "#B254A5", "FALSE" = "#D5D5D5"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  coord_cartesian(xlim = range(SOD1_side_chains_all_int$abundance_score),
                  ylim = range(SOD1_side_chains_all_int$binding_score)) +
  labs(x = "Normalized abundance score", y = "Normalized binding score", color = "")+ 
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 33),
        legend.position = "top") +
  guides(color = guide_legend(
    override.aes = list(size = 13)))
p_loess_interface_del


p_combined_interface <- ggarrange(p_loess_interface_subs + theme(axis.title = element_blank()),
                                  p_loess_interface_ins + theme(axis.title.x = element_blank()),
                                  p_loess_interface_del + theme(axis.title.y = element_blank()),
                                  nrow = 3, align = "v", common.legend = TRUE)
p_combined_interface


ggsave(p_loess_interface_subs, path = path_fig3, file = "interface_correlation_subs.tiff", width = 8, height = 7, dpi = 300)
ggsave(p_loess_interface_ins, path = path_fig3, file = "interface_correlation_ins.tiff", width = 8, height = 7, dpi = 300)
ggsave(p_loess_interface_del, path = path_fig3, file = "interface_correlation_del.tiff", width = 8, height = 7, dpi = 300)

ggsave(p_combined_interface, path = path_fig3, file = "interface_correlation_subs_indels.tiff", width = 9, height = 10, dpi = 300)


#2. Second, I will use the FDR categories of the outliers to see which mutations have WT-like/increase abundance and decreased binding
#calculate the median residual by position of all interface residues
heatmap_res_median <- SOD1_side_chains_all_int %>%
  group_by(mutation_type,Pos) %>%
  summarise(residual_median = median(residuals, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(row = "row")


heatmap_res_median$Pos <- factor(
  heatmap_res_median$Pos,
  levels = as.character(sort(as.numeric(unique(heatmap_res_median$Pos))))
)
heatmap_res_median$mutation_type <- factor(heatmap_res_median$mutation_type,
                                           levels = c("subs", "ins", "del"))

p_heatmap_residual_interface<- ggplot(data = heatmap_res_median, aes(x = Pos, y = row, fill=residual_median)) +
  geom_tile(size = 0.2, color = "black") +
  scale_fill_gradient2(
    low = "#B99C33",
    mid = "grey95",
    high = "#AC7299",
    midpoint = 0,
    limits = c(min(heatmap_res_median$residual_median), 
               max(heatmap_res_median$residual_median))) +
  facet_wrap(~mutation_type, nrow = 3, labeller = as_labeller(c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions"))) +
  theme_classic() +
  labs(x = "Dimer interface position", y = "", fill = "Median residual") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        legend.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size=18),
        legend.position = "right"
  )
p_heatmap_residual_interface

ggsave(p_heatmap_residual_interface, path = path_fig4, file = "dimer_interface_residual_heatmap.tiff", width  = 16, height = 4, dpi = 300)


write.csv(SOD1_interface, file = "3D_interface_median_residual.csv")


#3. Now I want to plot the dimer interface positions in decrease order of residual
SOD1_subsitutions <- SOD1_subsitutions %>% 
  mutate(abs_residual = abs(residuals))

medians <- SOD1_subsitutions %>% 
  group_by(Pos) %>% 
  summarise(median_residual = abs(median(abs_residual, na.rm = TRUE)))

# 2. Crear factor ordenado en SOD1_subsitutions
SOD1_subsitutions <- SOD1_subsitutions %>%
  mutate(Pos = factor(Pos, levels = medians %>% arrange(desc(median_residual)) %>% pull(Pos)))

# 3. Graficar todos los puntos individuales
p_interface_res_order <- ggplot(SOD1_subsitutions, aes(x = Pos, y = abs_residual)) + 
  geom_violin(width = 1.2, fill = "lightgrey") +
  geom_point(data = medians, aes(x = Pos, y = median_residual), color = "black", size = 4) +
  theme_classic() +
  geom_hline(yintercept = median_residual_int_subs, linetype = "dashed", size = 1) +
  labs(x = "Dimer interface position", y = "|Binding residual|") +
  theme(axis.title = element_text(size = 33),
        axis.text.y = element_text(size = 28),
        axis.text.x = element_text(size = 24, angle = 90, hjust = 1))

p_interface_res_order

ggsave(p_interface_res_order, path = path_fig3, file = "importance_interface_positions.tiff", width = 12, height = 8, dpi = 300)

#### Define residuals intensity effects ####
#I will determine the magnitude if the residual effect by comparing the residual median with the residual mean and sd of the dimer interface

#1. Calculate the absolute median per position
median_pos_int <- SOD1_side_chains_all_int %>%
  group_by(mutation_type, Pos) %>%
  summarise(
    residual_median = median(residuals, na.rm = TRUE),
    residual_median_abs = median(abs(residuals), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(row = "row")


median_pos_int$Pos <- factor(
  median_pos_int$Pos,
  levels = as.character(sort(as.numeric(unique(median_pos_int$Pos))))
)
median_pos_int$mutation_type <- factor(median_pos_int$mutation_type,
                                       levels = c("subs", "ins", "del"))

# p_heatmap_residual_interface<- ggplot(data = median_pos_int, aes(x = Pos, y = row, fill=residual_median_abs)) +
#   geom_tile(size = 0.2, color = "black") +
#   scale_fill_gradient2(
#     low = "#B99C33",
#     mid = "grey95",
#     high = "#AC7299",
#     midpoint = 0,
#     limits = c(min(median_pos_int$residual_median_abs), 
#                max(median_pos_int$residual_median_abs))) +
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


#2. Calculate the mean and standard deviation of residuals medians of the dimer interface
mean_res_int <- mean(median_pos_int$residual_median_abs)
sd_res_int <- sd(median_pos_int$residual_median_abs)

median_pos_int <- median_pos_int %>%
  mutate(residual_effect = case_when(
    residual_median > mean_res_int + sd_res_int ~ "high residual",
    residual_median < mean_res_int - sd_res_int ~ "moderate residual",
    TRUE ~ "low residual"  
  ))
median_pos_int_abs$residual_effect <- factor(median_pos_int_abs$residual_effect,
                                             levels = c("high residual", "moderate residual", "low residual"))


median_pos_int <- median_pos_int %>%
  mutate(
    residual_intensity = case_when(
      residual_median_abs > mean_res_int + sd_res_int ~ "high",
      residual_median_abs < mean_res_int - sd_res_int ~ "low",
      TRUE ~ "medium"
    ),
    residual_direction = case_when(
      residual_median > 0 ~ "positive",
      residual_median < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    residual_category = paste(residual_intensity, residual_direction, sep = "_")
  )


median_pos_int$residual_category <- factor(median_pos_int$residual_category,
                                           levels = c("high_positive", "medium_positive", "low_positive",
                                                      "low_negative", "medium_negative", "high_negative"))

#3. Plot residual categories
p_res_categories <- ggplot(median_pos_int, aes(x = Pos, y = row, fill=residual_category)) +
  geom_tile(size = 0.2, color = "black") +
  scale_fill_manual(values = c("high_positive" = "#7C4D79","high_negative" = "#AD9024","low_positive" = "#DCBDC9", "low_negative" = "#E2D5C3",
                               "medium_positive" = "#AF7CA1", "medium_negative" = "#C6A844")) +
  facet_wrap(~mutation_type, nrow = 3, labeller = as_labeller(c("subs" = "Substitutions", "ins" = "Insertions", "del" = "Deletions"))) +
  theme_classic() +
  labs(x = "Dimer interface position", y = "", fill = "Residual category") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        legend.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.text = element_text(size=18),
        legend.position = "right"
  )
p_res_categories

ggsave(p_res_categories, path = path_fig4, file = "dimer_interface_residual_categories.tiff", width  = 16, height = 4, dpi = 300)

p_residual_interface <- ggarrange(p_heatmap_residual_interface + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                                  p_res_categories,
                                  nrow = 2, align = "hv", common.legend = FALSE)
p_residual_interface

ggsave(p_residual_interface, path = path_fig4, file = "dimer_interface_residual_categories.tiff", width  = 16, height = 8, dpi = 300)

#save to color 3D structure (only by residual)
SOD1_side_chains_all_int_3d <- SOD1_side_chains_all_int %>%
  select(ID, WT_AA, Mut_AA, Pos, residuals, mutation_type)

write.csv(SOD1_side_chains_all_int_3d, file = "3D_interface_median_residual.csv", row.names = FALSE)


#### Define the role of dimer interlace positions in dimerization ####
#by defining residual categories, I could determine which are the dimer interface positions that are less influenced by abundance
#now, I want to determine the importance of each of them

#1. Heatmap for top dimer interface positions
top_interface_subs <- SOD1_side_chains_all_int %>% 
  filter(mutation_type == "subs") %>% 
  filter(Pos %in% c(112, 113, 114, 149))

#add syn
positions <- 1:153
syn.df <- data.frame(
  "WT_AA" = peptide_seq,
  "Mut_AA" = peptide_seq,
  "Pos" = positions,
  "binding_score" = 0,
  "ID" = "syn"
)



heatmap_top_int_subs<-rbind(top_interface_subs[,c("WT_AA", "Mut_AA", "Pos", "binding_score", "ID")], syn.df)


heatmap_top_int_subs$label<-""
heatmap_top_int_subs[heatmap_top_int_subs$ID=="syn",]$label<-"WT"
heatmap_top_int_subs$Pos <- as.numeric(heatmap_top_int_subs$Pos)

heatmap_top_int_subs <- heatmap_top_int_subs %>%
  complete(Pos = 1:153, Mut_AA = vectorAA, fill = list(binding_score = NA, WT_AA = NA, ID = NA, label = ""))


min<-min(heatmap_top_int_subs$binding_score, na.rm = TRUE)
max<-max(heatmap_top_int_subs$binding_score, na.rm = TRUE)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

heatmap_top_int_subs <- heatmap_top_int_subs %>%
  filter(Pos %in% c(112, 113, 114, 149))

heatmap_top_int_subs <- heatmap_top_int_subs %>%
  mutate(Pos = factor(Pos, levels = c(112, 113, 114, 149)))


p_heatmap_top_int_subs<-ggplot(heatmap_top_int_subs)+
  geom_tile(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), fill=binding_score), color='white', size=0.2)+
  scale_x_discrete(labels = peptide_seq_pos[c(112, 113, 114, 149)])  +
  theme_minimal()+
  theme()+
  labs(x="SOD1 top dimer interface position", y="Mutant amino acid", fill="Binding score")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14)) +
  geom_text(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey") +
  guides(fill = guide_colorbar(
    barwidth = 1,   
    barheight = 5,
    direction = "vertical",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"))
p_heatmap_top_int_subs


top_interface_ins <- SOD1_side_chains_all_int %>% 
  filter(mutation_type == "ins") %>% 
  filter(Pos %in% c(112, 114, 149))

#add syn
positions <- 1:153
syn.df <- data.frame(
  "WT_AA" = peptide_seq,
  "Mut_AA" = peptide_seq,
  "Pos" = positions,
  "binding_score" = 0,
  "ID" = "syn"
)



heatmap_top_int_ins<-rbind(top_interface_ins[,c("WT_AA", "Mut_AA", "Pos", "binding_score", "ID")], syn.df)


heatmap_top_int_ins$label<-""
heatmap_top_int_ins[heatmap_top_int_ins$ID=="syn",]$label<-"WT"
heatmap_top_int_ins$Pos <- as.numeric(heatmap_top_int_ins$Pos)

heatmap_top_int_ins <- heatmap_top_int_ins %>%
  complete(Pos = 1:153, Mut_AA = vectorAA, fill = list(binding_score = NA, WT_AA = NA, ID = NA, label = ""))


min<-min(heatmap_top_int_ins$binding_score, na.rm = TRUE)
max<-max(heatmap_top_int_ins$binding_score, na.rm = TRUE)
neg_colors <- max(round((-min/(-min+max)*100)-0.5), 1)
pos_colors <- max(round((max/(-min+max)*100)-0.5), 1)


cols <- c(colorRampPalette(c("darkorange", "grey95"))(neg_colors),"grey95",
          colorRampPalette(c("grey95", "darkgreen"))(pos_colors))

heatmap_top_int_ins <- heatmap_top_int_ins %>%
  filter(Pos %in% c(112, 114, 149))

heatmap_top_int_ins <- heatmap_top_int_ins %>%
  mutate(Pos = factor(Pos, levels = c(112, 114, 149)))


p_heatmap_top_int_ins<-ggplot(heatmap_top_int_ins)+
  geom_tile(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), fill=binding_score), color='white', size=0.2)+
  scale_x_discrete(labels = peptide_seq_pos[c(112, 114, 149)])  +
  theme_minimal()+
  theme()+
  labs(x="SOD1 top dimer interface position", y="Mutant amino acid", fill="Binding score")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14)) +
  geom_text(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey") +
  guides(fill = guide_colorbar(
    barwidth = 1,   
    barheight = 5,
    direction = "vertical",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"))
p_heatmap_top_int_ins














#### Residual vs minimal distance #### 
#load atom distance data
min_distance <- read.csv("min_atom_distance_SOD1.csv")
min_distance <- rename(min_distance, Pos = resno_A)


########## Substitutions allosteric sites ###
#determine type or residual (< or >)
SOD1_subsitutions$Pos <- as.numeric(SOD1_subsitutions$Pos)
# SOD1_subsitutions <- SOD1_subsitutions %>%
#   mutate(residual_type = case_when(
#     residuals > 0 ~ "positive",
#     residuals < 0 ~ "negative"
#   ))

SOD1_subsitutions_dist <- left_join(SOD1_subsitutions, min_distance, by = "Pos")

#per position
#calculate de median residual in interface residues
SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>%
  group_by(Pos) %>%
  mutate(median_residual_abs = median(abs(residuals), na.rm = TRUE),
         median_residual = median(residuals),
         residual_type = case_when(
           median_residual > 0 ~ "positive",
           median_residual < 0 ~ "negative")) %>%
  select(ID, median_residual_abs, median_residual, Pos, interface, min_dist_to_B, location, residual_type, s_location)



SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>%
  distinct(Pos, .keep_all = TRUE)
SOD1_subsitutions_dist <- na.omit(SOD1_subsitutions_dist)

interface_subs <- SOD1_subsitutions %>% 
  filter(interface == TRUE)

p_min_dist_pos_subs <- ggplot(SOD1_subsitutions_dist, aes(x = min_dist_to_B, y = median_residual_abs)) +
  geom_point(aes(color = interface, shape = residual_type), size = 6, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(aes(label = Pos, color = interface), size = 6) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 33),
        axis.title = element_text(size = 36),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 6.477411, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = 0.11, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_pos_subs

#per mutation
SOD_final_subs_dist$residuals <- abs(SOD_final_subs_dist$residuals)
p_min_dist_mut_subs <- ggplot(SOD_final_subs_dist, aes(x = min_dist_to_B, y = residuals)) +
  geom_point(aes(color = interface, shape = residual_type), size = 4, alpha = 0.5) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  #geom_text_repel(aes(label = Pos, color = interface), size = 6) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 6.477411, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_residual_int_subs, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_mut_subs

high_res_int_count <- SOD_final_subs_dist %>%
  filter(interface == TRUE) %>% 
  group_by(Pos) %>%
  summarise(
    n_high_res = sum(abs(residuals) > median_residual_int_subs, na.rm = TRUE),
    n_total = sum(!is.na(residuals)),
    frac_alosteric = n_high_res / n_total)

high_res_int_count <- high_res_int_count %>%
  mutate(Pos = as.numeric(Pos)) %>%
  arrange(Pos)

p_high_res_mut_pos_subs <- ggplot(high_res_int_count, aes(x = factor(Pos, levels = sort(Pos)), y = n_high_res)) +
  geom_col(width = 0.8, fill = "lightgrey", color = "black") +
  theme_classic() +
  labs(x = "Dimer interface position", y = "Number of mutations\nwith residual > median residual interface")
p_high_res_mut_pos_subs

########## Insertions allosteric sites ###
#determine type or residual (< or >)
SOD_final_ins$Pos <- as.numeric(SOD_final_ins$Pos)
SOD_final_ins <- SOD_final_ins %>%
  mutate(residual_type = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  ))

SOD_final_ins_dist <- left_join(SOD_final_ins, min_distance, by = "Pos")

#per position
SOD_final_ins_dist_pos <- SOD_final_ins_dist %>%
  group_by(Pos) %>%
  mutate(median_residual = abs(median(residuals, na.rm = TRUE))) %>%
  select(ID, median_residual, Pos, interface, min_dist_to_B, location, residual_type)


SOD_final_ins_dist_pos <- SOD_final_ins_dist_pos %>%
  distinct(Pos, .keep_all = TRUE)
SOD_final_ins_dist_pos <- na.omit(SOD_final_ins_dist_pos)

#calculate de median residual in interface residues
interface_ins <- SOD_final_ins %>% 
  filter(interface == TRUE)

median_residual_int_ins <- median(abs(interface_ins$residuals))


p_min_dist_pos_ins <- ggplot(SOD_final_ins_dist_pos, aes(x = min_dist_to_B, y = median_residual)) +
  geom_point(aes(color = interface, shape = residual_type), size = 4, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(aes(label = Pos, color = interface), size = 6) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 7.290945, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_residual_int_ins, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_pos_ins


#per mutation
SOD_final_ins_dist$residuals <- abs(SOD_final_ins_dist$residuals)
SOD_final_ins_dist <- na.omit(SOD_final_ins_dist)

p_min_dist_mut_ins <- ggplot(SOD_final_ins_dist, aes(x = min_dist_to_B, y = residuals)) +
  geom_point(aes(color = interface, shape = residual_type), size = 4, alpha = 0.5) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 7.290945, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_residual_int_ins, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_mut_ins


########## Deletions allosteric sites ###
#determine type or residual (< or >)
SOD_final_del$Pos <- as.numeric(SOD_final_del$Pos)
SOD_final_del <- SOD_final_del %>%
  mutate(residual_type = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  ))

SOD_final_del_dist <- left_join(SOD_final_del, min_distance, by = "Pos")

#per position
SOD1_subsitutions_dist_pos <- SOD_final_del_dist %>%
  group_by(Pos) %>%
  mutate(median_residual = abs(median(residuals, na.rm = TRUE))) %>%
  select(ID, median_residual, Pos, interface, min_dist_to_B, location, residual_type)


SOD_final_del_dist_pos <- SOD_final_del_dist_pos %>%
  distinct(Pos, .keep_all = TRUE)
SOD_final_del_dist_pos <- na.omit(SOD_final_del_dist_pos)

#calculate de median residual in interface residues
interface_del <- SOD_final_del %>% 
  filter(interface == TRUE)

median_residual_int_del <- median(abs(interface_del$residuals))


p_min_dist_pos_del <- ggplot(SOD_final_del_dist_pos, aes(x = min_dist_to_B, y = median_residual)) +
  geom_point(aes(color = interface, shape = residual_type), size = 4, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(aes(label = Pos, color = interface), size = 6) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 7.290945, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_residual_int_del, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_pos_del

#per mutation
SOD_final_del_dist$residuals <- abs(SOD_final_del_dist$residuals)
SOD_final_del_dist <- na.omit(SOD_final_del_dist)

p_min_dist_mut_del <- ggplot(SOD_final_del_dist, aes(x = min_dist_to_B, y = residuals)) +
  geom_point(aes(color = interface, shape = residual_type), size = 4, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 7.290945, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = median_residual_int_del, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_mut_del

ggsave(p_min_dist_pos_subs, path = path_fig4, file = "allosteric_sites_subs.tiff", width = 10, height = 10, dpi = 300)
ggsave(p_min_dist_pos_ins, path = path_fig4, file = "allosteric_sites_ins.tiff", width = 10, height = 10, dpi = 300)
ggsave(p_min_dist_pos_del, path = path_fig4, file = "allosteric_sites_del.tiff", width = 10, height = 10, dpi = 300)

ggsave(p_min_dist_mut_subs, path = path_fig4, file = "allosteric_mut_subs.tiff", width = 10, height = 10, dpi = 300)
ggsave(p_min_dist_mut_ins, path = path_fig4, file = "allosteric_mut_ins.tiff", width = 10, height = 10, dpi = 300)
ggsave(p_min_dist_mut_del, path = path_fig4, file = "allosteric_mut_del.tiff", width = 10, height = 10, dpi = 300)

#### Clustering and analysis of allosteric sites ####
SOD1_subsitutions_dist <- SOD1_subsitutions_dist %>% 
  mutate(allosteric_cluster = case_when(
    Pos %in% c(104, 105, 106) ~ "Cluster 1",
    Pos %in% c(117,145,119) ~ "Cluster 2",
    Pos %in% c(84) ~ "Allosteric connector",
    TRUE ~ "rest"
  ))



p_min_dist_pos_subs_clusters <- ggplot(SOD1_subsitutions_dist, aes(x = min_dist_to_B, y = median_residual)) +
  geom_point(aes(color = allosteric_cluster, shape = residual_type), size = 4, alpha = 1) +
  scale_color_manual(values = c("Cluster 1" = "#7C4D79", "Cluster 2" = "#C884B1", "Allosteric connector" = "#EBB5DA", "rest" = "lightgrey")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(aes(label = Pos, color = interface), size = 6) +
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median binding residual|") +
  theme(axis.text = element_text(size = 33),
        axis.title = element_text(size = 36),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 28),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 6.477411, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = 0.11, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 9)),
    shape = guide_legend(override.aes = list(size = 9))
  )

p_min_dist_pos_subs_clusters
ggsave(p_min_dist_pos_subs_clusters, path = path_fig4, file = "allosteric_clusters_subs_points.tiff", width = 10, height = 10, dpi = 300)


SOD1_subsitutions <- left_join(SOD1_subsitutions, min_distance, by = "Pos")


SOD1_subsitutions<- SOD1_subsitutions %>% 
  mutate(allosteric_cluster = case_when(
    Pos %in% c(104, 105, 106) ~ "Cluster 1",
    Pos %in% c(117,145,119) ~ "Cluster 2",
    Pos %in% c(84) ~ "Allosteric connector",
    TRUE ~ "rest"
  )) %>% 
  mutate(allosteric_mutation = case_when(
    residuals > 0.11 & min_dist_to_B > 7.290945 ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(residuals = abs(residuals))


p_dist_clusters <- ggplot(SOD1_subsitutions, aes(x = min_dist_to_B, y = residuals, color = allosteric_cluster)) +
  geom_point(size = 3, alpha = 0.8) + 
  theme_classic() +
  scale_color_manual(values = c("Cluster 1" = "#7C4D79", "Cluster 2" = "#E6B9D9", 
                                "Cluster 3" = "#B99C33", "rest" = "lightgrey")) +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Binding residual|", color = "Allosteric cluster") +
  geom_vline(xintercept = 6.477411, linetype = "dashed", color = "black", size = 1) +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 33),
        legend.text = element_text(size = 25),
        legend.position = "top")
p_dist_clusters

ggsave(p_dist_clusters, path = path_fig4, file = "p_distance_allosteric_clusters.tiff", width = 9, height = 7, dpi = 300)



allosteric_clusters <- SOD1_subsitutions_dist %>% 
  filter(allosteric_cluster != "rest")

allosteric_clusters <- allosteric_clusters %>% 
  group_by(allosteric_cluster, residual_type) %>% 
  summarise(median_res_clstr = median(median_residual)) %>% 
  mutate(allosteric_cluster = factor(allosteric_cluster,
                                     levels = c("Cluster 1","Cluster 2",
                                                "Allosteric connector")))


p_allosteric_clusters <- ggplot(allosteric_clusters, aes(x = allosteric_cluster, y = median_res_clstr, fill = allosteric_cluster)) +
  geom_col(color = "black", linewidth = 1, width = 0.4) +
  scale_fill_manual(values = c("Cluster 1" = "#7C4D79", "Cluster 2" = "#C884B1", "Allosteric connector" = "#EBB5DA")) +
  scale_x_discrete(labels = c(
    "Cluster 1" = "Allosteric\ncluster 1",
    "Cluster 2" = "Allosteric\ncluster 2",
    "Allosteric connector" = "Allosteric\nconnector"
  )) +
  geom_hline(yintercept = 0.11, linetype = "dashed", size = 1) +
  labs(y = "Median residual\nby allosteric cluster", x = "", fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 28),
        legend.position = "none")
p_allosteric_clusters

ggsave(p_allosteric_clusters, path = path_fig4, file = "median_res_allosteric_clusters.tiff", width = 11, height = 8, dpi = 300)


#calculate number of allosteric and non-allosteric mutations in each cluster
count_allosteric_mutations <- SOD1_subsitutions %>%
  count(allosteric_cluster, allosteric_mutation) %>%
  filter(allosteric_cluster != "rest") %>%
  group_by(allosteric_cluster) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()
count_allosteric_mutations$allosteric_cluster <- factor(count_allosteric_mutations$allosteric_cluster,
                                                        levels = c("Cluster 1", "Cluster 2", "Allosteric connector"))


p_mut_clusters <- ggplot(count_allosteric_mutations,aes(x = allosteric_cluster, y = prop, fill = allosteric_mutation)) +
  geom_col(color = "black",linewidth = 1,width = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("TRUE" = "#2F2F2F", "FALSE" = "grey"),
                    labels = c("TRUE" = "Allosteric mutations","FALSE" = "Non-allosteric mutations")) +
  scale_x_discrete(labels = c(
    "Cluster 1" = "Allosteric\ncluster 1",
    "Cluster 2" = "Allosteric\ncluster 2",
    "Allosteric connector" = "Allosteric\nconnector"
  )) +
  labs(x = "", y = "Percentage\nof mutations", fill = "") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 33),
    legend.text = element_text(size = 30),
    legend.position = "none")
p_mut_clusters

p_allostery_bar <- ggarrange(p_allosteric_clusters, 
                             p_mut_clusters,
                             ncol = 2, common.legend = F, align = "h")

p_allostery_bar

ggsave(p_allostery_bar, path = path_fig4, file = "allosteric_clusters_bar.tiff", width = 16, height = 8, dpi = 300)


#allosteric decay
#see the distance of each secondary structure to the dimer interface
SOD_final_subs_dist <- SOD_final_subs_dist %>% 
  mutate(
    s_location = recode(s_location,
                        "Cu binding\nresidues" = "Electrostatic\nloop",
                        "Zn binding\nresidues" = "Zn binding\nloop")) %>% 
  mutate(s_location= factor(s_location, levels = c("N term", "β1", "loop 1", "β2", "loop 2", "β3", "loop 3",
                                                   "β4", "Zn binding\nloop", "β5", "loop 5", "β6", "loop 6",
                                                   "β7", "Electrostatic\nloop", "β8", "loop 8")))

dist_interface <- SOD_final_subs_dist %>% 
  group_by(s_location) %>% 
  mutate(median_dist_pos = median(min_dist_to_B),
         median_residual = median(abs(residuals))) %>% 
  ungroup()

dist_interface <- dist_interface %>%
  mutate(s_location = factor(s_location, 
                             levels = dist_interface %>%
                               group_by(s_location) %>%
                               summarize(median_dist_pos = median(min_dist_to_B, na.rm = TRUE)) %>%
                               arrange(median_dist_pos) %>%
                               pull(s_location)))

dist_interface <- dist_interface %>% 
  mutate(abs_residual = abs(residuals))

region_summary <- dist_interface %>%
  group_by(s_location) %>%
  summarize(median_res = median(abs(residuals), na.rm = TRUE),
            median_dist = median(min_dist_to_B, na.rm = TRUE))

ggplot(dist_interface, aes(x = s_location, y = abs(residuals))) +
  geom_violin(width = 0.8, fill = "lightgrey") +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  geom_point(data = region_summary, aes(x = s_location, y = median_res), color = "red", size = 2) +
  geom_line(data = region_summary, aes(x = s_location, y = median_res, group = 1), color = "red") +
  theme_classic() +
  labs(x = "Region", y = "Magnitude of binding residual") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_dist_to_interface <- ggplot(dist_interface,aes(x = s_location, y = abs_residual, fill = median_residual)) +
  geom_violin(size = 1, width = 2, position = position_dodge(width = 1)) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  scale_fill_gradient2(low = "#2A5676", mid = "white",  high = "#C4D8F3", midpoint = 0.05) +
  theme_classic() +
  labs(x = "", y = "Binding residual", fill = "Minimal dist. to\ndimer interface") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 28),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "top")
p_dist_to_interface

ggsave(p_dist_to_interface, path = path_fig4, file = "sec_struc_dist_to_interface.tiff", width = 15, height = 8, dpi = 300)

p_dist_res <- ggplot(dist_interface, aes(x = min_dist_to_B, y = residuals)) +
  geom_point()
p_dist_res


all_allosteric_sites <- SOD_final_subs_dist %>% 
  filter(Pos %in% c(71,80,82,83,84,86,93,97,79,85,87,95,104,106,108,117,119,138,145))

all_allosteric_sites <- all_allosteric_sites %>% 
  group_by(s_location) %>% 
  mutate(median_dist_pos = median(min_dist_to_B),
         median_residual = abs(median(residuals))) %>% 
  ungroup()

all_allosteric_sites <- all_allosteric_sites %>%
  group_by(s_location) %>%
  mutate(median_location = median(min_dist_to_B, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(s_location = factor(s_location,
                             levels = names(sort(tapply(min_dist_to_B, s_location, median),
                                                 decreasing = TRUE))))


p_decay <- ggplot(all_allosteric_sites, aes(x = s_location, y = min_dist_to_B, fill = median_residual)) +
  geom_violin(drop = FALSE) +
  scale_fill_gradient2(low = "#2A5676", mid = "white",  high = "#C4D8F3", midpoint = 8) +
  theme_classic()
p_decay


##### Allosteric pockets ####
#load pockets prediction by CASTp
pockets_pred <- read.csv("SOD1_pockets_fpocket.csv", sep = ";")

#filter by allosteric positions
SOD_final_subs_dist_pos <- SOD_final_subs_dist_pos %>% 
  mutate(allosteric_pos = case_when(
    median_residual > median_residual_int_subs & interface == FALSE ~ TRUE,
    TRUE ~ FALSE
  ))

allosteric_pos_subs <- SOD_final_subs_dist_pos %>% 
  filter(allosteric_pos == TRUE)

#join allosteric sites with allosteric pockets data frames
pockets_assignment <- inner_join(pockets_pred, SOD_final_subs_dist, by = "Pos")

heatmap_pockets_subs <- pockets_assignment %>%
  group_by(Pocket) %>%
  mutate(Pos_fac = factor(Pos, levels = sort(unique(Pos)))) %>%
  ungroup()




min<-min(heatmap_pockets_subs$binding_score, na.rm = TRUE)
max<-max(heatmap_pockets_subs$binding_score, na.rm = TRUE)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))


min_res<-min(heatmap_pockets_subs$residuals, na.rm = TRUE)
max_res<-max(heatmap_pockets_subs$residuals, na.rm = TRUE)
cols_res <- c(colorRampPalette(c( "#B99C33", "grey95"))((-min_res/(-min_res+max_res)*100)-0.5), colorRampPalette("grey95")(1),
              colorRampPalette(c("grey95",  "#AC7299"), bias=1)((max_res/(-min_res+max_res)*100)-0.5))


p_heatmap_pockets_subs<-ggplot(heatmap_pockets_subs) +
  geom_tile(aes(x = factor(Pos, levels = sort(unique(Pos))), 
                y = factor(Mut_AA, levels = rev(vectorAA)), fill = binding_score),
            color = 'white', size = 0.2) +
  facet_wrap(~Pocket, scales = "free_x") +
  theme_minimal()+
  theme()+
  labs(x="SOD1 position", y="Mutant amino acid", fill="Binding score")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14)) +
  #geom_text(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey") +
  guides(fill = guide_colorbar(
    barwidth = 1,   
    barheight = 5,
    direction = "vertical",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"))
p_heatmap_pockets_subs

p_heatmap_pockets_subs_res<-ggplot(heatmap_pockets_subs) +
  geom_tile(aes(x = factor(Pos, levels = sort(unique(Pos))), 
                y = factor(Mut_AA, levels = rev(vectorAA)), fill = residuals),
            color = 'white', size = 0.2) +
  facet_wrap(~Pocket, ncol = 5, scales = "free_x", labeller = as_labeller(function(x) paste("Pocket", x))) +
  theme_minimal()+
  theme()+
  labs(x="SOD1 WT AA position", y="Mutant amino acid", fill="Binding residual")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 22)) +
  #geom_text(aes(Pos,factor(Mut_AA, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols_res, limits=c(min_res,max_res), na.value = "grey") +
  guides(fill = guide_colorbar(
    barwidth = 1,   
    barheight = 5,
    direction = "vertical",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"))
p_heatmap_pockets_subs_res

ggsave(p_heatmap_pockets_subs_res, path = path_fig4, file = "pockets_heatmaps_residuals.tiff", width = 15, height = 18, dpi = 300)



#calculate the median of binding and residual per position on each pocket
median_pockets <- pockets_assignment %>%
  group_by(Pocket, Pos) %>%
  summarise(
    residual_median = median(residuals, na.rm = TRUE),
    residual_median_abs = median(abs(residuals), na.rm = TRUE),
    binding_median = median(binding_score, na.rm = TRUE),
    binding_median_abs = median(abs(binding_score), na.rm = TRUE)
  ) %>%
  mutate(Pos_fac = factor(Pos, levels = sort(unique(Pos)))) %>%
  ungroup() %>%
  mutate(row = "row")


median_pockets$Pos <- factor(
  median_pockets$Pos,
  levels = as.character(sort(as.numeric(unique(median_pockets$Pos))))
)

p_median_res_pockets <- ggplot(median_pockets, aes(x = factor(Pos, levels = sort(unique(Pos))), y = residual_median_abs)) +
  geom_col(fill = "darkgrey", color = "black") +
  facet_wrap(~Pocket,  ncol = 5, scales = "free_x", labeller = as_labeller(function(x) paste("Pocket", x))) +
  geom_hline(yintercept = median_residual_int_subs, linetype = "dashed", size = 1) +
  labs(x="SOD1 WT AA position", y="Mutant amino acid")+
  theme_classic()
p_median_res_pockets

ggsave(p_median_res_pockets, path = path_fig4, file = "pockets_median_residuals.tiff", width = 10, height = 8, dpi = 300)



#check if allosteric clusters are in any allosteric pocket
cluster_pocket <- inner_join(allosteric_clusters, pockets_pred, by = "Pos")



#get outliers to define allosteric threshold
#IQR to detect outliers
Q1 <- quantile(SOD1_insertions_dist_pos$median_residual, 0.25, na.rm = TRUE)
Q3 <- quantile(SOD1_insertions_dist_pos$median_residual, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

SOD1_insertions_dist_pos$outlier_decay_iqr <- SOD1_insertions_dist_pos$median_residual < lower_bound | 
  SOD1_insertions_dist_pos$median_residual > upper_bound


#plot residual vs min distance
p_min_dist_pos_ins <- ggplot(SOD1_insertions_dist_pos, aes(x = min_dist_to_B, y = median_residual)) +
  geom_point(aes(color = interface, shape = residual_type), size = 3, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(aes(label = Pos, color = interface), size = 4) +  # evita solapamiento
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median abundance-binding residual|") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 7.290945, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = 0.1667026, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    shape = guide_legend(override.aes = list(size = 6))
  )

p_min_dist_pos_ins

p_min_dist_pos_ins <- ggMarginal(
  p_min_dist_pos_ins,
  type = "density",     # densidad en lugar de histograma
  fill = "grey70",      # color de relleno
  alpha = 0.6
)

p_min_dist_pos_ins


#same for deletions
#join with dms data
SOD_final_del$Pos <- as.numeric(SOD_final_del$Pos)
SOD_final_del <- SOD_final_del %>%
  mutate(residual_type = case_when(
    residuals > 0 ~ "positive",
    residuals < 0 ~ "negative"
  ))

SOD1_deletions_dist <- left_join(SOD_final_del, min_distance, by = "Pos")

#per position
SOD1_deletions_dist_pos <- SOD1_deletions_dist %>%
  group_by(Pos) %>%
  mutate(residuals = abs(residuals)) %>%
  select(ID, residuals, Pos, interface, min_dist_to_B, location, residual_type)


SOD1_deletions_dist_pos <- SOD1_deletions_dist_pos %>%
  distinct(Pos, .keep_all = TRUE)
SOD1_deletions_dist_pos <- na.omit(SOD1_deletions_dist_pos)

#get outliers to define allosteric threshold
#IQR to detect outliers
Q1 <- quantile(SOD1_deletions_dist_pos$residuals, 0.25, na.rm = TRUE)
Q3 <- quantile(SOD1_deletions_dist_pos$residuals, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

SOD1_deletions_dist_pos$outlier_decay_iqr <- SOD1_deletions_dist_pos$residuals < lower_bound | 
  SOD1_deletions_dist_pos$residuals > upper_bound


#plot residual vs min distance
p_min_dist_pos_del <- ggplot(SOD1_deletions_dist_pos, aes(x = min_dist_to_B, y = residuals)) +
  geom_point(aes(color = interface, shape = residual_type), size = 3, alpha = 1) +
  scale_color_manual(
    values = c("TRUE" = "#2A5676", "FALSE" = "#9E3D22"),
    labels = c("TRUE" = "Interface residue", "FALSE" = "Non-interface residue")) +
  scale_shape_manual(
    values = c("positive" = 16, "negative" = 17),
    labels = c("positive" = "Residual > 0", "negative" = "Residual < 0")
  ) +
  geom_text_repel(aes(label = Pos, color = interface), size = 4) +  # evita solapamiento
  theme_classic() +
  labs(x = "Minimal distance to dimer interface (Å)", y = "|Median abundance-binding residual|") +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.title = element_blank(),
        legend.title.align = 0.5,
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.box = "vertical") +
  geom_vline(xintercept = 7.290945, linetype = "dashed", color = "darkgrey", size = 1) +
  geom_hline(yintercept = 0.2547288, linetype = "dashed", color = "darkgrey", size = 1) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    shape = guide_legend(override.aes = list(size = 6))
  )

p_min_dist_pos_del

p_min_dist_pos <- ggMarginal(
  p_min_dist_pos,
  type = "density",     # densidad en lugar de histograma
  fill = "grey70",      # color de relleno
  alpha = 0.6
)

p_min_dist_pos



#### Find outlier mutations####
#model with all the dataset (subs, ins an del)
#for substitutions
SOD1_subsitutions_interface <- SOD1_subsitutions %>%
  filter(side_chain == "dimer interface")

focal_points <- data.frame(
  abundance_score = c(-1,0),
  binding_score = c(-1,0),
  weights = c(1e04, 1e04)
)

# Añade tus datos reales con peso normal = 1
SOD1_subsitutions_interface <- SOD1_subsitutions_interface %>%
  mutate(weights = 1) %>%
  bind_rows(focal_points)

# Ajusta el modelo con los pesos correctos
model <- loess(binding_score ~ abundance_score,
               data = SOD1_subsitutions_interface,
               span = 0.5,
               weights = SOD1_subsitutions_interface$weights)

# Predicciones del modelo loess
pred <- data.frame(
  abundance_score = seq(min(SOD1_subsitutions_interface$abundance_score),
                        max(SOD1_subsitutions_interface$abundance_score),
                        length.out = 200)
)
pred$binding_pred <- predict(model, newdata = pred)

SOD1_subsitutions_interface$fitted <- fitted(model)
SOD1_subsitutions_interface$residuals <- SOD1_subsitutions_interface$binding_score - SOD1_subsitutions_interface$fitted

max_resid <- max(SOD1_subsitutions_interface, na.rm = TRUE)

desvio_estandar <- sd(SOD1_subsitutions_interface$residuals)
SOD1_subsitutions_interface <- SOD1_subsitutions_interface %>%
  mutate(outlier = case_when(
    abs(residuals) > (2 * desvio_estandar) & binding_score > -0.85 ~ TRUE,
    TRUE ~ FALSE
  ))

outliers <- SOD1_subsitutions_interface %>% filter(outlier)




p_subs_interface <- ggplot(SOD1_subsitutions_interface, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = outlier), size = 5, shape=16) +
  scale_color_manual(values = c("TRUE" = "#B254A5", "FALSE" = "#D5D5D5"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.5) +
  coord_cartesian(xlim = range(SOD1_subsitutions_interface$abundance_score),
                  ylim = range(SOD1_subsitutions_interface$binding_score)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(x = "Normalized abundance score", y = "Normalized binding score", color = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 33),
        legend.position = "top") +
  guides(color = guide_legend(
    override.aes = list(size = 13)))

p_subs_interface

ggsave(p_subs_interface, path = path_fig3, file = "interface_correlation_subs.tiff", width = 8, height = 7, dpi = 300)


#for insertions
SOD1_insertions_interface <- SOD1_insertions %>%
  filter(side_chain == "dimer interface")

focal_points <- data.frame(
  abundance_score = c(-1,0),
  binding_score = c(-1,0),
  weights = c(1e04, 1e04)
)

SOD1_insertions_interface <- SOD1_insertions_interface %>%
  mutate(weights = 1) %>%
  bind_rows(focal_points)

model <- loess(binding_score ~ abundance_score,
               data = SOD1_insertions_interface,
               span = 0.5,
               weights = SOD1_insertions_interface$weights)

pred <- data.frame(
  abundance_score = seq(min(SOD1_insertions_interface$abundance_score),
                        max(SOD1_insertions_interface$abundance_score),
                        length.out = 200)
)
pred$binding_pred <- predict(model, newdata = pred)

SOD1_insertions_interface$fitted <- fitted(model)
SOD1_insertions_interface$residuals <- SOD1_insertions_interface$binding_score - SOD1_insertions_interface$fitted

max_resid <- max(SOD1_insertions_interface, na.rm = TRUE)

desvio_estandar <- sd(SOD1_insertions_interface$residuals)
SOD1_insertions_interface <- SOD1_insertions_interface %>%
  mutate(outlier = case_when(
    abs(residuals) > (2 * desvio_estandar) & binding_score > -0.85 ~ TRUE,
    TRUE ~ FALSE
  ))

outliers <- SOD1_insertions_interface %>% filter(outlier)

p_ins_interface <- ggplot(SOD1_insertions_interface, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = outlier), size = 5, shape=16) +
  scale_color_manual(values = c("TRUE" = "#B254A5", "FALSE" = "#D5D5D5"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.5) +
  coord_cartesian(xlim = range(SOD1_subsitutions_interface$abundance_score),
                  ylim = range(SOD1_subsitutions_interface$binding_score)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(x = "Normalized abundance score", y = "Normalized binding score", color = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 33),
        legend.position = "top") +
  guides(color = guide_legend(
    override.aes = list(size = 13)))

p_ins_interface

ggsave(p_subs_interface, path = path_fig3, file = "interface_correlation_ins.tiff", width = 8, height = 7, dpi = 300)


#for deletions
SOD1_deletions_interface <- SOD1_deletions %>%
  filter(side_chain == "dimer interface")

focal_points <- data.frame(
  abundance_score = c(-1,0),
  binding_score = c(-1,0),
  weights = c(1e04, 1e04)
)

SOD1_deletions_interface <- SOD1_deletions_interface %>%
  mutate(weights = 1) %>%
  bind_rows(focal_points)

model <- loess(binding_score ~ abundance_score,
               data = SOD1_deletions_interface,
               span = 0.5,
               weights = SOD1_deletions_interface$weights)

pred <- data.frame(
  abundance_score = seq(min(SOD1_deletions_interface$abundance_score),
                        max(SOD1_deletions_interface$abundance_score),
                        length.out = 200)
)
pred$binding_pred <- predict(model, newdata = pred)

SOD1_deletions_interface$fitted <- fitted(model)
SOD1_deletions_interface$residuals <- SOD1_deletions_interface$binding_score - SOD1_deletions_interface$fitted

max_resid <- max(SOD1_deletions_interface, na.rm = TRUE)

desvio_estandar <- sd(SOD1_deletions_interface$residuals)
SOD1_deletions_interface$outlier <- abs(SOD1_deletions_interface$residuals) > (2 * desvio_estandar)
outliers <- SOD1_deletions_interface %>% filter(outlier)

p_del_interface <- ggplot(SOD1_deletions_interface, aes(x = abundance_score, y = binding_score)) +
  geom_point(aes(color = outlier), size = 5, shape=16) +
  scale_color_manual(values = c("TRUE" = "#B254A5", "FALSE" = "#D5D5D5"),
                     labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.5) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(x = "Normalized abundance score", y = "Normalized binding score", color = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 33),
        legend.position = "top") +
  guides(color = guide_legend(
    override.aes = list(size = 13)))

p_del_interface

ggsave(p_del_interface, path = path_fig3, file = "interface_correlation_del.tiff", width = 8, height = 7, dpi = 300)

#combine the three plots
x_range <- range(
  c(SOD1_subsitutions_interface$abundance_score,  
    SOD1_insertions_interface$abundance_score,
    SOD1_deletions_interface$abundance_score),
  na.rm = TRUE
)

y_range <- range(
  c(SOD1_subsitutions_interface$binding_score,  
    SOD1_insertions_interface$binding_score,
    SOD1_deletions_interface$binding_score),
  na.rm = TRUE
)

p_subs_interface <- p_subs_interface +
  coord_cartesian(xlim = x_range, ylim = y_range)

p_ins_interface <- p_ins_interface +
  coord_cartesian(xlim = x_range, ylim = y_range)

p_del_interface <- p_del_interface +
  coord_cartesian(xlim = x_range, ylim = y_range)

ggsave(p_subs_interface, path = path_fig3, file = "interface_correlation_subs.tiff", width = 8, height = 7, dpi = 300)
ggsave(p_ins_interface, path = path_fig3, file = "interface_correlation_ins.tiff", width = 8, height = 7, dpi = 300)
ggsave(p_del_interface, path = path_fig3, file = "interface_correlation_del.tiff", width = 8, height = 7, dpi = 300)


p_combined_interface <- ggarrange(p_subs_interface + theme(axis.title = element_blank()),
                                  p_ins_interface + theme(axis.title.x = element_blank()),
                                  p_del_interface + theme(axis.title.y = element_blank()),
                                  nrow = 3, align = "v", common.legend = TRUE)

p_combined_interface

ggsave(p_combined_interface, path = path_fig3, file = "interface_correlation_subs_indels.tiff", width = 8, height = 10, dpi = 300)

medianas_por_pos <- SOD1_subsitutions_interface %>%
  group_by(Pos) %>%
  summarise(mediana_binding = median(binding_score, na.rm = TRUE))

SOD1_subsitutions_interface <- SOD1_subsitutions_interface  %>%
  mutate(Pos = factor(Pos, levels = sort(unique(as.numeric(Pos)))),
         category_fdr_bpca = factor(category_fdr_bpca, 
                                    levels = c("FS_inc", "WT-like", "FS_dec", "Stop-like", "FS_dec_stop")),
         category_fdr_bpca = recode(category_fdr_bpca, "FS_dec_stop" = "Stop-like"))



p_interface_median <- ggplot(SOD1_subsitutions_interface, aes(x = Pos, y = binding_score)) +
  geom_violin(scale = "width", trim = TRUE, size = 1.5) +
  geom_jitter(width = 0.1, aes(color = category_fdr_bpca), alpha = 1, size = 3) +
  scale_color_manual(values = c("FS_inc" = "darkgreen", 
                                "WT-like" = "lightgrey", 
                                "FS_dec" = "#F69541", 
                                "Stop-like" = "#E44244"),
                     labels = c("FS_inc" = "Increased", 
                                "WT-like" = "WT-like", 
                                "FS_dec" = "Decreased", 
                                "Stop-like" = "Stop-like")) +
  labs(x = "Dimer interface position", y = "Binding score", color = "Binding FDR = 0.1") +
  geom_crossbar(data = medianas_por_pos,
                aes(x = factor(Pos), y = mediana_binding,
                    ymin = mediana_binding, ymax = mediana_binding),
                color = "black", width = 0.5, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 28),
        legend.position = "top"
  )


p_interface_median


medianas_por_pos_ins <- SOD1_insertions_interface %>%
  group_by(Pos) %>%
  summarise(mediana_binding = median(binding_score, na.rm = TRUE))

SOD1_insertions_interface <- SOD1_insertions_interface  %>%
  mutate(Pos = factor(Pos, levels = sort(unique(as.numeric(Pos)))),
         category_fdr_bpca = factor(category_fdr_bpca, 
                                    levels = c("FS_inc", "WT-like", "FS_dec", "Stop-like", "FS_dec_stop")),
         category_fdr_bpca = recode(category_fdr_bpca, "FS_dec_stop" = "Stop-like"))



p_interface_ins_median <- ggplot(SOD1_insertions_interface, aes(x = Pos, y = binding_score)) +
  geom_violin(scale = "width", trim = TRUE, size = 1.5) +
  geom_jitter(width = 0.1, aes(color = category_fdr_bpca), alpha = 1, size = 3) +
  scale_color_manual(values = c("FS_inc" = "darkgreen", 
                                "WT-like" = "lightgrey", 
                                "FS_dec" = "#F69541", 
                                "Stop-like" = "#E44244"),
                     labels = c("FS_inc" = "Increased", 
                                "WT-like" = "WT-like", 
                                "FS_dec" = "Decreased", 
                                "Stop-like" = "Stop-like")) +
  labs(x = "Dimer interface position", y = "Binding score", color = "Binding FDR = 0.1") +
  geom_crossbar(data = medianas_por_pos_ins,
                aes(x = factor(Pos), y = mediana_binding,
                    ymin = mediana_binding, ymax = mediana_binding),
                color = "black", width = 0.5, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 28),
        legend.position = "top"
  )


p_interface_ins_median

SOD1_deletions_interface <- na.omit(SOD1_deletions_interface)

SOD1_deletions_interface <- SOD1_deletions_interface  %>%
  mutate(Pos = factor(Pos, levels = sort(unique(as.numeric(Pos)))),
         category_fdr_bpca = factor(category_fdr_bpca, 
                                    levels = c("FS_inc", "WT-like", "FS_dec", "Stop-like", "FS_dec_stop")),
         category_fdr_bpca = recode(category_fdr_bpca, "FS_dec_stop" = "Stop-like"))

p_interface_deletions <- ggplot(SOD1_deletions_interface, aes(x = Pos,  y = binding_score)) +
  geom_col(aes(fill = category_fdr_bpca), linewidth = 1, color = "black", width = 0.5) +
  scale_fill_manual(values = c("FS_inc" = "darkgreen", 
                               "WT-like" = "lightgrey", 
                               "FS_dec" = "#F69541", 
                               "Stop-like" = "#E44244"),
                    labels = c("FS_inc" = "Increased", 
                               "WT-like" = "WT-like", 
                               "FS_dec" = "Decreased", 
                               "Stop-like" = "Stop-like")) +
  labs(x = "Dimer interface position", y = "Binding score", fill = "Binding FDR = 0.1") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_classic() +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 33),
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 28),
        legend.position = "none"
  )
p_interface_deletions


p_interface_median <- p_interface_median +
  coord_cartesian( ylim = y_range)

p_interface_ins_median <- p_interface_ins_median +
  coord_cartesian(ylim = y_range)

p_interface_deletions <- p_interface_deletions +
  coord_cartesian(ylim = y_range)


p_interface_all <- ggarrange(p_interface_median + theme(axis.text.x = element_blank(), axis.title= element_blank()),
                             p_interface_ins_median + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                             p_interface_deletions + theme(axis.title.y = element_blank()),
                             nrow = 3, align = "v", common.legend = TRUE)
p_interface_all

ggsave(p_interface_all, path = path_fig3, file = "interface_violin_all.tiff", width = 24, height = 15, dpi = 300)

SOD1_subsitutions_interface <- na.omit(SOD1_subsitutions_interface)

SOD1_subsitutions_interface_enrich <- SOD1_subsitutions_interface %>%
  group_by(Pos, outlier) %>%
  summarise(n_mut = n(), .groups = "drop_last") %>%
  mutate(total_mut = sum(n_mut)) %>%
  ungroup() %>%
  mutate(enrichment = n_mut / total_mut) %>%
  mutate(Pos = factor(Pos, levels = sort(unique(as.numeric(Pos)))))



p_enrich_subs <- ggplot(SOD1_subsitutions_interface_enrich,
                        aes(x = Pos, y = enrichment, fill = outlier)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", linewidth = 1.2, width = 0.7) +
  labs(x = "Dimer interface position", y = "Fraction of mutations") +
  scale_fill_manual(values = c("TRUE" = "#B254A5", "FALSE" = "lightgrey"),
                    labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 33),
    axis.text = element_text(size = 30),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.position = "top"
  )
p_enrich_subs

ggsave(p_number_mut_subs, path = path_fig3, file = "interface_correlation_subs_nmut.tiff", width = 9, height = 7, dpi = 300)


SOD1_insertions_interface <- na.omit(SOD1_insertions_interface)

SOD1_insertions_interface_enrich <- SOD1_insertions_interface %>%
  group_by(Pos, outlier) %>%
  summarise(n_mut = n(), .groups = "drop_last") %>%
  mutate(total_mut = sum(n_mut)) %>%
  ungroup() %>%
  mutate(enrichment = n_mut / total_mut) %>%
  mutate(Pos = factor(Pos, levels = sort(unique(as.numeric(Pos)))))



p_enrich_ins <- ggplot(SOD1_insertions_interface_enrich,
                       aes(x = Pos, y = enrichment, fill = outlier)) +
  geom_col(position = position_dodge(width = 0.8),
           color = "black", linewidth = 1.2, width = 0.7) +
  labs(x = "Dimer interface position", y = "Fraction of mutations") +
  scale_fill_manual(values = c("TRUE" = "#B254A5", "FALSE" = "lightgrey"),
                    labels = c("TRUE" = "Outlier", "FALSE" = "Non-outlier")) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 33),
    axis.text = element_text(size = 30),
    legend.title = element_blank(),
    legend.text = element_text(size = 30),
    legend.position = "top"
  )
p_enrich_ins

ggsave(p_number_mut_ins, path = path_fig3, file = "interface_correlation_ins_nmut.tiff", width = 9, height = 7, dpi = 300)


outliers_del<- SOD1_deletions_interface %>%
  filter(outlier == "TRUE" & binding_score > -0.8) %>%
  group_by(Pos) %>%   
  summarize(n_mut = n(),
            mean_resid = mean(residuals),
            max_resid = max(residuals),
            min_resid = min(residuals)) %>%
  mutate(Pos = factor(Pos, levels = sort(unique(as.numeric(Pos)))))



p_number_mut_del <- ggplot(outliers_del, aes(x = Pos, y = n_mut)) +
  geom_col(fill = "#B254A5", color = "black", width = 0.3, linewidth  = 1.2) +
  labs(x = "Dimer interface position", y = "Number of mutations") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30))
p_number_mut_del

ggsave(p_number_mut_del, path = path_fig3, file = "interface_correlation_del_nmut.tiff", width = 9, height = 7, dpi = 300)



SOD1_subsitutions_interface <- na.omit(SOD1_subsitutions_interface)
SOD1_subsitutions_interface <- SOD1_subsitutions_interface %>%
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
    WT_AA %in% "P" ~ "Proline"
  )) %>%
  mutate(residual_type = case_when(residuals < 0 ~ "Negative residual",
                                   residuals > 0 ~ "Positive residual"))

outlier_int_subs <- SOD1_subsitutions_interface %>%
  filter(outlier == "TRUE")

no_outlier_int_subs <- SOD1_subsitutions_interface %>%
  filter(outlier == "FALSE")


prop_df <- SOD1_subsitutions_interface %>%
  group_by(outlier, Property_wt) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(outlier) %>%
  mutate(prop = n / sum(n))

p_int_subs_aatype <- ggplot(prop_df, aes(x = Property_wt, y = prop)) +
  geom_col() +
  facet_wrap(~outlier) +
  labs(y = "Proporción", x = "Tipo de mutación") +
  theme_classic()
p_int_subs_aatype


p_int_subs_aatype <- ggplot(outlier_int_subs, aes(x = Property_wt)) +
  geom_bar() +
  facet_wrap(~residual_type)
p_int_subs_aatype




