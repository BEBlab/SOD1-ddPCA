#load required packages
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(patchwork)

#define path
path_fig2 <- "C:/Users/tquiroga/OneDrive - IBEC/Projects/From owncloud/SOD1 binding PCA/SOD1 manuscript/Plots/Figure 2"
#path_fig2 <- "C:/Users/USUARIO/OneDrive/OneDrive - IBEC/Projects/From owncloud/SOD1 binding PCA/SOD1 manuscript/Plots/Figure 2"


#load aPCA and bPCA data
SOD_final <- read.csv("SOD1_abundance_binding_scores.csv")

SOD_final <- SOD_final %>%
  mutate(mutation_type = recode(mutation_type, "syn" = "subs"
  )) %>%
  mutate(WT_AA = str_extract(ID, "^[A-za-z]"),
         Pos = str_extract(ID, "\\d+"),
         Mut = str_extract(ID, "[A-za-z]$")) 


#Abundance heatmaps ####
##Single substitutions ####
#filter to work with substitutions and stops
SOD_final_subs <- SOD_final %>%
  filter(mutation_type %in% c("subs", "stop"))

SOD_final_subs <- SOD_final_subs %>%
  mutate(Mut = case_when(
    mutation_type == "stop" ~ "stop",
    TRUE ~ Mut
  )) %>%
  mutate(Pos = as.numeric(str_extract(ID, "\\d+")))

#add amino acid side-chain information
predictions <- read.csv("SOD1_rSASA.csv", sep = ";")
predictions$rSASA.monomer <- readr::parse_number(predictions$rSASA.monomer)
predictions$rSASA.dimer <- readr::parse_number(predictions$rSASA.dimer)
predictions$Pos <- readr::parse_number(predictions$Pos)

predictions <- predictions %>%
  mutate(
    side_chain_monomer = case_when(
      rSASA.monomer <=25 ~ "core",
      rSASA.monomer >25 ~ "surface")) %>% 
  mutate(
    side_chain_dimer = case_when(
      rSASA.dimer <=25 ~ "core",
      rSASA.dimer > 25 ~ "surface"
    )
  ) %>%
  mutate(side_chain_monomer = case_when(
    Pos %in% c(71, 80, 83) ~ "Zn binding",
    Pos %in% c(46, 48, 120) ~ "Cu binding",
    Pos == 63 ~ "Zn and Cu binding",
    TRUE ~ side_chain_monomer)) %>%
  mutate(side_chain_dimer = case_when(
    Pos %in% c(71, 80, 83) ~ "Zn binding",
    Pos %in% c(46, 48, 120) ~ "Cu binding",
    Pos == 63 ~ "Zn and Cu binding",
    Pos %in% c(5,7,50,51,52,53,54,113,114,148,150,151,152,153) ~ "dimer interface",
    TRUE ~ side_chain_dimer))
predictions$row <- "row"



#prepare data for heatmap
peptide_seq<-'ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "stop")

#add syn
positions <- 1:153
syn.df <- data.frame(
  "WT_AA" = peptide_seq,
  "Mut" = peptide_seq,
  "Pos" = positions,
  "abundance_score" = 0,
  "ID" = "syn"
)
heatmap_df<-rbind(SOD_final_subs[,c("WT_AA", "Mut", "Pos", "abundance_score", "ID")], syn.df)
heatmap_df <- heatmap_df %>%
  mutate(Mut = if_else(is.na(Mut), "stop", Mut))
heatmap_df$label<-""
heatmap_df[heatmap_df$ID=="syn",]$label<-"WT"
heatmap_df <- heatmap_df %>%
  complete(Pos = 1:153, Mut = vectorAA, fill = list(abundance_score = NA, WT_AA = NA, ID = NA, label = ""))

min<-min(heatmap_df$abundance_score, na.rm = TRUE)
max<-max(heatmap_df$abundance_score, na.rm = TRUE)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

p_heatmap<-ggplot(heatmap_df)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=abundance_score), color='white', size=0.2)+
  scale_x_continuous(breaks=seq(1,153), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  theme()+
  labs(x="SOD1 WT amino acid position", y="Mutant amino acid", fill="Abundance\nscore")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text(aes(Pos,factor(Mut, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-1, 0, 1), labels = c(-1, 0, 1)) +
  guides(fill = guide_colorbar(
    barwidth = 5,   
    barheight = 1,
    direction = "horizontal",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"
  )) 
p_heatmap

#violin plot
heatmap_df <- na.omit(heatmap_df)
heatmap_df <- heatmap_df %>%
  filter(Mut != "stop")
subsp_median_df<-as.data.frame(heatmap_df %>% group_by(Pos) %>% dplyr::summarise(median_p=median(abundance_score)))
heatmap_df<-left_join(heatmap_df, subsp_median_df)

p_violin_p <- ggplot(heatmap_df, aes(x = as.factor(Pos), y = abundance_score)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_violin(scale = "width", size = 0.2, aes(fill = median_p)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.2) +
  scale_x_discrete(breaks = seq(1, 153), labels = peptide_seq_pos, expand = c(0, 0)) +
  theme_classic() +
  labs(x = "SOD1 WT amino acid position", y = "Median\nabundance score", fill = "Median FS") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(0, 0, 0.5, 0), 'cm'),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = cols, limits = c(min, max), breaks = c(-4, -2, 0, 2, 4)) +
  theme(legend.position = "none")
p_violin_p

#add side-chain tiles
p_side_chain <- ggplot(predictions, aes(x = factor(Pos), y = row, fill = side_chain_monomer)) +
   geom_tile(color = "white", height = 0.1, size = 0.2) +
   scale_fill_manual(
     values = c(
       "core" = "#26456E",
       "surface" = "#C4D8F3",
       "dimer interface" = "#41B7C4",
       "Zn binding" = "#BB173A",
       "Cu binding" = "#DECC61",
       "Zn and Cu binding" = "#B254A5")
   ) +
   scale_y_discrete(expand = c(0,0)) + 
   theme(
     legend.position = "none",
     axis.text = element_blank(),
     axis.title = element_blank(),
     axis.ticks = element_blank(),
     plot.margin = unit(c(0,0,0.5,0), 'cm')
   )
p_side_chain

p_heatmap <- p_heatmap + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

#combine heatmap with side-chain tiles 
p_heatmap_side_chain <- ggarrange(p_side_chain + theme(axis.title = element_blank(), axis.text = element_blank()),
                                   p_heatmap,
                                   p_violin_p,
                                   heights = c(0.05, 1, 0.2), nrow = 3, align = "v")
p_heatmap_side_chain

ggsave(p_heatmap_side_chain, path = path_fig2, file="SOD1_aPCA_full_heatmap_violin_subs.tiff",width=30, height=10, dpi = 600)

## Single insertions ####
#filter to work with insertions
SOD_final_ins <- SOD_final %>%
  filter(mutation_type %in% c("ins", "stop"))

#add position 1. We have not mutated it but we show it as NA
aa_vector <- sort(unique(SOD_final_ins$Mut))  
missing_pos1 <- data.frame(
  ID = paste0("ins1", aa_vector),
  Pos = 1,
  Mut = aa_vector,
  abundance_score = NA
)

SOD_final_ins$Pos <- as.numeric(SOD_final_ins$Pos)
SOD_final_ins <- bind_rows(SOD_final_ins, missing_pos1)

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "stop")

#add syn
positions <- 1:153  
SOD_final_ins <- SOD_final_ins %>%
  mutate(Pos = as.numeric(Pos)) %>%
  mutate(WT_AA = peptide_seq[Pos],
         Mut = case_when(
           mutation_type == "stop" ~ "stop",
           TRUE ~ Mut
         )) 
heatmap_df_apca_ins<-SOD_final_ins

heatmap_df_apca_ins <- heatmap_df_apca_ins %>%
  complete(Pos = 1:153, Mut = vectorAA, fill = list(abundance_score = NA, WT_AA = NA, ID = NA, label = ""))

min<-min(heatmap_df_apca_ins$abundance_score, na.rm = TRUE)
max<-max(heatmap_df_apca_ins$abundance_score, na.rm = TRUE)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

p_heatmap_ins_apca<-ggplot(heatmap_df_apca_ins)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=abundance_score), color='white', size=0.2)+
  scale_x_continuous(breaks=positions, labels = positions, expand=c(0,0))+
  theme_minimal()+
  theme()+
  labs(x="Position of inserted position", y="Inserted amino acid", fill="Abundance score")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-1, 0, 1), labels = c(-1, 0, 1)) +
  guides(fill = guide_colorbar(
    barwidth = 5,   
    barheight = 1,
    direction = "horizontal",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"
  )) 

p_heatmap_ins_apca

#violinplot
heatmap_df_apca_ins <- heatmap_df_apca_ins %>%
  filter(Mut != "stop")

inssp_median_df_bpca <- heatmap_df_apca_ins %>%
  group_by(Pos) %>%
  summarise(median_p = median(abundance_score, na.rm = TRUE))

all_positions <- 1:max(heatmap_df_apca_ins$Pos)
heatmap_df_apca_ins <- heatmap_df_apca_ins %>%
  right_join(data.frame(Pos = all_positions), by = "Pos") %>%
  left_join(inssp_median_df_bpca, by = "Pos")  


p_violin_ins_apca <- ggplot(heatmap_df_apca_ins, aes(x = as.factor(Pos), y = abundance_score)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_violin(scale = "width", size = 0.2, aes(fill = median_p)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.2) +
  scale_x_discrete(breaks = positions, labels = positions, expand = c(0, 0)) +
  theme_classic() +
  labs(x = "Position of inserted amino acid", y = "Median\nabundance score", fill = "Median FS") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(0, 0, 0.5, 0), 'cm'),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = cols, limits = c(min, max), breaks = c(-4, -2, 0, 2, 4)) +
  theme(legend.position = "none")

p_violin_ins_apca


#combine insertions heatmap, violin and side-chain tiles
p_heatmap_ins_side_chain <- ggarrange(p_side_chain + theme(axis.title = element_blank(), axis.text = element_blank()),
                                      p_heatmap_ins_apca + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                                      p_violin_ins_apca,
                                      heights = c(0.05, 1, 0.2), nrow = 3, align = "v")

p_heatmap_ins_side_chain
ggsave(p_heatmap_ins_side_chain,, path = path_fig2, file="SOD1_aPCA_full_heatmap_violin_ins.tiff",width=30, height=10, dpi = 600)


## Single deletions ####
#filter to wotk with deletions
SOD_final_del <- SOD_final %>%
  filter(mutation_type=="del")


SOD_final_del <- SOD_final_del %>%
  separate_rows(ID, sep = ";") %>%
  mutate(
    del_pos = str_extract(ID, "\\d+"),
    del_aa = str_extract(ID, "^[A-Z]")
  )

del_apca <- SOD_final_del %>%
  mutate(assay = "apca",
         fitness = abundance_score,
         sigma = abundance_sigma)

del_bpca <- SOD_final_del %>%
  mutate(assay = "bpca",
         fitness = binding_score,
         sigma = binding_sigma)

SOD_final_del <- rbind(del_apca, del_bpca)


#prepare data for plot
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

heatmap_df_del<-(SOD_final_del)

min<-min(heatmap_df_del$abundance_score)
max<-max(heatmap_df_del$abundance_score)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

heatmap_df <- heatmap_df_del %>%
  distinct(ID, .keep_all = TRUE)

p_single_deletion_apca <- ggplot(heatmap_df_del, aes(x = factor(del_pos, levels = 1:153), y = abundance_score)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_errorbar(aes(ymin = abundance_score - 1.96 * abundance_sigma, ymax = abundance_score + 1.96 * abundance_sigma),
                width = 0, size = 0.1) +
  geom_point(data = heatmap_df_del,
             aes(fill = abundance_score), size = 6, shape = 21, stroke = 1.2) +
  labs(y = "Abundance score", x = "Position of deleted amino acid", fill = "Abundance/binding\nscore") +
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

p_single_deletion_apca

#combine deletions plot with side-chain tiles
p_del_apca <- ggarrange(p_side_chain + theme(axis.title = element_blank(), axis.text = element_blank()),
                        p_single_deletion_apca,
                        heights = c(0.1, 1), nrow = 2, align = "v")

p_del_apca
ggsave(p_del_apca , path = path_fig2, file="SOD1_aPCA_full_map_deletions.tiff",width=30, height=5, dpi = 600)


#Heterodimerization heatmaps ####
##Single substitutions ####
#filter to work with substitutions and stops
peptide_seq<-'ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "stop")

#add syn
positions <- 1:153
syn.df <- data.frame(
  "WT_AA" = peptide_seq,
  "Mut" = peptide_seq,
  "Pos" = positions,
  "binding_score" = 0,
  "ID" = "syn"
)
heatmap_binding <- rbind(SOD_final_subs[,c("WT_AA", "Mut", "Pos", "binding_score", "ID")], syn.df)
heatmap_binding <- heatmap_binding %>%
  mutate(Mut = if_else(is.na(Mut), "stop", Mut))
heatmap_binding$label<-""
heatmap_binding[heatmap_binding$ID=="syn",]$label<-"WT"
heatmap_binding <- heatmap_binding %>%
  complete(Pos = 1:153, Mut = vectorAA, fill = list(binding_score = NA, WT_AA = NA, ID = NA, label = ""))

min<-min(heatmap_binding$binding_score, na.rm = TRUE)
max<-max(heatmap_binding$binding_score, na.rm = TRUE)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

p_heatmap_binding<-ggplot(heatmap_binding)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=binding_score), color='white', size=0.2)+
  scale_x_continuous(breaks=seq(1,153), labels = peptide_seq_pos, expand=c(0,0))+
  theme_minimal()+
  theme()+
  labs(x="SOD1 WT amino acid position", y="Mutant amino acid", fill="Heterodim\nscore")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  geom_text(aes(Pos,factor(Mut, levels=rev(vectorAA)), label=label), color="black", size=2)+
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-1, 0, 1), labels = c(-1, 0, 1)) +
  guides(fill = guide_colorbar(
    barwidth = 5,   
    barheight = 1,
    direction = "horizontal",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"
  )) 
p_heatmap_binding

#violin plot
heatmap_binding <- na.omit(heatmap_binding)
heatmap_binding <- heatmap_binding %>%
  filter(Mut != "stop")
subsp_median_df<-as.data.frame(heatmap_binding %>% group_by(Pos) %>% dplyr::summarise(median_p=median(binding_score)))
heatmap_binding<-left_join(heatmap_binding, subsp_median_df)

p_violin_p <- ggplot(heatmap_binding, aes(x = as.factor(Pos), y = binding_score)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_violin(scale = "width", size = 0.2, aes(fill = median_p)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.2) +
  scale_x_discrete(breaks = seq(1, 153), labels = peptide_seq_pos, expand = c(0, 0)) +
  theme_classic() +
  labs(x = "SOD1 WT amino acid position", y = "Median\nhet. score", fill = "Median FS") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(0, 0, 0.5, 0), 'cm'),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = cols, limits = c(min, max), breaks = c(-4, -2, 0, 2, 4)) +
  theme(legend.position = "none")
p_violin_p

#add side-chain tiles
p_side_chain_b <- ggplot(predictions, aes(x = factor(Pos), y = row, fill = side_chain_dimer)) +
  geom_tile(color = "white", height = 0.1, size = 0.2) +
  scale_fill_manual(
    values = c(
      "core" = "#26456E",
      "surface" = "#C4D8F3",
      "dimer interface" = "#41B7C4",
      "Zn binding" = "#BB173A",
      "Cu binding" = "#DECC61",
      "Zn and Cu binding" = "#B254A5")) +
  scale_y_discrete(expand = c(0,0)) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0,0,0.5,0), 'cm')
  ) +
  labs(fill = "Side chain type")
p_side_chain_b


#combine heatmap with side-chain tiles 
p_heatmap_side_chain <- ggarrange(p_side_chain_b + theme(axis.title = element_blank(), axis.text = element_blank()),
                                  p_heatmap_binding + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                                  p_violin_p,
                                  heights = c(0.05, 1, 0.2), nrow = 3, align = "v")
p_heatmap_side_chain

ggsave(p_heatmap_side_chain, path = path_fig2, file="SOD1_bPCA_full_heatmap_violin_subs.tiff",width=30, height=10, dpi = 600)

## Single insertions ####
peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "stop")

heatmap_df_bpca_ins<-SOD_final_ins

heatmap_df_bpca_ins <- heatmap_df_bpca_ins %>%
  complete(Pos = 1:153, Mut = vectorAA, fill = list(binding_score = NA, WT_AA = NA, ID = NA, label = ""))

min<-min(heatmap_df_bpca_ins$binding_score, na.rm = TRUE)
max<-max(heatmap_df_bpca_ins$binding_score, na.rm = TRUE)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

p_heatmap_ins_bpca<-ggplot(heatmap_df_bpca_ins)+
  geom_tile(aes(Pos,factor(Mut, levels=rev(vectorAA)), fill=binding_score), color='white', size=0.2)+
  scale_x_continuous(breaks=positions, labels = positions, expand=c(0,0))+
  theme_minimal()+
  theme()+
  labs(x="Position of inserted position", y="Inserted amino acid", fill="Heterodim score")+
  theme(legend.title = element_text(size=18),
        legend.title.align = 0.5,
        legend.text = element_text(size=12), 
        axis.title = element_text(size = 22, face="bold"),
        axis.text.y = element_text(size=14),
        plot.margin = unit(c(0,0,0.5,0), 'cm'),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours=cols, limits=c(min,max), na.value = "grey",breaks=c(-1, 0, 1), labels = c(-1, 0, 1)) +
  guides(fill = guide_colorbar(
    barwidth = 5,   
    barheight = 1,
    direction = "horizontal",
    title.position = "top",
    ticks = TRUE, ticks.colour = "black", frame.colour = "white"
  )) 

p_heatmap_ins_bpca

#violinplot
heatmap_df_bpca_ins <- heatmap_df_bpca_ins %>%
  filter(Mut != "stop")

inssp_median_df_bpca <- heatmap_df_bpca_ins %>%
  group_by(Pos) %>%
  summarise(median_p = median(binding_score, na.rm = TRUE))

all_positions <- 1:max(heatmap_df_bpca_ins$Pos)
heatmap_df_bpca_ins <- heatmap_df_bpca_ins %>%
  right_join(data.frame(Pos = all_positions), by = "Pos") %>%
  left_join(inssp_median_df_bpca, by = "Pos")  


p_violin_ins_bpca <- ggplot(heatmap_df_bpca_ins, aes(x = as.factor(Pos), y = binding_score)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_violin(scale = "width", size = 0.2, aes(fill = median_p)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.2) +
  scale_x_discrete(breaks = positions, labels = positions, expand = c(0, 0)) +
  theme_classic() +
  labs(x = "Position of inserted amino acid", y = "Median\nheterodim. score", fill = "Median FS") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 14),
        plot.margin = unit(c(0, 0, 0.5, 0), 'cm'),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = cols, limits = c(min, max), breaks = c(-4, -2, 0, 2, 4)) +
  theme(legend.position = "none")

p_violin_ins_bpca


#combine insertions heatmap, violin and side-chain tiles
p_heatmap_ins_side_chain <- ggarrange(p_side_chain_b + theme(axis.title = element_blank(), axis.text = element_blank()),
                                      p_heatmap_ins_bpca + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                                      p_violin_ins_bpca,
                                      heights = c(0.05, 1, 0.2), nrow = 3, align = "v")

p_heatmap_ins_side_chain
ggsave(p_heatmap_ins_side_chain,, path = path_fig2, file="SOD1_bPCA_full_heatmap_violin_ins.tiff",width=30, height=10, dpi = 600)
    

## Single deletions ####
#prepare data for plot
peptide_seq<-c(strsplit(peptide_seq, '')[[1]])

peptide_seq_pos<-c()
for (n in seq_along(peptide_seq)){
  peptide_seq_pos<-c(peptide_seq_pos, paste0(peptide_seq[n], 0+n))
}

vectorAA <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P")

heatmap_df_del<-(SOD_final_del)

min<-min(heatmap_df_del$binding_score)
max<-max(heatmap_df_del$binding_score)
cols <- c(colorRampPalette(c( "darkorange", "grey95"))((-min/(-min+max)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkgreen"), bias=1)((max/(-min+max)*100)-0.5))

heatmap_df <- heatmap_df_del %>%
  distinct(ID, .keep_all = TRUE)

p_single_deletion_bpca <- ggplot(heatmap_df_del, aes(x = factor(del_pos, levels = 1:153), y = binding_score)) +
  geom_hline(yintercept = 0, size = 0.1) +
  geom_errorbar(aes(ymin = binding_score - 1.96 * binding_sigma, ymax = binding_score + 1.96 * binding_sigma),
                width = 0, size = 0.1) +
  geom_point(data = heatmap_df_del,
             aes(fill = binding_score), size = 6, shape = 21, stroke = 1.2) +
  labs(y = "Heterodimerization score", x = "Position of deleted amino acid", fill = "") +
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

p_single_deletion_bpca

#combine deletions plot with side-chain tiles
p_del_bpca <- ggarrange(p_side_chain_b + theme(axis.title = element_blank(), axis.text = element_blank()),
                        p_single_deletion_bpca,
                        heights = c(0.1, 1), nrow = 2, align = "v")

p_del_bpca
ggsave(p_del_bpca , path = path_fig2, file="SOD1_bPCA_full_map_deletions.tiff",width=30, height=5, dpi = 600)