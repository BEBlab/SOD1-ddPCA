#load required packages
library(dplyr)
library(stringr)
library(tidyr)
library(ggsignif)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(mclust)


#load aPCA and bPCA data
apca <- read.csv("SOD_final_dataset_aPCA_new.csv")
bpca <- read.csv("SOD_final_dataset_bpca.csv")


#adjust the data
apca <- rename(apca, fitness = ascore)
bpca <- rename(bpca, fitness = bscore)


apca <- select(apca, ID, fitness, sigma, type)
bpca <- select(bpca, ID, fitness, sigma, type)

apca$ID <- gsub("-", "", apca$ID)
bpca$ID <- gsub("-", "", bpca$ID)




apca <- apca %>%
  mutate(stop = str_detect(ID, "\\*"),
         WT_AA = str_extract(ID, "^[A-Za-z]"),
         Mut_AA = str_extract(ID, "[A-Za-z]$")) %>%
  mutate(WT = 
           case_when(WT_AA ==Mut_AA ~ "syn",
                     WT_AA !=Mut_AA ~ "subs"))
  
bpca <- bpca %>%
  mutate(stop = str_detect(ID, "\\*"),
         WT_AA = str_extract(ID, "^[A-Za-z]"),
         Mut_AA = str_extract(ID, "[A-Za-z]$")) %>%
  mutate(WT = 
           case_when(WT_AA ==Mut_AA ~ "syn",
                     WT_AA !=Mut_AA ~ "subs"))



#see stops distributions
apca$assay <- "aPCA"
bpca$assay <- "bPCA"

ddpca <- rbind(apca, bpca)
ddpca_stops <- ddpca %>%
  filter(stop == TRUE) %>%
  mutate(Mut = "stop")

ddpca_syn <- ddpca %>%
  filter(WT == "syn") %>%
  mutate(Mut = "syn")


ddpca <- rbind(ddpca_stops, ddpca_syn)

p_hist <- ggplot(ddpca, aes(x = fitness)) +
  geom_histogram(binwidth = 0.2, position="identity", aes(fill=assay), alpha=0.5) +
  facet_grid(Mut~.)

p_hist

p_hist_stop_sigma <- ggplot(ddpca, aes(x = sigma)) +
  geom_histogram(binwidth = 0.2, position="identity", aes(fill=assay), alpha=0.5)

p_hist_stop_sigma 



p_violin_stop <- ggplot(ddpca, aes(x = assay, y = fitness)) +
  geom_violin()

p_violin_stop



#Normalize data with median of stops####
stops_apca <- apca %>%
  filter(stop == TRUE)

stops_bpca <- bpca %>%
  filter(stop == TRUE)


median_stops_apca <- median(stops_apca$fitness, na.rm = FALSE)
median_stops_bpca <- median(stops_bpca$fitness, na.rm = FALSE)

#now I divide each fitness with its respective stop median
apca <- apca %>%
  mutate(norm_fitness = fitness / abs(median_stops_apca),
         norm_sigma = sigma / abs(median_stops_apca)) #to simplify, I can normalize the sigma with the stops fitness median (no need to use sigma median)

bpca <- bpca %>%
  mutate(norm_fitness = fitness / abs(median_stops_bpca),
         norm_sigma = sigma / abs(median_stops_bpca))

ddpca_norm <- rbind(apca, bpca)

#plot fitness distributions
p_hist_fitness <- ggplot(ddpca_norm, aes(x=fitness))+
  geom_histogram(binwidth = 0.2, position="stack", aes(fill=assay), alpha=0.5, color = "black")+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
  scale_fill_manual(values=c("#26456E", "#E9AF1D"))+
  labs(x="Fitness score", y="Counts")+
  theme(legend.position = c(0.2,0.9),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18),
        strip.text = element_text(size=14),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist_fitness

#plot normalized scores distribution
p_hist_fitness_norm <- ggplot(ddpca_norm, aes(x=norm_fitness))+
  geom_histogram(binwidth = 0.05, position="stack", alpha=0.5, aes(fill = assay), color = "black")+ 
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
  scale_fill_manual(values=c("#26456E", "#E9AF1D"))+
  labs(x="Normalized fitness score", y="Counts")+
  theme(legend.position = "none",
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18),
        strip.text = element_text(size=14),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist_fitness_norm

p_both <- ggarrange(p_hist_fitness,
                    p_hist_fitness_norm + theme(axis.title.y = element_blank()),
                    ncol = 2)

p_both

ggsave(p_both, file = "ddpca_fitness_dist.jpg", width = 15, height = 7)


#see how the distribution as violin plots
p_vio_fitness <- ggplot(ddpca_norm, aes(x=assay, y=fitness))+
  geom_violin(scale = "width", trim = TRUE, size = 0.8) +
  geom_jitter(data = subset(ddpca_norm, stop == TRUE), aes(x=assay, y=fitness), color = "red", fill = "red", shape = 21, size = 3, width = 0.1,alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  labs(x="", y="Fitness")+
  theme_classic(base_size = 14, base_line_size = 0.2) +
  theme(legend.position = "none",
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18))
p_vio_fitness

p_vio_fitness_norm <- ggplot(ddpca_norm, aes(x=assay, y=norm_fitness))+
  geom_violin(scale = "width", trim = TRUE, size = 0.8) +
  geom_jitter(data = subset(ddpca_norm, stop == TRUE), aes(x=assay, y=norm_fitness), color = "red", fill = "red", shape = 21, size = 3, width = 0.1,alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  labs(x="", y="Normalized fitness")+
  theme_classic(base_size = 14, base_line_size = 0.2) +
  theme(legend.position = "none",
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18))
p_vio_fitness_norm

p_comb_vio <- ggarrange(p_vio_fitness,
                        p_vio_fitness_norm,
                        ncol = 2, align = "v")

p_comb_vio

ggsave(p_comb_vio, file = "p_violin_fitness_dist_comparison.jpg", width = 12, height = 6)


#Sigma distribution####
p_hist_sigma <- ggplot(ddpca_norm, aes(x=sigma))+
  geom_histogram(binwidth = 0.1, position="stack", aes(fill=assay), alpha=0.5, color = "black")+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
  scale_fill_manual(values=c("#26456E", "#E9AF1D"))+
  labs(x="Sigma", y="Counts")+
  theme(legend.position = c(0.2,0.9),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18),
        strip.text = element_text(size=14),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist_sigma

p_hist_sigma_norm <- ggplot(ddpca_norm, aes(x=norm_sigma))+
  geom_histogram(binwidth = 0.03, position="stack", aes(fill=assay), alpha=0.5, color = "black")+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
  scale_fill_manual(values=c("#26456E", "#E9AF1D"))+
  labs(x="Normalized sigma", y="Counts")+
  theme(legend.position = "none",
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18),
        strip.text = element_text(size=14),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist_sigma_norm

p_both_sigma <- ggarrange(p_hist_sigma,
                    p_hist_sigma_norm + theme(axis.title.y = element_blank()),
                    ncol = 2)

p_both_sigma

ggsave(p_both_sigma, file = "ddpca_sigma_dist.jpg", width = 15, height = 7)

#same for stops and synonymous
syn_ddpca <- ddpca_norm %>%
  filter(WT == "syn")

stops_ddpca <- ddpca_norm %>%
  filter(stop == TRUE)

syn_ddpca$Mut <- "synonymous"
stops_ddpca$Mut <- "stop"

dist <- rbind(syn_ddpca, stops_ddpca)

p_hist_dist <- ggplot(dist, aes(x=fitness))+
  geom_histogram(binwidth = 0.2, position="stack", aes(fill=assay), alpha=0.5, color = "black")+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=1)+
  theme_bw()+
  facet_grid(Mut ~ .)+
  scale_fill_manual(values=c("#26456E", "#E9AF1D"))+
  labs(x="Fitness score", y="Counts")+
  theme(legend.position = c(0.2,0.9),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18),
        strip.text = element_text(size=20),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist_dist


p_hist_dist_norm <- ggplot(dist, aes(x=norm_fitness))+
  geom_histogram(binwidth = 0.2, position="stack", aes(fill=assay), alpha=0.5, color = "black")+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=-1), color="black", linetype="dashed", size=1)+
  theme_bw()+
  facet_grid(Mut ~ .)+
  scale_fill_manual(values=c("#26456E", "#E9AF1D"))+
  labs(x="Normalized fitness score", y="Counts")+
  theme(legend.position = "none",
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 24),
        axis.text = element_text(size=18),
        strip.text = element_text(size=20),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist_dist_norm


p_both_dist <- ggarrange(p_hist_dist,
                         p_hist_dist_norm + theme(axis.title.y = element_blank()),
                    ncol = 2)

p_both_dist

ggsave(p_both_dist, file = "ddpca_fitness_dist_syn_stops.jpg", width = 15, height = 7, dpi = 300)


#Define the synonymous threshold for abundance and binding, respectively
syn_apca <- syn_ddpca %>%
  filter(assay == "aPCA")

apca_norm <- ddpca_norm %>%
  filter(assay == "aPCA")

syn_threshold_apca <- quantile(syn_apca$norm_fitness, probs = 0.05, na.rm = TRUE, type = 7)

xlims <- range(apca_norm$norm_fitness, na.rm = TRUE)

p_syn_dens_apca <- ggplot(data = subset(apca_norm, type == "syn"), aes(x = norm_fitness)) +
  geom_density(position = "identity", fill = "#26456E", size = 1.2) +
  scale_x_continuous(limits = range(apca$norm_fitness, na.rm = TRUE)) +
  geom_vline(xintercept = syn_threshold_apca, color = "black", size = 1.5, linetype = "dashed") +
  labs(x = "Normalized score", y = "Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30))
p_syn_dens_apca


syn_bpca <- syn_ddpca %>%
  filter(assay == "bPCA")

bpca_norm <- ddpca_norm %>%
  filter(assay == "bPCA")

syn_threshold_bpca <- quantile(syn_bpca$norm_fitness, probs = 0.05, na.rm = TRUE, type = 7)

xlims <- range(bpca_norm$norm_fitness, na.rm = TRUE)

p_syn_dens_bpca <- ggplot(data = subset(bpca_norm, type == "syn"), aes(x = norm_fitness)) +
  geom_density(position = "identity", fill = "#E9AF1D", size = 1.2) +
  scale_x_continuous(limits = range(apca$norm_fitness, na.rm = TRUE)) +
  geom_vline(xintercept = syn_threshold_bpca, color = "black", size = 1.5, linetype = "dashed") +
  labs(x = "Normalized score", y = "Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30))
p_syn_dens_bpca

p_syn_ddpca_thres <- ggarrange(p_syn_dens_apca + theme(axis.title.x = element_blank()),
                               p_syn_dens_bpca,
                               nrow = 2, align = "v")
p_syn_ddpca_thres

ggsave(p_syn_ddpca_thres, file = "syn_ddpca_thres.tiff", width = 10, height = 8, dpi = 300)

#Define the stops threshold. This will serve to define extreme low abundance/binding variants (stop-like variants)
stop_apca <- stops_ddpca %>%
  filter(assay == "aPCA")


stop_threshold_apca <- quantile(stop_apca$norm_fitness, probs = 0.95, na.rm = TRUE, type = 7)

xlims <- range(apca_norm$norm_fitness, na.rm = TRUE)

p_stop_dens_apca <- ggplot(data = subset(apca_norm, stop == TRUE), aes(x = norm_fitness)) +
  geom_density(position = "identity", fill = "#26456E", size = 1.2) +
  scale_x_continuous(limits = range(apca$norm_fitness, na.rm = TRUE)) +
  geom_vline(xintercept = stop_threshold_apca, color = "black", size = 1.5, linetype = "dashed") +
  labs(x = "Normalized score", y = "Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30))
p_stop_dens_apca


stop_bpca <- stops_ddpca %>%
  filter(assay == "bPCA")

stop_threshold_bpca <- quantile(stop_bpca$norm_fitness, probs = 0.95, na.rm = TRUE, type = 7)

xlims <- range(bpca_norm$norm_fitness, na.rm = TRUE)

p_stop_dens_bpca <- ggplot(data = subset(bpca_norm, stop == TRUE), aes(x = norm_fitness)) +
  geom_density(position = "identity", fill = "#E9AF1D", size = 1.2) +
  scale_x_continuous(limits = range(bpca$norm_fitness, na.rm = TRUE)) +
  geom_vline(xintercept = stop_threshold_bpca, color = "black", size = 1.5, linetype = "dashed") +
  labs(x = "Normalized score", y = "Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 33),
        axis.text = element_text(size = 30))
p_stop_dens_bpca

p_stop_ddpca_thres <- ggarrange(p_stop_dens_apca + theme(axis.title.x = element_blank()),
                               p_stop_dens_bpca,
                               nrow = 2, align = "v")
p_stop_ddpca_thres

ggsave(p_stop_ddpca_thres, file = "stop_ddpca_thres.tiff", width = 10, height = 8, dpi = 300)


#Split into substitutions, insertions and deletions
subs_apca <- apca %>%
  filter(type %in% c("subs", "WT", "syn", "stop")) %>%
  mutate(ins_pos = "",
         ins_aa = "",
         del_pos = "",
         del_aa = "")

subs_bpca <-bpca %>%
  filter(type %in% c("subs", "WT", "syn", "stop")) %>%
  mutate(ins_pos = "",
         ins_aa = "",
         del_pos = "",
         del_aa = "")

insertions_apca <- apca %>%
  filter(type == "ins") %>%
  separate_rows(ID, sep = ";") %>%
  mutate(
    ins_pos = str_extract(ID, "\\d+"),
    ins_aa = str_extract(ID, "[A-Z]$"),
    del_pos = "",
    del_aa = ""
  )


insertions_bpca <- bpca %>%
  filter(type == "ins") %>%
  separate_rows(ID, sep = ";") %>%
  mutate(
    ins_pos = str_extract(ID, "\\d+"),
    ins_aa = str_extract(ID, "[A-Z]$"),
    del_pos = "",
    del_aa = ""
  )


deletions_apca <- apca %>%
  filter(type == "del") %>%
  separate_rows(ID, sep = ";") %>%
  mutate(
    del_pos = str_extract(ID, "\\d+"),
    del_aa = str_extract(ID, "^[A-Z]"),
    ins_pos = "",
    ins_aa = "",
    WT_AA = "",
    Mut_AA = ""
  )

deletions_bpca <- bpca %>%
  filter(type== "del") %>%
  separate_rows(ID, sep = ";") %>%
  mutate(
    del_pos = str_extract(ID, "\\d+"),
    del_aa = str_extract(ID, "^[A-Z]"),
    ins_pos = "",
    ins_aa = "",
    WT_AA = "",
    Mut_AA = ""
  )


apca <- rbind(subs_apca, insertions_apca, deletions_apca)
apca <- apca %>% 
  filter(type != "WT") %>% 
  mutate(type = factor(type, levels = c("subs", "ins", "del", "syn", "stop")))
bpca <- rbind(subs_bpca, insertions_bpca, deletions_bpca)
bpca <- bpca %>% 
  filter(type != "WT") %>% 
  mutate(type = factor(type, levels = c("subs", "ins", "del", "syn", "stop")))


#plot normalized abundance scores distribution for all mutations
p_vio_apca <- ggplot(apca, aes(x = type, y = norm_fitness)) +
  geom_violin(aes(color = type), scale = "width", trim = TRUE, size = 1, alpha = 0.5) +
  geom_jitter(aes(color = type, alpha = type),width = 0.1, size = 3) +
  scale_alpha_manual(values = c("subs" = 0.1, "ins" = 0.1, "del" = 1, "syn" = 1, "stop" = 1)) +
  scale_color_manual(values = c("subs" = "#0B2130", "ins" = "#1E6085", "del" = "#A1C8DC", "syn" = "grey", "stop" = "#49525E")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.7) +
  geom_hline(yintercept = - 0.96, size = 0.7) +
  scale_x_discrete(labels = c("subs" = "Substitutions\n(n=2889)", "ins" = "Insertions\n(n=2758)", "del" = "Deletions\n(n=153)", "syn" = "Synonymous\n(n=87)", "stop" = "Premature\nstops (n=42)")) +
  labs(x = "", y = "Abundance score", color = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.position = "none")
p_vio_apca
ggsave(p_vio_apca, file = "vio_apca_mutations.tiff", width = 12, height = 7)

#plot normalized binding scores distribution for all mutations
p_vio_bpca <- ggplot(bpca, aes(x = type, y = norm_fitness)) +
  geom_violin(aes(color = type), scale = "width", trim = TRUE, size = 1, alpha = 0.5) +
  geom_jitter(aes(color = type, alpha = type),width = 0.1, size = 3) +
  scale_alpha_manual(values = c("subs" = 0.1, "ins" = 0.1, "del" = 1, "syn" = 1, "stop" = 1)) +  scale_color_manual(values = c("subs" = "#0B2130", "ins" = "#1E6085", "del" = "#A1C8DC", "syn" = "grey", "stop" = "#49525E")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.7) +
  geom_hline(yintercept = - 0.96, size = 0.7) +
  scale_x_discrete(labels = c("subs" = "Substitutions\n(n=2898)", "ins" = "Insertions\n(n=2758)", "del" = "Deletions\n(n=153)", "syn" = "Synonymous\n(n=88)", "stop" = "Premature\nstops (n=42)")) +
  labs(x = "", y = "Heterodimerization score", color = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 25),
        legend.position = "none")
p_vio_bpca
ggsave(p_vio_bpca, file = "vio_bpca_mutations.tiff", width = 12, height = 7)



#divide into substitutions, insertions and deletions
apca <- apca %>%
  mutate(type = recode(type, "syn" = "subs")) %>%
  filter(type %in% c("subs", "stop", "ins", "del"))

bpca <- bpca %>%
  mutate(type = recode(type, "syn" = "subs")) %>%
  filter(type %in% c("subs", "stop", "ins", "del"))

#Data categorization with FDR####
##Abundance####
# FDR=0.1 correction and assignment into categories
#First, I test FDR against synonymous
median_syn_apca <- median(syn_apca$norm_fitness)
apca$category <- NA

apca$norm_sigma<-replace_na(apca$norm_sigma, 0)
apca$zscore<-apca$norm_fitness/apca$norm_sigma
apca$p.adjust<-p.adjust(2*pnorm(-abs(apca$zscore)), method = "BH")

apca$sig_10_norm<-FALSE
apca[apca$p.adjust<0.1,]$sig_10_norm<-TRUE


apca$category[apca$norm_fitness < syn_threshold_apca &
                apca$p.adjust < 0.1] <- "low abundance"

apca$category[is.na(apca$category)] <- "WT-like abundance"

apca$category[apca$norm_fitness > median_syn_apca &
                apca$p.adjust < 0.1] <- "high abundance"


#Second, I test FDR against stops
apca$zscore_stop <- (apca$norm_fitness - stop_threshold_apca) / apca$norm_sigma
apca$p.adjust_stop <- p.adjust(2*pnorm(-abs(apca$zscore_stop)), method = "BH")

apca$sig_10_stop <- FALSE
apca[apca$p.adjust_stop < 0.1, ]$sig_10_stop <- TRUE

idx_fsdec <- apca$category == "low abundance"

apca$category[idx_fsdec & apca$norm_fitness > stop_threshold_apca &
                apca$p.adjust < 0.1] <- "low abundance"

apca$category[idx_fsdec & apca$norm_fitness < stop_threshold_apca &
                apca$p.adjust < 0.1] <- "stop-like abundance"


apca$category <- factor(apca$category, levels = c("high abundance","WT-like abundance","low abundance", "stop-like abundance"))

p_cat <- ggplot(data = apca, aes(x = norm_fitness)) +
  geom_histogram(bins = 30,aes(fill = category), color = "black", size = 1) +
  scale_fill_manual(values = c("high abundance" = "darkgreen", "WT-like abundance" = "lightgrey",
                    "low abundance" = "darkorange", "stop-like abundance" = "#DD4636")) +
  facet_wrap(~category, nrow = 6) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = syn_threshold_apca, linetype = "dashed") +
  geom_vline(xintercept = stop_threshold_apca) +
  labs(x = "Normalized abundance score", y = "Number of mutations", fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 22),
        strip.text = element_text(size = 20),
        legend.position = "none")
p_cat


##Binding####
#First, I test FDR against synonymous
median_syn_bpca <- median(syn_bpca$norm_fitness)

bpca$category <- NA

bpca$norm_sigma<-replace_na(bpca$norm_sigma, 0)
bpca$zscore<-bpca$norm_fitness/bpca$norm_sigma
bpca$p.adjust<-p.adjust(2*pnorm(-abs(bpca$zscore)), method = "BH")

bpca$sig_10_norm<-FALSE
bpca[bpca$p.adjust<0.1,]$sig_10_norm<-TRUE


bpca$category[bpca$norm_fitness < syn_threshold_bpca &
                bpca$p.adjust < 0.1] <- "low heterodimerization"

bpca$category[is.na(bpca$category)] <- "WT-like heterodimerization"

bpca$category[bpca$norm_fitness > median_syn_bpca &
                bpca$p.adjust < 0.1] <- "high heterodimerization"


#Second, I test FDR against stops
bpca$zscore_stop <- (bpca$norm_fitness - stop_threshold_bpca) / bpca$norm_sigma
bpca$p.adjust_stop <- p.adjust(2*pnorm(-abs(bpca$zscore_stop)), method = "BH")

bpca$sig_10_stop <- FALSE
bpca[bpca$p.adjust_stop < 0.1, ]$sig_10_stop <- TRUE

idx_fsdec <- bpca$category == "low heterodimerization"

bpca$category[idx_fsdec & bpca$norm_fitness > stop_threshold_bpca &
                bpca$p.adjust < 0.1] <- "low heterodimerization"

bpca$category[idx_fsdec & bpca$norm_fitness < stop_threshold_bpca &
                bpca$p.adjust < 0.1] <- "stop-like heterodimerization"



bpca$category <- factor(bpca$category, levels = c("high heterodimerization","WT-like heterodimerization", "low heterodimerization", "stop-like heterodimerization"))

p_cat_bpca <- ggplot(data = bpca, aes(x = norm_fitness)) +
  geom_histogram(bins = 30,aes(fill = category), color = "black", size = 1) +
  scale_fill_manual(values = c("high heterodimerization" = "darkgreen", "WT-like heterodimerization" = "lightgrey",
                               "low heterodimerization" = "darkorange", "stop-like heterodimerization" = "#DD4636")) +
  facet_wrap(~category, nrow = 6) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = syn_threshold_bpca, linetype = "dashed") +
  geom_vline(xintercept = stop_threshold_bpca) +
  labs(x = "Normalized heterodimerization score", y = "Number of mutations", fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 22),
        strip.text = element_text(size = 20),
        legend.position = "none")
p_cat_bpca

ggsave(p_cat, path = path_fig1, file = "categories_apca.tiff", width = 10, height = 10, dpi = 300)
ggsave(p_cat_bpca, path = path_fig1, file = "categories_bpca.tiff", width = 10, height = 10, dpi = 300)





ddPCA <- inner_join(apca, bpca, by = "ID")

ddPCA <-ddPCA %>%
  rename(abundance_score = norm_fitness.x, abundance_sigma = norm_sigma.x, binding_score = norm_fitness.y, binding_sigma = norm_sigma.y, mutation_type = type.x, category_fdr_apca = category.x, category_fdr_bpca = category.y) %>%
  select(ID, abundance_score, abundance_sigma, binding_score, binding_sigma, mutation_type, category_fdr_apca, category_fdr_bpca) %>%
  filter(mutation_type != "WT" & ID != "WT")

ddPCA <- ddPCA %>%
  mutate(category_fdr_apca = recode(category_fdr_apca,
                                    "high abundance" = "high",
                                    "WT-like abundance" = "WT-like",
                                    "low abundance" = "low",
                                    "stop-like abundance" = "stop-like")) %>% 
  mutate(category_fdr_bpca = recode(category_fdr_bpca,
                                    "high heterodimerization" = "high",
                                    "WT-like heterodimerization" = "WT-like",
                                    "low heterodimerization" = "low",
                                    "stop-like heterodimerization" = "stop-like"))

write.csv(ddPCA, file = "SOD1_abundance_binding_scores.csv")

