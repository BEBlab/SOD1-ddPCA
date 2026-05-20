#load required packages
library(tidyverse)
library(ggpubr)


#load dimsum output and library design file. This will allow me to find in the dimsum output my designed variants
load("SOD2_final_fitness_replicates.RData")

lib2_design <- read.delim("SOD1_lib2.tsv")

#by doing this, since synonymous mutations will lead to the same aa_seq, they will be replicated in the df. To avoid working with them, I will remove them from this df and add them again from the syn df
all_variants <- inner_join(all_variants, lib2_design, by= "aa_seq")

all_variants<-rename(all_variants, ascore = fitness, ascore1 = fitness1_uncorr, ascore2 = fitness2_uncorr,
                     ascore3 = fitness3_uncorr)

all_variants <- all_variants %>%
  mutate(type =case_when(
    grepl("TRUE", STOP) ~ "stop",
    TRUE ~type))

#remove "wrong" syn mutations
all_variants$ID <- gsub("-", "", all_variants$ID)

all_variants <- all_variants %>%
  mutate(WT_AA = str_extract(ID, "^[A-Za-z]"),
         Mut_AA = str_extract(ID, "[A-Za-z*]$"), 
         syn = WT_AA == Mut_AA) %>%
  filter(is.na(syn) | !syn)

all_variants <- select(all_variants, -WT_AA, -Mut_AA, -syn, -nt_seq.x)


#now I can add the correct synonymous mutations to the all variants df
#select designed synonymous mutations
lib2_design$nt_seq <- tolower(lib2_design$nt_seq)
syn_designed <- inner_join(synonymous, lib2_design, by = "nt_seq")
syn_designed <- syn_designed %>%
  mutate(type = recode(type, 
                       "subs" = "syn",
                       "WT" = "syn")) %>%
  rename(ascore = fitness, ascore1 = fitness1_uncorr, ascore2 = fitness2_uncorr, ascore3 = fitness3_uncorr)

all_variants <- full_join(all_variants, syn_designed)

all_variants<-select(all_variants, aa_seq, starts_with('ascore') | starts_with('count'), ID, type, sigma, mean_count)
synonymous<-rename(synonymous, ascore = fitness, ascore1 = fitness1_uncorr, ascore2 = fitness2_uncorr, 
                   ascore3 = fitness3_uncorr)
                   
silent<-synonymous
all_variants$Pos<-""
all_variants$Mut <-""
all_variants$WT_AA <- ""

all_variants <- all_variants %>%
  mutate(Pos=as.numeric(str_extract(ID, "\\d+"))) %>%
  mutate(WT_AA=str_extract(ID, "^[A-Za-z]")) %>%
  mutate(Mut=str_extract(ID, "[A-Za-z]$"))

all_variants$ID <- gsub("-","", all_variants$ID)

all_variants <- all_variants %>%
  mutate(
    type = case_when(
      str_sub(ID, 1, 1) == str_sub(ID, -1, -1) ~ "syn",  
      str_detect(ID, "\\*") ~ "stop",  
      TRUE ~ type))

#the synonymous of this library were desgined to be by 2 nt changes when possible, and if not, a second synoymous was added (=2 nt changes in 2 different codons)

##Synonymous scores distribution####
#check the synonymous fitness distribution for nmut codons =2
mean_syn_1codon<-weighted.mean(silent[silent$Nmut_codons==2,]$ascore, silent[silent$Nmut_codons==2,]$sigma^-2, na.rm = T)

silent$ascore_c<-as.numeric(paste(silent$ascore+(-mean_syn_1codon)))
silent$ID<-"silent"
silent$Mut<-"silent"

p_hist<-ggplot(silent, aes(x=ascore_c))+
  geom_histogram(binwidth = 0.2, position="identity", alpha=0.5)+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  theme_bw()+
  scale_color_manual(values=c("blue", "red"))+
  labs(x="Abundance score", y="Counts")+
  theme(legend.position = c(0.1,0.8),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))
p_hist

ggsave(p_hist, file="SOD1_syn_fitness_distribution_nmutcodon_2.jpg")

##Center mutation scores to synonymous mean####
all_variants$ascore_c<-as.numeric(paste(as.numeric(all_variants$ascore)+(-mean_syn_1codon)))
all_variants$ascore1_c<-as.numeric(paste(as.numeric(all_variants$ascore1)+(-mean_syn_1codon)))
all_variants$ascore2_c<-as.numeric(paste(as.numeric(all_variants$ascore2)+(-mean_syn_1codon)))
all_variants$ascore3_c<-as.numeric(paste(as.numeric(all_variants$ascore3)+(-mean_syn_1codon)))

#FDR=0.1 correction and assignment into categories
all_variants$zscore<-all_variants$ascore_c/all_variants$sigma
all_variants$p.adjust<-p.adjust(2*pnorm(-abs(all_variants$zscore)), method = "BH")

all_variants$sig_10<-FALSE
all_variants[all_variants$p.adjust<0.1,]$sig_10<-TRUE

all_variants$category_10<-"WT-like"
all_variants[all_variants$sig_10==T & all_variants$ascore_c<0,]$category_10<-"FS_dec"
all_variants[all_variants$sig_10==T & all_variants$ascore_c>0,]$category_10<-"FS_inc"

#with stops
all_variants_stop<-all_variants

#remove stops
all_variants<-all_variants[all_variants$Mut!="*",]

#mean input and ouput reads counts
all_variants_stop <- all_variants_stop %>%
  rowwise() %>%
  mutate(input_mean_count=mean(c(count_e1_s0, count_e2_s0, count_e3_s0)),
         output_mean_count=mean(c(count_e1_s1, count_e2_s1, count_e3_s1)))

##Normalize sigma to interquartile range####
#first, I define interquartile range == fitness range
summary(all_variants_stop$ascore_c)
iqr<-IQR(all_variants_stop$ascore_c)

first_to_wt<-summary(all_variants_stop$ascore_c)[[2]]

p_iqr<-ggplot(all_variants_stop, aes(x=ascore_c))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  labs(x="abundance score", title="abundance range")+
  
  geom_vline(xintercept=0)+
  annotate("text", label="WT binidng= 0", x=0.1, y=110, color="black")+
  
  geom_vline(xintercept=summary(all_variants_stop$ascore_c)[[2]], color="red")+
  annotate("text", label=paste0("1st Qu= ", round(summary(all_variants_stop$ascore_c)[[2]], 2)), 
           x=summary(all_variants_stop$ascore_c)[[2]]-0.8, y=140, color="red")+
  
  geom_vline(xintercept=summary(all_variants_stop$ascore_c)[[5]], color="red")+
  annotate("text", label=paste0("3rd Qu= ", round(summary(all_variants_stop$ascore_c)[[5]], 2)), 
           x=summary(all_variants_stop$ascore_c)[[5]]+0.8, y=140, color="red")+
  
  geom_vline(xintercept=summary(all_variants_stop$ascore_c)[[3]], color="blue")+
  annotate("text", label=paste0("median= ",round(summary(all_variants_stop$ascore_c)[[3]], 2) ), 
           x=summary(all_variants_stop$ascore_c)[[3]]+0.1, y=100, color="blue")
p_iqr

ggsave(p_iqr, file="p_iqr_SOD1_allvariants_.jpg", width = 5, height = 3)

#then, I check the sigma distribution
p_sigma_dist<-ggplot(all_variants_stop, aes(x=sigma))+
  geom_histogram(color="black", fill="grey", bins=100)+
  theme_bw()+
  scale_x_continuous(limits = c(0,3))
p_sigma_dist

ggsave(p_sigma_dist, file="p_sigma_dist_sod1_allvariants.jpg", width = 5, height = 3)

#and finally, I normalise the sigmas to the fitness range (==IQR)
all_variants_stop$sigma_norm_iqr<-""
fitness_range_iqr=abs(IQR(all_variants_stop$ascore_c))

#or if fitness range is 1rst to WT fitness
all_variants_stop$sigma_norm_first_toWT<-""
fitness_range_first_toWT=abs(summary(all_variants_stop$ascore_c)[[2]])

all_variants_stop$sigma_norm_iqr<-all_variants_stop$sigma / fitness_range_iqr
all_variants_stop$sigma_norm_first_toWT<-all_variants_stop$sigma / fitness_range_first_toWT

all_variants_stop$sigma_norm_iqr<-as.numeric(all_variants_stop$sigma_norm_iqr)
all_variants_stop$sigma_norm_first_toWT<-as.numeric(all_variants_stop$sigma_norm_first_toWT)

#distribution of sigma normalized to IQR
p1<-ggplot(all_variants_stop, aes(x=sigma_norm_iqr))+
  geom_histogram(bins=200, aes(fill=factor(category_10, levels=c("FS_dec", "FS_inc", "WT-like"))))+
  scale_fill_manual(values=c("#FA9B4F", "#5CA259", "#F2F2F2"))+
  facet_wrap(~category_10, ncol=1)+
  labs(title="sigma normalised to fitness range (1st Qu to 3r Qu)", fill="Category_FDR")+
  scale_x_continuous(limits = c(0,1))+
  geom_vline(xintercept = 0.3)+
  theme_bw()+
  theme(legend.title=element_text(size=14, face = "bold"),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        plot.margin = unit(c(0,0,0,0), 'cm'))
p1

ggsave(p1, file="p_sigma_normalised_SOD1_allvariants.jpg",  width = 8, height = 5)

## exclude those variants that are over 30% of the fitness range
all_variants_stop$low_sigma<-FALSE
all_variants_stop[all_variants_stop$sigma_norm_iqr<=0.30,]$low_sigma<-TRUE

# Categories based on sigma normalized to fitness range
all_variants_stop$category_10_sigma<-"unclassified"

all_variants_stop[all_variants_stop$category_10=="FS_inc"& all_variants_stop$low_sigma==T,]$category_10_sigma<-"FS_inc"
all_variants_stop[all_variants_stop$category_10=="FS_dec" & all_variants_stop$low_sigma==T,]$category_10_sigma<-"FS_dec"
all_variants_stop[all_variants_stop$category_10=="WT-like" & all_variants_stop$low_sigma==T,]$category_10_sigma<-"WT-like"

save(all_variants_stop, file="all_variants_SOD1_stops_newfilter_L2_aPCA.RData")
save(silent, all_variants, all_variants_stop, file="Fscore_df_SOD1_all_variants_identity_L2_aPCA.RData")



#Save Lib1 scores table####
apca_l2 <- select(all_variants_stop, aa_seq, ID, type, count_e1_s0, count_e2_s0, count_e3_s0, count_e1_s1, count_e2_s1, count_e3_s1, input_mean_count, output_mean_count, ascore1_c, ascore2_c, ascore3_c, ascore_c, sigma, sigma_norm_iqr, low_sigma)
apca_l2$library <- 2

#add aa_seq to synonymous
wt_seq <- apca_l2 %>%
  filter(type == "WT") %>%
  pull(aa_seq) %>%
  unique()

apca_l2 <- apca_l2 %>%
  mutate(
    aa_seq = ifelse(type == "syn", wt_seq[1], aa_seq)
  )


write.csv(apca_l2, file = "aPCA_L2_table.csv", row.names = FALSE)

