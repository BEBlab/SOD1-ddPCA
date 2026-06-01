#This script uses the fitness scores in the overlap regions between libraries to transform fitness values and merge the 3 libraries

#load required packages
library(tidyverse)
library(ggpubr)
library(ggplot2)


#I will merge libraries per separate. First, I will merge library 1 and 2, normalizing to lib2, and then library 2 and 3, normalizing to lib2. I use lib2 as normalization as it is the one that shares overlap with lib1 and lib3
#this data is already centered to synonymous


# LIBRARY 1 - LIBRARY 2####
## Load L1 data#### 
# First get all the variants that are trustable (low_sigma or sig_10)
load("Fscore_df_SOD1_all_variants_identity.RData")
all_variants_stops_trustable_SOD1<-all_variants_stop[all_variants_stop$low_sigma == T | all_variants_stop$sig_10 == T,]
ID_trustables_SOD1<-all_variants_stops_trustable_SOD1$ID

all_SOD1<-all_variants_stops_trustable_SOD1
all_SOD1_stop <- all_variants_stop
silent_SOD1<-silent
silent_SOD1$Region<-'1'
silent_SOD1$Mut<-'silent'
all_SOD1_stop$Region <- '1'

all_SOD1 <- as.data.frame(all_SOD1)

all_SOD1 <- all_SOD1 %>% mutate(mean_out=rowMeans(select(all_SOD1, ends_with('_s1'))))
all_SOD1<-select(all_SOD1, aa_seq, starts_with('bscore'), mean_out, ID, mean_count, Pos, Mut, WT_AA, sigma_norm_iqr, type)
all_SOD1 <- rename(all_SOD1, sigma = sigma_norm_iqr)


## Load L2 data####
# First get all the variants that are trustable (low_sigma or sig_10)
load("Fscore_df_SOD2_all_variants_identity.RData")
all_variants_stops_trustable_SOD2<-all_variants_stop[all_variants_stop$low_sigma == T| all_variants_stop$sig_10 == T, ]
ID_trustables_SOD2<-all_variants_stops_trustable_SOD2$ID
all_SOD2<-all_variants_stops_trustable_SOD2
all_SOD2_stop <- all_variants_stop
silent_SOD2<-silent
silent_SOD2$Region<-'2'
silent_SOD2$Mut<-'silent'
all_SOD2_stop$Region <- '2'

all_SOD2 <- as.data.frame(all_SOD2)

all_SOD2 <- all_SOD2 %>% mutate(mean_out=rowMeans(select(all_SOD2, ends_with('_s1'))))
all_SOD2<-select(all_SOD2, aa_seq, starts_with('bscore'), mean_out, ID, mean_count, Pos, Mut, WT_AA, sigma_norm_iqr, type)
all_SOD2 <- rename(all_SOD2, sigma = sigma_norm_iqr)

##Load L3 data####
# First get all the variants that are trustable (low_sigma or sig_10)
load("Fscore_df_SOD3_all_variants_identity.RData")
all_variants_stops_trustable_SOD3<-all_variants_stop[all_variants_stop$low_sigma == T| all_variants_stop$sig_10 == T,]
ID_trustables_SOD3<-all_variants_stops_trustable_SOD3$ID


all_SOD3<-all_variants_stops_trustable_SOD3
all_SOD3_stop <- all_variants_stop
silent_SOD3<-silent
silent_SOD3$Region<-'3'
silent_SOD3$Mut<-'silent'
all_SOD3_stop$Region <- '3'
all_SOD3 <- as.data.frame(all_SOD3)
all_SOD3 <- all_SOD3 %>% mutate(mean_out=rowMeans(select(all_SOD3, ends_with('_s1'))))
all_SOD3<-select(all_SOD3, aa_seq, starts_with('bscore'), mean_out, ID, mean_count, Pos, Mut, WT_AA, sigma_norm_iqr, type)
all_SOD3 <- rename(all_SOD3, sigma = sigma_norm_iqr)

#I will work with L1-L2 and L2-L3 per separate
SOD12<-inner_join(all_SOD1,all_SOD2, by="ID")
SOD23<-inner_join(all_SOD2,all_SOD3, by= "ID")


SOD12<-rename(SOD12, mean_count_1 = mean_count.x, mean_count_2 = mean_count.y, 
              sigma = sigma.x, sigma_2 = sigma.y, WT_AA = WT_AA.x, Mut = Mut.x, Pos = Pos.x,
              mean_out_1 = mean_out.x, mean_out_2 = mean_out.y)
SOD12<-select(SOD12, WT_AA, Mut, Pos, ID,
              mean_count_1, sigma, mean_count_2, sigma_2, bscore_c.x, bscore_c.y, mean_out_1, mean_out_2, type.x, type.y)

SOD12<-rename(SOD12, bscore_1_c = bscore_c.x, bscore_2_c = bscore_c.y)


SOD23<-rename(SOD23, mean_count_1 = mean_count.x, mean_count_2 = mean_count.y, 
              sigma = sigma.x, sigma_2 = sigma.y, WT_AA = WT_AA.x, Mut = Mut.x, Pos = Pos.x,
              mean_out_1 = mean_out.x, mean_out_2 = mean_out.y)

SOD23<-select(SOD23, WT_AA, Mut, Pos, ID,
              mean_count_1, sigma, mean_count_2, sigma_2, bscore_c.x, bscore_c.y, mean_out_1, mean_out_2, type.x, type.y)


##Synonymous and Stops distribution#### 
stops_SOD1 <- all_SOD1_stop %>%
  filter(type=="stop")
stops_SOD2 <- all_SOD2_stop %>%
  filter(type=="stop")
stops_SOD3 <- all_SOD3_stop %>%
  filter(type=="stop")

stops_SOD1$Mut <- "stop"
stops_SOD2$Mut <- "stop"
stops_SOD3$Mut <- "stop"


silentSOD123<-rbind(select(silent_SOD1, bscore_c, Mut, Region, Nmut_codons, sigma), select(silent_SOD2, bscore_c, Mut, Region, Nmut_codons, sigma), select(silent_SOD3, bscore_c, Mut, Region, Nmut_codons, sigma))
stops_SOD123<-rbind(select(stops_SOD1, bscore_c, Mut, Region), select(stops_SOD2, bscore_c, Mut, Region),select(stops_SOD3, bscore_c, Mut, Region))




dist<-rbind(silentSOD123[,c('bscore_c', 'Mut', 'Region')], stops_SOD123[,c('bscore_c', 'Mut', 'Region')])

p_hist<-ggplot(dist, aes(x=bscore_c))+
  geom_histogram(binwidth = 0.2, position="identity", aes(fill=Region), alpha=0.5)+
  #geom_histogram(binwidth = 0.2, position="identity", alpha=0.8, color=NA, aes(fill=Region))+
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=0.5)+
  facet_grid(Mut ~ .)+
  theme_bw()+
  scale_color_manual(values=c("blue", "red"))+
  #scale_color_manual(values=c("green", "black"))+
  labs(x="Heterodimerization score", y="Counts")+
  theme(legend.position = c(0.1,0.8),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5,"cm"),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        axis.line = element_line(color='black'),
        axis.title  = element_text(size = 16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=14),
        plot.margin = margin(t=1, r=1, b=1, l=1, unit = 'cm'))+
  xlim(-4,1)
p_hist


ggsave(p_hist, file="SOD123_overlapp_hist_identity_bpca.jpg", width=8, height=4)

### Linear Regression of SOD1-SOD2 ####
# Check correlation between the centered nscores calculated in each library
corr<-cor.test(SOD12$bscore_1_c, SOD12$bscore_2_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

# Linear Regression 
# x = SOD2 and y = SOD1 Then i can predict the values of SOD1 adjusted to SOD2
LR_SOD12<-lm(bscore_1_c~bscore_2_c, data = SOD12)
summary(LR_SOD12)

#Check the linear model
# Residuals vs Fitted Plot 
# not have any obvious distinct pattern. While it is slightly curved, it has equally spread residuals around the horizontal line without a distinct pattern. This is a good indication it is not a non-linear relationship.
plot(LR_SOD12, which=1, col=c("blue")) 
# Q-Q Plot
# Residuals should be normally distributed and the Q-Q Plot will show this. If residuals follow close to a straight line on this plot, it is a good indication they are normally distributed.
plot(LR_SOD12, which=2, col=c("red"))  
# Scale-Location Plot
# that the residuals have equal variance along the regression line. It is also called the Spread-Location plot.
plot(LR_SOD12, which=3, col=c("blue"))  
# Residuals vs Leverage
# An influential case will appear in the top right or bottom left of the chart inside a red line which marks Cook’s Distance. 
plot(LR_SOD12, which=5, col=c("blue")) 

ID<-SOD12$ID
y<-SOD12$bscore_1_c
x<-SOD12$bscore_2_c
SD2<-2*sd(resid(LR_SOD12))

LR_SOD12_df<-data.frame(x, y)
LR_SOD12_df$residuals_abs<-abs(LR_SOD12$residuals)
LR_SOD12_df$outliers<-FALSE
LR_SOD12_df[LR_SOD12_df$residuals_abs>SD2,]$outliers<-TRUE
LR_SOD12_df$ID<-ID
LR_SOD12_df[LR_SOD12_df$outliers == FALSE,]$ID<-''

p_corr<-ggplot(LR_SOD12_df, aes(x=x, y=y))+
  geom_point(aes(color=outliers))+
  theme_bw()+
  labs(x="Overlap in library 2", y="Overlap in library 1")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 16),
          axis.text = element_text(size=14))+
  annotate("text", x = -0.5, y = 0.5, label = paste0("R=", round(R, 2)), size=4)+
  annotate("text", x = -0.5, y = 0.4, label = paste0("p=",format(p, digits = 2, scientific = T)), size=4)+
  scale_color_manual(values=c('grey', 'red'))+
  geom_text(aes(label=ID), hjust = 0, nudge_x = 0.1, size=2)+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")
p_corr

ggsave(p_corr, file="SOD1_SOD2_overlapp_corr_bpca.jpg", width=8, height=4)

#### SOD1-SOD2 data transformation ####
##With this LR I can transform the values of SOD1
all_SOD1 <- all_SOD1 %>%
  rename(bscore_2_c = bscore_c)  

all_SOD1$bscore_2_c <- predict(LR_SOD12, newdata = all_SOD1[, c("aa_seq","bscore_2_c")])

#now SOD1 has its nscore_c.y normalized to SOD2. I will modify SOD2 to then join them
all_SOD2$bscore_1_c<-NA
all_SOD2$sigma_2<-all_SOD2$sigma
all_SOD2$sigma<-NA
all_SOD2$mean_out_2<-all_SOD2$mean_out
all_SOD2$mean_out_1<-NA

#now I modify SOD1 to keep only the values that come from SOD1
all_SOD1$bscore_1_c<-NA
all_SOD1$sigma_2<-NA
all_SOD1$sigma<-all_SOD1$sigma
all_SOD1$mean_out_1<-all_SOD1$mean_out
all_SOD1$mean_out_2<-NA

#now I join both SOD1 and SOD2, with the normalized data
SOD12<-inner_join(all_SOD1, all_SOD2, by="ID")

SOD12<-rename(SOD12, mean_count_1 = mean_count.x, mean_count_2 = mean_count.y, 
                     sigma = sigma.x, sigma_2 = sigma.y, WT_AA = WT_AA.x, Mut = Mut.x, Pos = Pos.x,
                     mean_out_1 = mean_out.x, mean_out_2 = mean_out.y, 
                     bscore_1_c = bscore_2_c, bscore_2_c = bscore_c, aa_seq = aa_seq.x)

SOD12<-select(SOD12, WT_AA, Mut, Pos, ID, aa_seq, 
                     mean_count_1, sigma, mean_count_2, sigma_2, bscore_1_c, bscore_2_c, mean_out_1, mean_out_2, type.x)


SOD12<-SOD12[, c('aa_seq', 'WT_AA', 'Mut', 'Pos', 'ID', 'sigma', 'sigma_2',
                               'bscore_1_c', 'bscore_2_c', 'mean_out_1', 'mean_out_2', 'type.x')]

all_SOD2<-rename(all_SOD2, bscore_2_c=bscore_c)
all_SOD1<-all_SOD1[, c('aa_seq', 'WT_AA', 'Mut', 'Pos', 'ID', 'sigma', 'sigma_2',
                                   'bscore_1_c', 'bscore_2_c', 'mean_out_1', 'mean_out_2', 'type')]
all_SOD2<-all_SOD2[, c('aa_seq', 'WT_AA', 'Mut', 'Pos', 'ID', 'sigma', 'sigma_2',
                                   'bscore_1_c', 'bscore_2_c', 'mean_out_1', 'mean_out_2', 'type')]


all_SOD2<-all_SOD2[!all_SOD2$ID %in% all_SOD1$ID,]

SOD12_final<-rbind(all_SOD1, all_SOD2)
#this df is ready to be merged with SOD23_final, that I prepare from here



#merging SOD1 libraries 2 and 3####
# SOD1 LIBRARY 2
# First get all the variants that are trustable (low_sigma or sig_10)
load("Fscore_df_SOD2_all_variants_identity.RData")
all_variants_stops_trustable_SOD2<-all_variants_stop[all_variants_stop$low_sigma == T | all_variants_stop$sig_10 == T,]
ID_trustables_SOD2<-all_variants_stops_trustable_SOD2$ID
all_SOD2<-all_variants_stops_trustable_SOD2
all_SOD2_stop <- all_variants_stop


silent_SOD2<-silent
silent_SOD2$Region<-'2'
silent_SOD2$Mut<-'silent'
all_SOD2_stop$Region <- '2'

all_SOD2 <- as.data.frame(all_SOD2)

all_SOD2 <- all_SOD2 %>% mutate(mean_out=rowMeans(select(all_SOD2, ends_with('_s1'))))
all_SOD2<-select(all_SOD2, aa_seq, starts_with('bscore'), mean_out, ID, mean_count, Pos, Mut, WT_AA, sigma_norm_iqr, type)
all_SOD2 <- rename(all_SOD2, sigma = sigma_norm_iqr)

#SOD1 LIBRARY 3
# First get all the variants that are trustable (low_sigma or sig_10)
load("Fscore_df_SOD3_all_variants_identity.RData")
all_variants_stops_trustable_SOD3<-all_variants_stop[all_variants_stop$low_sigma == T | all_variants_stop$sig_10 == T,]
ID_trustables_SOD3<-all_variants_stops_trustable_SOD3$ID


all_SOD3<-all_variants_stops_trustable_SOD3
all_SOD3_stop <- all_variants_stop
silent_SOD3<-silent
silent_SOD3$Region<-'3'
silent_SOD3$Mut<-'silent'
all_SOD3_stop$Region <- '3'

all_SOD3 <- as.data.frame(all_SOD3)

all_SOD3 <- all_SOD3 %>% mutate(mean_out=rowMeans(select(all_SOD3, ends_with('_s1'))))
all_SOD3<-select(all_SOD3, aa_seq, starts_with('bscore'), mean_out, ID, mean_count, Pos, Mut, WT_AA, sigma_norm_iqr, type)
all_SOD3 <- rename(all_SOD3, sigma = sigma_norm_iqr)


#I will work with SOD1-SOD2 and SOD2-SOD3 per separate
SOD23<-inner_join(all_SOD2,all_SOD3, by= "ID")

SOD23<-rename(SOD23, mean_count_1 = mean_count.x, mean_count_2 = mean_count.y, 
              sigma = sigma.x, sigma_2 = sigma.y, WT_AA = WT_AA.x, Mut = Mut.x, Pos = Pos.x,
              mean_out_1 = mean_out.x, mean_out_2 = mean_out.y)

SOD23<-select(SOD23, WT_AA, Mut, Pos, ID,
              mean_count_1, sigma, mean_count_2, sigma_2, bscore_c.x, bscore_c.y, mean_out_1, mean_out_2, type.x, type.y)
SOD23<-rename(SOD23, bscore_2_c = bscore_c.x, bscore_3_c = bscore_c.y)



#### Linear Regression of SOD2-SOD3 ####
corr<-cor.test(SOD23$bscore_2_c, SOD23$bscore_3_c, use="complete.obs")
R<-corr$estimate
p<-corr$p.value

# Linear Regression 
# x = SOD3 and y = SOD2 Then i can predict the values of SOD1 adjusted to SOD2
LR_SOD23<-lm(bscore_3_c~bscore_2_c, data = SOD23)
summary(LR_SOD23)

#Check the linear model
# Residuals vs Fitted Plot 
# not have any obvious distinct pattern. While it is slightly curved, it has equally spread residuals around the horizontal line without a distinct pattern. This is a good indication it is not a non-linear relationship.
plot(LR_SOD23, which=1, col=c("blue")) 
# Q-Q Plot
# Residuals should be normally distributed and the Q-Q Plot will show this. If residuals follow close to a straight line on this plot, it is a good indication they are normally distributed.
plot(LR_SOD23, which=2, col=c("red"))  
# Scale-Location Plot
# that the residuals have equal variance along the regression line. It is also called the Spread-Location plot.
plot(LR_SOD23, which=3, col=c("blue"))  
# Residuals vs Leverage
# An influential case will appear in the top right or bottom left of the chart inside a red line which marks Cook’s Distance. 
plot(LR_SOD23, which=5, col=c("blue")) 
#observations 54 and 70 since to be highlly influential, I will remove them

ID<-SOD23$ID
y<-SOD23$bscore_3_c
x<-SOD23$bscore_2_c
SD2<-2*sd(resid(LR_SOD23))

LR_SOD23_df<-data.frame(x, y)
LR_SOD23_df$residuals_abs<-abs(LR_SOD23$residuals)
LR_SOD23_df$outliers<-FALSE
LR_SOD23_df[LR_SOD23_df$residuals_abs>SD2,]$outliers<-TRUE
LR_SOD23_df$ID<-ID
LR_SOD23_df[LR_SOD23_df$outliers == FALSE,]$ID<-''

p_corr<-ggplot(LR_SOD23_df, aes(x=x, y=y))+
  geom_point(aes(color=outliers))+
  theme_bw()+
  labs(x="Overlap in library 2", y="Overlap in library 3")+
  theme(  panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color='black'),
          axis.title  = element_text(size = 16),
          axis.text = element_text(size=14))+
  annotate("text", x = -0.5, y = 0.5, label = paste0("R=", round(R, 2)), size=4)+
  annotate("text", x = -0.5, y = 0.4, label = paste0("p=",format(p, digits = 2, scientific = T)), size=4)+
  scale_color_manual(values=c('grey', 'red'))+
  geom_text(aes(label=ID), hjust = 0, nudge_x = 0.1, size=2)+
  geom_smooth(method = "lm", linetype = 2, size=1, se = F, color = "black")
p_corr

ggsave(p_corr, file="SOD2_SOD3_overlapp_corr_identity_bpca.jpg", width=8, height=4)



#### SOD2-SOD3 data transformation ####

###################
## With this LR I can transform the values of SOD3
all_SOD3 <- all_SOD3 %>%
  rename(bscore_2_c = bscore_c)  


all_SOD3$bscore_2_c <- predict(LR_SOD23, newdata = all_SOD3[, c("aa_seq","bscore_2_c")])

#now SOD1 has its nscore_c.y normalized to SOD2. I will modify SOD2 to then join them
all_SOD2$bscore_3_c<-NA
all_SOD2$sigma_2<-all_SOD2$sigma
all_SOD2$sigma<-NA
all_SOD2$mean_out_2<-all_SOD2$mean_out
all_SOD2$mean_out_1<-NA

#now I modify SOD1 to keep only the values that come from SOD1
all_SOD3$bscore_3_c<-NA
all_SOD3$sigma_2<-NA
all_SOD3$sigma<-all_SOD3$sigma
all_SOD3$mean_out_1<-all_SOD3$mean_out
all_SOD3$mean_out_2<-NA


#now I join both SOD1 and SOD2, with the normalized data
SOD23<-inner_join(all_SOD2, all_SOD3, by="ID")

SOD23<-rename(SOD23, mean_count_1 = mean_count.x, mean_count_2 = mean_count.y, 
              sigma = sigma.x, sigma_2 = sigma.y, WT_AA = WT_AA.x, Mut = Mut.x, Pos = Pos.x,
              mean_out_1 = mean_out.x, mean_out_2 = mean_out.y, 
              bscore_2_c = bscore_c, bscore_3_c = bscore_2_c, aa_seq = aa_seq.x)

SOD23<-select(SOD23, WT_AA, Mut, Pos, ID, aa_seq, 
              mean_count_1, sigma, mean_count_2, sigma_2, bscore_2_c, bscore_3_c, mean_out_1, mean_out_2, type.x)


SOD23<-SOD23[, c('aa_seq', 'WT_AA', 'Mut', 'Pos', 'ID', 'sigma', 'sigma_2',
                 'bscore_2_c', 'bscore_3_c', 'mean_out_1', 'mean_out_2', 'type.x')]

all_SOD2<-rename(all_SOD2, bscore_2_c=bscore_c)

all_SOD2<-all_SOD2[, c('aa_seq', 'WT_AA', 'Mut', 'Pos', 'ID', 'sigma', 'sigma_2',
                       'bscore_2_c', 'bscore_3_c', 'mean_out_1', 'mean_out_2', 'type')]
all_SOD3<-all_SOD3[, c('aa_seq', 'WT_AA', 'Mut', 'Pos', 'ID', 'sigma', 'sigma_2',
                       'bscore_2_c', 'bscore_3_c', 'mean_out_1', 'mean_out_2', 'type')]


all_SOD3<-all_SOD3[!all_SOD3$ID %in% all_SOD2$ID,]


SOD23_final<-rbind(all_SOD2, all_SOD3)

SOD23_final <- subset(SOD23_final, Pos != 38)

SOD12_final <- SOD12_final %>%
  mutate(sigma_final=coalesce(sigma, sigma_2)) %>%
  select(-sigma, -sigma_2)

SOD23_final <- SOD23_final %>%
  mutate(sigma_final=coalesce(sigma, sigma_2)) %>%
  select(-sigma, -sigma_2)




SOD_final <- full_join(SOD12_final, SOD23_final, by = c("ID", "Pos", "aa_seq"))
SOD_final <- SOD_final %>%
  mutate(bscore = coalesce(bscore_2_c.x, bscore_2_c.y)) %>%  
  select(-bscore_2_c.x, -bscore_2_c.y)  %>%
  mutate(type= coalesce(type.x,type.y)) %>%
  select(-type.x, -type.y) %>%
  mutate(WT_AA = coalesce(WT_AA.x, WT_AA.y)) %>%  
  select(-WT_AA.x, -WT_AA.y)  %>%
  mutate(Mut = coalesce(Mut.x, Mut.y)) %>% 
  select(-Mut.x, -Mut.y) %>%
  mutate(sigma=coalesce(sigma_final.x,sigma_final.y)) %>%
  select(-sigma_final.x, sigma_final.y)


SOD_final <- SOD_final %>%
  mutate(overlapp_pos = Pos %in% c(51,52,53))

SOD_final_overlap <- SOD_final %>%
  filter(overlapp_pos == TRUE)


SOD_final_overlap <- SOD_final_overlap %>%
  group_by(ID) %>%
  summarise(
    across(-bscore, first),
    bscore = mean(bscore, na.rm = TRUE),
    .groups = "drop"
  )

# SOD_final_overlap <- SOD_final_overlap %>%
#   group_by(ID) %>%
#   summarise(bscore = mean(bscore, na.rm = TRUE))

SOD_final <- SOD_final %>%
  filter(overlapp_pos == FALSE)

SOD_final <- bind_rows(SOD_final, SOD_final_overlap)

SOD_final_df <- select(SOD_final, ID, bscore, type, sigma)


write.csv(SOD_final, file="SOD_final_dataset_bPCA.csv")

