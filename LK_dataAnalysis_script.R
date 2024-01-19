setwd("/Users/rmaqsood/Documents/pre_post_covid/virome/")
library(RColorBrewer)
library(gplots)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)
library(vegan)
library(tidyverse)
library(tidyr)
library(nlme)
library(phyloseq)
library(ape)
library(corrplot)
library(taxonomizr)
#library(VirFinder)
library("qvalue")
library(glmnet)
library(Rcpp)
library(decontam)

###################Virome Infants vs Mother###################
setwd("/Users/rmaqsood/Documents/pre_post_covid/virome/")
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
#person_code significance alpha and richness
set.seed(100)
allrichness<-lme(Richness~  person_code+ pp.time +covidstatus2+ everAntibiotics, random=~1|studyID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ person_code+ pp.time +covidstatus2+ everAntibiotics, random=~1|studyID, data=m.1)
summary(allalphaDiversity)

# AllCompareRichness<-wilcox.test(Richness ~ person_code, data=m.1) 
# AllCompareRichness
# AllCompareAlpha<-wilcox.test(AlphaDiversity ~ person_code, data=m.1) 
# AllCompareAlpha

#plots
loessPlotRichness <- ggplot(m.1, aes(x = pp.time, y = Richness, color = person_code, fill = person_code)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$pp.time), max(m.1$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="Virome Richness", title = "Mothers vs infants" )
loessPlotRichness

loessPlotShannon <- ggplot(m.1, aes(x = pp.time, y = AlphaDiversity, color = person_code, fill = person_code)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$pp.time), max(m.1$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="Virome Alpha Diversity", title = "Mothers vs infants" )
loessPlotShannon

#betadiverstiy
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
d.2<-d.1[,colnames(d.1) %in% m.2$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetavirome_MotherInfant_042423.csv")

#PCoA
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC3 <- as.matrix(res$vectors[,3])
PC1
PC2
data<-read.delim("PC_virome_AllMotherInfant_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
data$PC3 <- PC3
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-.7, .7) + ylim(-.7, .7) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#personCode
p <- ggplot(data, aes(x=PC1, y=PC3, color=as.factor(person_code))) + 
  geom_point() + xlim(-.7, .7) + ylim(-.7, .7) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Infant","Mother"),
                     values = c("green", "salmon")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = m.2$Sample, # sample ID
                   id = m.2$ptnum, # subject ID 
                   tp = m.2$pp.time,
                   abx = m.2$everAntibiotics,
                   hivstatus=m.2$hivStatus,
                   covidstatus2 = m.2$covidstatus2,
                   neverPostive = m.2$NeverPositive, 
                   person_code=m.2$person_code,
                   studyID=m.2$studyID) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetavirome_MotherInfant_042423.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   person_code + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
#melt format
library(reshape2)
df <- melt(as.matrix(beta_dist))
p    <- t(apply(df[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])
p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))
df   <- df[-c(rmv1,rmv2),] 
write.csv(df,"mother_infantPair_WeightedBeta_v2.csv")

#Maaslin2
setwd('/Users/rmaqsood/Documents/pre_post_covid/virome/')
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
#infant
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
d.2<-d.1[,colnames(d.1) %in% infant.virome$Sample]
write.table(infant.virome,"infantvirome_metadata.txt")
write.table(d.2,"infantvirome_Relabun.txt")
#mother
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
d.3<-d.1[,colnames(d.1) %in% mother.virome$Sample]
write.table(mother.virome,"mothervirome_metadata.txt")
write.table(d.3,"mothervirome_Relabun.txt")
#all
library(Maaslin2)
input_data <- read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("Analysis_virome_030723.txt", row.names = 1)
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_MotherInfant_time_011824', transform = "NONE",
  fixed_effects = c('person_code'),
  random_effects = c('studyID'),
  normalization = 'NONE',
  standardize = TRUE)
#separate maaslin2
library(Maaslin2)
input_data <- read.delim("infantvirome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("infantvirome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_Infant_weaning_111323', transform = "NONE",
  fixed_effects = c('TimeSinceWeaning','covidstatus2','NeverPositive','everAntibiotics'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)

library(Maaslin2)
input_data <- read.delim("infantvirome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("infantvirome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_Infant_pptime_111323', transform = "NONE",
  fixed_effects = c('covidstatus2','NeverPositive','everAntibiotics','pp.time'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)

library(Maaslin2)
input_data <- read.delim("infantvirome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("infantvirome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_Infant_pptime_111323', transform = "NONE",
  fixed_effects = c('pp.time','covidstatus2','NeverPositive','everAntibiotics'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)
library(Maaslin2)
input_data <- read.delim("mothervirome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("mothervirome_metadata.txt", row.names = 2,sep = "")

fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_Mother_011924', transform = "NONE",
  fixed_effects = c('covidstatus2','NeverPositive','pp.time','hivStatus','everAntibiotics'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)
#maasline2plots
library(ggplot2)
library(RColorBrewer)
data<-read.delim("maaslin2_motherInfant_personcode.txt")
ggplot(data, aes(x = person_code, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

data<-read.delim("maaslin2_infantOnly_pp.time.txt")
ggplot(data, aes(x = pp.time, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

data<-read.delim("maaslin2_infantOnly_weaning.txt")
ggplot(data, aes(x = weaning, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

data<-read.delim("maaslin2_infantOnly_covidstatus2.txt")
ggplot(data, aes(x = covidstatus2, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

###################Bacteriome Infants vs Mother###################
setwd("/Users/rmaqsood/Documents/pre_post_covid/bacteria/")
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
#person_code significance
set.seed(100)
allrichness<-lme(Richness~  person_code+ pp.time +covidstatus2+ everAntibiotics, random=~1|studyID, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ person_code+ pp.time +covidstatus2+ everAntibiotics, random=~1|studyID, data=m.1)
summary(allalphaDiversity)

#beta diversity
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
d.2<-d.1[,colnames(d.1) %in% m.1$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaBacteriomeMotherInfants_033023.csv")

#PCoA
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_MotherInfantBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive","NeverPositive"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#personCode
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(person_code))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Infant","Mother"),
                     values = c("green", "salmon")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#PERMANOVA
sub_info <- cbind( smp = m.2$Sample, # sample ID
                   id = m.2$ptnum, # subject ID 
                   tp = m.2$pp.time,
                   abx = m.2$everAntibiotics,
                   hivstatus=m.2$hivStatus,
                   covidstatus2 = m.2$covidstatus2,
                   neverPostive = m.2$NeverPositive, 
                   person_code=m.2$person_code, 
                   studyID=m.2$studyID) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaBacteriomeMotherInfants_033023.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
sub_clust$studyID <- as.factor(sub_clust$studyID)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$studyID))
adonis2(wudm ~   person_code + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#Maaslin2
setwd('/Users/rmaqsood/Documents/pre_post_covid/bacteria/')
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance_033023.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
#Infant
infant.Bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
d.2<-d.1[,colnames(d.1) %in% infant.Bacteriome$Sample]
write.table(infant.Bacteriome,"infantBacteriome_metadata.txt")
write.table(d.2,"infantBacteriome_Relabun.txt")
#mother
mother.Bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.3<-d.1[,colnames(d.1) %in% mother.Bacteriome$Sample]
write.table(mother.Bacteriome,"motherBacteriome_metadata.txt")
write.table(d.3,"motherBacteriome_Relabun.txt")
#all
library(Maaslin2)
input_data <- read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance_033023.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("Analysis_bacterial_031423.txt", row.names = 1)
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_MotherInfant_time_011824', transform = "NONE",
  fixed_effects = c('person_code'),
  random_effects = c('studyID'),
  normalization = 'NONE',
  standardize = TRUE)
#separate maaslin2
library(Maaslin2)
input_data <- read.delim("infantBacteriome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("infantBacteriome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_InfantInteraction_011924', transform = "NONE",
  fixed_effects = c('covidstatus2','NeverPositive','pp.time','everAntibiotics'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)
library(Maaslin2)
input_data <- read.delim("infantBacteriome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("infantBacteriome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_Infant_weaning_111323', transform = "NONE",
  fixed_effects = c('timeSinceWeaning','covidstatus2','NeverPositive','everAntibiotics'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)

library(Maaslin2)
input_data <- read.delim("motherBacteriome_Relabun.txt", row.names = 1,sep = "")
input_data2<-t(input_data)
input_metadata <- read.delim("motherBacteriome_metadata.txt", row.names = 2,sep = "")
fit_data <- Maaslin2(
  input_data2, input_metadata, 'maaslinOutput_MotherInteraction_011924', transform = "NONE",
  fixed_effects = c('covidstatus2','NeverPositive','pp.time','hivStatus','everAntibiotics'),
  random_effects = c('ptnum'),
  normalization = 'NONE',
  standardize = TRUE)

#maaslin2 mother+infant plots
library(ggplot2)
library(RColorBrewer)
data<-read.delim("MaaslinMotherInfant_RelAbun_personCode.txt")
ggplot(data, aes(x = person_code, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()
data<-read.delim("MaaslinMotherInfant_RelAbun_pp.time.txt")
ggplot(data, aes(x = Time, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

#maaslin2 infant plots
data<-read.delim("MaaslinInfant_RelAbun_covidstatus2.txt")
ggplot(data, aes(x = covidstatus2, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()
data<-read.delim("MaaslinInfant_RelAbun_abx.txt")
ggplot(data, aes(x = Abx, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()
data<-read.delim("MaaslinInfant_RelAbun_pp.time.txt")
ggplot(data, aes(x = pp.time, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()
data<-read.delim("MaaslinInfant_RelAbun_weaning.txt")
ggplot(data, aes(x = weaning, y = feature.1, size = Prevalence, color = Mean)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 25, 50, 75, 100),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic()

#Kmeans
library(factoextra)
library(cluster)
#calculate optimal number of clusters (sum of squares)
# WUDM<-read.csv("weightedBetaBacteriomeMotherInfants_033023.csv",row.names=1,header = T)
# a<-fviz_nbclust(WUDM, pam, method = "wss")
# a
# b<-fviz_nbclust(WUDM, pam, method = "silhouette")
# b
# #number of clusters vs. gap statistic
# gap_stat <- clusGap(WUDM,
#                     FUN = pam,
#                     K.max = 20, #max clusters to consider
#                     B = 50) #total bootstrapped iterations
# fviz_gap_stat(gap_stat)
#kmeans clusters
set.seed(123)
km.res2<-kmeans(WUDM, 4, iter.max = 10, nstart = 25)
km.res2
set.seed(123)
km.res3<-kmeans(WUDM, 5, iter.max = 10, nstart = 25)
km.res3
set.seed(123)
km.res4<-kmeans(WUDM, 6, iter.max = 10, nstart = 25)
km.res4
set.seed(123)
km.res5<-kmeans(WUDM, 7, iter.max = 10, nstart = 25)
km.res5
set.seed(123)
km.res6<-kmeans(WUDM, 8, iter.max = 10, nstart = 25)
km.res6

#multinomial
setwd("/Users/rmaqsood/Documents/pre_post_covid/bacteria/")
library(mclogit)
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.factor(m.1$ptnum)
m.1$covid_status<-as.factor(m.1$covid_status)
m.1$NeverPositive<-as.factor(m.1$NeverPositive)
m.1$covidstatus2<-as.factor(m.1$covidstatus2)
m.1$hivStatus<-as.factor(m.1$hivStatus)
m.1$everAntibiotics<-as.factor(m.1$everAntibiotics)
m.1$person_code<-as.factor(m.1$person_code)
m.1$kmeans.all<-as.factor(m.1$kmeans.motherInfant)
Subset1<-subset(m.1,kmeans.all =='2'|kmeans.all =='3'|kmeans.all =='4'|kmeans.all =='5')
Subset2<-subset(m.1,kmeans.all =='3'|kmeans.all =='4'|kmeans.all =='5')
Subset3<-subset(m.1,kmeans.all =='4'|kmeans.all =='5')
set.seed(123)
(comm.mblogit <- mblogit(kmeans.all~person_code+NeverPositive+covidstatus2+pp.time+hivStatus+everAntibiotics,  data = m.1, random=~1|studyID))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(kmeans.all~person_code+NeverPositive+covidstatus2+pp.time+hivStatus+everAntibiotics,  data = Subset1, random=~1|ptnum))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(kmeans.all~person_code+NeverPositive+covidstatus2+pp.time+hivStatus+everAntibiotics,  data = Subset2, random=~1|ptnum))
summary(comm.mblogit.subset2)
set.seed(123)
(comm.mblogit.subset3 <- mblogit(kmeans.all~person_code+NeverPositive+covidstatus2+pp.time+hivStatus+everAntibiotics,  data = Subset3, random=~1|ptnum))
summary(comm.mblogit.subset3)


#################HIV*SARS-CoV-2 interaction in mothers#################
#virome
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
summary(m.1)
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
#lme
set.seed(100)
Mother.richness<-lme(Richness~  (covidstatus2*hivStatus)+ pp.time + everAntibiotics, random=~1|ptnum, data=mother.virome)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ (covidstatus2*hivStatus)+ pp.time + everAntibiotics, random=~1|ptnum, data=mother.virome)
summary(Mother.alphaDiversity)
#beta
sub_info <- cbind( smp = mother.virome$Sample, # sample ID
                   id = mother.virome$ptnum, # subject ID 
                   tp = mother.virome$pp.time,
                   abx = mother.virome$everAntibiotics,
                   hivstatus=mother.virome$hivStatus,
                   covidstatus2 = mother.virome$covidstatus2,
                   neverPostive = mother.virome$NeverPositive, 
                   person_code=mother.virome$person_code,
                   studyID=mother.virome$studyID) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetavirome_MotherInfant_042423.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~ covidstatus2*hivstatus + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#bacteriome
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
#lme
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
set.seed(100)
Mother.richness<-lme(Richness~  (covidstatus2*hivStatus)+ pp.time + everAntibiotics, random=~1|ptnum, data=mother.bacteriome)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ (covidstatus2*hivStatus)+ pp.time + everAntibiotics, random=~1|ptnum, data=mother.bacteriome)
summary(Mother.alphaDiversity)
#beta
sub_info <- cbind( smp = mother.bacteriome$Sample, # sample ID
                   id = mother.bacteriome$ptnum, # subject ID 
                   tp = mother.bacteriome$pp.time,
                   abx = mother.bacteriome$everAntibiotics,
                   hivstatus=mother.bacteriome$hivStatus,
                   covidstatus2 = mother.bacteriome$covidstatus2,
                   neverPostive = mother.bacteriome$NeverPositive) # time point
sub_info<-as.data.frame(sub_info, id, tp)
WUDM<-read.csv("weightedBetaBacteriomeAllWomen_Hypothesis2_012723.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~   hivstatus*covidstatus2 + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#################Virome trajectory change#################
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
summary(m.1)
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
#mothers
#interaction
set.seed(100)
Mother.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time + hivStatus+everAntibiotics + (covidstatus2*pp.time), random=~1|ptnum, data=mother.virome)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time +hivStatus+ everAntibiotics + (covidstatus2*pp.time), random=~1|ptnum, data=mother.virome)
summary(Mother.alphaDiversity)
#No interaction
set.seed(100)
Mother.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time + hivStatus+everAntibiotics, random=~1|ptnum, data=mother.virome)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time +hivStatus+ everAntibiotics, random=~1|ptnum, data=mother.virome)
summary(Mother.alphaDiversity)
#plots
loessPlotShannon <- ggplot(mother.virome, aes(x = pp.time, y = AlphaDiversity, color = NeverPositive, fill = NeverPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(mother.virome$pp.time), max(mother.virome$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="virome Alpha Diversity", title = "All women " )
loessPlotShannon

loessPlotRichness <- ggplot(mother.virome, aes(x = pp.time, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(mother.virome$pp.time), max(mother.virome$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="virome Richness", title = "All women" )
loessPlotRichness

#timeSinceCovid ever-SARS-CoV-2 only
#interaction timesincecovid*sars-cov-2-infection
prePostCovid<-mother.virome[ which(mother.virome$NeverPositive=="Positive"),]
set.seed(100)
mother.richness<-lme(Richness~  (covidstatus2*timeSinceCovid)+timeSinceCovid + covidstatus2 + hivStatus+everAntibiotics , random=~1|ptnum, data=prePostCovid)
summary(mother.richness)
set.seed(100)
mother.alphaDiversity<-lme(AlphaDiversity~ (covidstatus2*timeSinceCovid)+timeSinceCovid + covidstatus2 +hivStatus+ everAntibiotics, random=~1|ptnum, data=prePostCovid)
summary(mother.alphaDiversity)
#No interaction
set.seed(100)
mother.richness<-lme(Richness~  pp.time+timeSinceCovid + covidstatus2 + hivStatus+everAntibiotics , random=~1|ptnum, data=prePostCovid)
summary(mother.richness)
set.seed(100)
mother.alphaDiversity<-lme(AlphaDiversity~ pp.time+timeSinceCovid + covidstatus2 +hivStatus+ everAntibiotics, random=~1|ptnum, data=prePostCovid)
summary(mother.alphaDiversity)
#plots
loessPlotShannon <- ggplot(prePostCovid, aes(x = timeSinceCovid, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(prePostCovid$timeSinceCovid), max(prePostCovid$timeSinceCovid), by =1 ),1)) + labs (x= "Time since covid",y="virome Alpha Diversity", title = "Evercovid women " )
loessPlotShannon

loessPlotRichness <- ggplot(prePostCovid, aes(x = timeSinceCovid, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(prePostCovid$timeSinceCovid), max(prePostCovid$timeSinceCovid), by =1 ),1)) + labs (x= "Time since covid",y="virome Richness", title = "Evercovid women" )
loessPlotRichness

#Mother beta diversity
#unweighted women 
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
d.2<-d.1[,colnames(d.1) %in% mother.virome$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"unweightedBetaviromeAllWomen_Hypothesis2_041923.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC3 <- as.matrix(res$vectors[,3])
PC1
PC2
data<-read.delim("PC_virome_allWomen_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
data$PC3 <- PC3
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#hivStatus
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(hivStatus))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.4, .4) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("skyblue","darkorchid4")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = mother.virome$Sample, # sample ID
                   id = mother.virome$ptnum, # subject ID 
                   tp = mother.virome$pp.time,
                   abx = mother.virome$everAntibiotics,
                   hivstatus=mother.virome$hivStatus,
                   covidstatus2 = mother.virome$covidstatus2,
                   neverPostive = mother.virome$NeverPositive) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("unweightedBetaviromeAllWomen_Hypothesis2_041923.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~ tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#weighted all women
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
d.2<-d.1[,colnames(d.1) %in% mother.virome$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaviromeAllWomen_Hypothesis2_041923.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_virome_allWomen_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.4, .4) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#covidstatus2
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covidstatus2))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.4, .4) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#hivStatus
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(hivStatus))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.4, .4) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("skyblue","darkorchid4")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.4, .4) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = mother.virome$Sample, # sample ID
                   id = mother.virome$ptnum, # subject ID 
                   tp = mother.virome$pp.time,
                   abx = mother.virome$everAntibiotics,
                   hivstatus=mother.virome$hivStatus,
                   covidstatus2 = mother.virome$covidstatus2,
                   neverPostive = mother.virome$NeverPositive) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaviromeAllWomen_Hypothesis2_041923.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#evercovid
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
prePostCovid<-mother.virome[ which(mother.virome$NeverPositive=="Positive"),]
prePostCovid$timeSinceCovid<-prePostCovid[, 13] + 10
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
d.2<-d.1[,colnames(d.1) %in% prePostCovid$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaViromeprePostCovidWomen_Hypothesis2_051823.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllWomenBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = prePostCovid$Sample, # sample ID
                   id = prePostCovid$ptnum, # subject ID 
                   tp = prePostCovid$pp.time,
                   abx = prePostCovid$everAntibiotics,
                   hivstatus=prePostCovid$hivStatus,
                   covidstatus2 = prePostCovid$covidstatus2,
                   neverPostive = prePostCovid$NeverPositive, 
                   timeCovid=prePostCovid$timeSinceCovid) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaViromeprePostCovidWomen_Hypothesis2_051823.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~   timeCovid*covidstatus2 + tp +covidstatus2 + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  timeCovid +covidstatus2  + tp+abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#all infant
#interaction infant age
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time +everAntibiotics + (covidstatus2*pp.time), random=~1|ptnum, data=infant.virome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time + everAntibiotics + (covidstatus2*pp.time), random=~1|ptnum, data=infant.virome)
summary(infant.alphaDiversity)
#interaction since weaning
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +TimeSinceWeaning +everAntibiotics + (covidstatus2*TimeSinceWeaning), random=~1|ptnum, data=infant.virome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +TimeSinceWeaning + everAntibiotics + (covidstatus2*TimeSinceWeaning), random=~1|ptnum, data=infant.virome)
summary(infant.alphaDiversity)
#No interaction infant age
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time +everAntibiotics, random=~1|ptnum, data=infant.virome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time + everAntibiotics, random=~1|ptnum, data=infant.virome)
summary(infant.alphaDiversity)
#No interaction since weaning
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +TimeSinceWeaning +everAntibiotics, random=~1|ptnum, data=infant.virome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +TimeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.virome)
summary(infant.alphaDiversity)

#plots pp.time 
loessPlotShannon <- ggplot(infant.virome, aes(x = pp.time, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.virome$pp.time), max(infant.virome$pp.time), by =1 ),1)) + labs (x= "Month of life",y="virome Alpha Diversity", title = "All infants " )
loessPlotShannon

loessPlotRichness <- ggplot(infant.virome, aes(x = pp.time, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +  
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.virome$pp.time), max(infant.virome$pp.time), by =1 ),1)) + labs (x= "Month of life",y="virome Richness", title = "All infants" )
loessPlotRichness

#plots weaning 
loessPlotShannon <- ggplot(infant.virome, aes(x = TimeSinceWeaning, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.virome$TimeSinceWeaning), max(infant.virome$TimeSinceWeaning), by =1 ),1)) + labs (x= "Time since weaning",y="virome Alpha Diversity", title = "All infants " )
loessPlotShannon

loessPlotRichness <- ggplot(infant.virome, aes(x = TimeSinceWeaning, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.virome$TimeSinceWeaning), max(infant.virome$TimeSinceWeaning), by =1 ),1)) + labs (x= "Time since weaning",y="virome Richness", title = "All infants" )
loessPlotRichness

#unweighted Infant beta diversity
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
d.2<-d.1[,colnames(d.1) %in% infant.virome$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"unweightedBetaviromeInfantsHypothesis2_042023.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_virome_infants_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants by month of life")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#time since weaning 
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants by time since weaning")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.virome$Sample, # sample ID
                   id = infant.virome$ptnum, # subject ID 
                   firstPos = infant.virome$firstPositive, # trt, unt
                   tp = infant.virome$pp.time,
                   covidstatus2 = infant.virome$covidstatus2,
                   neverPostive = infant.virome$NeverPositive,
                   abx = infant.virome$everAntibiotics,
                   weaning=infant.virome$TimeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("unweightedBetaviromeInfantsHypothesis2_042023.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   covidstatus2*weaning + covidstatus2 +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  weaning +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")

#weighted infants betaDiversity
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
d.2<-d.1[,colnames(d.1) %in% infant.virome$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaviromeInfantsHypothesis2_042023.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_virome_infants_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.5, .5) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.5, .5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#weaning
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-.6, .6) + ylim(-.5, .5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.virome$Sample, # sample ID
                   id = infant.virome$ptnum, # subject ID 
                   firstPos = infant.virome$firstPositive, # trt, unt
                   tp = infant.virome$pp.time,
                   covidstatus2 = infant.virome$covidstatus2,
                   neverPostive = infant.virome$NeverPositive,
                   abx = infant.virome$everAntibiotics,
                   weaning=infant.virome$TimeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaviromeInfantsHypothesis2_042023.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   covidstatus2*weaning + weaning +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  weaning +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")


#################Bacteriome trajectory change#################
################all women
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
m.1$studyID<-as.character(m.1$studyID)
#lme
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
infant.bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
#mothers
#interaction
set.seed(100)
Mother.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time + hivStatus+everAntibiotics + (covidstatus2*pp.time), random=~1|ptnum, data=mother.bacteriome)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time +hivStatus+ everAntibiotics + (covidstatus2*pp.time), random=~1|ptnum, data=mother.bacteriome)
summary(Mother.alphaDiversity)
#No interaction
set.seed(100)
Mother.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time + hivStatus+everAntibiotics , random=~1|ptnum, data=mother.bacteriome)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time +hivStatus+ everAntibiotics, random=~1|ptnum, data=mother.bacteriome)
summary(Mother.alphaDiversity)

#plots
loessPlotShannon <- ggplot(mother.bacteriome, aes(x = pp.time, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(mother.bacteriome$pp.time), max(mother.bacteriome$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="Bacteriome Alpha Diversity", title = "All women " )
loessPlotShannon

loessPlotRichness <- ggplot(mother.bacteriome, aes(x = pp.time, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(mother.bacteriome$pp.time), max(mother.bacteriome$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="Bacteriome Richness", title = "All women" )
loessPlotRichness

#timeSinceCovid ever-SARS-CoV-2 only
#interaction timesincecovid*sars-cov-2-infection
#timeSinceCovid
prePostCovid<-mother.bacteriome[ which(mother.bacteriome$NeverPositive=="Positive"),]
set.seed(100)
Mother.richness<-lme(Richness~  (covidstatus2*timeSinceCovid)+timeSinceCovid + covidstatus2 + hivStatus+everAntibiotics , random=~1|ptnum, data=prePostCovid)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ (covidstatus2*timeSinceCovid)+timeSinceCovid + covidstatus2 +hivStatus+ everAntibiotics, random=~1|ptnum, data=prePostCovid)
summary(Mother.alphaDiversity)

set.seed(100)
Mother.richness<-lme(Richness~timeSinceCovid + covidstatus2 + hivStatus+everAntibiotics , random=~1|ptnum, data=prePostCovid)
summary(Mother.richness)

precovid<-prePostCovid[ which(prePostCovid$covidstatus2=="Negative"),]
postcovid<-prePostCovid[ which(prePostCovid$covidstatus2=="Positive"),]
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ timeSinceCovid +pp.time+hivStatus+ everAntibiotics, random=~1|ptnum, data=precovid)
summary(Mother.alphaDiversity)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ timeSinceCovid +pp.time+hivStatus+ everAntibiotics, random=~1|ptnum, data=postcovid)
summary(Mother.alphaDiversity)

#plots
loessPlotShannon <- ggplot(prePostCovid, aes(x = timeSinceCovid, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(prePostCovid$timeSinceCovid), max(prePostCovid$timeSinceCovid), by =1 ),1)) + labs (x= "Time since covid",y="virome Alpha Diversity", title = "Evercovid women " )
loessPlotShannon

loessPlotRichness <- ggplot(prePostCovid, aes(x = timeSinceCovid, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(prePostCovid$timeSinceCovid), max(prePostCovid$timeSinceCovid), by =1 ),1)) + labs (x= "Time since covid",y="virome Richness", title = "Evercovid women" )
loessPlotRichness

#unweighted all women
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.Bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.2<-d.1[,colnames(d.1) %in% mother.bacteriome$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"unweightedBetaBacteriomeAllWomen_Hypothesis2_012723.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllWomenBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANONA
sub_info <- cbind( smp = mother.bacteriome$Sample, # sample ID
                   id = mother.bacteriome$ptnum, # subject ID 
                   tp = mother.bacteriome$pp.time,
                   abx = mother.bacteriome$everAntibiotics,
                   hivstatus=mother.bacteriome$hivStatus,
                   covidstatus2 = mother.bacteriome$covidstatus2,
                   neverPostive = mother.bacteriome$NeverPositive) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("unweightedBetaBacteriomeAllWomen_Hypothesis2_012723.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~ tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#weighted all women
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.Bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.2<-d.1[,colnames(d.1) %in% mother.bacteriome$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaBacteriomeAllWomen_Hypothesis2_012723.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllWomenBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#abx
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Abx))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Yes","No"),
                     values = c("hotpink","black")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = mother.bacteriome$Sample, # sample ID
                   id = mother.bacteriome$ptnum, # subject ID 
                   tp = mother.bacteriome$pp.time,
                   abx = mother.bacteriome$everAntibiotics,
                   hivstatus=mother.bacteriome$hivStatus,
                   covidstatus2 = mother.bacteriome$covidstatus2,
                   neverPostive = mother.bacteriome$NeverPositive) # time point
sub_info<-as.data.frame(sub_info, id, tp)
WUDM<-read.csv("weightedBetaBacteriomeAllWomen_Hypothesis2_012723.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  tp +covidstatus2 + neverPostive + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")

#evercovid women
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
prePostCovid<-mother.bacteriome[ which(mother.bacteriome$NeverPositive=="Positive"),]
#prePostCovid$timeSinceCovid<-prePostCovid[, 13] + 10
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
d.2<-d.1[,colnames(d.1) %in% prePostCovid$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaBacteriomeprePostCovidWomen_Hypothesis2_051623.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllWomenBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = prePostCovid$Sample, # sample ID
                   id = prePostCovid$ptnum, # subject ID 
                   tp = prePostCovid$pp.time,
                   abx = prePostCovid$everAntibiotics,
                   hivstatus=prePostCovid$hivStatus,
                   covidstatus2 = prePostCovid$covidstatus2,
                   neverPostive = prePostCovid$NeverPositive, 
                   timeCovid=prePostCovid$timeSinceCovid) # time point
sub_info<-as.data.frame(sub_info, id, tp)
WUDM<-read.csv("weightedBetaBacteriomeprePostCovidWomen_Hypothesis2_051623.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~   timeCovid*covidstatus2 + tp +covidstatus2 + abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  timeCovid +covidstatus2  + tp+abx+hivstatus, data = sub_clust, permutations = perm, by = "margin")


################all infant
#interaction infant age
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time + everAntibiotics + (covidstatus2*pp.time) , random=~1|ptnum, data=infant.bacteriome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time + everAntibiotics + (covidstatus2*pp.time) , random=~1|ptnum, data=infant.bacteriome)
summary(infant.alphaDiversity)
#interaction since weaning
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +timeSinceWeaning + everAntibiotics + (covidstatus2*timeSinceWeaning), random=~1|ptnum, data=infant.bacteriome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +timeSinceWeaning + everAntibiotics + (covidstatus2*timeSinceWeaning) , random=~1|ptnum, data=infant.bacteriome)
summary(infant.alphaDiversity)
#No interaction infant age
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +pp.time + everAntibiotics, random=~1|ptnum, data=infant.bacteriome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +pp.time + everAntibiotics, random=~1|ptnum, data=infant.bacteriome)
summary(infant.alphaDiversity)
#No interaction since weaning
set.seed(100)
infant.richness<-lme(Richness~ NeverPositive + covidstatus2 +timeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.bacteriome)
summary(infant.richness)
set.seed(100)
infant.alphaDiversity<-lme(AlphaDiversity~ NeverPositive + covidstatus2 +timeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.bacteriome)
summary(infant.alphaDiversity)

#plots pp.time 
loessPlotShannon <- ggplot(infant.bacteriome, aes(x = pp.time, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.bacteriome$pp.time), max(infant.bacteriome$pp.time), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Alpha Diversity", title = "Infants" )
loessPlotShannon

loessPlotRichness <- ggplot(infant.bacteriome, aes(x = pp.time, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.bacteriome$pp.time), max(infant.bacteriome$pp.time), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Richness", title = "Infants" )
loessPlotRichness

loessPlotRichness <- ggplot(infant.bacteriome, aes(x = pp.time, y = Richness, color = NeverPositive, fill = NeverPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.bacteriome$pp.time), max(infant.bacteriome$pp.time), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Richness", title = "Infants" )
loessPlotRichness

#plots weaning 
loessPlotShannon <- ggplot(infant.bacteriome, aes(x = timeSinceWeaning, y = AlphaDiversity, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.bacteriome$timeSinceWeaning), max(infant.bacteriome$timeSinceWeaning), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Alpha Diversity", title = "Infants" )
loessPlotShannon

loessPlotRichness <- ggplot(infant.bacteriome, aes(x = timeSinceWeaning, y = Richness, color = covid_status, fill = covid_status)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.bacteriome$timeSinceWeaning), max(infant.bacteriome$timeSinceWeaning), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Richness", title = "Infants" )
loessPlotRichness

#beta diversity
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.Bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
d.2<-d.1[,colnames(d.1) %in% infant.Bacteriome$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"BetaBacteriomeInfantsHypothesis2_020623_unweighted.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllInfantsBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#time since weaning 
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.Bacteriome$Sample, # sample ID
                   id = infant.Bacteriome$ptnum, # subject ID 
                   firstPos = infant.Bacteriome$firstPositive, # trt, unt
                   tp = infant.Bacteriome$pp.time,
                   covidstatus2 = infant.Bacteriome$covidstatus2,
                   neverPostive = infant.Bacteriome$NeverPositive,
                   abx = infant.Bacteriome$everAntibiotics,
                   exclsuveBF=infant.Bacteriome$exclusiveBF,
                   weaning=infant.Bacteriome$timeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("BetaBacteriomeInfantsHypothesis2_020623_unweighted.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   covidstatus2*weaning + covidstatus2 +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  weaning +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")

#weighted infants betaDiversity
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.Bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
d.2<-d.1[,colnames(d.1) %in% infant.Bacteriome$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"BetaBacteriomeInfantsHypothesis2_020623_weighted.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllInfantsBacteriome_hypothesis2.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#covid_status
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(covid_status))) + 
  geom_point() + xlim(-1.5, 1.5) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("neverpositive","precovid","postcovid"),
                     values = c("darkgreen", "lightpink","deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.Bacteriome$Sample, # sample ID
                   id = infant.Bacteriome$ptnum, # subject ID 
                   firstPos = infant.Bacteriome$firstPositive, # trt, unt
                   tp = infant.Bacteriome$pp.time,
                   covidstatus2 = infant.Bacteriome$covidstatus2,
                   neverPostive = infant.Bacteriome$NeverPositive,
                   abx = infant.Bacteriome$everAntibiotics,
                   exclsuveBF=infant.Bacteriome$exclusiveBF, 
                   weaning=infant.Bacteriome$timeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("BetaBacteriomeInfantsHypothesis2_020623_weighted.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~   tp*covidstatus2 + tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~   covidstatus2*weaning + weaning +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  tp +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(100)
adonis2(wudm ~  weaning +covidstatus2 + neverPostive + abx, data = sub_clust, permutations = perm, by = "margin")


#################Virome immediate change, first positive#################
################all women
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
firstPosOnly<-mother.virome[ which(mother.virome$firstPositive!='Positive.2'),]
firstPosOnly<-firstPosOnly[ which(firstPosOnly$firstPositive!='Positive.3'),]
#lme
set.seed(100)
Mother.richness<-lme(Richness~ firstPositive + pp.time + everAntibiotics+ hivStatus, random=~1|ptnum, data=firstPosOnly)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ firstPositive + pp.time + everAntibiotics + hivStatus, random=~1|ptnum, data=firstPosOnly)
summary(Mother.alphaDiversity)
#plots
loessPlotShannon <- ggplot(firstPosOnly, aes(x = pp.time, y = AlphaDiversity, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(firstPosOnly$pp.time), max(firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="virome Alpha Diversity", title = "All women " )
loessPlotShannon

loessPlotRichness <- ggplot(firstPosOnly, aes(x = pp.time, y = Richness, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(firstPosOnly$pp.time), max(firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="virome Richness", title = "All women" )
loessPlotRichness

#unweighted beta diversity
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
firstPosOnly<-mother.virome[ which(mother.virome$firstPositive!='Positive.2'),]
firstPosOnly<-firstPosOnly[ which(firstPosOnly$firstPositive!='Positive.3'),]
d.2<-d.1[,colnames(d.1) %in% firstPosOnly$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"unweightedBetaViromeHypothesis1_041823.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_virome_allWomen_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = firstPosOnly$Sample, # sample ID
                   id = firstPosOnly$ptnum, # subject ID 
                   firstPos = firstPosOnly$firstPositive, # trt, unt
                   tp = firstPosOnly$pp.time,
                   abx = firstPosOnly$everAntibiotics,
                   hivstatus=firstPosOnly$hivStatus,
                   exclsuveBF=firstPosOnly$exclusiveBF) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("unweightedBetaViromeHypothesis1_041823.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~  tp + firstPos + abx + hivstatus, data = sub_clust, permutations = perm, by = "margin")

#weighted all women betaDiversity
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.virome<-m.1[ which(m.1$person_code=='Mother'),]
firstPosOnly<-mother.virome[ which(mother.virome$firstPositive!='Positive.2'),]
firstPosOnly<-firstPosOnly[ which(firstPosOnly$firstPositive!='Positive.3'),]
d.2<-d.1[,colnames(d.1) %in% firstPosOnly$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaViromeHypothesis1_041823.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC3 <- as.matrix(res$vectors[,3])
PC1
PC2
data<-read.delim("PC_virome_allWomen_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-.3, .3) + ylim(-.5, .5) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-.3, .3) + ylim(-.5, .5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#hivStatus
p <- ggplot(data, aes(x=PC2, y=PC3, color=as.factor(hivStatus))) + 
  geom_point() + xlim(-.5, .5) + ylim(-.5, .5) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("skyblue","darkorchid4")) +
  theme(aspect.ratio=1)+labs(title="All women by covid status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = firstPosOnly$Sample, # sample ID
                   id = firstPosOnly$ptnum, # subject ID 
                   firstPos = firstPosOnly$firstPositive, # trt, unt
                   tp = firstPosOnly$pp.time,
                   abx = firstPosOnly$everAntibiotics,
                   hivstatus=firstPosOnly$hivStatus,
                   exclsuveBF=firstPosOnly$exclusiveBF) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaViromeHypothesis1_041823.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
adonis2(wudm ~  tp + firstPos + abx + hivstatus, data = sub_clust, permutations = perm, by = "margin")

################All infant
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
infant.firstPosOnly<-infant.virome[ which(infant.virome$firstPositive!='Positive.2'),]
infant.firstPosOnly<-infant.firstPosOnly[ which(infant.firstPosOnly$firstPositive!='Positive.3'),]
#lme
set.seed(100)
Infant.richness<-lme(Richness~ firstPositive + pp.time + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.richness)
set.seed(100)
Infant.alphaDiversity<-lme(AlphaDiversity~ firstPositive + pp.time + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.alphaDiversity)
#Infants timesinceWeaing
set.seed(100)
Infant.richness<-lme(Richness~ firstPositive + TimeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.richness)
set.seed(100)
Infant.alphaDiversity<-lme(AlphaDiversity~ firstPositive + TimeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.alphaDiversity)

#plots
loessPlotShannon <- ggplot(infant.firstPosOnly, aes(x = pp.time, y = AlphaDiversity, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.firstPosOnly$pp.time), max(infant.firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Month of life",y="virome Alpha Diversity", title = "Infants" )
loessPlotShannon

loessPlotRichness <- ggplot(infant.firstPosOnly, aes(x = pp.time, y = Richness, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.firstPosOnly$pp.time), max(infant.firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Month of life",y="virome Richness", title = "Infants" )
loessPlotRichness
#timesinceweaning
loessPlotShannon <- ggplot(infant.firstPosOnly, aes(x = TimeSinceWeaning, y = AlphaDiversity, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.firstPosOnly$TimeSinceWeaning), max(infant.firstPosOnly$TimeSinceWeaning), by =1 ),1)) + labs (x= "Time since weaning",y="virome Alpha Diversity", title = "Infants" )
loessPlotShannon

loessPlotRichness <- ggplot(infant.firstPosOnly, aes(x = TimeSinceWeaning, y = Richness, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.firstPosOnly$TimeSinceWeaning), max(infant.firstPosOnly$TimeSinceWeaning), by =1 ),1)) + labs (x= "Time since weaning",y="virome Richness", title = "Infants" )
loessPlotRichness

#unweighted infants betaDiversity
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
infant.firstPosOnly<-infant.virome[ which(infant.virome$firstPositive!='Positive.2'),]
infant.firstPosOnly<-infant.firstPosOnly[ which(infant.firstPosOnly$firstPositive!='Positive.3'),]
d.2<-d.1[,colnames(d.1) %in% infant.firstPosOnly$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"unweightedBetaViromeInfantsHypothesis1_041923.csv")

res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_virome_infants_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#time since weaning
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative by time since weaning")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.firstPosOnly$Sample, # sample ID
                   id = infant.firstPosOnly$ptnum, # subject ID 
                   firstPos = infant.firstPosOnly$firstPositive, # trt, unt
                   tp = infant.firstPosOnly$pp.time,
                   abx = infant.firstPosOnly$everAntibiotics,
                   exclsuveBF=infant.firstPosOnly$exclusiveBF, 
                   weaning=infant.firstPosOnly$TimeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("unweightedBetaViromeInfantsHypothesis1_041923.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(123)
adonis2(wudm ~  tp + firstPos + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(123)
adonis2(wudm ~  weaning +firstPos + abx, data = sub_clust, permutations = perm, by = "margin")

#weighted infants betaDiversity
m.1 <- read.delim("Analysis_virome_030723.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePost_CleanNormalizedRPKM_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.virome<-m.1[ which(m.1$person_code=='Infant'),]
infant.firstPosOnly<-infant.virome[ which(infant.virome$firstPositive!='Positive.2'),]
infant.firstPosOnly<-infant.firstPosOnly[ which(infant.firstPosOnly$firstPositive!='Positive.3'),]
d.2<-d.1[,colnames(d.1) %in% infant.firstPosOnly$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaViromeInfantsHypothesis1_041923.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_virome_infants_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1, 1) + ylim(-.5, .5) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-.5, .5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#timeSinceWeaning
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1, 1) + ylim(-.5, .5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.firstPosOnly$Sample, # sample ID
                   id = infant.firstPosOnly$ptnum, # subject ID 
                   firstPos = infant.firstPosOnly$firstPositive, # trt, unt
                   tp = infant.firstPosOnly$pp.time,
                   abx = infant.firstPosOnly$everAntibiotics,
                   exclsuveBF=infant.firstPosOnly$exclusiveBF, 
                   weaning=infant.firstPosOnly$TimeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaViromeInfantsHypothesis1_041923.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(123)
adonis2(wudm ~  tp + firstPos + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(123)
adonis2(wudm ~  weaning +firstPos + abx, data = sub_clust, permutations = perm, by = "margin")

#################Bacteriome immediate change, first positive#################
################all women
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
firstPosOnly<-mother.bacteriome[ which(mother.bacteriome$firstPositive!='Positive.2'),]
firstPosOnly<-firstPosOnly[ which(firstPosOnly$firstPositive!='Positive.3'),]
#lme
set.seed(100)
Mother.richness<-lme(Richness~ firstPositive + pp.time + everAntibiotics+ hivStatus, random=~1|ptnum, data=firstPosOnly)
summary(Mother.richness)
set.seed(100)
Mother.alphaDiversity<-lme(AlphaDiversity~ firstPositive + pp.time + everAntibiotics + hivStatus, random=~1|ptnum, data=firstPosOnly)
summary(Mother.alphaDiversity)
#plots
loessPlotShannon <- ggplot(firstPosOnly, aes(x = pp.time, y = AlphaDiversity, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(firstPosOnly$pp.time), max(firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="Bacteriome Alpha Diversity", title = "All women " )
loessPlotShannon

loessPlotRichness <- ggplot(firstPosOnly, aes(x = pp.time, y = Richness, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(firstPosOnly$pp.time), max(firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Post-partum month",y="Bacteriome Richness", title = "All women" )
loessPlotRichness

#unweighted WLHIV
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
firstPosOnly<-mother.bacteriome[ which(mother.bacteriome$firstPositive!='Positive.2'),]
firstPosOnly<-firstPosOnly[ which(firstPosOnly$firstPositive!='Positive.3'),]
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.Bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.2<-d.1[,colnames(d.1) %in% firstPosOnly$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"unweightedBetaBacteriomeAllWomen_Hypothesis1_012723.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllWomenBacteriome_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = firstPosOnly$Sample, # sample ID
                   id = firstPosOnly$ptnum, # subject ID 
                   firstPos = firstPosOnly$firstPositive, # trt, unt
                   tp = firstPosOnly$pp.time,
                   abx = firstPosOnly$everAntibiotics,
                   hivstatus=firstPosOnly$hivStatus,
                   exclsuveBF=firstPosOnly$exclusiveBF) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("unweightedBetaBacteriomeAllWomen_Hypothesis1_012723.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~  tp + firstPos + abx + hivstatus, data = sub_clust, permutations = perm, by = "margin")

#weighted all women
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
mother.bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
firstPosOnly<-mother.bacteriome[ which(mother.bacteriome$firstPositive!='Positive.2'),]
firstPosOnly<-firstPosOnly[ which(firstPosOnly$firstPositive!='Positive.3'),]
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
mother.Bacteriome<-m.1[ which(m.1$person_code=='Mother'),]
d.2<-d.1[,colnames(d.1) %in% firstPosOnly$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBetaBacteriomeAllWomen_Hypothesis1_012723.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllWomenBacteriome_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1, 1) + ylim(-.5, .5) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-.5, .5) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All women first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#Abx
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(Abx))) + 
  geom_point() + xlim(-1, 1) + ylim(-.5, .5) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Yes","No"),
                     values = c("deeppink","black")) +
  theme(aspect.ratio=1)+labs(title="All women by HIV status")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = firstPosOnly$Sample, # sample ID
                   id = firstPosOnly$ptnum, # subject ID 
                   firstPos = firstPosOnly$firstPositive, # trt, unt
                   tp = firstPosOnly$pp.time,
                   abx = firstPosOnly$everAntibiotics,
                   hivstatus=firstPosOnly$hivStatus,
                   exclsuveBF=firstPosOnly$exclusiveBF) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("weightedBetaBacteriomeAllWomen_Hypothesis1_012723.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(100)
adonis2(wudm ~  tp + firstPos + abx + hivstatus, data = sub_clust, permutations = perm, by = "margin")


################all infants
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
infant.bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
infant.firstPosOnly<-infant.bacteriome[ which(infant.bacteriome$firstPositive!='Positive.2'),]
infant.firstPosOnly<-infant.firstPosOnly[ which(infant.firstPosOnly$firstPositive!='Positive.3'),]
#pp.time
set.seed(100)
Infant.richness<-lme(Richness~ firstPositive + pp.time + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.richness)
set.seed(100)
Infant.alphaDiversity<-lme(AlphaDiversity~ firstPositive + pp.time + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.alphaDiversity)
#weaning
set.seed(100)
Infant.richness<-lme(Richness~ firstPositive + timeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.richness)
set.seed(100)
Infant.alphaDiversity<-lme(AlphaDiversity~ firstPositive + timeSinceWeaning + everAntibiotics, random=~1|ptnum, data=infant.firstPosOnly)
summary(Infant.alphaDiversity)

loessPlotShannon <- ggplot(infant.firstPosOnly, aes(x = pp.time, y = AlphaDiversity, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.firstPosOnly$pp.time), max(infant.firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Alpha Diversity", title = "Infants" )
loessPlotShannon

loessPlotRichness <- ggplot(infant.firstPosOnly, aes(x = pp.time, y = Richness, color = firstPositive, fill = firstPositive)) +
  geom_point() +
  stat_smooth(method="lm", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(infant.firstPosOnly$pp.time), max(infant.firstPosOnly$pp.time), by =1 ),1)) + labs (x= "Month of life",y="Bacteriome Richness", title = "Infants" )
loessPlotRichness

#unweighted infants betaDiversity
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.Bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
infant.firstPosOnly<-infant.Bacteriome[ which(infant.Bacteriome$firstPositive!='Positive.2'),]
infant.firstPosOnly<-infant.firstPosOnly[ which(infant.firstPosOnly$firstPositive!='Positive.3'),]
d.2<-d.1[,colnames(d.1) %in% infant.firstPosOnly$Sample]
d.2[d.2>0] <-1
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"BetaBacteriomeInfantsHypothesis1_012723_unweighted.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllInfantsBacteriome_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#ExclusiveBF
p <- ggplot(data, aes(x=PC1, y=PC2, color=(exclusiveBF))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Yes","No"),
                     values = c("navajowhite2", "azure4")) +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative by breast feeding")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#time since weaning
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1, 1) + ylim(-1, 1) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative by time since weaning")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.firstPosOnly$Sample, # sample ID
                   id = infant.firstPosOnly$ptnum, # subject ID 
                   firstPos = infant.firstPosOnly$firstPositive, # trt, unt
                   tp = infant.firstPosOnly$pp.time,
                   abx = infant.firstPosOnly$everAntibiotics,
                   exclsuveBF=infant.firstPosOnly$exclusiveBF, 
                   weaning=infant.firstPosOnly$timeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("BetaBacteriomeInfantsHypothesis1_012723_unweighted.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(123)
adonis2(wudm ~  tp + firstPos + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(123)
adonis2(wudm ~  weaning +firstPos + abx, data = sub_clust, permutations = perm, by = "margin")

#weighted infants betaDiversity
m.1 <- read.delim("Analysis_bacterial_031423.txt", header=TRUE)
m.1$ptnum<-as.character(m.1$ptnum)
d.1<-read.delim("LKPrePostCovid_bacteriaData_cleanNormalized_relativeAbundance.txt", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
m.2<-m.1[!m.1$Richness=="#N/A",]
infant.Bacteriome<-m.1[ which(m.1$person_code=='Infant'),]
infant.firstPosOnly<-infant.Bacteriome[ which(infant.Bacteriome$firstPositive!='Positive.2'),]
infant.firstPosOnly<-infant.firstPosOnly[ which(infant.firstPosOnly$firstPositive!='Positive.3'),]
d.2<-d.1[,colnames(d.1) %in% infant.firstPosOnly$Sample]
beta_dist <- vegdist(t(d.2),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"BetaBacteriomeInfantsHypothesis1_012723_weighted.csv")
res <- pcoa(distMatrix)
PC1 <- as.matrix(res$vectors[,1])
PC2 <- as.matrix(res$vectors[,2])
PC1
PC2
data<-read.delim("PC_AllInfantsBacteriome_hypothesis1.txt", header=TRUE, row.names=1) #Use data table in betaDiversity.2 tab in supplementary excel file
data$PC1 <- PC1
data$PC2 <- PC2
summary(data)
#firstPositive
p <- ggplot(data, aes(x=PC1, y=PC2, color=as.factor(firstPositive))) + 
  geom_point() + xlim(-1.5, 1.2) + ylim(-1.2, 1.2) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c("springgreen4", "deeppink")) +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#pp.time
p <- ggplot(data, aes(x=PC1, y=PC2, color=(pp.time))) + 
  geom_point() + xlim(-1.5, 1.2) + ylim(-1.2, 1.2) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#timeSinceWeaning
p <- ggplot(data, aes(x=PC1, y=PC2, color=(timeSinceWeaning))) + 
  geom_point() + xlim(-1.5, 1.2) + ylim(-1.2, 1.2) +
  scale_color_gradient (low="yellow", high="red") +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative over time")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#ExclusiveBF
p <- ggplot(data, aes(x=PC1, y=PC2, color=(exclusiveBF))) + 
  geom_point() + xlim(-1.5, 1.2) + ylim(-1.2, 1.2) +
  stat_ellipse() +
  scale_color_manual(breaks = c("Yes","No"),
                     values = c("navajowhite2", "azure4")) +
  theme(aspect.ratio=1)+labs(title="All infants first postive vs all Negative by breast feeding")
p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#PERMANOVA
sub_info <- cbind( smp = infant.firstPosOnly$Sample, # sample ID
                   id = infant.firstPosOnly$ptnum, # subject ID 
                   firstPos = infant.firstPosOnly$firstPositive, # trt, unt
                   tp = infant.firstPosOnly$pp.time,
                   abx = infant.firstPosOnly$everAntibiotics,
                   exclsuveBF=infant.firstPosOnly$exclusiveBF,
                   weaning=infant.firstPosOnly$timeSinceWeaning) # time point
sub_info<-as.data.frame(sub_info)
WUDM<-read.csv("BetaBacteriomeInfantsHypothesis1_012723_weighted.csv", header=T, row.names = 1)
x <- setdiff(sub_info$smp, rownames(WUDM))
sub_info[sub_info$smp %in% x, 'smp'] <- NA
sub_wide <- reshape(sub_info, timevar = "tp", idvar = "id", v.names = "smp", direction = "wide")
sub_clust <- sub_info[which(sub_info$smp %in% rownames(WUDM)), ]
sub_clust$id <- as.factor(sub_clust$id)
sub_clust$tp <- as.factor(sub_clust$tp)
ord <- match(sub_clust$smp, rownames(WUDM))
WUDM <- WUDM[ord, ord]
wudm <- as.dist(WUDM)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = sub_clust$id))
set.seed(123)
adonis2(wudm ~  tp + firstPos + abx, data = sub_clust, permutations = perm, by = "margin")
set.seed(123)
adonis2(wudm ~  weaning +firstPos + abx, data = sub_clust, permutations = perm, by = "margin")













