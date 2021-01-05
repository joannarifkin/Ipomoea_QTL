library(dplyr)
library(stringr)


#Import data ############
setwd("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Writing/Manuscript/Final_data_for_upload/")
#setwd("D:/Dropbox/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Writing/Manuscript/Final_data_for_upload/")

flowerdata=read.csv("Clean_raw_floral_phenotypes.csv",header=TRUE, stringsAsFactors = F) 

phenologydata=read.csv("Clean_raw_phenology.csv",header=TRUE, stringsAsFactors = F)

pollendata_imageJ=read.csv("ImageJ_raw_4_3_15.csv",header=TRUE, stringsAsFactors = F)
  
pollendata_Coulter<-read.csv("Coulter_pollen.csv",header=TRUE, stringsAsFactors = F)

phenology2013<-read.csv("2013_growout_phenology_clean.csv",header=TRUE, stringsAsFactors = F)
flowers2013<-read.csv("2013_growout_floral_clean.csv",header=TRUE, stringsAsFactors = F)



#Convert measurements in flower data ############

flowerdata$L_over_W<-flowerdata$Cor_Length/flowerdata$Cor_Width
flowerdata$LV_over_L<-flowerdata$Cor_tissue_Length/flowerdata$Cor_Length
flowerdata$Nectar.volume.uL<-pi*0.2*0.2*flowerdata$Nectar_volume_mm
flowerdata$Date_measured<-as.Date(flowerdata$Date_measured, "%m/%d/%Y")
flowerdata$Sterility[is.na(flowerdata$Sterility)]<-0

phenologydata$Germinated<-as.Date(phenologydata$Germinated, "%m/%d/%Y")
phenologydata$Planted<-as.Date(phenologydata$Planted, "%m/%d/%Y")
phenologydata$Leaf.1.open<-as.Date(phenologydata$Leaf.1.open, "%m/%d/%Y")
phenologydata$Leaf.2.open<-as.Date(phenologydata$Leaf.2.open, "%m/%d/%Y")
phenologydata$Leaf.3.open<-as.Date(phenologydata$Leaf.3.open, "%m/%d/%Y")

pollendata_imageJ$Date<-as.Date(pollendata_imageJ$Date, "%m/%d")

flowers2013$L_over_W<-flowers2013$Cor..Length/flowers2013$Cor..Width
flowers2013$Nectar.volume.uL<-pi*0.2*0.2*flowers2013$Nectar.volume

########################## Inspect raw data ##########

colnames(flowerdata)

par(mfrow=c(3,3))
hist(flowerdata$Date_measured, breaks="weeks")
#View(flowerdata)

hist(flowerdata$Color_numeric)
hist(flowerdata$Cyme_length)
hist(flowerdata$flowers_per_day)
hist(flowerdata$Flowers_buds_on_cymepeduncle)
hist(flowerdata$Cor_Length)
hist(flowerdata$Cor_tissue_Length)
hist(flowerdata$Cor_Width)
hist(flowerdata$Anther_stigma_position)
hist(flowerdata$Nectar_volume_mm)
hist(flowerdata$Nectar.volume.uL)
hist(flowerdata$Style_Length)
hist(flowerdata$Sterility)

colnames(phenologydata)
hist(phenologydata$Germinated, breaks="days")
hist(phenologydata$Leaf.1.open, breaks="days")
hist(phenologydata$Leaf.2.open, breaks="days")
hist(phenologydata$Leaf.3.open, breaks="days")
hist(phenologydata$Height_mm_day_21)
hist(phenologydata$Internode_1_day_21)
hist(phenologydata$Internode_2_day_21)
hist(phenologydata$Internode_3_day_21)
# 5/19 / 2020 flower and phenology data look clean

colnames(pollendata_Coulter)
hist(pollendata_Coulter$avg.number)
hist(pollendata_Coulter$avg.size)
hist(pollendata_Coulter$samples)
#View(pollendata_Coulter)

colnames(pollendata_imageJ)
hist(pollendata_imageJ$Mean_aliquots)
# 5/ 20 / 2020 pollen data look clean

par(mfrow=c(3,3))
hist(flowers2013$Opening)
hist(flowers2013$Color)
hist(flowers2013$Peduncle.length)
hist(flowers2013$Flowers.buds.on.peduncle)
hist(flowers2013$internode.below)
hist(flowers2013$Cor..Length)
hist(flowers2013$Cor..Width)
hist(flowers2013$Anther.stigma.position)
hist(flowers2013$Stigma.Length)
hist(flowers2013$Nectar.volume.uL)


hist(phenology2013$Germination.Day)
hist(phenology2013$Cotyledon.day)
hist(phenology2013$Height_mm_day_21)
hist(phenology2013$Internode_1_day_21)
hist(phenology2013$Internode_2_day_21)
hist(phenology2013$Internode_3_day_21)
hist(phenology2013$Number_leaves_day_23)
hist(phenology2013$Leaf.1.open)
hist(phenology2013$Leaf.2.open)
hist(phenology2013$Leaf.3.open)
hist(phenology2013$Flower.1.Open)


###### Average measurements by individual #######


colnames(flowerdata)

flowerdata<-flowerdata%>%select(., ID, ID_numeric, Date_measured,  Color_numeric, Cyme_length, flowers_per_day, Flowers_buds_on_cymepeduncle, Cor_Length, Cor_tissue_Length, Cor_Width, Anther_stigma_position,  Nectar_volume_mm, Nectar.volume.uL, Style_Length, Sterility, L_over_W,  LV_over_L)
select(flowerdata, ID_numeric)

means_flowerdata<-flowerdata %>% group_by(ID) %>% add_tally() %>% summarise_all(., mean, na.rm=T)
#View(means_flowerdata) 

#View(avg_flowerdata)
colnames(means_flowerdata)[3]<-"Flower_date"
colnames(means_flowerdata)[18]<-"n_flowers"


avg_pollendata_imageJ<-pollendata_imageJ %>% group_by(ID) %>% summarise_all(., mean)
avg_pollendata_imageJ$Pollen_per_flower<-(avg_pollendata_imageJ$Mean_aliquots*5)
#View(avg_pollendata_imageJ)

pollendata_Coulter
avg_pollendata_Coulter<-pollendata_Coulter %>% group_by(ID) %>% summarise_all(., mean)
avg_pollendata_Coulter



colnames(flowers2013)

flowers2013<-flowers2013%>%select(., Species, ID, Opening, Color, Peduncle.length, Flowers.buds.on.peduncle, internode.below, Cor..Length, Cor..Width, Anther.stigma.position, Stigma.Length, Nectar.volume, L_over_W, Nectar.volume.uL)

means_flowers_2013<-flowers2013 %>% group_by(ID, Species) %>% add_tally() %>% summarise_all(., mean,na.rm=T)
warnings()
#View(means_flowerdata) 

phenology2013<-phenology2013%>%select(., Species, ID, Germination.Day, Height_mm_day_21, Internode_1_day_21, Internode_2_day_21, Internode_3_day_21, Number_leaves_day_23, Leaf.1.open, Leaf.2.open, Leaf.3.open, Flower.1.Open)


##### Join dataframes together ##### 

phenotypes <- left_join(phenologydata, means_flowerdata, by="ID")
phenotypes <- left_join(phenotypes, avg_pollendata_imageJ, by="ID")
phenotypes <- left_join(phenotypes, avg_pollendata_Coulter, by="ID")
#View(phenotypes)

phenotypes2013<-left_join(phenology2013, means_flowers_2013, by="ID")
#View(phenotypes2013)
View(means_flowers_2013)

##### Add transformed phenotype columns #######
colnames(phenotypes)
phenotypes$TwoOrMore<-as.numeric(phenotypes$n_flowers>=2)
phenotypes$ThreeOrMore<-as.numeric(phenotypes$n_flowers>=3)
phenotypes$Flower_date<-as.Date(phenotypes$Flower_date)
phenotypes$Germinated<-as.Date(phenotypes$Germinated)
phenotypes$Day_germ<-as.numeric(phenotypes$Flower_date-phenotypes$Germinated)
phenotypes$LnNectaruLPlusOne<-log(phenotypes$Nectar.volume.uL+1)
phenotypes$more_than_one_flower_per_day<-as.numeric(phenotypes$flowers_per_day>1)
phenotypes$flowers_per_day_greater_than_one<-case_when(
  (phenotypes$flowers_per_day>1)~phenotypes$flowers_per_day,
    (!phenotypes$flowers_per_day>1) ~ NA_real_)
phenotypes$more_than_one_flower_on_cyme<-as.numeric(phenotypes$Flowers_buds_on_cymepeduncle>1)
phenotypes$flowers_per_cyme_greater_than_one<-case_when(
  (phenotypes$Flowers_buds_on_cymepeduncle>1)~phenotypes$Flowers_buds_on_cymepeduncle,
  (!phenotypes$Flowers_buds_on_cymepeduncle>1) ~ NA_real_)
phenotypes$LN_Nectar.volume<-log(phenotypes$Nectar_volume_mm+1)
phenotypes$LN_Nectar.volume.uL<-log(phenotypes$Nectar.volume.uL+1)
phenotypes$Leaf.1.adjusted<-as.numeric(as.Date(phenotypes$Leaf.1.open)-as.Date(phenotypes$Germinated))
phenotypes$Leaf.2.adjusted<-as.numeric(as.Date(phenotypes$Leaf.2.open)-as.Date(phenotypes$Germinated))
phenotypes$Leaf.3.adjusted<-as.numeric(as.Date(phenotypes$Leaf.3.open)-as.Date(phenotypes$Germinated))
phenotypes$neg_recip_L1_adj<-(-1/as.numeric(phenotypes$Leaf.1.adjusted))
phenotypes$neg_recip_L2_adj<-(-1/as.numeric(phenotypes$Leaf.2.adjusted))
phenotypes$neg_recip_L3_adj<-(-1/as.numeric(phenotypes$Leaf.3.adjusted))
phenotypes$n_flowers<-case_when(
  (is.na(phenotypes$Leaf.1.adjusted)==F & is.na(phenotypes$n_flowers == T)) ~ 0,
  (is.na(phenotypes$Leaf.1.adjusted)==F & phenotypes$n_flowers>0) ~ phenotypes$n_flowers,
  (is.na(phenotypes$Leaf.1.adjusted)==T ~ NA_real_))
phenotypes$Ever_flowered<-case_when(
  (is.na(phenotypes$Leaf.1.adjusted)==F)~as.numeric(phenotypes$n_flowers>0),
  TRUE ~ NA_real_)

######## Export phenotype files ##########



#write.csv(phenotypes, quote=F, paste("All_phenotype_means_",Sys.Date(),".csv",sep = ""))
#write.csv(phenotypes2013, quote=F, paste("2013_phenotype_means_",Sys.Date(),".csv",sep = ""))

#phenotypes<-read.csv("All_phenotype_means.csv", stringsAsFactors = F)
#View(phenotypes)

##### Export F2 data for QTL ##### 

#phenotypes<-read.csv("All_phenotype_means.csv")
phenotypes<-read.csv("All_phenotype_means_2020-05-25.csv")
F2means<-subset(phenotypes, phenotypes$species=="F2")

F2means$ID<-str_pad(F2means$ID, 3, "left", "0")
F2means$ID<-str_c(F2means$species, F2means$ID, sep="_") #Change ID names to match map data 

colnames(F2means)
F2means_QTL<-select(F2means, c(ID, Leaf.1.adjusted, Leaf.2.adjusted, Leaf.3.adjusted, neg_recip_L1_adj, neg_recip_L2_adj, neg_recip_L3_adj, Internode_1_day_21, Internode_2_day_21, Internode_3_day_21, Height_mm_day_21, Number_leaves_day_21, Day_germ, Ever_flowered, more_than_one_flower_per_day, flowers_per_day_greater_than_one, more_than_one_flower_on_cyme, flowers_per_cyme_greater_than_one, TwoOrMore, ThreeOrMore, Color_numeric, Cyme_length, flowers_per_day, Flowers_buds_on_cymepeduncle, Cor_Length, Cor_tissue_Length, Cor_Width, L_over_W, LV_over_L, Style_Length, Anther_stigma_position, Nectar.volume.uL, LnNectaruLPlusOne, n_flowers, Pollen_per_flower, avg.number, avg.size, Sterility))

#mutate(F2means, )
#paste species "_" id

head(F2means)
#str(F2means)
#write.csv(F2means_QTL, quote=F, row.names=F, paste("F2_means_for_QTL_analysis_",Sys.Date(),".csv",sep = ""))
#write.csv(F2means_QTL, quote=F, row.names=F, paste("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/F2_means_for_QTL_analysis_",Sys.Date(),".csv",sep = ""))
#F2means<-read.csv("F2_means_for_QTL_analysis.csv")
getwd()



######## Summaries by species #####

str(phenotypes)
species_means_SDs<-phenotypes %>% group_by(species) %>% summarise_all(funs(mean,  sd), na.rm = TRUE)
warnings()
#View(species_means_SDs)


write.csv(species_means_SDs, quote=F, paste("Summary_data_by_species_",Sys.Date(),".csv", sep = ""))


species_means_SDs_2013<-phenotypes2013 %>% group_by(Species.x) %>% summarise_all(funs(mean,  sd), na.rm = TRUE)
#write.csv(species_means_SDs_2013, quote=F, paste("2013_Summary_data_by_species_", Sys.Date(), ".csv"))



str(phenotypes2013)

###################### Variation among individuals and number flowers sampled ############

Five_sampled<-means_flowerdata[ID_freq_index5,] #Have to do this on MEANS, not on the original data
Five_sampled
names(Five_sampled)
summary(Five_sampled)

Three_sampled<-means_flowerdata[ID_freq_index3,] #Have to do this on MEANS, not on the original data
Three_sampled
names(Three_sampled)
summary(Three_sampled)

Two_sampled<-means_flowerdata[ID_freq_index2,] #Have to do this on MEANS, not on the original data
Two_sampled
names(Two_sampled)
summary(Two_sampled)


#############Ratio of purple to white flowers##############
if(!require(dplyr)){install.packages("dplyr")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(grid)){install.packages("grid")}
if(!require(pwr)){install.packages("pwr")}

table(phenotypes$Color_numeric)
observed = c(283, 73)        # observed frequencies
expected = c(0.75, 0.25)      # expected proportions

chisq.test(x = observed,
           p = expected)


#data:  observed
#X-squared = 3.8352, df = 1, p-value = 0.05019

#The overall ratio of purple:white did not differ from 3:1

table(subset(phenotypes, phenotypes$n_flowers>4)$Color_numeric)

observed = c(136, 27)        # observed frequencies
expected = c(0.75, 0.25)      # expected proportions

chisq.test(x = observed,
           p = expected)
#The ratio did differ in plants that produced at least 5 flowers, suggesting underlying phenological differences
#X-squared = 6.1861, df = 1, p-value = 0.01288


136/(136+27)
27/(136+27)


table(phenotypes$Sterility)
348+11

observed = c(348, 11)        # observed frequencies
expected = c((1-1/15), (1/15))      # expected proportions

chisq.test(x = observed,
           p = expected)
#> 11/359
#[1] 0.03064067
#> 1/16
#[1] 0.0625
#Less than 1/16 had sterile morphology
#X-squared = 7.4883, df = 1, p-value = 0.00621




############  Compare phenotypes of F2 individuals with more vs. less flowers ###################
#F2s all means with Ns

#phenotypes<-read.csv("All_phenotype_means.csv")

#F2means<-subset(phenotypes, phenotypes$species=="F2")
head(F2means)

#Comparing more-flowering to less-flowering plants

bartlett.test(F2means$Day_germ~F2means$ThreeOrMore,data=F2means)
t.test(F2means$Day_germ~F2means$ThreeOrMore,data=F2means,var.equal=FALSE, verbose=T)

bartlett.test(F2means$Nectar.volume.uL~F2means$ThreeOrMore,data=F2means)
t.test(F2means$Nectar.volume.uL~F2means$ThreeOrMore,data=F2means,var.equal=FALSE, verbose=T)

bartlett.test(F2means$Cyme_length~F2means$ThreeOrMore,data=F2means)
t.test(F2means$Cyme_length~F2means$ThreeOrMore,data=F2means,var.equal=FALSE, verbose=T)

bartlett.test(F2means$Flowers_buds_on_cymepeduncle~F2means$ThreeOrMore,data=F2means)
t.test(F2means$Flowers_buds_on_cymepeduncle~F2means$ThreeOrMore,data=F2means,var.equal=FALSE, verbose=T)

bartlett.test(F2means$flowers_per_day~F2means$ThreeOrMore,data=F2means)
t.test(F2means$flowers_per_day~F2means$ThreeOrMore,data=F2means,var.equal=FALSE, verbose=T)


bartlett.test(F2means$Cor_Length~F2means$ThreeOrMore,data=F2means)
t.test(F2means$Cor_Length~F2means$ThreeOrMore,data=F2means,var.equal=FALSE, verbose=T)

################## F2 means vs. midparent #########################


F2means
head(F2means)
#attach(F2means)

observed    = F2means$Cyme_length
theoretical = 18.39728507
t.test(observed,
       mu = theoretical,
       conf.int = 0.95, na.rm=T)



observed    = Flowers_buds_on_cymepeduncle
theoretical = 1.28280543
t.test(observed,
       mu = theoretical,
       conf.int = 0.95, na.rm=T)

observed    = Cor_Length
theoretical = 24.75995475
t.test(observed,
       mu = theoretical,
       conf.int = 0.95, na.rm=T)

observed    = Cor_Length_variant
theoretical = 28.33557692
t.test(observed,
       mu = theoretical,
       conf.int = 0.95, na.rm=T)

observed    = Cor._Width
theoretical = 19.98755656
t.test(observed,
       mu = theoretical,
       conf.int = 0.95, na.rm=T)

observed    = L_over_W
theoretical = 1.271971556
t.test(observed,
       mu = theoretical,
       conf.int = 0.95, na.rm=T)

observed    = LV_over_L
  theoretical = 1.141553871
  t.test(observed,
         mu = theoretical,
         conf.int = 0.95, na.rm=T)

  observed    = Style_Length
  theoretical = 14.34253394
  t.test(observed,
         mu = theoretical,
         conf.int = 0.95, na.rm=T)

  observed    = Anther_stigma_position
  observed  = subset(F2means, F2means$Sterility==0)$Anther_stigma_position
  sd(subset(F2means, F2means$Sterility==0)$Anther_stigma_position)
  theoretical = -0.053167421
  t.test(observed,
         mu = theoretical,
         conf.int = 0.95, na.rm=T)
  
  
  observed    = Nectar.volume.uL
  theoretical = 0.655299268
  t.test(observed,
         mu = theoretical,
         conf.int = 0.95, na.rm=T)

  
  observed    = Leaf.1.adjusted
  theoretical =11.37536232
  t.test(observed,
         mu = theoretical,
         conf.int = 0.95, na.rm=T)

  observed    = Leaf.2.adjusted
    theoretical = 13.42149758
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
  observed    = Leaf.3.adjusted
    theoretical = 15.77222222
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
  observed    = Height_mm_day_21
    theoretical = 383.9426457
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
  observed    = Internode_1_day_21
    theoretical = 10.83186864
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  

  observed    = Internode_2_day_21
    theoretical = 12.37257169
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
  observed    = Internode_3_day_21
    theoretical = 45.06822387
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
    
    observed    = Number_leaves_day_21
    theoretical = 5.965309898
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
    
    observed    = Pollen_per_flower
    theoretical = 698.125
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
 
    
    observed    = avg.number
    theoretical = 423
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    
    observed    = avg.size
    theoretical = 75.415
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
  
 ##F1 vs. midparent means (2013 growout)   ################
    View(phenology2013)
    colnames(phenotypes2013)
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Height_mm_day_21
    theoretical = 413.9051724
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Internode_1_day_21
    theoretical = 9.681034483
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Internode_2_day_21
    theoretical = 14.61206897
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Internode_3_day_21
    theoretical = 52.93965517
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Number_leaves_day_23
    theoretical = 7.974137931
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    
    colnames(phenotypes2013)
    
    
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Leaf.1.open
    theoretical = 13.19889163
        t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
    
    observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Leaf.2.open
    theoretical = 15.51631773
        t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
          observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Leaf.3.open
    theoretical = 18.37315271
      t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
    
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Opening
      theoretical = 86.74310516
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Peduncle.length
      theoretical = 8.707328386
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Flowers.buds.on.peduncle
      theoretical = 1.203756957
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      
    
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$internode.below
      theoretical = 51.10844988
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Cor..Length
      theoretical = 22.65847724
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Cor..Width
      theoretical = 19.10105965
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
    
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Anther.stigma.position
      theoretical = -0.323001701
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Stigma.Length
      theoretical = 12.15438776
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Nectar.volume
      theoretical = 0.440456865
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$L_over_W
      theoretical = 1.21999658
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      
      observed    = subset(phenotypes2013,phenotypes2013$Species.x=="F1")$Nectar.volume.uL
      theoretical = 0.055349442
      t.test(observed,
             mu = theoretical,
             conf.int = 0.95, na.rm=T)
      

###########F1 vs.  midparent floral traits (flowers from clone) #######

    #cyme length
         observed    = c(15,14.7,16.7,13.7,17.2,19)
    theoretical = 18.39728507
    t.test(observed,
           mu = theoretical,
           conf.int = 0.95, na.rm=T)
  
    #corolla length
    observed    = c(24.5,24.6,24.6,24.2,24.1,25.3)
       theoretical = 24.75995475
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
 
    #corolla tissue length
       observed    = c(28.1,28.5,29,28.3,26.6,28.8)
       theoretical = 28.33557692
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
      
       
       #corolla width
       observed    = c(20.7,21.4,21.3,20.5,18.3,20.6)
       theoretical = 19.98755656
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
       
       
       #Corolla shape (length / width)
       observed    = c(1.183574879,1.14953271,1.154929577,1.180487805,1.316939891,1.22815534)
       theoretical = 1.271971556
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
       
       #Corolla shape (tissue length / length)
       
       observed    = c(1.146938776,1.158536585,1.178861789,1.169421488,1.10373444,1.138339921)
       theoretical = 1.141553871

         t.test(observed,
                mu = theoretical,
                conf.int = 0.95, na.rm=T)
       
         
         #Nectar volume uL
       observed    = c(0.452389342,0.376991118,0.540353936,0.251327412,0,0.753982237)
       theoretical = 0.655299268
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
       
       #Style length
       observed    = c(15.3,15.1,15.8,15.6,15.6,16.1)
       theoretical = 14.34253394
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
       
       #F1 vs.  midparent early growth traits (individuals)
       phenologydata
       View(phenologydata)
       #Height
       colnames(phenotypes)
       observed    = subset(phenotypes, phenotypes$species=="F1")$Height_mm_day_21
       theoretical = 383.9426457
       t.test(observed,
              mu = theoretical,
              conf.int = 0.95, na.rm=T)
       
        
       "Internode_1_day_21"                "Internode_2_day_21"               
       [13] "Internode_3_day_21"                "Number_leaves_day_21"     
       "Leaf.1.adjusted"                  
       [51] "Leaf.2.adjusted"                   "Leaf.3.adjusted"        
         


############# Correlations ##################

rcorrAllSpearman<-rcorr(as.matrix(F2means,type="spearman"))
RvalsSp<-rcorrAllSpearman$r
PvalsSp<-rcorrAllSpearman$P
#write.csv(PvalsSp, "All_P_spearman.csv")
#write.csv(RvalsSp, "All_R_spearman.csv")


############# Correlation plots #############
    
#Import F2 means
attach(F2means)
library(GGally)
library(ggcorrplot)
library(ggplot2)
library(reshape2)
library(corrplot)
library(dplyr)   
library(psych)

colnames(F2means)

### Subset traits for correlation matrix #####
F2_corr_subset<-as.matrix(
  select(
    F2means,c(
      Cor_Length, L_over_W, LV_over_L, Cor_tissue_Length, Cor_Width, Style_Length, Cyme_length, Flowers_buds_on_cymepeduncle, Height_mm_day_21, Internode_1_day_21, Internode_2_day_21, Internode_3_day_21, Number_leaves_day_21, Day_germ, n_flowers, flowers_per_day, Leaf.1.adjusted, Leaf.2.adjusted, Leaf.3.adjusted, Nectar.volume.uL, avg.number, Pollen_per_flower, avg.size
      )))

colnames(F2_corr_subset)<-c("CL", "LW", "TLL", "TL", "CW", "SL", "CY", "FC", "H", "INT1", "INT2", "INT3", "NL", "FD", "NF", "FPD", "L1", "L2", "L3", "NVL", "P1", "P2", "PS")


##### Plot with ggpairs  #####


warnings()

plot<-ggpairs(as.data.frame(F2_corr_subset),use = "pairwise.complete.obs",  lower = list(continuous = wrap("points", alpha=0.3, size=0.5), combo = "facethist", discrete = "facetbar", na = "na"),
        upper = list(continuous = wrap("cor", combo = "box", discrete = "facetbar", na = "na", size=2)),label_round=0,
)+labs(title="Correlations between traits")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=10), axis.text = element_text(size=4))
plot
ggsave("Correlation_plot_1-4-21.pdf",width = 12, height = 9, units = "in")
ggsave("Correlation_plot_1-4-21.png",width = 12, height = 9, units = "in")
