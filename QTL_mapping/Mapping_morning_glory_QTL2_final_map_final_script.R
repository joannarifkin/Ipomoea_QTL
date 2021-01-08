#install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
#update.packages(ask=FALSE, repos="http://rqtl.org/qtl2cran")
#install.packages("qtl")
#install.packages("iplots")
#install.packages("devtools")
#install_github("jtlovell/qtlTools")
library(qtl2)
library(qtl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
library(tibble)
library(iplots)
library(plotly)
library(ASMap)
library(devtools)
library(qtlTools)
#sessionInfo()
#?cbind.scan1
#Import data from file ------------
ipomoea<- read_cross2 ("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/ipomoea.yaml")
setwd("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/")
#ipomoea<- read_cross2 ("ipomoea.yaml")

###############Check phenotype distributions ############
as.list(colnames(ipomoea$pheno))

par(mfrow=c(2,2))
for (i in c(1:37))
     {
    print(i)
      hist(ipomoea$pheno[,i], main = colnames(ipomoea$pheno)[i])
            }

for (i in c(1:37))
{
  qqnorm(ipomoea$pheno[,i], main = colnames(ipomoea$pheno)[i])
}

hist(ipomoea$pheno[,6])

# Leaf.1.adjusted # skewed, transformed to negative reciprocal
# Leaf.2.adjusted # skewed, transformed to negative reciprocal
# Leaf.3.adjusted # skewed, transformed to negative reciprocal
# neg_recip_L1_adj # normal
# neg_recip_L2_adj # normal
# neg_recip_L3_adj # normal
# Internode_1_day_21 # normal
# Internode_2_day_21 # skewed, but other internodes are normal
# Internode_3_day_21 # normal
# Height_mm_day_21 # somewhat skewed
# Number_leaves_day_21 # somewhat skewed
# Day_germ # normal
# Ever_flowered # binary
# more_than_one_flower_per_day # binary
# flowers_per_day_greater_than_one # skewed
# more_than_one_flower_on_cyme # binary
# flowers_per_cyme_greater_than_one # skewed
# TwoOrMore # binary
# ThreeOrMore # binary
# Color_numeric # binary
# Cyme_length # normal
# flowers_per_day # skewed
# Flowers_buds_on_cymepeduncle # skewed
# Cor_Length # normal
# Cor_tissue_Length # normal
# Cor_Width # normal
# L_over_W # normal
# LV_over_L # normal
# Style_Length # normal
# Anther_stigma_position # narrow?
# Nectar.volume.uL # skewed
# LnNectaruLPlusOne # still skewed
# n_flowers # skewed
# Pollen_per_flower # normalish
# avg.number # normalish
# avg.size # normalish
# Sterility # binary


############### Make R/QTL1 input file #########
setwd("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/QTL1_input_files")
map<-read.csv("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/ipomoea_gmap.csv", stringsAsFactors = F, header = T)
genotypes<-read.csv("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/ipomoea_geno_F2_only.csv", stringsAsFactors = F, header = T)
phenotypes<-read.csv("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Final_map/F2_means_for_QTL_analysis_2020-05-25.csv", stringsAsFactors = F, header = T)
colnames(map)
colnames(genotypes)
rownames(phenotypes)
map_w_genotypes<-left_join(map, genotypes, by=c(marker="CHR"))
map_w_genotypes[map_w_genotypes==21]<-12
colnames(map_w_genotypes)[2:3]<-""
View(map_w_genotypes)
write.csv(map_w_genotypes, "ipomoea_gen_rot.csv", quote = F, row.names = F)
head(phenotypes)
rownames(phenotypes)<-phenotypes$ID
transposed_phenotypes<-as.data.frame(t(as.matrix(phenotypes)))
transposed_phenotypes <- transposed_phenotypes[c(2:38,1),]
#colnames(transposed_phenotypes)<-NULL
#IDs<-as.list(transposed_phenotypes[38,])

#transposed_phenotypes[38,]
#View(transposed_phenotypes)
transposed_phenotypes$marker<-rownames(transposed_phenotypes)
transposed_phenotypes <- transposed_phenotypes[,c(509,1:508)]

#write.csv(transposed_phenotypes, "ipomoea_phe_rot.csv", row.names = F, quote=F)
#ipomoea1<-read.cross("csvsr", ".",  "ipomoea_gen_rot.csv","ipomoea_phe_rot.csv",  alleles = c("1", "2"))
ipomoea1<-read.cross("csvsr", "E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/QTL1_input_files",  "ipomoea_gen_rot.csv","ipomoea_phe_rot.csv",  genotypes = c("11", "12", "22", "1-", "2-"))

########## Inspect map ############
help("est.rf")

class(ipomoea$gmap)
plot.map(ipomoea$gmap)
plotRF(ipomoea1)
plotRF(ipomoea1,chr = 12)
plotRF(ipomoea1,chr = 14)
plotRF(ipomoea1,chr = 15)
#plotRF(listeria.a4)
help(plotRF)
############## Data diagnostics #################
"est_rf_f2"

percent_missing <- n_missing(ipomoea, "ind", "prop")*100
#View(percent_missing)
#No missing data because of Lep-Map's phased output

cg <- compare_geno(ipomoea, cores=0)
summary(cg)
#No pairs with proportion matches above 0.9.

pr <- calc_genoprob(ipomoea, error_prob=0.002, map_function="c-f", cores=0)
m <- maxmarg(pr, minprob=0.5, cores=0)
nxo <- count_xo(m, cores=0)
totxo <- rowSums(nxo)
plot(totxo)
sum(as.numeric(totxo>70))
#10 individuals have more than 70 crossovers, and should probably be removed. 
too_many_xo<-names(subset(totxo, totxo>70))
okay_individuals<-subset(ind_ids(ipomoea),(!ind_ids(ipomoea) %in% too_many_xo))
ipomoea_filtered<-subset(ipomoea, ind=okay_individuals)
#Remove 10 individuals with many double crossovers

e <- calc_errorlod(ipomoea_filtered, pr, cores=0)
e <- do.call("cbind", e)
errors_ind <- rowSums(e>2)/n_typed(ipomoea_filtered)*100
#No unlikely genotypes
ipomoea<-ipomoea_filtered

############# Segregation distortion ###########
ipomoea2 <-convert2bcsft(ipomoea1, F.gen = 2, BC.gen=0,estimate.map = FALSE)

profileMark(ipomoea2, stat.type = c("seg.dist", "prop", "dxo", "recomb"),
               layout = c(1, 5), type = "l")
ipomoea3 <- pullCross(ipomoea2, type = "seg.distortion", pars =
                      list(seg.thresh = "bonf"))

distortion_table<-ipomoea3$seg.distortion$table
getwd()
write.csv(distortion_table, "Bonferroni_distorted_markers.csv")

profileMark(ipomoea2, stat.type = c("seg.dist", "prop", "recomb"),
            layout = c(1, 7), type = "l")

profileMark(ipomoea2, stat.type = c("seg.dist", "prop", "recomb"),
            layout = c(1, 7), type = "l", chr = 3)

profileMark(ipomoea2, stat.type = c("seg.dist", "prop", "recomb"),
            layout = c(1, 7), type = "l", chr = 10)
profileMark(ipomoea2, stat.type = c("seg.dist", "prop", "recomb"),
            layout = c(1, 7), type = "l", chr = 11)
profileMark(ipomoea2, stat.type = c("seg.dist", "prop", "recomb"),
            layout = c(1, 7), type = "l", chr = 15)

###Impute genotypes --------------
map <- insert_pseudomarkers(ipomoea$gmap, step=1)
#pull.map(ipomoea$gmap)
pr <- calc_genoprob(ipomoea, map, error_prob=0.002)
#pr <- calc_genoprob(ipomoea, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)

#version()
#operm <- scan1perm(pr[,"1"], ipomoea$pheno,  n_perm=10)
#colnames(ipomoea$pheno)
#warnings()

#Permutation test to identify significance thresholds ---------------------
operm <- scan1perm(pr, ipomoea$pheno, n_perm=1000)
operm_fig <- scan1perm(pr, ipomoea$pheno, n_perm=100)

colnames(ipomoea$pheno)
#write.table(operm, "Permutation_test_6-3-20.txt")
operm_summary<-summary(operm, alpha=c(0.2, 0.1, 0.05,0.01))
#write.table(operm_summary,"Summary_permutation_test_6-3-20.txt")
operm_bin <- scan1perm(pr, ipomoea$pheno[,c(13:14, 16, 18:19)], model="binary",
                       n_perm=1000) #adding bintol and eta_max parameters
operm_bin_color <- scan1perm(pr, ipomoea$pheno[,c(20)], model="binary",
                       n_perm=1000,bintol = 1e-4, eta_max=15) #adding bintol and eta_max 
operm_bin_color_summary<-summary(operm_bin_color, alpha=c(0.2, 0.1, 0.05,0.01))
#write.table(operm_bin_color_summary,"Summary_bin_permutation_test_color_6-8-20-bin_tol_1e4_etamax_15.txt")

operm_bin_sterility <- scan1perm(pr, ipomoea$pheno[,c(37)], model="binary",
                             n_perm=1000,bintol = 1e-4, eta_max=15) #adding bintol and eta_max 
operm_bin_sterility_summary<-summary(operm_bin_sterility, alpha=c(0.2, 0.1, 0.05,0.01))
#write.table(operm_bin_sterility_summary,"Summary_bin_permutation_test_sterility_6-8-20-bin_tol_1e4_etamax_15.txt")

#write.table(operm_bin, "Binary_permutation_test_6-4-20_did_not_converge.txt")
#write.table(operm_bin, "Binary_permutation_test_6-4-20-bin_tol_1e4_etamax_18_still_not_converged_fewer_warnings.txt")
#write.table(operm_bin, "Binary_permutation_test_6-4-20-bin_tol_1e4_etamax_16_still_not_converged_one_warning.txt")
#write.table(operm_bin, "Binary_permutation_test_6-4-20-bin_tol_1e4_etamax_15_still_not_converged_one_warning.txt")
operm_bin<-read.table( "Binary_permutation_test_6-4-20-bin_tol_1e4_etamax_15_still_not_converged_one_warning.txt")
operm_bin<-as.scanoneperm(operm_bin)
?summary.scanoneperm
operm_bin_summary<-summary.scanoneperm(operm_bin, alpha=c(0.2, 0.1, 0.05,0.01))
#write.table(operm_bin_summary,"Summary_bin_permutation_test_6-4-20-bin_tol_1e4_etamax_15_did_not_converge.txt")
#write.table(operm_bin_summary,"Summary_bin_permutation_test_6-4-20-bin_tol_1e4_etamax_16_did_not_converge.txt")
#write.table(operm_bin_summary,"Summary_bin_permutation_test_6-4-20-bin_tol_1e4_etamax_18_did_not_converge.txt")
#write.table(operm_bin_summary,"Summary_bin_permutation_test_6-4-20_did_not_converge.txt")
#write.table(operm_bin_summary,"Summary_bin_permutation_test_5-28-20_bintol_1e3_did_not_converge.txt")
warnings()

colnames(ipomoea$pheno)
operm_bin_ever_flowered <- scan1perm(pr, ipomoea$pheno[,c(13)], model="binary",
                       n_perm=1000,bintol = 1e-4, eta_max=15) #adding bintol and eta_max 


setwd("D:/Dropbox/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Chromosome_specific_permutation_tests/")

for (i in 1:15) {
  operm <- scan1perm(pr[,i], ipomoea$pheno,  n_perm=1000)
  summary_operm<-summary(operm)
  write.table(operm, paste(i,"_permutations_6_8_20.txt"))
  write.csv(summary_operm, paste(i,"_permutation_summary_6_8_20.txt"))
  operm_binary <- scan1perm(pr[,i], ipomoea$pheno, model="binary", n_perm=1000)
  summary_operm_binary<-summary(operm_binary)
  write.table(operm_binary, paste(i,"binary_permutations_6_8_20.txt"))
  write.csv(summary_operm_binary, paste(i,"binary_permutation_summary_6_8_20.txt"))
}
setwd("E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/Chromosome_specific_permutation_tests/")
library(qtl2)




warnings()

#write.table(operm,"Permutation_test_8-14-19.txt")
write.table(operm,"Permutation_test_10-3-19.txt")
#write.table(operm_bin, "Binary_permutation_test_8-14-19.txt")
#write.table(operm_bin, "Binary_permutation_test_8-20-19.txt")
write.table(operm_bin, "Binary_permutation_test_10-3-19.txt")
operm_summary<-summary(operm, alpha=c(0.2, 0.1, 0.05,0.01))
#write.table(operm_summary,"Summary_permutation_test_8-14-19.txt")











getwd()
summary(operm_bin, alpha=c(0.2,0.1, 0.05,0.01))

colnames(ipomoea$pheno)
operm_SFS <- scan1perm(pr, ipomoea$pheno[,c(31)], model="binary",
                       n_perm=1000)

operm_more1flower <- scan1perm(pr, ipomoea$pheno[,c(39,41)], model="binary",
                       n_perm=1000)

operm_bin <- scan1perm(pr, ipomoea$pheno[,c(2:4,7,31)], model="binary",
                       n_perm=1000)
operm<-read.table("Permutation_test_5-22-18.txt")
operm_bin <- scan1perm(pr, ipomoea$pheno[,c(2:4, 7, 31, 39,41)], model="binary",
                       n_perm=1000)


#Scan for QTL---------------------
#Normal traits		--------------	
#7 Internode_1_day_21 # normal
#8 Internode_2_day_21 # skewed, but other internodes are normal
#9 Internode_3_day_21 # normal
#10 Height_mm_day_21 # somewhat skewed
#11 Number_leaves_day_21 # somewhat skewed
#12 Day_germ # normal
#15 flowers_per_day_greater_than_one # skewed
#17 flowers_per_cyme_greater_than_one # skewed
#21 Cyme_length # normal
#24 Cor_Length # normal
#25 Cor_tissue_Length # normal
#26 Cor_Width # normal
#27 L_over_W # normal
#28 LV_over_L # normal
#29 Style_Length # normal
#31 Nectar.volume.uL # skewed
#33 n_flowers # skewed
#34 Pollen_per_flower # normalish
#35 avg.number # normalish
#36 avg.size # normalish
colnames(ipomoea$pheno)
out <- scan1(pr, ipomoea$pheno[,c(7:9,10:12, 15,17, 21,24:29, 31,33:36 )])
#out <- scan1(pr, ipomoea$pheno[,c(5,6)])
#QTLs_normal<-find_peaks(out, map, threshold=2, drop = 1.5, expand2markers = T)
############# Find and export peaks #######################
QTLs_normal<-find_peaks(out, map, threshold=2.4, drop = 1.5, expand2markers = T) #set relative to lowest possible chromosome specific value
QTLs_normal$pos_marker<-find_marker(map, QTLs_normal$chr, QTLs_normal$pos) #Expanded to nearest marker
QTLs_normal$pos_no_pseudo_marker<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$pos) #using ipomoea$gmap rather than map allows me to find the marker rather htan the pseudomarker
QTLs_normal$pos_ci_lo<-find_marker(map, QTLs_normal$chr, QTLs_normal$ci_lo) #Expanded to nearest marker 
QTLs_normal$pos_ci_hi<-find_marker(map, QTLs_normal$chr, QTLs_normal$ci_hi) #Expanded to nearest marker



#QTLs_normal$pos_marker<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$pos) #using ipomoea$gmap rather than map allows me to find the marker rather htan the pseudomarker
#QTLs_normal$pos_ci_lo<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$ci_lo)
#QTLs_normal$pos_ci_hi<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$ci_hi)
find_marker(map, c3.loc156)
out <- scan1(pr, ipomoea$pheno)
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out) # overall maximum LOD score
plot(out, map, lodcolumn=6, col="slateblue", ylim=c(0, ymx*1.02))
find_peaks(out, map, threshold=4, drop=1.5)
find_peaks(out, map, threshold=2, drop=1.5)
help(find_peaks)

#QTL effects exported as table with loop ---------------
?fit1
#out_fit1 <- fit1(pr[[7]][,,152.281], ipomoea$pheno[,"Height_mm_day_21"], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
#out_fit1$lod
#out_fit1$coef
#out_fit1$SE
#out_fit1$coef
#mean(c(418.6686,448.7646,473.4724))
#colnames(ipomoea$pheno)
#?scan1coef 
QTLs_normal[,12:24]<-NA
par(mfrow=c(2,2))
for (row in 1:nrow(QTLs_normal)) {
  chr_var<-as.numeric(QTLs_normal$chr[row])
  pos_var<-QTLs_normal$pos_marker[row]
  pos_cm<-QTLs_normal$pos[row]
  pheno_var<-QTLs_normal$lodcolumn[row]
  print(row)
  print(chr_var)
  print(class(chr_var))
  print(pos_var)
  print(class(pos_var))
  print(pheno_var)
  print(class(pheno_var))
  mu_a_d_effects<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef
  allele_effects<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var])$coef
  mu_a_d_effects_no_zero<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)),zerosum=FALSE)$coef
  allele_effects_no_zero<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var],zerosum=FALSE)$coef
  print(mu_a_d_effects)
  print(allele_effects)
  QTLs_normal[row,12:14]<-mu_a_d_effects
  QTLs_normal[row,15:18]<-allele_effects
  QTLs_normal[row,19:21]<-mu_a_d_effects_no_zero
  QTLs_normal[row,22:24]<-allele_effects_no_zero
  g <- maxmarg(pr, map, chr=chr_var, pos=pos_cm, return_char=TRUE)
  plot_pxg(g, ipomoea$pheno[,paste(pheno_var)], ylab=pheno_var, main=c(chr_var, pos_cm))
}
View(QTLs_normal)
warnings()
colnames(QTLs_normal)
colnames(QTLs_normal)[12:24]<-c("mu_zero","a_zero","d_zero","11_zero","12_zero","22_zero","intercept_zero","mu_non_zero","a_non_zero","d_non_zero","11_non_zero","12_non_zero","22_non_zero")
head(QTLs_normal)
#write.csv(QTLs_normal, "QTL_effects_normal_6_10-20.csv")
#write.csv(QTLs_normal, "QTL_effects_normal_6_8-20.csv")
#write.csv(QTLs_normal, "QTL_effects_normal_8_15.csv")
#write.csv(QTLs_normal, "QTL_effects_normal_no_pseudomarkers_8_15.csv")
#write.csv(QTLs_normal, "LowLOD_QTL_effects_normal_no_pseudomarkers_8_21.csv")
#write.csv(QTLs_normal, "QTL_effects_normal_all_8_26.csv")
#write.csv(QTLs_normal, "QTL_effects_normal_all_10_2.csv")

View(QTLs_normal)

#effects<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef
#effects<-fit1(pr[[1]][,,0], ipomoea$pheno[,"Height_mm_day_21"], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef

#effects<-fit1(pr[[13]][,,138.102], ipomoea$pheno[,"Pollen.count_Stephen"], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef
#out_fit1 <- fit1(pr[[7]][,,152.281], ipomoea$pheno[,"Height_mm_day_21"], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef
out_fit1$

#QTLs_normal[,11:13]<- fit1(pr[[as.numeric(QTLs_normal$chr[4])]][,,QTLs_normal$pos], ipomoea$pheno[,QTLs_normal$lodcolumn], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coeff


#QTLs_normal %>% mutate (mu = fit1(pr[[as.numeric(chr)]][,,QTLs_normal$pos], ipomoea$pheno[,QTLs_normal$lodcolumn], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coeff[1])


#QTL effects



?scan1coef
rownames(c2eff)
c2eff <- scan1coef(pr[,"5"], ipomoea$pheno[,"Color"], model="binary")
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["5"], columns=1:3, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
max(c2eff[,3])
View(c2eff)
c2effB <- scan1coef(pr[,"5"], ipomoea$pheno[,"Color"],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB, map["5"], columns=2:3, col=col)
last_coef <- unclass(c2effB)[nrow(c2effB),2:3] # last two coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
max(c2effB[,1])
max(c2effB[,2])
max(c2effB[,3])

View(c2effB)

#Plot QTLs for normal traits#########################################################

colnames(out)
plot
warnings()



plot(out, map, lodcolumn=1, col="slateblue", main="Internode_1_day_21", ylim=c(0,6))
abline(h=c(4.097476013,4.880234629), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=2, col="slateblue", main="Internode_2_day_21", ylim=c(0,6))
abline(h=c(4.124270924,4.819989357), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=3, col="slateblue", main="Internode_3_day_21", ylim=c(0,6))
abline(h=c(4.097356447,4.764062094), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=4, col="slateblue", main="Height_mm_day_21", ylim=c(0,6))
abline(h=c(4.17494575, 4.947938713 ), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=5, col="slateblue", main="Number_leaves_day_21", ylim=c(0,7))
abline(h=c(4.045449907,5.310292574), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=6, col="slateblue", main="Day_germ",ylim=c(0,6))
abline(h=c(4.074185078, 4.809182532), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=7, col="slateblue", main="flowers_per_day_greater_than_one", ylim=c(0,6))
abline(h=c(4.169153453,     5.442333329), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=8, col="slateblue", main="flowers_on_cyme_greater_than_one", ylim=c(0,6))
abline(h=c(4.215169924,     5.185238944), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=9, col="slateblue", main="Peduncle_length",ylim=c(0,5))
abline(h=c(4.062029123,    4.786526072), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=10, col="slateblue", main="Cor_Length")
abline(h=c(4.042632469,4.71485117), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=11, col="slateblue", main="Cor_tissue_Length")
abline(h=c(4.207653076,4.791886693), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=12, col="slateblue", main="Cor_Width",ylim=c(0,5))
abline(h=c(4.132319085,4.897233709), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=13, col="slateblue", main="L_over_W",ylim=c(0,9))
abline(h=c( 4.055946277,4.829795325), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=14, col="slateblue", main="LV_over_L",ylim=c(0,6))
abline(h=c(4.031208111, 4.94928622), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=15, col="slateblue", main="Style_Length")
abline(h=c(4.106584701,4.715630803), col=c("blue", "red"), lty=c(3,3))


plot(out, map, lodcolumn=18, col="slateblue", main="Pollen.count_Stephen", ylim=c(0,5))
abline(h=c(4.115537144,4.910052868), col=c("blue", "red"), lty=c(3,3))

plot(out, map, lodcolumn=19, col="slateblue", main="Pollen_Count_Chang_Lab", ylim=c(0,6))
abline(h=c(4.400759479,5.26497719), col=c("blue", "red"), lty=c(3,3))


plot(out, map, lodcolumn=20, col="slateblue", main="Pollen_diameter_Chang_lab", ylim=c(0,6))
abline(h=c(4.104771254,5.315614386), col=c("blue", "red"), lty=c(3,3))



#Transformed traits##########################################

colnames(ipomoea$pheno)



#	flowers_per_day 	 skewed and stepwise. Log transforming and then splitting into two variables 	 Greater than / equal to 1 flower, and then a quantitative variable for all other #s of flowers. 
#		Flowers_buds_peduncle 	 skewed and stepwise. Log transforming and then splitting into two variables 	 Greater than / equal to 1 flower, and then a quantitative variable for all other #s of flowers. 		
#		Nectar.volume.uL 	 skewed, log +1 transformed		

#(31 Nectar_volume.uL    )
#32 LN_Nectar.volume.uL    
#3 neg_recip_L1_adj  #  
#4 neg_recip_L2_adj  #  
#5 neg_recip_L3_adj  #  
colnames(ipomoea$pheno)

out_transformed <- scan1(pr, ipomoea$pheno[,c(4:6, 31:32)])

QTLs_transformed<-find_peaks(out_transformed, map, threshold=2.4, drop = 1.5, expand2markers = T)
QTLs_transformed<-find_peaks(out_transformed, map, threshold=2, drop = 1.5, expand2markers = T)
?find_peaks
str(out)
?scan1
?find_marker
find_marker()
str(QTLs_normal)
QTLs_transformed$pos_marker<-find_marker(map, QTLs_transformed$chr, QTLs_transformed$pos) #Expanded to nearest marker
QTLs_transformed$pos_no_pseudo_marker<-find_marker(ipomoea$gmap, QTLs_transformed$chr, QTLs_transformed$pos) #using ipomoea$gmap rather than map allows me to find the marker rather htan the pseudomarker
QTLs_transformed$pos_ci_lo<-find_marker(map, QTLs_transformed$chr, QTLs_transformed$ci_lo) #Expanded to nearest marker 
QTLs_transformed$pos_ci_hi<-find_marker(map, QTLs_transformed$chr, QTLs_transformed$ci_hi) #Expanded to nearest marker
#QTLs_normal$pos_marker<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$pos) #using 
#QTLs_transformed$pos_marker<-find_marker(ipomoea$gmap, QTLs_transformed$chr, QTLs_transformed$pos) #using ipomoea$gmap rather than map allows me to find the marker rather htan the pseudomarker
#QTLs_transformed$pos_ci_lo<-find_marker(ipomoea$gmap, QTLs_transformed$chr, QTLs_transformed$ci_lo)
#QTLs_transformed$pos_ci_hi<-find_marker(ipomoea$gmap, QTLs_transformed$chr, QTLs_transformed$ci_hi)
par(mfrow=c(2,2))
QTLs_transformed[,12:24]<-NA
for (row in 1:nrow(QTLs_transformed)) {
  chr_var<-as.numeric(QTLs_transformed$chr[row])
  pos_var<-QTLs_transformed$pos_marker[row]
  pheno_var<-QTLs_transformed$lodcolumn[row]
  pos_cm<-QTLs_transformed$pos[row]
  pheno_var<-QTLs_transformed$lodcolumn[row]
  print(row)
  print(chr_var)
  print(class(chr_var))
  print(pos_var)
  print(class(pos_var))
  print(pheno_var)
  print(class(pheno_var))
  mu_a_d_effects<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef
  allele_effects<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var])$coef
  mu_a_d_effects_no_zero<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var], contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)),zerosum=FALSE)$coef
  allele_effects_no_zero<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var],zerosum=FALSE)$coef
  print(mu_a_d_effects)
  print(allele_effects)
  QTLs_transformed[row,12:14]<-mu_a_d_effects
  QTLs_transformed[row,15:18]<-allele_effects
  QTLs_transformed[row,19:21]<-mu_a_d_effects_no_zero
  QTLs_transformed[row,22:24]<-allele_effects_no_zero
  g <- maxmarg(pr, map, chr=chr_var, pos=pos_cm, return_char=TRUE)
  plot_pxg(g, ipomoea$pheno[,paste(pheno_var)], ylab=pheno_var, main=c(chr_var, pos_cm))
}
warnings()
colnames(QTLs_transformed)
colnames(QTLs_transformed)[12:24]<-c("mu_zero","a_zero","d_zero","11_zero","12_zero","22_zero","intercept_zero","mu_non_zero","a_non_zero","d_non_zero","11_non_zero","12_non_zero","22_non_zero")

#write.csv(QTLs_transformed, "QTL_effects_transformed_6_10_20.csv")
#write.csv(QTLs_transformed, "QTL_effects_transformed_6_8_20.csv")

#Plot QTLs for transformed traits#########################################################

colnames(out_transformed)


plot(out_transformed, map, lodcolumn=1, col="slateblue", main="neg_recip_L1_adj")
abline(h=c(4.120504444,4.999967039), col=c("blue", "red"), lty=c(3,3))

plot(out_transformed, map, lodcolumn=2, col="slateblue", main="neg_recip_L2_adj")
abline(h=c(4.197412924,4.957310873), col=c("blue", "red"), lty=c(3,3))

plot(out_transformed, map, lodcolumn=3, col="slateblue", main="neg_recip_L3_adj")
abline(h=c(4.121141975, 4.818345351), col=c("blue", "red"), lty=c(3,3))


plot(out_transformed, map, lodcolumn=4, col="slateblue", main="Nectar_volume_uL", ylim=c(0,9))
abline(h=c(4.146476854,4.838683699), col=c("blue", "red"), lty=c(3,3))

plot(out_transformed, map, lodcolumn=5, col="slateblue", main="LN_Nectar.volume.uL", ylim=c(0,8))
abline(h=c(4.197040973,4.810052233), col=c("blue", "red"), lty=c(3,3))




#Binary traits#########################################################
#18 TwoOrMore  #  
#19 ThreeOrMore  #  
#13 Ever_flowered   binary  
#20 Color   binary  
#37 SFS   binary  
#14 more_than_one_flower_per_day  #  
#16 more_than_one_flower_on_cyme  #  
colnames(ipomoea$pheno)
out_binary <- scan1(pr, ipomoea$pheno[,c(18,19,13,14,16)], model="binary")

out_binary <- scan1(pr, ipomoea$pheno[,c(20)], model="binary", eta_max=22) #color
out_binary <- scan1(pr, ipomoea$pheno[,c(37)], model="binary", eta_max=26) #sterility
#Color converges up to 22
#Sterility converges up to 26 

plot(out_binary, map, lodcolumn=1, col="slateblue", main="color")
plot(out_binary, map, lodcolumn=1, col="slateblue", main="sterility")
plot(out_binary, map, lodcolumn=1, col="slateblue", main="sterility", chr=10)
plot(out_binary, map, lodcolumn=1, col="slateblue", main="sterility", chr=14)
?scan1coef
c2eff <- scan1coef(pr[,"6"], ipomoea$pheno[,c(20)])
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["6"], columns=1:3, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
c2effB <- scan1coef(pr[,"6"], ipomoea$pheno[,c(20)],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
plot(c2effB, map["6"], columns=2:3, col=col)
last_coef <- unclass(c2effB)[nrow(c2effB),2:3] # last two coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])

g <- maxmarg(pr, map, chr=6, pos=56.84, return_char=TRUE)
plot_pxg(g, ipomoea$pheno[,c(20)], ylab="Color phenotype")

g <- maxmarg(pr, map, chr=7, pos=176, return_char=TRUE)
plot_pxg(g, ipomoea$pheno[,c(20)], ylab="Color phenotype")

g <- maxmarg(pr, map, chr=10, pos=50, return_char=TRUE)
plot_pxg(g, ipomoea$pheno[,c(37)], ylab="sterility phenotype")

g <- maxmarg(pr, map, chr=14, pos=73, return_char=TRUE)
plot_pxg(g, ipomoea$pheno[,c(37)], ylab="sterility phenotype")



colnames(ipomoea$pheno)

QTLs_binary<-find_peaks(out_binary, map, threshold=2.4, drop = 1.5, expand2markers = T) #Threshhold set based on lowest possible chromosome-specific value
#QTLs_binary<-find_peaks(out_binary, map, threshold=2, drop = 1.5, expand2markers = T)

QTLs_binary$pos_marker<-find_marker(map, QTLs_binary$chr, QTLs_binary$pos) #Expanded to nearest marker
QTLs_binary$pos_no_pseudo_marker<-find_marker(map, QTLs_binary$chr, QTLs_binary$pos) #Expanded to nearest marker
QTLs_binary$pos_ci_lo<-find_marker(map, QTLs_binary$chr, QTLs_binary$ci_lo) #Expanded to nearest marker 
QTLs_binary$pos_ci_hi<-find_marker(map, QTLs_binary$chr, QTLs_binary$ci_hi) #Expanded to nearest marker

#QTLs_binary$pos_ci_lo<-find_marker(map, QTLs_binary$chr, QTLs_binary$ci_lo) #Expanded to #nearest marker 
#QTLs_binary$pos_ci_hi<-find_marker(map, QTLs_binary$chr, QTLs_binary$ci_hi) #Expanded to nearest marker
#QTLs_normal$pos_marker<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$pos) #using 
#QTLs_normal$pos_marker<-find_marker(ipomoea$gmap, QTLs_normal$chr, QTLs_normal$pos) #using 
#QTLs_binary$pos_marker<-find_marker(ipomoea$gmap, QTLs_binary$chr, QTLs_binary$pos) #using ipomoea$gmap rather than map allows me to find the marker rather htan the pseudomarker
#QTLs_binary$pos_ci_lo<-find_marker(ipomoea$gmap, QTLs_binary$chr, QTLs_binary$ci_lo)
#QTLs_binary$pos_ci_hi<-find_marker(ipomoea$gmap, QTLs_binary$chr, QTLs_binary$ci_hi)
?fit1
#rm(QTLs_binary)
par(mfrow=c(2,2))
QTLs_binary[,12:24]<-NA
for (row in 1:nrow(QTLs_binary)) {
  chr_var<-as.numeric(QTLs_binary$chr[row])
  pos_var<-QTLs_binary$pos_marker[row]
  pos_cm<-QTLs_binary$pos[row]
  pheno_var<-QTLs_binary$lodcolumn[row]
  #pheno_var<-"Sterility"
  #pheno_var<-"Color"
  print(row)
  print(chr_var)
  print(class(chr_var))
  print(pos_var)
  print(class(pos_var))
  print(pheno_var)
  print(class(pheno_var))
  mu_a_d_effects<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var],model="binary", contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))$coef
  allele_effects<-fit1(pr[[chr_var]][,,pos_var], model="binary", ipomoea$pheno[,pheno_var])$coef
  mu_a_d_effects_no_zero<-fit1(pr[[chr_var]][,,pos_var], ipomoea$pheno[,pheno_var],model="binary", contrasts = cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)),zerosum=FALSE)$coef
  allele_effects_no_zero<-fit1(pr[[chr_var]][,,pos_var], model="binary", ipomoea$pheno[,pheno_var],zerosum=FALSE)$coef
  print(mu_a_d_effects)
  print(allele_effects)
  QTLs_binary[row,12:14]<-mu_a_d_effects
  QTLs_binary[row,15:18]<-allele_effects
  QTLs_binary[row,19:21]<-mu_a_d_effects_no_zero
  QTLs_binary[row,22:24]<-allele_effects_no_zero
  g <- maxmarg(pr, map, chr=chr_var, pos=pos_cm, return_char=TRUE)
  plot_pxg(g, ipomoea$pheno[,paste(pheno_var)], ylab=pheno_var, main=c(chr_var, pos_cm))
}

colnames(QTLs_binary)
colnames(QTLs_binary)[12:24]<-c("mu_zero","a_zero","d_zero","11_zero","12_zero","22_zero","intercept_zero","mu_non_zero","a_non_zero","d_non_zero","11_non_zero","12_non_zero","22_non_zero")
#write.csv(QTLs_binary, "QTL_effects_binary_6_10_2020.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_6_9_2020.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_color_6_9_2020.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_sterility_6_9_2020.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_8_19_no_pseudomarkers.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_8_21.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_8_21_no_pseudomarkers.csv")
#write.csv(QTLs_binary, "Low_LOD_QTL_effects_binary_8_21.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_all_8_26.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_all_9_11.csv")
#write.csv(QTLs_binary, "QTL_effects_binary_all_10_4.csv")
getwd()
#Plot QTLs for binary traits#########################################################

colnames(out_binary)
plot
warnings()
plot(out_binary, map, lodcolumn=1, col="slateblue", main="TwoOrMore", ylim=c(0,6))
abline(h=c(4.037339185,4.89345557), col=c("blue", "red"), lty=c(3,3))

plot(out_binary, map, lodcolumn=2, col="slateblue",main="ThreeOrMore", ylim=c(0,5))
abline(h=c(4.090567039, 4.807677067), col=c("blue", "red"), lty=c(3,3))

plot(out_binary, map, lodcolumn=3, col="slateblue",main="EverFlowered", ylim=c(0,5))
abline(h=c(4.15719203, 4.840112207), col=c("blue", "red"), lty=c(3,3))

plot(out_binary, map, lodcolumn=1, col="slateblue",main="Color")
abline(h=c(4.1919993371545,4.83771145819556), col=c("blue", "red"), lty=c(3,3))
plot(out_binary, map, lodcolumn=4, col="slateblue",main="Color", ylim=c(0,10))

plot(out_binary, map, lodcolumn=1, col="slateblue",main="SFS")
abline(h=c(3.98026271621382,4.50090755414526), col=c("blue", "red"), lty=c(3,3))

plot(out_binary, map, lodcolumn=4, col="slateblue",main="MoreThanOneFlowerPerDay", ylim=c(0,6))
abline(h=c(4.082564563,4.798497777), col=c("blue", "red"), lty=c(3,3))

plot(out_binary, map, lodcolumn=5, col="slateblue",main="MoreThanOneFlowerPerCyme", ylim=c(0,8))
abline(h=c(4.071676355, 4.921293393), col=c("blue", "red"), lty=c(3,3))
     
3.62531891,4.329288607
3.66131693,4.478948444

#Plot map #################
ipomoea1<-read.cross("csvsr", "E:/Users/Joanna/Dropbox/Professional/Duke/Rausher_lab/ChapterThreeMapping_Cross_Experiment/Analysis/Reanalysis2018_2019/2019/QTL2/QTL1_input_files",  "ipomoea_gen_rot.csv","ipomoea_phe_rot.csv",  genotypes = c("11", "12", "22", "1-", "2-"))


#normqtls<-read.csv("QTL_effects_normal_6_10-20.csv")
getwd()
plot_qtls<-read.csv("All_QTLs_for_plotting_final.csv") #Including all peaks, including ns
plot_GWS_qtls<-subset(plot_qtls, plot_qtls$lod>3.1)

#genome wide only = LOD greater than 3.119986721
#all = LOD greater than 2.707544048

par(xpd=TRUE)
par(mar=rep(5,4))
par(mar=c(2.5,5,2.5,6))
segmentsOnMap(ipomoea1, phe = plot_qtls$lod_abbrev, chr=plot_qtls$chr, l = plot_qtls$ci_lo, h=plot_qtls$ci_hi, palette=rainbow, legendPosition = "topright", leg.inset = c(-.05,0.01), legendCex=.5)
#View(plot_qtls)

par(xpd=TRUE)
par(mar=rep(5,4))
par(mar=c(2.5,5,2.5,6))
segmentsOnMap(ipomoea1, phe = plot_GWS_qtls$lod_abbrev, chr=plot_GWS_qtls$chr, l = plot_GWS_qtls$ci_lo, h=plot_GWS_qtls$ci_hi, palette=rainbow, legendPosition = "topright", leg.inset = c(-0.05,0.01), legendCex=.5)
#View(plot_qtls)

