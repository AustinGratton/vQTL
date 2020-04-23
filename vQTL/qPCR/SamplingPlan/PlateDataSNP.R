### Sample plan data cleaning ###
library(stringr)
library(tidyverse)
library(dplyr)

setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/SamplingPlan")

#############################################
##### formatting the Sampling Plan data #####
#############################################
dat = read.csv(file = "../2016_Clayton/Field Book (2016) - Clayton - Sampling Plan_TIDIED.csv", header = TRUE)


#filter out the recongizable MO###
dat = dat %>% filter(str_detect(dat$Genotype, "Mo") == TRUE)
#create the breedtype category
BreedType = ifelse(substr(dat$Genotype, 1,1)=="M", "Inbred", "Hybrid")
dat = cbind(BreedType, dat)


snpFull = read.csv(file = "../IBM94markerset08seq.csv", header = TRUE)
snp = snpFull[,-(1:5)]

#############################################
### matching the Mo### to their snp Values ##
#############################################
library(data.table)
#function takes the last n characters of a string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
#creating datasets exclusively containing the last three numbers 
dat.num =cbind(dat, as.integer(substrRight(as.character(dat$Genotype), 3)))
#produces warning sign about NA's for data that is not availible for invalid Mo###
snp.num = as.integer(substrRight(as.character(colnames(snp)), 3))

snp.new = data.frame(lapply(snp,as.character),stringsAsFactors=FALSE)
snpMatch = rbind(snp.num,snp.new)
snpMatch = transpose(snpMatch)
# creating matching Genotype columns to merge the data
#need = c(1,5,7,8,9,10,16) #column numbers of needed variables
#dat.num = dat.num[,need]
colnames(snpMatch) = "GenotypeNum"
#datnames =names(dat.num)
colnames(dat.num)[17] = "GenotypeNum"
dat2 = merge(dat.num,snpMatch, by.x = "GenotypeNum", by.y = "GenotypeNum",all.MoNum = all)
dat2 = dat2[order(dat2$Genotype, decreasing = FALSE),]
dat2 = dat2[,-c(1, 4)] #omit GenotypeNum
dat2[1:20, 1:20]
#write.csv(dat2, "SamplingPlan_dat2.csv")

########################################
######### Including plate data #########
########################################

colnames(dat2)[4] = "sampleID.exp"

# ### by month framing ####
# 
# # 2018_11 plate data #
# ### the file in plate will come from the qPCR output including stress ratios ###
# plate_11 = read.csv(file = "../2018_11/2018_11_withStress.csv", header = TRUE)
# colnames(plate_11)[2] = "sampleID"
# full_11 = merge(plate_11, dat2, by = "sampleID")
# 
# # 2018_6 plate data #
# plate_6 = read.csv(file = , header = TRUE)
# colnames(plate_6)[2] = "sampleID"
# full_6 = merge(plate_6, dat2, by = "sampleID")
# 
# # 2018_8 plate data #
# plate_8 = read.csv(file = , header = TRUE)
# colnames(plate_8)[2] = "sampleID"
# full_8 = merge(plate_8, dat2, by = "sampleID")
# 
# #create one big dataframe containing all plate months
# #full = rbind(full_11, full_6, full_8)
# 
# #### Heirarchical Model Framing ####
HPlate = read.csv(file = "../Heirarchical/Hierarchical_exp_data_stress.csv")


#split isolate the June month since the sampleID is availible for it
June = HPlate %>% filter(str_detect(HPlate$month, "june") == TRUE)
colnames(June)[2] = "Barcode"
June = merge(June, dat2, by = "Barcode")

#August and November Months
AugNov = HPlate %>% filter(str_detect(HPlate$month, "aug")|str_detect(HPlate$month, "nov"))
AugNov = merge(AugNov, dat2, by = "sampleID.exp")

#making sure that non of the genotype information was matched to more than on obs
bad =intersect(June$Barcode, AugNov$Barcode)
length(bad)
baddoup = AugNov[match(bad, AugNov$Barcode),1:20]
repJune = June[c(21,32),] # we want to replace the repeated values in June with N/A
June[c(21,32), ] = NA

full = rbind(June, AugNov)
#moving stress to the first column
head(names(full), 20)
full = cbind(full$stress, full[,-12])
colnames(full)[1] = "stress"
#####Adding marker location and chromosome#####

addmarker <- function(full, plate){
  aux = matrix(snpFull$incre_new, nrow= 1)
  aux = rbind(snpFull$Chromosome, aux)
  zeros = dim(full)[2]-dim(aux)[2]
  fillnames = names(full)[1:zeros]
  other = as.data.frame(matrix(rep(0,2*zeros), nrow = 2)) #repeat the number of 0 as the number of variables
  aux = cbind(other,aux)
  colnames(aux) = rep("",dim(aux)[2])
  colnames(full) = rep("",dim(aux)[2])
  dat3 = rbind(aux,full)
  colnames(dat3) = c(fillnames, as.character(snpFull$markername))
  return(dat3)
  
}

#make sure to delete extra 0's and also delete '('
full = addmarker(full)
#unnecessary columns for vQTL input, we only need: Stress, breedtype, genotype, barcode
notneed = c(3:7,9:12,14:17,19:26)
full = full[,-notneed]
write.csv(full, file = "../Heirarchical/vqtlinput.csv" ,row.names = FALSE)


#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####
head(names(full), 30) #BreedType = index 13

#needed blank space for vqtl
blanksp = full[1:2,]

FullInb = full %>% filter(str_detect(full$BreedType, "Inbred") == TRUE)
FullHyb = full %>% filter(str_detect(full$BreedType, "Hybrid") == TRUE)
write.csv(rbind(blanksp,FullInb), file = "../Heirarchical/FullInbvqtlinput.csv" ,row.names = FALSE)
write.csv(rbind(blanksp,FullHyb), file = "../Heirarchical/FullHybvqtlinput.csv" ,row.names = FALSE)

