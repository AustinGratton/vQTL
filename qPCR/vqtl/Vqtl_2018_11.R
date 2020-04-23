
### qtl anlysis on 2018_11 plate data ###
library("qtl")
library("vqtl")
library("dplyr")
library("stringr")
#library(beepr)
setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
#now with the Stress_11 dataset
Stress_11 <-read.cross(file = "./2018_11/vqtlinput_11.csv")
#####small sample of snpHeight
#####Stress_11 = read.cross(file = "sampSnpHeight.csv")
Stress_11 <- drop.nullmarkers(Stress_11)

##### CORTY code #####
Stress_11$pheno$BreedType = factor(Stress_11$pheno$BreedType)

#Stress_11$pheno$DumBreed = ifelse(Stress_11$pheno$BreedType=="Inbred", 0, 1)
# ggplot of Stress_11 type and height
# library(tidyverse)
# Stress_11$pheno %>%
#   ggplot(mapping = aes(x = BreedType, y = stress)) +
#   geom_jitter(width = 0.2)

Stress_11 <- calc.genoprob(Stress_11)

###### additive model #####
# outv <- scanonevar(cross = Stress_11,
#                    mean.formula = stress ~ BreedType + mean.QTL.add + mean.QTL.dom,
#                    var.formula = ~ BreedType + var.QTL.add + var.QTL.dom,
#                    return.covar.effects = TRUE)
# #beep()
# 
# #outv$result %>% glimpse()
# 
# write.csv(outv$result, file = "AdditiveModelStress_11_Output.csv")


##### interactive model ####

# eliminate all '(' in the ID in the stress files #
outv <- scanonevar(cross = Stress_11,
                   mean.formula = stress ~ BreedType * (mean.QTL.add + mean.QTL.dom),
                   var.formula = ~ BreedType * (var.QTL.add + var.QTL.dom),
                   return.covar.effects = TRUE)

outv$result %>% glimpse()
write.csv(outv$result, file = "InteractiveModelStress_11_Output.csv")
