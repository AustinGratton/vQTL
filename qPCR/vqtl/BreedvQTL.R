#install.packages("vqtl")
#install.packages("qtl")
library("qtl")
library("vqtl")
library("dplyr")
library("stringr")
#library(beepr)
setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
#now with the breed dataset
breed <-read.cross(file = "snpHeight.csv")
#####small sample of snpHeight
#####breed = read.cross(file = "sampSnpHeight.csv")
breed <- drop.nullmarkers(breed)

##### CORTY code #####
breed$pheno$BreedType = factor(breed$pheno$BreedType)

# ggplot of breed type and height
# library(tidyverse)
# breed$pheno %>%
#   ggplot(mapping = aes(x = BreedType, y = Height)) +
#   geom_jitter(width = 0.2)

breed <- calc.genoprob(breed)


# outv <- scanonevar(cross = breed,
#                    mean.formula = Height ~ BreedType + mean.QTL.add + mean.QTL.dom,
#                    var.formula = ~ BreedType + var.QTL.add + var.QTL.dom,
#                    return.covar.effects = TRUE)
# #beep()
# 
# #outv$result %>% glimpse()
# 
# write.csv(outv$result, file = "AdditiveModelBreed_Output.csv")

outv <- scanonevar(cross = breed,
                   mean.formula = Height ~ BreedType * (mean.QTL.add + mean.QTL.dom),
                   var.formula = ~ BreedType * (var.QTL.add + var.QTL.dom),
                   return.covar.effects = TRUE)

outv$result %>% glimpse()
write.csv(outv$result, file = "InteractiveModelBreed_Output.csv")

