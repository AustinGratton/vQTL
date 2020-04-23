#install.packages("vqtl")
#install.packages("qtl")
library("qtl")
library("vqtl")
library("dplyr")
library("stringr")
#library(beepr)
setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/vqtl/")

stressprod = read.cross(file = "../Heirarchical/vqtlinput.csv")
stressprod <- drop.nullmarkers(stressprod)

##### CORTY code #####
stressprod$pheno$BreedType = factor(stressprod$pheno$BreedType)
stressprod$pheno$month = factor(stressprod$pheno$month)
# ggplot of stress type and height
# library(tidyverse)
# stressprod$pheno %>%
#   ggplot(mapping = aes(x = BreedType, y = stress)) +
#   geom_jitter(width = 0.2)

stressprod <- calc.genoprob(stressprod)


# outv <- scanonevar(cross = stressprod,
#                    mean.formula = stress ~ BreedType + mean.QTL.add + mean.QTL.dom,
#                    var.formula = ~ BreedType + var.QTL.add + var.QTL.dom,
#                    return.covar.effects = TRUE)
# #beep()
# 
# #outv$result %>% glimpse()
# 
# write.csv(outv$result, file = "AdditiveModelstress_Output.csv")

outv <- scanonevar(cross = stressprod,
                   mean.formula = stress ~ BreedType * (mean.QTL.add + mean.QTL.dom),
                   var.formula = ~ BreedType * (var.QTL.add + var.QTL.dom),
                   return.covar.effects = TRUE)

outv$result %>% glimpse()
write.csv(outv$result, file = "InteractiveModelstress_Output.csv")

### don't need .dom ###