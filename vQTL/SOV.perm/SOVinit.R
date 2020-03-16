library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

setwd("/work/04902/azg5169/stampede2/vQTL/vQTL/SOV.perm")
vQTLsub <- read.cross(file = "../ManchingStressData_Covar.csv" )

vQTLsub <- drop.nullmarkers(vQTLsub)
vQTLsub <- calc.genoprob(vQTLsub)


##INTERACTIVE
#with ManchingStressData_covar.csv
colnames(vQTLsub$pheno)[1] = "Height"

intRealOneVar <- scanonevar(cross = vQTLsub, mean.formula = Height ~ Env*(mean.QTL.add + mean.QTL.dom), var.formula = ~ Env*(var.QTL.add + var.QTL.dom), return.covar.effects = TRUE)

write_rds(intRealOneVar, "InteractiveResult.rds")
