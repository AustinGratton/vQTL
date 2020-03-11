library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

vQTLsub <- read.cross(file = "../ManchingStressData_Covar.csv" )

vQTLsub <- drop.nullmarkers(vQTLsub)
vQTLsub <- calc.genoprob(vQTLsub)


##INTERACTIVE
#with ManchingStressData_covar.csv
intRealOneVar <- scanonevar(cross = vQTLsub, 
                            mean.formula = ï..Height ~ Env*(mean.QTL.add + mean.QTL.dom),
                            var.formula = ~ Env*(var.QTL.add + var.QTL.dom),
                            return.covar.effects = TRUE)
table(intRealOneVar$result$mvQTL.asymp.p <= .05)

write_rds(intRealOneVar, "InteractiveResult.rds")
