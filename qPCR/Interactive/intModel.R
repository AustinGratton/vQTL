library(gridExtra)
library(qtl)
library(vqtl)
library(ggplot2)

test_hyb <-read.cross(file = "qPCR/FullHybvqtlinput.csv")

test_inbr <-read.cross(file = "qPCR/FullInbvqtlinput.csv")

test_full <- read.cross(file = "qPCR/Fullvqtlinput.csv")

test_full <- calc.genoprob(test_full, error.prob = .001)
full_so <- scanone(cross = test_full, pheno.col = 'stress')
full_sov <- scanonevar(cross = test_full, 
                     mean.formula = stress ~ BreedType + mean.QTL.add + mean.QTL.dom,
                     var.formula = ~ BreedType + var.QTL.add + var.QTL.dom,
                     return.covar.effects = TRUE)

SOV.perm = full_sov

SOV.perm$result$loc.name[SOV.perm$result$mQTL.asymp.p <= .0001]
SOV.perm$result$loc.name[SOV.perm$result$vQTL.asymp.p <= .0001]
SOV.perm$result$loc.name[SOV.perm$result$mvQTL.asymp.p <= .0001]


identical(test_hyb$pheno$Genotype, test_inbr$pheno$Genotype)


in_p1 <- scanone(cross = test_full, pheno.col = 'stress')
inv_p1 <- scanonevar(cross = test_full,
                     mean.formula = stress ~ mean.QTL.add,
                     var.formula = ~ var.QTL.add)
