library(gridExtra)
library(qtl)
library(vqtl)
library(ggplot2)

test_hyb <-read.cross(file = "../qPCR/FullHybvqtlinput.csv")

test_inbr <-read.cross(file = "../qPCR/FullInbvqtlinput.csv")

df = read.csv(file = "../qPCR/FullInbvqtlinput.csv")
#calc.genoprob(test_hyb)
#?calc.genoprob
#head(as.numeric(as.factor(test_cross$pheno$Barcode)))
# do scans in units of LOD scores

test_hyb <- calc.genoprob(test_hyb, error.prob = .001)
hy_p1 <- scanone(cross = test_hyb, pheno.col = 'stress')
hyv_p1 <- scanonevar(cross = test_hyb, 
                     mean.formula = stress ~ mean.QTL.add, 
                     var.formula = ~ var.QTL.add, 
                     return.covar.effects = TRUE)



test_inbr <- calc.genoprob(test_inbr, error.prob = .001)
in_p1 <- scanone(cross = test_inbr, pheno.col = 'stress')
inv_p1 <- scanonevar(cross = test_inbr,
										 mean.formula = stress ~ mean.QTL.add,
										 var.formula = ~ var.QTL.add)


# so_p3 <- scanone(cross = test_cross, pheno.col = 'phenotype3', addcovar = test_cross$pheno$sex)
# sov_p3 <- scanonevar(cross = test_cross,
# 										 mean.formula = phenotype3 ~ sex + mean.QTL.add + mean.QTL.dom,
# 										 var.formula = ~sex + var.QTL.add + var.QTL.dom)
# 
# 
# so_p4 <- scanone(cross = test_cross, pheno.col = 'phenotype4', addcovar = test_cross$pheno$sex)
# sov_p4 <- scanonevar(cross = test_cross,
# 										 mean.formula = phenotype4 ~ sex + mean.QTL.add + mean.QTL.dom,
# 										 var.formula = ~sex + var.QTL.add + var.QTL.dom)
# 
# saveRDS(object = list(so_p1 = so_p1, sov_p1 = sov_p1,
# 											so_p2 = so_p2, sov_p2 = sov_p2,
# 											so_p3 = so_p3, sov_p3 = sov_p3,
# 											so_p3 = so_p4, sov_p4 = sov_p4),
# 				file = 'saves/LOD_scans.RDS')
# 


# make scans for fig 1 -- LOD score scans
ymax <- 6



chr1_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(1))
chr2_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(2))
chr3_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(3))
chr4_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(4))
chr5_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(5))
chr6_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(6))
chr7_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(7))
chr8_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(8))
chr9_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(9))
chr10_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax, chr = c(10))

table(as.character(test_hyb$geno$`1`$data))

test_hyb$geno$`1`$data[is.na(test_hyb$geno$`1`$data)] <- 2
test_hyb$geno$`2`$data[is.na(test_hyb$geno$`2`$data)] <- 2
test_hyb$geno$`3`$data[is.na(test_hyb$geno$`3`$data)] <- 2
test_hyb$geno$`4`$data[is.na(test_hyb$geno$`4`$data)] <- 2
test_hyb$geno$`5`$data[is.na(test_hyb$geno$`5`$data)] <- 2
test_hyb$geno$`6`$data[is.na(test_hyb$geno$`6`$data)] <- 2
test_hyb$geno$`7`$data[is.na(test_hyb$geno$`7`$data)] <- 2
test_hyb$geno$`8`$data[is.na(test_hyb$geno$`8`$data)] <- 2
test_hyb$geno$`9`$data[is.na(test_hyb$geno$`9`$data)] <- 2
test_hyb$geno$`10`$data[is.na(test_hyb$geno$`10`$data)] <- 2


mean_var_plot_model_based(test_hyb, 'stress', 'gpm27', genotype.names = c('A','-','B'))

hyv_p1$result$loc.name[hyv_p1$result$mvQTL.asymp.p <= .05]

table(hyv_p1$result$mQTL.asymp.p <= .05)
table(hyv_p1$result$mvQTL.asymp.p <= .05)
table(1-hyv_p1$result$vQTL.asymp.p <= .05)

hyv_p1$result$loc.name[hyv_p1$result$mQTL.asymp.p <= .05]


table(test_inbr$geno$`1`$data)
test_inbr$geno$`1`$data[is.na(test_inbr$geno$`1`$data)] <- NA
test_inbr$geno$`2`$data[is.na(test_inbr$geno$`2`$data)] <- 2
test_inbr$geno$`3`$data[is.na(test_inbr$geno$`3`$data)] <- 2
test_inbr$geno$`4`$data[is.na(test_inbr$geno$`4`$data)] <- 2
test_inbr$geno$`5`$data[is.na(test_inbr$geno$`5`$data)] <- 2
test_inbr$geno$`6`$data[is.na(test_inbr$geno$`6`$data)] <- 2
test_inbr$geno$`7`$data[is.na(test_inbr$geno$`7`$data)] <- 2
test_inbr$geno$`8`$data[is.na(test_inbr$geno$`8`$data)] <- 2
test_inbr$geno$`9`$data[is.na(test_inbr$geno$`9`$data)] <- 2
test_inbr$geno$`10`$data[is.na(test_inbr$geno$`10`$data)] <- 2

table(inv_p1$result$mQTL.asymp.p <= .05)
table(inv_p1$result$mvQTL.asymp.p <= .05)
table(inv_p1$result$vQTL.asymp.p <= .05)

sig.loc = inv_p1$result$loc.name[inv_p1$result$vQTL.asymp.p <= .05]


sig.plot1 = mean_var_plot_model_based(test_hyb, 'stress', sig.loc[1], genotype.names = c('A','B'))
sig.plot2 = mean_var_plot_model_based(test_hyb, 'stress', sig.loc[2], genotype.names = c('A','B'))
sig.plot3 = mean_var_plot_model_based(test_hyb, 'stress', sig.loc[3], genotype.names = c('A','B'))
sig.plot4 = mean_var_plot_model_based(test_hyb, 'stress', 'gpm603', genotype.names = c('A','B'))
sig.plot5 = mean_var_plot_model_based(test_hyb, 'stress', sig.loc[5], genotype.names = c('A','B'))

require(gridExtra)
grid.arrange(sig.plot1, sig.plot2, sig.plot3, sig.plot4, sig.plot5, ncol = 2)

test_cross = readRDS("../qPCR/rcorty/saves/test_cross2.RDS")
mean_var_plot_model_based(test_cross, 'phenotype4', 'D3M6')

require(gridExtra)
grid.arrange(chr1_LOD_scan, chr2_LOD_scan, chr3_LOD_scan, chr4_LOD_scan, chr5_LOD_scan, chr6_LOD_scan, chr7_LOD_scan, chr8_LOD_scan, chr9_LOD_scan, chr10_LOD_scan, ncol=3)


plot(p2_LOD_scan <- plot(x = sov_p2, y = so_p2, ymax = ymax))



# plot(p3_LOD_scan <- plot(x = sov_p3, y = so_p3, ymax = ymax))
# plot(p4_LOD_scan <- plot(x = sov_p4, y = so_p4, ymax = ymax))
# 
# 
# 
# 
# 
# ggsave(plot = p1_LOD_scan, filename = 'images/LOD_scan_phen1.pdf', height = 2.5, width = 9)
# ggsave(plot = p2_LOD_scan, filename = 'images/LOD_scan_phen2.pdf', height = 2.5, width = 9)
# ggsave(plot = p3_LOD_scan, filename = 'images/LOD_scan_phen3.pdf', height = 2.5, width = 9)
# ggsave(plot = p4_LOD_scan, filename = 'images/LOD_scan_phen4.pdf', height = 2.5, width = 9)
