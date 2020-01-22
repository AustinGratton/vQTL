library(vqtl)
library(ggplot2)

test_hyb <-read.cross(file = "../qPCR/FullHybvqtlinput.csv")

test_inbr <-read.cross(file = "../qPCR/FullInbvqtlinput.csv")

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
so_p2 <- scanone(cross = test_inbr, pheno.col = 'stress')
sov_p2 <- scanonevar(cross = test_inbr,
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


plot(p1_LOD_scan <- plot(x = hyv_p1, y = hy_p1, ymax = ymax))
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
