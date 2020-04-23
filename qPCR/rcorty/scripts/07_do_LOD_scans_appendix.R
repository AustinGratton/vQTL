library(vqtl)
library(ggplot2)

# do scans in units of LOD scores
so_p1x <- scanone(cross = test_cross, pheno.col = 'phenotype1x', addcovar = test_cross$pheno$sex)
sov_p1x <- scanonevar(cross = test_cross,
											mean.formula = phenotype1x ~ sex + mean.QTL.add + mean.QTL.dom,
											var.formula = ~sex + var.QTL.add + var.QTL.dom)


so_p2x <- scanone(cross = test_cross, pheno.col = 'phenotype2x', addcovar = test_cross$pheno$sex)
sov_p2x <- scanonevar(cross = test_cross,
											mean.formula = phenotype2x ~ sex + mean.QTL.add + mean.QTL.dom,
											var.formula = ~ sex + var.QTL.add + var.QTL.dom)


so_p3x <- scanone(cross = test_cross, pheno.col = 'phenotype3x', addcovar = test_cross$pheno$sex)
sov_p3x <- scanonevar(cross = test_cross,
											mean.formula = phenotype3x ~ sex + mean.QTL.add + mean.QTL.dom,
											var.formula = ~sex + var.QTL.add + var.QTL.dom)


so_p4x <- scanone(cross = test_cross, pheno.col = 'phenotype4x', addcovar = test_cross$pheno$sex)
sov_p4x <- scanonevar(cross = test_cross,
											mean.formula = phenotype4x ~ sex + mean.QTL.add + mean.QTL.dom,
											var.formula = ~sex + var.QTL.add + var.QTL.dom)


saveRDS(object = list(so_p1x = so_p1x, sov_p1x = sov_p1x,
											so_p2x = so_p2x, sov_p2x = sov_p2x,
											so_p3x = so_p3x, sov_p3x = sov_p3x,
											so_p3x = so_p4x, sov_p4x = sov_p4x),
				file = 'saves/LOD_scans_appendix.RDS')



# make scans for fig 1 -- LOD score scans
ymax <- 6

plot(p1x_LOD_scan <- plot(x = sov_p1x, y = so_p1x, ymax = ymax))
plot(p2x_LOD_scan <- plot(x = sov_p2x, y = so_p2x, ymax = ymax))
plot(p3x_LOD_scan <- plot(x = sov_p3x, y = so_p3x, ymax = ymax))
plot(p4x_LOD_scan <- plot(x = sov_p4x, y = so_p4x, ymax = ymax))


ggsave(plot = p1x_LOD_scan, filename = 'images/LOD_scan_phen1x.pdf', height = 2.5, width = 9)
ggsave(plot = p2x_LOD_scan, filename = 'images/LOD_scan_phen2x.pdf', height = 2.5, width = 9)
ggsave(plot = p3x_LOD_scan, filename = 'images/LOD_scan_phen3x.pdf', height = 2.5, width = 9)
ggsave(plot = p4x_LOD_scan, filename = 'images/LOD_scan_phen4x.pdf', height = 2.5, width = 9)
