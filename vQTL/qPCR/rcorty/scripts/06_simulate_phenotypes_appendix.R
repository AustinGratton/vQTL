library(qtl)
library(vqtl)

set.seed(27599)

test_cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(100, 3), n.mar = 11, eq.spacing = TRUE, include.x = FALSE, anchor.tel = TRUE),
														 n.ind = 400,
														 type = 'f2')
test_cross$pheno$sex <- rep(x = c(0, 1), each = 200)
test_cross <- qtl::calc.genoprob(cross = test_cross, step = 2)

test_cross$pheno$phenotype1x <- rnorm(n = qtl::nind(test_cross))

test_cross$pheno$phenotype2x <- rnorm(n = qtl::nind(test_cross),
																			mean = 0.26*(test_cross$geno$`1`$data[,6] - 2),
																			sd = exp(test_cross$pheno$sex - 0.5))

test_cross$pheno$phenotype3x <- rnorm(n = qtl::nind(test_cross),
																			sd = exp(0.21*(test_cross$geno$`2`$data[,6] - 2) + (test_cross$pheno$sex - 0.5)))

test_cross$pheno$phenotype4x <- rnorm(n = qtl::nind(test_cross),
																			mean = 0.2*(test_cross$geno$`3`$data[,6] - 2),
																			sd = exp(0.15*(test_cross$geno$`3`$data[,6] - 2) + (test_cross$pheno$sex - 0.5)))

saveRDS(object = test_cross,
				file = 'saves/test_cross_appendix.RDS')
