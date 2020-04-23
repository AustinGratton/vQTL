library(qtl)


set.seed(27599)

test_cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(100, 3), n.mar = 11, eq.spacing = TRUE, include.x = FALSE, anchor.tel = TRUE),
														 n.ind = 400,
														 type = 'f2')
test_cross$pheno$sex <- rep(x = c(0, 1), each = 200)
test_cross <- qtl::calc.genoprob(cross = test_cross, step = 2)

test_cross$pheno$phenotype1 <- rnorm(n = qtl::nind(test_cross))

test_cross$pheno$phenotype2 <- rnorm(n = qtl::nind(test_cross),
																		 mean = 0.28*(test_cross$geno$`1`$data[,6] - 2))

test_cross$pheno$phenotype3 <- rnorm(n = qtl::nind(test_cross),
																		 sd = exp(0.23*(test_cross$geno$`2`$data[,6] - 2)))

test_cross$pheno$phenotype4 <- rnorm(n = qtl::nind(test_cross),
																		 mean = 0.24*(test_cross$geno$`3`$data[,6] - 2),
																		 sd = exp(0.16*(test_cross$geno$`3`$data[,6] - 2)))

saveRDS(object = test_cross,
				file = 'saves/test_cross.RDS')
