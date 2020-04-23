library(vqtl)
library(ggplot2)
library(evd)

# do scans in units of empirical p value

# make fig 2 -- empirircal p-value scans
num_perms <- 1000
num_cores <- 40
seed <- 27599

library(evd)

PGEV <- function(q, gev, ...) {

	if ((!'estimate' %in% names(gev))) {
		stop("argument 'gev' must have an element named 'estimate' (all gev objects do)")
	}
	if (length(gev$estimate) != 3) {
		stop("gev$estimate must have three elements")
	}

	return(pgev(q = q,
							loc = gev$estimate[1],
							scale = gev$estimate[2],
							shape = gev$estimate[3], ...))
}

sop_p1 <-  scanone(cross = test_cross, pheno.col = 'phenotype1', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p1 <- evd::fgev(x = sop_p1)
so_p1$empir.p <- PGEV(q = so_p1$lod, gev = evd_p1, lower.tail = FALSE)

sov_p1 <- scanonevar.perm(sov = sov_p1, n.perms = num_perms, random.seed = seed, n.cores = num_cores)


sop_p2 <-  scanone(cross = test_cross, pheno.col = 'phenotype2', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p2 <- evd::fgev(x = sop_p2)
so_p2$empir.p <- PGEV(q = so_p2$lod, gev = evd_p2, lower.tail = FALSE)

sov_p2 <- scanonevar.perm(sov = sov_p2, n.perms = num_perms, random.seed = seed, n.cores = num_cores)

sop_p3 <-  scanone(cross = test_cross, pheno.col = 'phenotype3', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p3 <- evd::fgev(x = sop_p3)
so_p3$empir.p <- PGEV(q = so_p3$lod, gev = evd_p3, lower.tail = FALSE)

sov_p3 <- scanonevar.perm(sov = sov_p3, n.perms = num_perms, random.seed = seed, n.cores = num_cores)


sop_p4 <-  scanone(cross = test_cross, pheno.col = 'phenotype4', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p4 <- evd::fgev(x = sop_p4)
so_p4$empir.p <- PGEV(q = so_p4$lod, gev = evd_p4, lower.tail = FALSE)

sov_p4 <- scanonevar.perm(sov = sov_p4, n.perms = num_perms, random.seed = seed, n.cores = num_cores)


saveRDS(object = list(so_p1, sov_p1,
											so_p2, sov_p2,
											so_p3, sov_p3,
											so_p4, sov_p4),
				file = 'saves/empir_p_scans.RDS')



# make plots for fig 2 -- empirircal p-value scans
ymax <- 3.5

plot(p1_empir_p_scan <- plot(x = sov_p1, y = so_p1, ymax = ymax))
plot(p2_empir_p_scan <- plot(x = sov_p2, y = so_p2, ymax = ymax))
plot(p3_empir_p_scan <- plot(x = sov_p3, y = so_p3, ymax = ymax))
plot(p4_empir_p_scan <- plot(x = sov_p4, y = so_p4, ymax = ymax))


ggplot2::ggsave(plot = p1_empir_p_scan, filename = 'images/empir_p_scan_phen1.pdf', height = 2.5, width = 9)
ggplot2::ggsave(plot = p2_empir_p_scan, filename = 'images/empir_p_scan_phen2.pdf', height = 2.5, width = 9)
ggplot2::ggsave(plot = p3_empir_p_scan, filename = 'images/empir_p_scan_phen3.pdf', height = 2.5, width = 9)
ggplot2::ggsave(plot = p4_empir_p_scan, filename = 'images/empir_p_scan_phen4.pdf', height = 2.5, width = 9)
