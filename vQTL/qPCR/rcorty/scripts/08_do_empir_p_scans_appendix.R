library(vqtl)
library(ggplot2)
library(evd)

# do scans in units of empirical p value
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

sop_p1x <-  scanone(cross = test_cross, pheno.col = 'phenotype1x', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p1x <- evd::fgev(x = sop_p1x)
so_p1x$empir.p <- PGEV(q = so_p1x$lod, gev = evd_p1x, lower.tail = FALSE)

sov_p1x <- scanonevar.perm(sov = sov_p1x, n.perms = num_perms, random.seed = seed, n.cores = num_cores)


sop_p2x <-  scanone(cross = test_cross, pheno.col = 'phenotype2x', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p2x <- evd::fgev(x = sop_p2x)
so_p2x$empir.p <- PGEV(q = so_p2x$lod, gev = evd_p2x, lower.tail = FALSE)

sov_p2x <- scanonevar.perm(sov = sov_p2x, n.perms = num_perms, random.seed = seed, n.cores = num_cores)

sop_p3x <-  scanone(cross = test_cross, pheno.col = 'phenotype3x', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p3x <- evd::fgev(x = sop_p3x)
so_p3x$empir.p <- PGEV(q = so_p3x$lod, gev = evd_p3x, lower.tail = FALSE)

sov_p3x <- scanonevar.perm(sov = sov_p3x, n.perms = num_perms, random.seed = seed, n.cores = num_cores)


sop_p4x <-  scanone(cross = test_cross, pheno.col = 'phenotype4x', addcovar = test_cross$pheno$sex,
									 n.perm = num_perms, n.cluster = num_cores)
evd_p4x <- evd::fgev(x = sop_p4x)
so_p4x$empir.p <- PGEV(q = so_p4x$lod, gev = evd_p4x, lower.tail = FALSE)

sov_p4x <- scanonevar.perm(sov = sov_p4x, n.perms = num_perms, random.seed = seed, n.cores = num_cores)


saveRDS(object = list(so_p1x, sov_p1x,
											so_p2x, sov_p2x,
											so_p3x, sov_p3x,
											so_p4x, sov_p4x),
				file = 'saves/empir_p_scans_appendix.RDS')



# make plots for fig 2 -- empirircal p-value scans
ymax <- 3.5

plot(p1x_empir_p_scan <- plot(x = sov_p1x, y = so_p1x, ymax = ymax))
plot(p2x_empir_p_scan <- plot(x = sov_p2x, y = so_p2x, ymax = ymax))
plot(p3x_empir_p_scan <- plot(x = sov_p3x, y = so_p3x, ymax = ymax))
plot(p4x_empir_p_scan <- plot(x = sov_p4x, y = so_p4x, ymax = ymax))


ggplot2::ggsave(plot = p1x_empir_p_scan, filename = 'images/empir_p_scan_phen1x.pdf', height = 2.5, width = 9)
ggplot2::ggsave(plot = p2x_empir_p_scan, filename = 'images/empir_p_scan_phen2x.pdf', height = 2.5, width = 9)
ggplot2::ggsave(plot = p3x_empir_p_scan, filename = 'images/empir_p_scan_phen3x.pdf', height = 2.5, width = 9)
ggplot2::ggsave(plot = p4x_empir_p_scan, filename = 'images/empir_p_scan_phen4x.pdf', height = 2.5, width = 9)
