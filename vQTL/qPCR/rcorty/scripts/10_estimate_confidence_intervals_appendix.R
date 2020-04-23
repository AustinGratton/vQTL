library(vqtl)

d <- readRDS(file = 'saves/LOD_scans_appendix.RDS')

# do bootstraps to estimate confidence intervals
bs_p2 <- scanonevar.boot(sov = d$sov_p2, n.resamples = 1000, chr = '1', qtl_type = 'mQTL', random.seed = 27599, n.cores = 40)

bs_p3 <- scanonevar.boot(sov = d$sov_p3, n.resamples = 1000, chr = '2', qtl_type = 'vQTL', random.seed = 27599, n.cores = 40)

bs_p4 <- scanonevar.boot(sov = d$sov_p4, n.resamples = 1000, chr = '3', qtl_type = 'mvQTL', random.seed = 27599, n.cores = 40)


saveRDS(object = list(bs_p2 = bs_p2,
											bs_p3 = bs_p3,
											bs_p4 = bs_p4),
				file = 'saves/bootstraps_appendix.RDS')

