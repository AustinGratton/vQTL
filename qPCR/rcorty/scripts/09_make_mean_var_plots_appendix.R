xlims <- c(-0.3, 0.5)
ylims <- c(0.6, 1.6)

plot(p1x <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype1x',
																		 focal.groups = 'D1M1',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))

plot(p2x <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype2x',
																		 focal.groups = 'D1M6',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))

plot(p3x <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype3x',
																		 focal.groups = 'D2M6',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))

plot(p4x <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype4x',
																		 focal.groups = 'D3M6',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))


ggplot2::ggsave(plot = p1x, filename = 'images/mean_var_plot_phen1x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p2x, filename = 'images/mean_var_plot_phen2x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p3x, filename = 'images/mean_var_plot_phen3x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p4x, filename = 'images/mean_var_plot_phen4x.pdf', height = 3, width = 4)
