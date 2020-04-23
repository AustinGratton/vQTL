xlims <- c(-0.3, 0.5)
ylims <- c(0.6, 1.6)

plot(p1 <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype1',
																		 focal.groups = 'D1M1',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))

plot(p2 <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype2',
																		 focal.groups = 'D1M6',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))

plot(p3 <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype3',
																		 focal.groups = 'D2M6',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))

plot(p4 <- mean_var_plot_model_based(cross = test_cross,
																		 phenotype.name = 'phenotype4',
																		 focal.groups = 'D3M6',
																		 point_size = 3,
																		 xlim = xlims,
																		 ylim = ylims))


ggplot2::ggsave(plot = p1, filename = 'images/mean_var_plot_phen1.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p2, filename = 'images/mean_var_plot_phen2.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p3, filename = 'images/mean_var_plot_phen3.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p4, filename = 'images/mean_var_plot_phen4.pdf', height = 3, width = 4)
