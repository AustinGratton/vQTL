library(qtl)
library(vqtl)
library(evd)
library(tidyverse)
library(RColorBrewer)

params <- expand.grid(num_inds = c(100, 200, 400, 800),
                      num_markers = c(100, 200, 400, 800),
                      t = NA)

for (i in 1:nrow(params)) {

  message('Staring row ', i)
  c <- sim.cross(map = sim.map(len = rep(100, 1), n.mar = params$num_markers[i], eq.spacing = TRUE), n.ind = params$num_inds[i])
  c <- calc.genoprob(cross = c, step = 0)

  c$pheno$sex <- sample(x = 0:1, size = nind(c), replace = TRUE)

  sov <- scanonevar(cross = c,
                    mean.formula = phenotype ~ sex + mean.QTL.add + mean.QTL.dom,
                    var.formula = ~ sex + var.QTL.add + var.QTL.dom)

  start <- Sys.time()
  scanonevar.perm(sov = sov, n.perms = 1000, n.cores = 32, random.seed = 27599)
  t <- difftime(time1 = Sys.time(), time2 = start, units = 's')

  params$t[i] <- round(x = t)

}

print(params)
saveRDS(object = params, file = 'saves/benchmark_sim_cross.RDS')

params <- readRDS(file = 'saves/benchmark_sim_cross.RDS')

params %>%
  ggplot(mapping = aes(x = factor(num_markers), y = factor(num_inds), fill = t/60)) +
  geom_tile(color = 'black') +
  geom_label(mapping = aes(label = paste(round(x = t/60, digits = 0), 'm')), fill = 'white') +
  scale_fill_continuous(low = 'white', high = 'black', name = 'time (min)') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(axis.title = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'line')),
        axis.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'line')),
        plot.margin = unit(x = c(1, 1, 1, 1)/2, units = 'lines'),
        legend.title.align = 0.2) +
  xlab('Number of Markers') +
  ylab('Number of Individuals')


params %>%
  ggplot(mapping = aes(x = num_markers, y = t/60, color = factor(num_inds))) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(name = 'population size',
                     values = brewer.pal(n = 6, name = 'Blues')[3:6]) +
  scale_x_continuous(name = 'number of markers') +
  scale_y_continuous(name = 'time (minutes)') +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.75),
        legend.background = element_rect())

ggsave(filename = 'images/benchmark_sim_cross.pdf', height = 4, width = 6)

# params %>%
#   ggplot(mapping = aes(x = num_inds, y = t/60, color = factor(num_markers))) +
#   geom_line() +
#   geom_point() +
#   scale_color_discrete(name = 'number of markers') +
#   scale_x_continuous(name = 'number of individuals') +
#   scale_y_continuous(name = 'time (minutes)') +
#   theme_minimal()


library(broom)

params %>%
  group_by(num_inds) %>%
  do(tidy(lm(formula = t ~ 0 + num_markers, data = .)))
