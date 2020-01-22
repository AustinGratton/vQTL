library(qtl)
library(vqtl)
library(evd)

c <- read.cross(format = 'csv',  file = 'data/Kumar2014.csv', genotypes = c('AA', 'AB', 'BB'))
c <- calc.genoprob(cross = c, step = 2)

sov <- scanonevar(cross = c,
                  mean.formula = Avg_Counts ~ sex + mean.QTL.add + mean.QTL.dom,
                  var.formula = ~ sex + var.QTL.add + var.QTL.dom)


# params <- expand.grid(num_perms = c(10, 20, 30, 40), num_cores = c(10, 20, 30, 40), t = NA)
params <- expand.grid(num_perms = c(10, 100, 500, 1000),
                      num_cores = c(4, 8, 16, 32),
                      t = NA)
params

for (i in 1:nrow(params)) {

  message('Starting on row ', i)

  t <- system.time(scanonevar.perm(sov = sov, n.perms = params$num_perms[i], n.cores = params$num_cores[i], random.seed = 27599), gcFirst = TRUE)[3]

  params$t[i] <- t

}

saveRDS(object = params, file = 'saves/benchmark_kumar.RDS')

library(tidyverse)

params <- readRDS(file = 'saves/benchmark_kumar.RDS')

params %>%
  ggplot(mapping = aes(x = factor(num_cores), y = factor(num_perms), fill = t/60)) +
  geom_tile(color = 'black') +
  geom_label(mapping = aes(label = paste(round(x = t/60, digits = 0), 'm')), fill = 'white') +
  scale_fill_continuous(low = 'white', high = 'black', name = 'time (m)') +
  theme_void() +
  theme(axis.title = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.margin = unit(x = c(1, 1, 1, 1)/2, units = 'lines'),
        legend.title.align = 0.2) +
  xlab('Number of cores') +
  ylab('Number of permutations')

params %>%
  ggplot(mapping = aes(x = num_cores, y = t, color = factor(num_perms))) +
  geom_line() +
  geom_point()

params %>%
  ggplot(mapping = aes(x = num_perms, y = t/60, color = factor(num_cores))) +
  geom_line() +
  geom_point() +
  scale_color_manual(name = 'number of cores',
                     values = brewer.pal(n = 6, name = 'Blues')[3:6]) +
  scale_y_continuous(name = 'time (minutes)') +
  scale_x_continuous(name = 'number of permutations') +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.75),
        legend.background = element_rect())

ggsave(filename = 'images/benchmark_kumar.pdf', height = 4, width = 6)


library(broom)

params %>%
  group_by(num_cores) %>%
  do(tidy(lm(formula = t ~ 0 + num_perms, data = .)))
