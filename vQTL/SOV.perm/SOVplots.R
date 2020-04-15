library(gridExtra)
library(qtl)
library(vqtl)
library(ggplot2)
library(tidyverse)

vQTLsub <- read.cross(file = "vQTL/ManchingStressData_Covar.csv" )

SOV.perm = read_rds("vQTL/SOV.perm/permResult50.rds")


table(SOV.perm$result$mQTL.asymp.p <= .0001)
table(SOV.perm$result$mvQTL.asymp.p <= .0001)
table(SOV.perm$result$vQTL.asymp.p <= .0001)

SOV.perm$result$loc.name[SOV.perm$result$mQTL.asymp.p <= .0001]
SOV.perm$result$loc.name[SOV.perm$result$vQTL.asymp.p <= .0001]
SOV.perm$result$loc.name[SOV.perm$result$mvQTL.asymp.p <= .0001]


sig.loc = SOV.perm$result$loc.name[SOV.perm$result$vQTL.asymp.p <= .0001]

colnames(vQTLsub$pheno)[1] = "Height"

sig.plot1 = mean_var_plot_model_based(vQTLsub, 'Height', sig.loc[1], genotype.names = c('A','B'))
sig.plot2 = mean_var_plot_model_based(vQTLsub, 'Height', sig.loc[2], genotype.names = c('A','B'))
sig.plot3 = mean_var_plot_model_based(vQTLsub, 'Height', sig.loc[3], genotype.names = c('A','B'))
sig.plot4 = mean_var_plot_model_based(vQTLsub, 'Height', sig.loc[4], genotype.names = c('A','B'))
sig.plot5 = mean_var_plot_model_based(vQTLsub, 'Height', sig.loc[5], genotype.names = c('A','B'))
sig.plot6 = mean_var_plot_model_based(vQTLsub, 'Height', sig.loc[6], genotype.names = c('A','B'))

require(gridExtra)
grid.arrange(sig.plot1, sig.plot2, sig.plot3, sig.plot4, sig.plot5, sig.plot6, ncol = 2)

all.loc = SOV.perm$result$loc.name

PDFpath = "D:/GitHub/vQTL/vQTL/SOV.perm/plots/allPlots.pdf"
pdf(file = PDFpath)

for (i in 1:length(all.loc)){
  rm(all.plot)
  all.plot = mean_var_plot_model_based(vQTLsub, 'Env', all.loc[i], genotype.names = c('A','B'))
  print(all.plot)
}
dev.off()





#######OTHER PLOTS

SOV.perm = read_rds("vQTL/SOV.perm/permResult50.rds")

SOV.res = SOV.perm$result

SOV.sig = SOV.res[which(SOV.perm$result$vQTL.asymp.p <= .0001),]
