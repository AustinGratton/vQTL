library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

#set directory and load scanonevar object
setwd("/work/04902/azg5169/stampede2/vQTL/vQTL/SOV.perm")
intResult = read_rds("InteractiveResult.rds")


#inputs: SOV object, number of permutations, random seed
#default inputs: n.cores = parallel::detectcores()-2
SOVperm = scanonevar.perm(intResult, n.perms = 100)

write_rds(SOVperm, "permResult100.rds")

