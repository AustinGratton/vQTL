library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)

setwd("/work/04902/azg5169/stampede2/vQTL/vQTL/SOV.perm")
intResult = read_rds("InteractiveResult.rds")

SOVperm = scanonevar.perm(intResult, n.perms = 2, random.seed = 3112020)

write_rds(SOVperm, "permResult.rds")