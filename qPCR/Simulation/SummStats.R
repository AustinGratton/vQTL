# Summary Statistics on full data set
library(readr)
library(tidyverse)

df = read_csv("Fullinbhyb.csv")

df = df[-2,]

View(df[1:10,1:10])

df2= df[-1,]
View(df2[1:10,1:10])

mkrs = df[2:364,6:3240]
sapply(df[2:364,6:3240],mean)
table(scan(mkrs))


library(gdata)
vector_data<- unmatrix(mkrs,byrow=T)
tbl = table(vector_data)[c(1,3,4)]
mkrsum = sum(tbl)
pcg = tbl/mkrsum
df2 = 
df %>%
  group_by(df2[,4])
