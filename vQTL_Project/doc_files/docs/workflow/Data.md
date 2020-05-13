---
title: "Data"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# Data Collected 

## qPCR Data
Not sure what this looks like

## vQTL Data
The [data](https://github.com/MRCopeland74/stapleton_lab/blob/master/vQTL/ManchingStressData_Covar.csv) we used for the vQTL data was all in one csv file known as the Manching Stress Product Data. This file conatained a phenotype which is the height of the corn crop. There are 8 columns of different environment combinations of low water, low nitrogen, or presence of a pathogen. There is an environment column numbering the different combinations from 1-8. These combinations are either 1 or 0.  Then there is 3235 columns of different gene names. These are either A, B or NA. The for the rows we have a row indicating the chromosones that the genes are on. There are 10 different chromosones. There is another row that indicates the distance the gene is on the chromosone. The next 6672 rows are different tests with varying gene combinations and varying envronmental combinations. The data is very similar however not the same for each environmental combination, and for some environmental combinations more tests have been done.
```{r, echo = FALSE, results = 'asis'}
 library(knitr)
 Maindata = read.csv("https://github.com/MRCopeland74/stapleton_lab/blob/master/vQTL/ManchingStressData_Covar.csv")
 top10 = Maindata[[1:10]][1:20]
 kable(top10)
 ```

