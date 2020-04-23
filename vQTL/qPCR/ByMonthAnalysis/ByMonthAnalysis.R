# By Month qpcr analysis model

library(MASS)
library(stringr)
## Ordinal Net package ##
library("ordinalNet")

#set working directory to ByMonthAnalysis folder
setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/ByMonthAnalysis")

#####################################
###### MONTH 1 (2018_6 / JUNE) ######
#####################################

#CALIBRATED DATA FRAME 
calib_data_6 = na.omit(read.csv("../2018_6/calib_2018_6.csv")[,-1])
calib_data_6$ztest1 = (calib_data_6$test1 - mean(calib_data_6$test1))/sd(calib_data_6$test1)
calib_data_6$zallP = (calib_data_6$allP - mean(calib_data_6$allP))/sd(calib_data_6$allP)
calib_data_6$month ='june'
#calib_data_6

plot(calib_data_6$allP,log(calib_data_6$startq), col = 'red', main = "June Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_6$test1,log(calib_data_6$startq), col = 'blue')
abline(lm(log(calib_data_6$startq)~calib_data_6$allP), col = 'red')
abline(lm(log(calib_data_6$startq)~calib_data_6$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#define ordinal model starq~zallP+ztest1
ordmod6 = ordinalNet(as.matrix(calib_data_6[,4:5]), as.factor(calib_data_6$startq))
summary(ordmod6)
coef(ordmod6, matrix=TRUE)
#kfold cv
set.seed(3)
ordfit6 = ordinalNetTune(as.matrix(calib_data_6[,4:5]), as.factor(calib_data_6$startq), family = "cumulative",
                        link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                        warn = FALSE, printProgress = FALSE)
head(ordfit6$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit6$loglik))
head(coef(ordfit6$fit, matrix = TRUE, whichLambda = bestLambdaIndex))

# MONTH 1 (2018_6 / JUNE) EXPERIMENTAL DATA FRAME 
exp_data_6 = na.omit(read.csv("../2018_6/exp_2018_6.csv")[,-c(1,5,6)])
exp_data_6$ztest1 = (exp_data_6$test1 - mean(exp_data_6$test1))/sd(exp_data_6$test1)
exp_data_6$zallP = (exp_data_6$allP - mean(exp_data_6$allP))/sd(exp_data_6$allP)
exp_data_6$month ='june'
#exp_data_6

#create a prob matrix using the model and month
probmat6 = predict(ordfit6$fit, as.matrix(exp_data_6[,4:5]))
probmat6[1:10,]

##### finding adjustment value and adjusted test 1 in calibrated #####
group = split.data.frame(calib_data_6, calib_data_6$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(k$allP, k$test1))
}

adj6 = as.data.frame(cbind(as.factor(unique(calib_data_6$startq)), unique(adjval)))
colnames(adj6) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
adj6$adj = (adj6$adj)*1000
# convert test1 and allp to femtograms
exp_data_6[,c(2,3)] = exp_data_6[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
# this creates a weighted average for the adjustment for each observation
exp_data_6$exp.adjust = probmat6%*%adj6$adj

# calculate the adjusted test 1
exp_data_6$exp.adjustTest1 = exp_data_6$test1.exp+exp_data_6$exp.adjust

# Create new column with stress product (VQTL input)
exp_data_6$stress = exp_data_6$allP.exp - exp_data_6$exp.adjustTest1
summary(exp_data_6$stress)
boxplot(exp_data_6$stress)
head(sort(exp_data_6$stress), 20)
### analyzing negative stress ###
negstress6 = na.omit(ifelse(exp_data_6$stress<0,exp_data_6$stress,NA))
hist(negstress6, col = "light green")
length(negstress6)/length(exp_data_6$stress) # 29.49% of stress is negative


######################## End of June ########################

#####################################
##### MONTH 2 (2018_8 / AUGUST) #####
#####################################

#CALIBRATED DATA FRAME
calib_data_8 = na.omit(read.csv("../2018_8/calib_2018_8.csv")[,-1])
calib_data_8$ztest1 = (calib_data_8$test1 - mean(calib_data_8$test1))/sd(calib_data_8$test1)
calib_data_8$zallP = (calib_data_8$allP - mean(calib_data_8$allP))/sd(calib_data_8$allP)
calib_data_8$month ='aug'
#calib_data_8

#define ordinal model starq~zallP+ztest1
ordmod8 = ordinalNet(as.matrix(calib_data_8[,4:5]), as.factor(calib_data_8$startq))
summary(ordmod8)
coef(ordmod8, matrix=TRUE)
#kfold cv
set.seed(3)
ordfit8 = ordinalNetTune(as.matrix(calib_data_8[,4:5]), as.factor(calib_data_8$startq), family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit8$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit8$loglik))
head(coef(ordfit8$fit, matrix = TRUE, whichLambda = bestLambdaIndex))

# MONTH 2 (2018_8 / AUGUST) EXPERIMENTAL DATA FRAME
exp_data_8 = na.omit(read.csv("../2018_8/exp_2018_8.csv")[,-1])
exp_data_8$ztest1 = (exp_data_8$test1 - mean(exp_data_8$test1))/sd(exp_data_8$test1)
exp_data_8$zallP = (exp_data_8$allP - mean(exp_data_8$allP))/sd(exp_data_8$allP)
exp_data_8$month ='aug'
exp_data_8$sampleID.exp = as.factor(exp_data_8$sampleID.exp)
#exp_data_8

#create a prob matrix using the model and month
probmat8 = predict(ordfit8$fit, as.matrix(exp_data_8[,4:5]))
probmat8[1:10,]

##### finding adjustment value and adjusted test 1 in calibrated #####
group = split.data.frame(calib_data_8, calib_data_8$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(k$allP, k$test1))
}

adj8 = as.data.frame(cbind(as.factor(unique(calib_data_8$startq)), unique(adjval)))
colnames(adj8) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
adj8$adj = (adj8$adj)*1000
# convert test1 and allp to femtograms
exp_data_8[,c(2,3)] = exp_data_8[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data_8$exp.adjust = probmat8%*%adj8$adj

# Create new column with stress product (VQTL input)
exp_data_8$exp.adjustTest1 = exp_data_8$test1.exp+exp_data_8$exp.adjust

# convert allP and adjusted test 1 to femtograms
exp_data_8$stress = exp_data_8$allP.exp - exp_data_8$exp.adjustTest1
summary(exp_data_8$stress)
boxplot(exp_data_8$stress)
head(sort(exp_data_8$stress), 20)
### analyzing negative stress ###
negstress8 = na.omit(ifelse(exp_data_8$stress<0,exp_data_8$stress,NA))
hist(negstress8, col = "light green")
length(negstress8)/length(exp_data_8$stress) # 13.17% of stress is negative

####################### End of August #######################

#####################################
#### MONTH 3 (2018_11 / NOVEMBER) ###
#####################################

#CALIBRATED DATA FRAME
calib_data_11 = na.omit(read.csv("../2018_11/calib_2018_11.csv")[,-1])
calib_data_11$ztest1 = (calib_data_11$test1 - mean(calib_data_11$test1))/sd(calib_data_11$test1)
calib_data_11$zallP = (calib_data_11$allP - mean(calib_data_11$allP))/sd(calib_data_11$allP)
calib_data_11$month ='nov'
#calib_data_11

#define ordinal model starq~zallP+ztest1
ordmod11 = ordinalNet(as.matrix(calib_data_11[,4:5]), as.factor(calib_data_11$startq))
summary(ordmod11)
coef(ordmod11, matrix=TRUE)
#kfold cv
set.seed(3)
ordfit11 = ordinalNetTune(as.matrix(calib_data_11[,4:5]), as.factor(calib_data_11$startq), family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit11$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit11$loglik))
head(coef(ordfit11$fit, matrix = TRUE, whichLambda = bestLambdaIndex))

# MONTH 3 (2018_11 / NOVEMBER) EXPERIMENTAL DATA FRAME
exp_data_11 = na.omit(read.csv("../2018_11/exp_2018_11.csv")[,-1])
exp_data_11$ztest1 = (exp_data_11$test1 - mean(exp_data_11$test1))/sd(exp_data_11$test1)
exp_data_11$zallP = (exp_data_11$allP - mean(exp_data_11$allP))/sd(exp_data_11$allP)
exp_data_11$month ='nov'
exp_data_11$sampleID.exp = as.factor(exp_data_11$sampleID.exp)
#exp_data_11

#create a prob matrix using the model and month
probmat11 = predict(ordfit11$fit, as.matrix(exp_data_11[,4:5]))
probmat11[1:10,]

##### finding adjustment value and adjusted test 1 in calibrated #####
group = split.data.frame(calib_data_11, calib_data_11$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(k$allP, k$test1))
}

adj11 = as.data.frame(cbind(as.factor(unique(calib_data_11$startq)), unique(adjval)))
colnames(adj11) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
adj11$adj = (adj11$adj)*1000
# convert test1 and allp to femtograms
exp_data_11[,c(2,3)] = exp_data_11[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data_11$exp.adjust = probmat11%*%adj11$adj

# Create new column with stress product (VQTL input)
exp_data_11$exp.adjustTest1 = exp_data_11$test1.exp+exp_data_11$exp.adjust

# convert allP and adjusted test 1 to femtograms
exp_data_11$stress = exp_data_11$allP.exp - exp_data_11$exp.adjustTest1
summary(exp_data_11$stress)
boxplot(exp_data_11$stress)
head(sort(exp_data_11$stress), 20)
### analyzing negative stress ###
negstress11 = na.omit(ifelse(exp_data_11$stress<0,exp_data_11$stress,NA))
hist(negstress11, col = "light green")
length(negstress11)/length(exp_data_11$stress) # 7.14% of stress is negative

###################### End of November ######################

allStress = c(exp_data_6$stress, exp_data_11$stress, exp_data_8$stress)
NegStress = c(negstress11, negstress6, negstress8)

length(NegStress)/length(allStress)
