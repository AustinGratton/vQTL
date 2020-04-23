### HEIARCHICAL MODEL ###

library(MASS)
library(stringr)
## Ordinal Net package ##
library("ordinalNet")

# set directory to the Heirarchical folder
setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/Heirarchical")

#### reading in and setting up calibrated data ####
# MONTH 1 (2018_6 / JUNE) CALIBRATED DATA FRAME 
calib_data_6 = na.omit(read.csv("../2018_6/calib_2018_6.csv")[,-1])
calib_data_6$ztest1 = (calib_data_6$test1 - mean(calib_data_6$test1))/sd(calib_data_6$test1)
calib_data_6$zallP = (calib_data_6$allP - mean(calib_data_6$allP))/sd(calib_data_6$allP)
calib_data_6$month ='june'
#calib_data_6

# MONTH 2 (2018_8 / AUGUST) CALIBRATED DATA FRAME
calib_data_8 = na.omit(read.csv("../2018_8/calib_2018_8.csv")[,-1])
calib_data_8$ztest1 = (calib_data_8$test1 - mean(calib_data_8$test1))/sd(calib_data_8$test1)
calib_data_8$zallP = (calib_data_8$allP - mean(calib_data_8$allP))/sd(calib_data_8$allP)
calib_data_8$month ='aug'
#calib_data_8

# MONTH 3 (2018_11 / NOVEMBER) CALIBRATED DATA FRAME
calib_data_11 = na.omit(read.csv("../2018_11/calib_2018_11.csv")[,-1])
calib_data_11$ztest1 = (calib_data_11$test1 - mean(calib_data_11$test1))/sd(calib_data_11$test1)
calib_data_11$zallP = (calib_data_11$allP - mean(calib_data_11$allP))/sd(calib_data_11$allP)
calib_data_11$month ='nov'
#calib_data_11

# Combined Calib d.f. for all months
calib_data = rbind(calib_data_6, calib_data_8, calib_data_11)

# Create dummy varible columns for each month
calib_data$june = ifelse(str_detect(calib_data[,6], "june"), 1, 0)
calib_data$aug = ifelse(str_detect(calib_data[,6], "aug"), 1, 0)
#calib_data = calib_data[,-6]

# create a subset only containing the startq, zvalues, and dummy months
calib_subset = calib_data[,c(1, 4,5,8,7)]
######

#### Plotting calibrated data ####
#graphing log starting quantity to cp values
plot(calib_data$allP,log(calib_data$startq), col = 'red', main = "Hierarchical Log Starting Quanitty vs. Cp Values", 
      ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data$test1,log(calib_data$startq), col = 'blue')
abline(lm(log(calib_data$startq)~calib_data$allP), col = 'red')
abline(lm(log(calib_data$startq)~calib_data$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)


#graphing log starting quantity to cp values by month
#June
plot(calib_data_6$allP,log(calib_data_6$startq), col = 'red', main = "June Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_6$test1,log(calib_data_6$startq), col = 'blue')
abline(lm(log(calib_data_6$startq)~calib_data_6$allP), col = 'red')
abline(lm(log(calib_data_6$startq)~calib_data_6$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#August
plot(calib_data_8$allP,log(calib_data_8$startq), col = 'red', main = "August Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_8$test1,log(calib_data_8$startq), col = 'blue')
abline(lm(log(calib_data_8$startq)~calib_data_8$allP), col = 'red')
abline(lm(log(calib_data_8$startq)~calib_data_8$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)

#November
plot(calib_data_11$allP,log(calib_data_11$startq), col = 'red', main = "November Log Starting Quanitty vs. Cp Values", 
     ylab = "log(Starting Quantity)", xlab = "Cp Values" , xlim = c(5, 25), ylim = c(-6, 1))
points(calib_data_11$test1,log(calib_data_11$startq), col = 'blue')
abline(lm(log(calib_data_11$startq)~calib_data_11$allP), col = 'red')
abline(lm(log(calib_data_11$startq)~calib_data_11$test1), col = 'blue')
legend('topright', legend=c("Test 1", "All Products"),
       col=c("blue", "red"), lty = 1, cex=0.8)
#####

hist(calib_data$allP, col = "red")
hist(calib_data$test1, col = "blue")

# ###### POLR models ######
# # Ordinal Logistic Regression Model 
# model = polr(as.factor(calib_subset$startq) ~ ., data=calib_subset, Hess = TRUE)
# #(summary(model))
# (ctable <- coef(summary(model)))
# ## calculate and store p values
# p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
# options(scipen=999)
# ## combined table
# (ctable <- cbind(ctable, "p value" = p))
# 
# ###### Ordinal Net models #######

#define ordinal model starq~zallP+ztest1+month
ordmod = ordinalNet(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq))
summary(ordmod)
coef(ordmod, matrix=TRUE)
#kfold cv
set.seed(3)
ordfit = ordinalNetTune(as.matrix(calib_subset[,2:5]), as.factor(calib_subset$startq), family = "cumulative",
                         link = "logit", parallelTerms = TRUE, nonparallelTerms = TRUE, 
                         warn = FALSE, printProgress = FALSE)
head(ordfit$loglik)
bestLambdaIndex = which.max(rowMeans(ordfit$loglik))
head(coef(ordfit$fit, matrix = TRUE, whichLambda = bestLambdaIndex))


###### Experimental Data #######

# MONTH 1 (2018_8 / AUGUST) EXPERIMENTAL DATA FRAME
exp_data_8 = na.omit(read.csv("../2018_8/exp_2018_8.csv")[,-1])
exp_data_8$ztest1 = (exp_data_8$test1 - mean(exp_data_8$test1))/sd(exp_data_8$test1)
exp_data_8$zallP = (exp_data_8$allP - mean(exp_data_8$allP))/sd(exp_data_8$allP)
exp_data_8$month ='aug'
exp_data_8$sampleID.exp = as.factor(exp_data_8$sampleID.exp)
#exp_data_8

# MONTH 2 (2018_6 / JUNE) EXPERIMENTAL DATA FRAME 
exp_data_6 = na.omit(read.csv("../2018_6/exp_2018_6.csv")[,-c(1,5,6)])
exp_data_6$ztest1 = (exp_data_6$test1 - mean(exp_data_6$test1))/sd(exp_data_6$test1)
exp_data_6$zallP = (exp_data_6$allP - mean(exp_data_6$allP))/sd(exp_data_6$allP)
exp_data_6$month ='june'
#exp_data_6

# MONTH 3 (2018_11 / NOVEMBER) EXPERIMENTAL DATA FRAME
exp_data_11 = na.omit(read.csv("../2018_11/exp_2018_11.csv")[,-1])
exp_data_11$ztest1 = (exp_data_11$test1 - mean(exp_data_11$test1))/sd(exp_data_11$test1)
exp_data_11$zallP = (exp_data_11$allP - mean(exp_data_11$allP))/sd(exp_data_11$allP)
exp_data_11$month ='nov'
exp_data_11$sampleID.exp = as.factor(exp_data_11$sampleID.exp)
#exp_data_11

# Combined exp d.f. for all months
exp_data = rbind(exp_data_6, exp_data_8, exp_data_11)

# Create dummy varible columns for each month
exp_data$june = ifelse(str_detect(exp_data[,6], "june"), 1, 0)
exp_data$aug = ifelse(str_detect(exp_data[,6], "aug"), 1, 0)


# Drop rows containing NA
exp_subset = exp_data[,c(1, 4, 5, 7, 8)]

#### calculating experimental starting quantity ####
probmat = predict(ordfit$fit, as.matrix(exp_subset[,2:5]))
probmat[1:10,]


##### finding adjustment value and adjusted test 1 in calibrated #####

group = split.data.frame(calib_data, calib_data$startq)

adj <- function(AllP, Test1){
  adjust = ave(AllP)-ave(Test1)
  return(adjust)
}

adjval = NULL
for (k in group){
  adjval = c(adjval,adj(k$allP, k$test1))
}


##### Creating a dataframe with the stq and adjustment values #########
calib_adj = as.data.frame(cbind(unique(calib_data$startq), unique(adjval)))
colnames(calib_adj) = c("startq", "adj")

#convert adjustment from picograms to femtograms (1:1000)
calib_adj$adj = (calib_adj$adj)*1000


# convert test1 and allp to femtograms
exp_data[,c(2,3)] = exp_data[,c(2,3)]*1000

# Apply probability matrix to the adjustment values using matrix multiplication 
exp_data$exp.adjust = probmat%*%calib_adj$adj

# Create new column with stress product (VQTL input)
exp_data$exp.adjustTest1 = exp_data$test1.exp+exp_data$exp.adjust

# convert allP and adjusted test 1 to femtograms
exp_data$stress = exp_data$allP.exp - exp_data$exp.adjustTest1

boxplot(exp_data$allP.exp, exp_data$test1.exp, exp_data$stress, main = "Boxplot of Experimental All Products, Test 1, and Stress",
        names =c("All Products", "Test 1", "Stress"), ylab = "Cp value", col =c("blue", "red", "green"))

hist(exp_data$stress, col = "light blue")
### analyzing negative stress ###
negstress = exp_data %>% filter(exp_data$stress<0)
hist(negstress$stress, col = "light green")
length(negstress$stress)/length(exp_data$stress)

# making negstress "NA" #
na.obs = match(negstress$stress,exp_data$stress) #returns list of incices
# we want to replace the negative values with N/A
exp_data[na.obs, ] = NA

### Write the exp data with the stress product as a new data frame for vqtl matching
write.csv(exp_data, "Hierarchical_exp_data_stress.csv")
