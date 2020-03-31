########################################################## 
############## QPCR PLATE & ADJUSTMENT MODEL #############
########################################################## 

library(tidyr)
#install.packages("pracma")
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
library(MASS)
#install.packages("glm.predict")
library(glm.predict)
#install.packages("Stack")
library(Stack)

### READ IN DERIVATIVE DATA ###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
#deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
#deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
#deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)
deriv_complete=read.csv(file = "../2018_6_1_qPCR_Output.csv", header=FALSE)

rm(test)

########################################################## 
################### Initial Data Framing #################
########################################################## 
deriv1 = deriv_complete
# Remove extra labels column 
deriv1 = deriv1[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv1 = as.data.frame(t(deriv1), header=FALSE)
# Remove blank column (4th)
#deriv1 = deriv1[,-5]
# Rename columns
colnames(deriv1)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv1['sampleID_NTC'] = grepl('NTC', deriv1$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv1$sampleID_NTC)
deriv1 = deriv1[-ntc,]
deriv1 = deriv1[,-c(6,7)]
# Indicate if sample is 'Plus' or 'Minus'
deriv1['sampleID_Minus'] = grepl('minus', deriv1$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv1$sampleID_Minus)
deriv1 = deriv1[-minus,]
deriv1 = deriv1[,-6]
deriv1$cpD1 = as.numeric(as.character(deriv1$cpD1))

### COMPLETED INITIAL DATA FRAMING ###

########################################################## 
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
deriv1 = deriv1 %>% filter(deriv1$cpD1 >= 10)
# Read in raw cycle data - may need to combine multiple files
cycle1 = read.csv(file = "../2018_6_1_plate.csv", header = FALSE)
# Create complete set of reaction data (derivative and cycle)
reaction = Stack(deriv_complete, cycle1)
# Remove repeat labeling
replace = reaction[7:10,]
reaction = reaction[-c(1:4, 7:10),]
reaction = Stack(replace, reaction)
# Transpose so column headers at top
reaction = as.data.frame(t(reaction))
reaction = reaction[,-c(6:7)]
# Replace column names with first row
colnames(reaction) <- as.character(unlist(reaction[1,]))
reaction = reaction[-1,]
colnames(reaction)[5] = "cpD1"
reaction$cpD1 = as.numeric(as.character(reaction$cpD1))
# Filter unusual observations (CT value less than 10)
unusual_obs_2018_6 = reaction %>% filter(reaction$cpD1 < 10)
# Write CSV file 
#write.csv(unusual_obs_2018_6, file="Unusual_Obs_2018_6.csv")

# ### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###


########################################################## 
################# Calibrated Data Framing ################
########################################################## 
#
library("rowr")
# Create/Write data frame for Calibrated values
calib_data = deriv1 %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data = calib_data[order(calib_data$starting_quantity),]

calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))


test1 = filter(calib_data, reaction_type=="test1")[,c(5,1)]
allP = filter(calib_data, reaction_type=="all_products")[,c(1,4:5)]


#Combine test1 and allP obs, with NA in blank cells
calib_data = as.data.frame(cbind.fill(allP, test1, fill = NA))
colnames(calib_data) = c("startq", 'allP', "test1" )
# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)

#Apply log scale to test1 and allP CT values
calib_data$allPln = log(calib_data$allP)
calib_data$test1ln = log(calib_data$test1)

write.csv(calib_data, file = "calib_2018_6.csv")


### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv1 %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find counts of each unique sampleID; for sample with a count not equal to 2, remove from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
# Remove invalid observations from data set
exp_data = exp_data[!exp_data$sampleID %in% countsne2$Var1,]

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID, convert to data frame
exp_data = as.data.frame(cbind(sampleID.exp, test1.exp, allP.exp))
exp_data$test1.exp = as.numeric(as.character(exp_data$test1.exp))
exp_data$allP.exp = as.numeric(as.character(exp_data$allP.exp))

# Apply log scale to test1 and allP CT values
exp_data$test1.exp.ln = log(exp_data$test1.exp)
exp_data$allP.exp.ln = log(exp_data$allP.exp)

write.csv(exp_data, file = "exp_2018_6.csv")
### COMPLETED EXPERIMENTAL DATA FRAME ###
