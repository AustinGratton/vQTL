########################################################## 
############## QPCR PLATE & ADJUSTMENT MODEL #############
########################################################## 

library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
library(MASS)
library(glm.predict)
library(Stack)
library(rowr)

### READ IN DERIVATIVE DATA ###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
deriv.1<-read.csv(file = "2018_11_1_qPCR_Output.csv", header=FALSE)
deriv.2<-read.csv(file = "2018_11_2_qPCR_Output.csv", header=FALSE)
deriv_complete=as.data.frame(cbind(deriv.1, deriv.2))

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)
#deriv=read.csv(file = "2018_06_01_plate_qPCR_output_2019_04_04.csv", header=FALSE)

########################################################## 
################### Initial Data Framing #################
########################################################## 
deriv = deriv_complete
# Remove extra column 
deriv = deriv[,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
# Rename columns
colnames(deriv)=c("plateID", "reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
### Removing NTC and gblock-Minus values ###
# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(6,7)]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
# IF "minus" RETURNS EMPTY VALUES, COMMENT OUT COMMAND BELOW
deriv = deriv[-minus,]
deriv = deriv[,-6]
# Remove two extra label rows from center of data frame
deriv['label.row'] = grepl('3', deriv$starting_quantity)
extra = which(deriv$label.row)
deriv = deriv[-extra,]
deriv = deriv[,-6]
deriv$cpD1 = as.numeric(as.character(deriv$cpD1))
### COMPLETED INITIAL DATA FRAMING ###


########################################################## 
############ Removing Ununsual Observations ##############
##########################################################

# Remove unusual observations from initial data frame (CT value less than 10)
deriv = deriv %>% filter(deriv$cpD1 >= 10)
# Read in raw cycle data - may need to combine multiple files
cycle1 = read.csv(file = "2018_11_1_plate.csv", header = FALSE)
cycle2 = read.csv(file = "2018_11_2_plate.csv", header = FALSE)
cycle = as.data.frame(cbind(cycle1, cycle2))
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
unusual_obs_2018_11 = reaction %>% filter(reaction$cpD1 < 10)
# Write CSV file 
#write.csv(unusual_obs_2018_11, file="Unusual_Obs_2018_11.csv")
### COMPLETED UNUSUAL OBSERVATIONS REMOVAL/REPORTING ###

########################################################## 
################# Calibrated Data Framing ################
##########################################################
# Create/Write data frame for Calibrated values
calib_data = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_data = calib_data[order(calib_data$starting_quantity),]

calib_data$starting_quantity = as.numeric(as.character(calib_data$starting_quantity))
calib_data$cpD1 = as.numeric(as.character(calib_data$cpD1))


test1 = filter(calib_data, reaction_type=="test1")[,5]
allP = filter(calib_data, reaction_type=="all_products")[,4:5]
# Combine test1 and allP obs, with NA in blank cells
calib_data = as.data.frame(cbind.fill(allP, test1, fill = NA))
colnames(calib_data) = c("startq", 'allP', "test1")

# Format starting quantity values as decimals, not scientific notation
calib_data$startq=as.factor(format(calib_data$startq, scientific=FALSE))
calib_data$startq=as.factor(calib_data$startq)
#write.csv(calib_data, file = "calib_2018_11.csv")
### COMPLETED CALIBRATED DATA FRAME ###


########################################################## 
############### Experimental Data Framing ################
########################################################## 

# Create/Write data frame for Experimental values
exp_data = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_data = exp_data[order(exp_data$starting_quantity),]
# # Remove first and last rows (unnecessary labeling)
# exp_data = exp_data[-1,]
# exp_data = exp_data[-nrow(exp_data),]
exp_data$cpD1 = as.numeric(as.character(exp_data$cpD1))
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
### Finding invalid observations ###
# Find invalid observations - Find counts of each unique sampleID; remove ones with count not equal to 2 from data frame
counts = as.data.frame(table(exp_data$sampleID))
countsne2 = as.data.frame(filter(counts, !counts$Freq==2))
countsne2$Var1 = as.numeric(as.character(countsne2$Var1)) 
# Remove observations with count not equal to 2 from data set
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
#write.csv(exp_data, file = "exp_2018_11.csv")
### COMPLETED EXPERIMENTAL DATA FRAME ###
