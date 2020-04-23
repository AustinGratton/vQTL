setwd("C:/Users/twili/Desktop/GIThub/Andrew/stapleton_lab/Stress_Splicing/2016_Clayton")
SamplingPlan=read.csv(file="Field Book (2016) - Clayton - Sampling Plan_RAW DATA.csv", header=TRUE)
sampNum = c(seq(1, dim(SamplingPlan)[1]))
SamplingPlan = cbind(sampNum, SamplingPlan)
SamplingPlanTidy=SamplingPlan[-c(505:516, 769:780, 949:960, 1129:1140, 1421:1424, 1813, 1838), ]

# Removed Rows:
#   -> 505-516: Plot 83, Notes: no plants
#   -> 769-780: Plot 127, Notes: no plants
#   -> 949-960: Plot 158, Notes: no plants
#   -> 1129-1140: Plot 187, Notes: no plants
#   -> 1421-1424: Plot 218 (Samples 1473-1476), Notes: samples not collected
#   -> 1813: Blank row in raw data
#   -> 1838: Variable/column labels for Plot 195
# Kept Rows: 
#   -> 1153-1164: Plot 190, Notes: send to ISU
#   -> 1177-1188: Plot 194, Notes: send to ISU
#   -> 1839-1938: Plot 195, Notes: "Chamber" notes in place of "Additional Notes"

write.csv(SamplingPlanTidy, file = "Field Book (2016) - Clayton - Sampling Plan_TIDIED.csv")
