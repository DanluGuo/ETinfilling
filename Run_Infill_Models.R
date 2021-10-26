# load example test dataset for missing morning, afternoon and evening
load("Testdata_MissingMorning.RData")
#load("Testdata_MissingAfternoon.RData")
#load("Testdata_MissingEvening.RData")

# examing the format of test data 
head(testET)
# column 'ET' is the ETa where gaps exist (to be infilled)

# load the functions for four ET infilling models
source("4_ETinfilling_Models.R")

# Executing each function to generate gap-free ET data
Model1(testET)
Model2(testET)
Model3(testET)
Model4(testET,trainET)
