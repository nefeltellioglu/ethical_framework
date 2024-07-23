# Script explaining some simple calculations used to obtain model parameters from 
# published sources

library(tidyverse)

# Infection hospitalization rates, by age
# https://www.cmaj.ca/content/195/42/E1427 
# British Columbia, sixth-seventh seroprevalence survey = delta/omicron

dat <- tibble(age = c("0-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "≥80"), 
              popn = c(142000, 152149, 321892, 478800 ,525414, 435537, 439590, 393225, 255419, 138469),
              #supp table 9
              IHR = c(0.13, 0.026, 0.015, 0.019, 0.041, 0.12, 0.31, 0.062, 3.27, 4.71))
              #6th/7th serosurvey Table 3

# under 70
sum(dat$IHR[1:8]*dat$popn[1:8])/sum(dat$popn[1:8])
# 0.09

# 70+
sum(dat$IHR[9:10]*dat$popn[9:10])/sum(dat$popn[9:10])
# 3.78


