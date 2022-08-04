setwd("~/GitHub/salmon_negbinmodel")

library(tidyverse)
library(rvest)
library(tabulizer)
library(lubridate)
library(rgdal)

source("functions.R")

data_OMR <- calc_OMR(dateStart = "1993-01-01",dateEnd = "2022-08-04", timing = "daily", extrap = T, proof = T)
#Regression of daily middle vs old river flow - adj R = 0.97
#data from 1992-12-19 to 2022-08-01

