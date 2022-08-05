setwd("~/GitHub/salmon_negbinmodel")

library(tidyverse)
library(rvest)
library(tabulizer)
library(lubridate)
library(rgdal)

source("functions.R")

###
#Pull OMR Data
data_OMR <- calc_OMR(dateStart = "1993-01-01",dateEnd = "2022-08-01", timing = "daily", extrap = T, proof = T)
#Regression of daily middle vs old river flow - adj R = 0.97
#data from 1992-12-19 to 2022-08-01

###
#Pull Dayflow Data
dayflow_1984_1996<- read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/cb04e626-9729-4105-af81-f6e5a37f116a/download/dayflow-results-1984-1996.csv")%>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))
dayflow_1997_2020<- read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/21c377fe-53b8-4bd6-9e1f-2025221be095/download/dayflow-results-1997-2020.csv")%>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))
dayflow_2021<-read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/83122ce7-e7f5-4ad1-b5e9-6a7032cb1117/download/dayflowcalculations2021.csv")

dayflow_combined <- bind_rows(dayflow_1984_1996,dayflow_1997_2020,dayflow_2021) %>% 
  # using the newer col names
  mutate(EXPORTS = coalesce(EXPORT, EXPORTS),
         DIVER = coalesce(DIVE, DIVER),
         EFFEC = coalesce(EFFECT, EFFEC),
         EFFDIV = coalesce(EFFD, EFFDIV), Month=month(Date)) %>% select(-c(Y, EXPORT, DIVE, EFFECT, EFFD))
remove(dayflow_1984_1996,dayflow_1997_2020,dayflow_2021)


