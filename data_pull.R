setwd("~/GitHub/salmon_negbinmodel")

library(tidyverse)
library(rvest)
library(tabulizer)
library(lubridate)
library(rgdal)
library(janitor)

source("functions.R")

#####
#Install "deltaFish" package if necessary
# Enable this universe
#options(repos = c(sbashevkin = 'https://sbashevkin.r-universe.dev',CRAN = 'https://cloud.r-project.org'))
## Install the package
#install.packages('deltafish')

library(deltafish)

###############
#Pull OMR Data

data_OMR <- calc_OMR(dateStart = "1993-01-01",dateEnd = "2022-08-01", timing = "daily", extrap = T, proof = T)
#Regression of daily middle vs old river flow - adj R = 0.97
#data from 1992-12-19 to 2022-08-01

#data_OMR summary
data_OMR_sum <- data_OMR %>% 
  filter(month %in% c(1,2,3,4,5,6,11,12)) %>% group_by(month, year) %>%summarise(omrFlowExtrap=mean(omrFlowExtrap, na.rm=T))

data_OMR_samplesize<-data_OMR_sum %>% mutate(samplesize=1) %>% group_by(year) %>% summarise(samplesize=sum(samplesize))
#Each month has data
remove(data_OMR_sum,data_OMR_samplesize)

#filter just relevant columns
data_OMR <- data_OMR %>% rename(Date=date) %>% select(Date, omrFlowExtrap)


###############
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
#filter just relevant columns
dayflow_combined <- dayflow_combined %>% select(Date, SAC, SJR, EXPORTS, OUT) %>% rename(Sac_flow=SAC, SanJoaquin_flow=SJR, Export=EXPORTS, DeltaOutflow=OUT)

###############
#Pull salvage data

salvage_data<-pull_salvage() 
str(salvage_data)
salvage_data_winter<- salvage_data %>% rename(Date='Sample Time',LAD_Race='LAD Race',LAD_Loss='LAD Loss',ExpandedSalvage='Expanded Salvage') %>%
  mutate(Date=as.Date(Date)) %>% filter(LAD_Race=="Winter") %>% group_by(Date) %>% 
  summarise(Winter_LAD_Loss=sum(LAD_Loss), Winter_ExpandedSalvage=sum(ExpandedSalvage))

salvage_data_spring<- salvage_data %>% rename(Date='Sample Time',LAD_Race='LAD Race',LAD_Loss='LAD Loss',ExpandedSalvage='Expanded Salvage') %>%
  mutate(Date=as.Date(Date)) %>% filter(LAD_Race=="Spring") %>% group_by(Date) %>% 
  summarise(Spring_LAD_Loss=sum(LAD_Loss), Spring_ExpandedSalvage=sum(ExpandedSalvage))


###############
#Pull monitoring data using 'deltafish' package
create_fish_db()
# open deltafish two data files
survey <- open_survey()
fish<- open_fish()

# filter for sources and taxa of interest
survey_DJFMP<- survey %>% 
  filter(Source == "DJFMP") 

fish_salmon <- fish %>% 
  filter(Taxa %in% c("Oncorhynchus tshawytscha"))

# do a join and collect the resulting data frame
# collect executes the sql query and gives you a table
sactrawl_data <- left_join(survey_DJFMP, fish_salmon) %>% 
  collect() %>% filter(Station %in% c("SR055M","SR055W","SR055X"))

# add length at date data
LAD<- read_csv("LAD.csv") %>% mutate(Date=as.Date(Date)) %>% 
  mutate(MonthDay = paste0(month(Date), "-", day(Date))) %>% rename(ForkLength = FL) %>% select(-Date)

# Summarize winter-run data
sactrawl_data_WR<- sactrawl_data %>% group_by(Station, Date, Datetime, SampleID) %>% rename(ForkLength=Length) %>%
  mutate(MonthDay = paste0(month(Date), "-", day(Date))) %>% 
  left_join(LAD) %>% filter(year(Date)>=1993) %>% 
  #Change to winter-run size only
  mutate(Count=ifelse(Race=="W",Count,0))%>% mutate(Count = replace_na(Count, 0)) %>%
  group_by(Station,SampleID,Date,Datetime) %>% summarise(Count=sum(Count)) %>% mutate(SampleSize=1) %>%
  group_by(Date) %>% summarise(Count=sum(Count),SampleSize=sum(SampleSize)) %>%
  rename(SacTrawlWR_Count=Count, SacTrawl_SampleSize=SampleSize)

# Summarize spring-run data
sactrawl_data_SR<- sactrawl_data %>% group_by(Station, Date, Datetime, SampleID) %>% rename(ForkLength=Length) %>%
  mutate(MonthDay = paste0(month(Date), "-", day(Date))) %>% 
  left_join(LAD) %>% filter(year(Date)>=1993) %>% 
  #Change to winter-run size only
  mutate(Count=ifelse(Race=="S",Count,0))%>% mutate(Count = replace_na(Count, 0)) %>%
  group_by(Station,SampleID,Date,Datetime) %>% summarise(Count=sum(Count)) %>% mutate(SampleSize=1) %>%
  group_by(Date) %>% summarise(Count=sum(Count),SampleSize=sum(SampleSize)) %>%
  rename(SacTrawlSR_Count=Count, SacTrawl_SampleSize=SampleSize)

#Combine
sactrawl_data<- full_join(sactrawl_data_WR,sactrawl_data_SR)
remove(sactrawl_data_WR,sactrawl_data_SR)


###############
###############
# Combine into a single data frame
final_df<- full_join(dayflow_combined,data_OMR) %>% filter(year(Date)>=1993) %>%
  left_join(salvage_data_winter) %>% left_join(salvage_data_spring) %>%
  mutate(Winter_LAD_Loss = replace_na(Winter_LAD_Loss, 0),
         Winter_ExpandedSalvage = replace_na(Winter_ExpandedSalvage, 0),
         Spring_LAD_Loss = replace_na(Spring_LAD_Loss, 0),
         Spring_ExpandedSalvage = replace_na(Spring_ExpandedSalvage, 0)) %>%
  left_join(sactrawl_data) %>% #Restrict to just 12-31-2020 based on Sac Trawl data
  filter(year(Date)<2021) %>%
  janitor::clean_names()
  
write.csv(final_df,file="negbinmodel_daily_dataset.csv", row.names = F)

###############
###############
# Monthly dataset
monthly_df <- final_df %>%
  mutate(year = year(date),
         month = month(date), 
         sac_trawl_wr_cpue = sac_trawl_wr_count/sac_trawl_sample_size,
         sac_trawl_sr_cpue = sac_trawl_sr_count/sac_trawl_sample_size) %>%
  select(-date) %>%
  group_by(year, month) %>%
  summarize(across(.cols = everything(), ~mean(.x,na.rm = TRUE)))%>%
  ungroup() %>%
  mutate(across(.cols = everything(), ~ifelse(is.nan(.x), NA, .x)))

write.csv(monthly_df, file = "negbinmodel_monthly_dataset.csv", row.names = F)

