
library("MASS")
library("MuMIn")
library("pscl")
library("tidyverse")
library("AICcmodavg")
library("knitr")
library("PerformanceAnalytics")
library("rsq")

setwd("~/GitHub/salmon_negbinmodel")

monthly_data <- read.csv("negbinmodel_monthly_dataset.csv")

str(monthly_data)

OMR_model<-lm(omr_flow_extrap~san_joaquin_flow+export,data=monthly_data)
summary(OMR_model)
#R^2 > 0.95

#Median value for Dec-Apr SJR flow
dec_apr<-monthly_data %>% filter(month %in% c(12,1,2,3,4))
median(dec_apr$san_joaquin_flow)
#2388.989

predict(OMR_model,newdata=data.frame(san_joaquin_flow=median(dec_apr$san_joaquin_flow),
                                     export=4700))
#-3494.425