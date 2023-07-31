setwd("~/GitHub/salmon_negbinmodel")

library("MASS")
library("MuMIn")
library("pscl")
library("tidyverse")
library("AICcmodavg")
library("lubridate")

monthly_data <- read.csv("negbinmodel_monthly_dataset.csv")


##################### Load steelhead salvage data from CDFW website
salvage_data_steelhead<-read.csv(file.path("original_data/SalvageExport_Steelhead.csv")) %>% mutate(Sample.Date=as.Date(Sample.Date)) %>% group_by(Sample.Date) %>% 
  summarise(Loss=sum(Loss)) %>%
  mutate(WY=ifelse(month(Sample.Date)>9,year(Sample.Date)+1,year(Sample.Date)),year=year(Sample.Date),month=month(Sample.Date)) %>%
  group_by(year,month) %>% summarise(steelhead_loss=sum(Loss))

monthly_data<- monthly_data %>% left_join(salvage_data_steelhead)


#Data for steelhead
data_STH<-monthly_data %>% filter(month %in% c(12,1,2,3,4,5,6)) %>%
  #standardize data by z-score and add Sac Trawl CPUE, also convert month to factor
  mutate(sac_flow_z = scale(sac_flow),
         san_joaquin_flow_z = scale(san_joaquin_flow),
         export_z = scale(export),
         delta_outflow_z = scale(delta_outflow),
         omr_flow_extrap_z = scale(omr_flow_extrap),
         month_factor=as.factor(month))

hist(data_STH$steelhead_loss)

#Remove missing loss data
data_STH <- data_STH %>% filter(!is.na(steelhead_loss))

###################DO VIF Analysis
str(data_STH)


sth_full_model<-glm.nb(steelhead_loss~ sac_flow_z + san_joaquin_flow_z + export_z + omr_flow_extrap_z, data=data_STH)
car::vif(sth_full_model)
#sac_flow_z san_joaquin_flow_z           export_z  omr_flow_extrap_z 
#3.137629          18.359772           8.053394          23.584154 

#Removed OMR flow due to excess collinearity (VIF > 10)

sth_full_model<-glm.nb(steelhead_loss~ sac_flow_z + san_joaquin_flow_z + export_z, data=data_STH)
car::vif(sth_full_model)
#sac_flow_z san_joaquin_flow_z           export_z 
#3.116307           2.860822           1.228729 

#We will go with 3 VIF threshold just as the other salmon

##############################
#Model selection for steelhead

Cand.set.STH <- list( )

Cand.set.STH[[1]] <-  glm.nb(steelhead_loss~NULL, data=data_STH)
Cand.set.STH[[2]] <-  glm.nb(steelhead_loss~month_factor, data=data_STH)
Cand.set.STH[[3]] <-  glm.nb(steelhead_loss~san_joaquin_flow_z, data=data_STH)
Cand.set.STH[[4]] <-  glm.nb(steelhead_loss~export_z, data=data_STH)
Cand.set.STH[[5]] <-  glm.nb(steelhead_loss~month_factor+san_joaquin_flow_z, data=data_STH)
Cand.set.STH[[6]] <-  glm.nb(steelhead_loss~month_factor+export_z, data=data_STH)
Cand.set.STH[[7]] <-  glm.nb(steelhead_loss~san_joaquin_flow_z+export_z, data=data_STH)
Cand.set.STH[[8]] <-  glm.nb(steelhead_loss~month_factor*san_joaquin_flow_z, data=data_STH)
Cand.set.STH[[9]] <-  glm.nb(steelhead_loss~month_factor*san_joaquin_flow_z+export_z, data=data_STH)
Cand.set.STH[[10]] <-  glm.nb(steelhead_loss~month_factor*export_z, data=data_STH)
Cand.set.STH[[11]] <-  glm.nb(steelhead_loss~month_factor*export_z+san_joaquin_flow_z, data=data_STH)
Cand.set.STH[[12]] <-  glm.nb(steelhead_loss~export_z*san_joaquin_flow_z, data=data_STH)
Cand.set.STH[[13]] <-  glm.nb(steelhead_loss~export_z*san_joaquin_flow_z+month_factor, data=data_STH)
Cand.set.STH[[14]] <-  glm.nb(steelhead_loss~export_z*san_joaquin_flow_z*month_factor, data=data_STH)


##create a vector of names to trace back models in set
Modnames <- paste("mod","Steelhead", 1:length(Cand.set.STH), sep = "_")

##generate AICc table
aictab(cand.set = Cand.set.STH, modnames = Modnames, sort = TRUE)

#Best model
summary(Cand.set.STH[[6]])
#2nd best model
summary(Cand.set.STH[[13]])

# Model coefficients
est <- cbind(Estimate = coef(Cand.set.STH[[6]]), confint(Cand.set.STH[[6]]))
exp(est)

# Model residuals

data_STH$residuals <- resid(Cand.set.STH[[6]])

hist(data_STH$residuals, breaks = 20)

# Looping prediction of best model

Pred.dat <- data.frame(Observed.dat = rep(NA, nrow(data_STH)), Pred.NegBin = rep(NA, nrow(data_STH)))
head(Pred.dat)

for(i in 1:nrow(data_STH)){
  Hat.dat <- data_STH[-i, ]
  CV.dat <- data_STH[i, ]
  Temp.mod1 <- glm.nb(steelhead_loss~month_factor+export_z, data=Hat.dat)
  Pred.dat[i, "Observed.dat"] <- CV.dat[, "steelhead_loss"]
  Pred.dat[i, "Pred.NegBin"] <- predict(Temp.mod1, newdata = CV.dat, type="response")
}

head(Pred.dat)

plot(Pred.dat$Observed.dat, Pred.dat$Pred.NegBin)

loocv.lm <- lm(log(Pred.NegBin+1)~log(Observed.dat+1), data=Pred.dat)
summary(loocv.lm)
#Multiple R-squared:  0.4756,	Adjusted R-squared:  0.4728 

plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.NegBin), xlim=c(0, 9), ylim=c(0, 9))
abline(loocv.lm, col="red")

###########Try 2nd best model
# Model residuals

data_STH$residuals <- resid(Cand.set.STH[[13]])

hist(data_STH$residuals, breaks = 20)

# Looping prediction of best model

Pred.dat <- data.frame(Observed.dat = rep(NA, nrow(data_STH)), Pred.NegBin = rep(NA, nrow(data_STH)))
head(Pred.dat)

for(i in 1:nrow(data_STH)){
  Hat.dat <- data_STH[-i, ]
  CV.dat <- data_STH[i, ]
  Temp.mod1 <- glm.nb(steelhead_loss~export_z*san_joaquin_flow_z+month_factor, data=Hat.dat)
  Pred.dat[i, "Observed.dat"] <- CV.dat[, "steelhead_loss"]
  Pred.dat[i, "Pred.NegBin"] <- predict(Temp.mod1, newdata = CV.dat, type="response")
}

head(Pred.dat)

plot(Pred.dat$Observed.dat, Pred.dat$Pred.NegBin)

loocv.lm <- lm(log(Pred.NegBin+1)~log(Observed.dat+1), data=Pred.dat)
summary(loocv.lm)
#Multiple R-squared:  0.4813,	Adjusted R-squared:  0.4786 

plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.NegBin), xlim=c(0, 9), ylim=c(0, 9))
abline(loocv.lm, col="red")

# Redo model with standard (original) values and ensure that results are essentially the same
model_sth_final<-glm.nb(steelhead_loss~month_factor+export, data=data_STH)

# Saving final spring-run model
saveRDS(model_sth_final, file = "model_steelhead_final.rds")