setwd("~/GitHub/salmon_negbinmodel")

library("MASS")
library("MuMIn")
library("pscl")
library("tidyverse")
library("AICcmodavg")

monthly_data <- read.csv("negbinmodel_monthly_dataset.csv")

str(monthly_data)
#Preliminary assessment for month to use
month_summary<- monthly_data %>% group_by(month) %>% summarise(mean_winter_lad_loss=mean(winter_lad_loss),sd_winter_lad_loss=sd(winter_lad_loss),
                                                               mean_spring_lad_loss=mean(spring_lad_loss),sd_spring_lad_loss=sd(spring_lad_loss))
#Based on the values on the summary table, we should use:
#December-April for winter-run
#March-June for Spring-run

#Data for winter-run
data_WR<-monthly_data %>% filter(month %in% c(12,1,2,3,4)) %>%
  #standardize data by z-score and add Sac Trawl CPUE, also convert month to factor
  mutate(sac_flow_z = scale(sac_flow),
         san_joaquin_flow_z = scale(san_joaquin_flow),
         export_z = scale(export),
         delta_outflow_z = scale(delta_outflow),
         omr_flow_extrap_z = scale(omr_flow_extrap),
         sac_trawl_wr_CPUE = sac_trawl_wr_count/sac_trawl_sample_size,
         sac_trawl_wr_CPUE_z = scale(sac_trawl_wr_CPUE),
         month_factor=as.factor(month))
hist(data_WR$winter_lad_loss)

#Data for spring-run
data_SR<-monthly_data %>% filter(month %in% c(3,4,5,6)) %>%
  #standardize data by z-score and add Sac Trawl CPUE
  mutate(sac_flow_z = scale(sac_flow),
         san_joaquin_flow_z = scale(san_joaquin_flow),
         export_z = scale(export),
         delta_outflow_z = scale(delta_outflow),
         omr_flow_extrap_z = scale(omr_flow_extrap),
         sac_trawl_sr_CPUE = sac_trawl_sr_count/sac_trawl_sample_size,
         sac_trawl_sr_CPUE_z = scale(sac_trawl_sr_CPUE),
         month_factor=as.factor(month)
  )

hist(data_SR$spring_lad_loss)

###################DO VIF Analysis
str(data_WR)

## Winter-run model VIF screening
wr_full_model<-glm.nb(winter_lad_loss~ sac_flow_z + san_joaquin_flow_z + export_z + omr_flow_extrap_z + sac_trawl_wr_CPUE_z, data=data_WR)
car::vif(wr_full_model)
#sac_flow_z  san_joaquin_flow_z            export_z   omr_flow_extrap_z sac_trawl_wr_CPUE_z 
#3.104473           26.271862           10.995967           33.578451            1.347934 

#Removed OMR flow due to excess collinearity (VIF > 10)

wr_full_model<-glm.nb(winter_lad_loss~ sac_flow_z + san_joaquin_flow_z + export_z + sac_trawl_wr_CPUE_z, data=data_WR)
car::vif(wr_full_model)
#sac_flow_z  san_joaquin_flow_z            export_z sac_trawl_wr_CPUE_z 
#3.081247            2.965193            1.335696            1.243651 

## Spring-run model VIF screening
sr_full_model<-glm.nb(spring_lad_loss~ sac_flow_z + san_joaquin_flow_z + export_z + omr_flow_extrap_z + sac_trawl_sr_CPUE_z, data=data_SR)
car::vif(sr_full_model)

#Removed OMR flow due to excess collinearity (VIF > 10)

sr_full_model<-glm.nb(spring_lad_loss~ sac_flow_z + san_joaquin_flow_z + export_z + sac_trawl_sr_CPUE_z, data=data_SR)
car::vif(sr_full_model)
#sac_flow_z  san_joaquin_flow_z            export_z sac_trawl_sr_CPUE_z 
#3.293633            3.078356            1.129226            1.092105 

#Plot correlation matrix just for winter-run
##
library("PerformanceAnalytics")

chart.Correlation(data_WR %>% select(winter_lad_loss, sac_flow_z, san_joaquin_flow_z, export_z, sac_trawl_wr_CPUE_z), histogram=TRUE, pch=19)
#0.77 correlation for Sac and SJ flows

#We will go with 3 VIF threshold instead then
wr_full_model<-glm.nb(winter_lad_loss~ san_joaquin_flow_z + export_z + sac_trawl_wr_CPUE_z, data=data_WR)
car::vif(wr_full_model)
#san_joaquin_flow_z   export_z sac_trawl_wr_CPUE_z 
#1.015942            1.174645            1.184954
sr_full_model<-glm.nb(spring_lad_loss~ san_joaquin_flow_z + export_z + sac_trawl_sr_CPUE_z, data=data_SR)
car::vif(sr_full_model)
#san_joaquin_flow_z   export_z  sac_trawl_sr_CPUE_z 
#1.100825   1.059054 1.070370 

# VIF values for both models look good ~1 after using VIF threshold of 3 instead

##############################
#Model selection for winter-run

Cand.set.WR <- list( )

Cand.set.WR[[1]] <-  glm.nb(winter_lad_loss~NULL, data=data_WR)
Cand.set.WR[[2]] <-  glm.nb(winter_lad_loss~month_factor, data=data_WR)
Cand.set.WR[[3]] <-  glm.nb(winter_lad_loss~san_joaquin_flow_z, data=data_WR)
Cand.set.WR[[4]] <-  glm.nb(winter_lad_loss~export_z, data=data_WR)
Cand.set.WR[[5]] <-  glm.nb(winter_lad_loss~sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[6]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z, data=data_WR)
Cand.set.WR[[7]] <-  glm.nb(winter_lad_loss~month_factor+export_z, data=data_WR)
Cand.set.WR[[8]] <-  glm.nb(winter_lad_loss~month_factor+sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[9]] <-  glm.nb(winter_lad_loss~san_joaquin_flow_z+export_z, data=data_WR)
Cand.set.WR[[10]] <-  glm.nb(winter_lad_loss~san_joaquin_flow_z+sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[11]] <-  glm.nb(winter_lad_loss~export_z+sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[12]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+export_z, data=data_WR)
Cand.set.WR[[13]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[14]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+sac_trawl_wr_CPUE_z+export_z, data=data_WR)
Cand.set.WR[[15]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z, data=data_WR)
Cand.set.WR[[16]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z+export_z, data=data_WR)
Cand.set.WR[[17]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z+sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[18]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z+sac_trawl_wr_CPUE_z+export_z, data=data_WR)
Cand.set.WR[[19]] <-  glm.nb(winter_lad_loss~month_factor*export_z, data=data_WR)
Cand.set.WR[[20]] <-  glm.nb(winter_lad_loss~month_factor*export_z+sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[21]] <-  glm.nb(winter_lad_loss~month_factor*export_z+san_joaquin_flow_z, data=data_WR)
Cand.set.WR[[22]] <-  glm.nb(winter_lad_loss~month_factor*export_z+sac_trawl_wr_CPUE_z+san_joaquin_flow_z, data=data_WR)
Cand.set.WR[[23]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_wr_CPUE_z, data=data_WR)
Cand.set.WR[[24]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_wr_CPUE_z+san_joaquin_flow_z, data=data_WR)
Cand.set.WR[[25]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_wr_CPUE_z+export_z, data=data_WR)
Cand.set.WR[[26]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_wr_CPUE_z+export_z+san_joaquin_flow_z, data=data_WR)


##create a vector of names to trace back models in set
Modnames <- paste("mod","WinterRun", 1:length(Cand.set.WR), sep = "_")

##generate AICc table
aictab(cand.set = Cand.set.WR, modnames = Modnames, sort = TRUE)

summary(Cand.set.WR[[14]])
summary(Cand.set.WR[[26]])

# Model coefficients
est <- cbind(Estimate = coef(Cand.set.WR[[14]]), confint(Cand.set.WR[[14]]))
exp(est)

# Model residuals
data_WR$residuals <- resid(Cand.set.WR[[14]])

hist(data_WR$residuals, breaks = 20)

# Looping prediction of best model

Pred.dat <- data.frame(Observed.dat = rep(NA, nrow(data_WR)), Pred.NegBin = rep(NA, nrow(data_WR)))
head(Pred.dat)

for(i in 1:nrow(data_WR)){
  Hat.dat <- data_WR[-i, ]
  CV.dat <- data_WR[i, ]
  Temp.mod1 <- glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+sac_trawl_wr_CPUE_z+export_z, data=Hat.dat)
  Pred.dat[i, "Observed.dat"] <- CV.dat[, "winter_lad_loss"]
  Pred.dat[i, "Pred.NegBin"] <- predict(Temp.mod1, newdata = CV.dat, type="response")
}

head(Pred.dat)

plot(Pred.dat$Observed.dat, Pred.dat$Pred.NegBin)

loocv.lm <- lm(log(Pred.NegBin+1)~log(Observed.dat+1), data=Pred.dat)
summary(loocv.lm)
#Multiple R-squared:  0.5108,	Adjusted R-squared:  0.5073 

plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.NegBin), xlim=c(0, 9), ylim=c(0, 9))
abline(loocv.lm, col="red")

# Redo model with standard (original) values and ensure that results are essentially the same
model_wr_final<-glm.nb(winter_lad_loss~month_factor+san_joaquin_flow+sac_trawl_wr_CPUE+export, data=data_WR)

##############################
#Model selection for spring-run

Cand.set.SR <- list( )

Cand.set.SR[[1]] <-  glm.nb(winter_lad_loss~NULL, data=data_SR)
Cand.set.SR[[2]] <-  glm.nb(winter_lad_loss~month_factor, data=data_SR)
Cand.set.SR[[3]] <-  glm.nb(winter_lad_loss~san_joaquin_flow_z, data=data_SR)
Cand.set.SR[[4]] <-  glm.nb(winter_lad_loss~export_z, data=data_SR)
Cand.set.SR[[5]] <-  glm.nb(winter_lad_loss~sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[6]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z, data=data_SR)
Cand.set.SR[[7]] <-  glm.nb(winter_lad_loss~month_factor+export_z, data=data_SR)
Cand.set.SR[[8]] <-  glm.nb(winter_lad_loss~month_factor+sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[9]] <-  glm.nb(winter_lad_loss~san_joaquin_flow_z+export_z, data=data_SR)
Cand.set.SR[[10]] <-  glm.nb(winter_lad_loss~san_joaquin_flow_z+sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[11]] <-  glm.nb(winter_lad_loss~export_z+sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[12]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+export_z, data=data_SR)
Cand.set.SR[[13]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[14]] <-  glm.nb(winter_lad_loss~month_factor+san_joaquin_flow_z+sac_trawl_sr_CPUE_z+export_z, data=data_SR)
Cand.set.SR[[15]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z, data=data_SR)
Cand.set.SR[[16]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z+export_z, data=data_SR)
Cand.set.SR[[17]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z+sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[18]] <-  glm.nb(winter_lad_loss~month_factor*san_joaquin_flow_z+sac_trawl_sr_CPUE_z+export_z, data=data_SR)
Cand.set.SR[[19]] <-  glm.nb(winter_lad_loss~month_factor*export_z, data=data_SR)
Cand.set.SR[[20]] <-  glm.nb(winter_lad_loss~month_factor*export_z+sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[21]] <-  glm.nb(winter_lad_loss~month_factor*export_z+san_joaquin_flow_z, data=data_SR)
Cand.set.SR[[22]] <-  glm.nb(winter_lad_loss~month_factor*export_z+sac_trawl_sr_CPUE_z+san_joaquin_flow_z, data=data_SR)
Cand.set.SR[[23]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_sr_CPUE_z, data=data_SR)
Cand.set.SR[[24]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_sr_CPUE_z+san_joaquin_flow_z, data=data_SR)
Cand.set.SR[[25]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_sr_CPUE_z+export_z, data=data_SR)
Cand.set.SR[[26]] <-  glm.nb(winter_lad_loss~month_factor*sac_trawl_sr_CPUE_z+export_z+san_joaquin_flow_z, data=data_SR)


##create a vector of names to trace back models in set
Modnames <- paste("mod","SpringRun", 1:length(Cand.set.SR), sep = "_")

##generate AICc table
aictab(cand.set = Cand.set.SR, modnames = Modnames, sort = TRUE)

#Cand.set.SR[[7]] <-  glm.nb(winter_lad_loss~month_factor+export_z, data=data_SR)
