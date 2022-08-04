
rm(list=ls())

library("MASS")
library("MuMIn")
library("pscl")
library("dplyr")
library("DHARMa")

setwd("~/R for ROC/NegBin_Salvage")


Salv.dat <- read.csv("Salvage_dat_3Mar2022.csv")
head(Salv.dat)


hist(Salv.dat$Salvage.Total, breaks = 100)


################# Fitting models to fall-run Salmon ###################

Fall.salv <- subset(Salv.dat, Run == "Fall")

head(Fall.salv)
nrow(Fall.salv)
hist(Fall.salv$Salvage.Total, breaks = 100)


Salv.mod <- glm(Salvage.Total ~ Month , family = "poisson", data = Fall.salv)
Salv.sim <- simulateResiduals(fittedModel = Salv.mod)
testZeroInflation(Salv.sim)


############### Fitting Zero Inflated Model to fall-run Salmon Data ####################

Mod.ZeroInf0 <- zeroinfl(Salvage.Total~Month | OMR, dist="negbin", data=Fall.salv)

############### Fitting NegBin Model to fall-run Salmon Data ####################

Mod.NegBin0 <- glm.nb(Salvage.Total~Month, data=Fall.salv)

########## Compare negbin vs zinfl #########
?vuong
vuong(Mod.ZeroInf0, Mod.NegBin0)


