
rm(list=ls())

library("MASS")
library("MuMIn")
library("pscl")
library("dplyr")

setwd("~/R for ROC/NegBin_Salvage")


Salv.dat <- read.csv("Salvage.csv")
head(Salv.dat)
Salv.dat$Year
max(Salv.dat$Year)

hist(Salv.dat$Salvage1, breaks = 100)
hist(Salv.dat$Salvage2)
plot(Salv.dat$Salvage1, Salv.dat$Salvage2)
unique(Salv.dat$Month)
unique(Salv.dat$Run)

Month.code <- data.frame(Month.num=seq(1, 12, 1), Month = c("January","February","March","April","May","June","July","August","September","October","November","December"))

Discharge.dat <- read.csv("Flow_Verona.csv")
Discharge.dat <- merge(Discharge.dat, Month.code, by="Month.num")
head(Discharge.dat)
Discharge.dat <- Discharge.dat[,c("Year", "Month", "Discharge.cfs")]

Salv.dat <- merge(Salv.dat, Discharge.dat, all.x=TRUE)
head(Salv.dat)

################# Fitting models to all Salmon ###################

head(Salv.dat)

Salv.dat$TPR[is.na(Salv.dat$TPR)] <- 0
Salv.dat$HRO[is.na(Salv.dat$HRO)] <- 0

Salv.dat2 <- transform(Salv.dat, Total.Exports = TPR+HRO, Salvage.Total = round(Salvage1+Salvage2, 0))

Salv.dat2$TPR
Salv.dat2$HRO
Salv.dat2$Total.Exports

head(Salv.dat2)

## Subsetting for OMR > -5000 ##
min(Salv.dat2$OMR, na.rm=TRUE)
Salv.dat2 <- subset(Salv.dat2, OMR >= -5000)
min(Salv.dat2$OMR, na.rm=TRUE)

## Subsetting for OMR > -3500 ##
min(Salv.dat2$OMR, na.rm=TRUE)
Salv.dat2 <- subset(Salv.dat2, OMR >= -3500)
min(Salv.dat2$OMR, na.rm=TRUE)

Salv.dat2 <- transform(Salv.dat2, Total.exports.z = scale(Total.Exports), OMR.z = scale(OMR), Discharge.z = scale(Discharge.cfs))
head(Salv.dat2)

Salv.dat3 <- Salv.dat2[, c("Year", "Month", "Run", "Salvage.Total", "Total.Exports", "OMR", "Total.exports.z", "OMR.z", "Discharge.z")]
Salv.dat3 <- na.omit(Salv.dat3)

#Salv.dat3 <- subset(Salv.dat3, Month %in% c("March", "April", "May"))
#Salv.dat3 <- subset(Salv.dat3, Month %in% c("February", "March", "April", "May", "June"))
Salv.WR <- subset(Salv.dat3, Month %in% c("December", "January", "February", "March", "April"))

nrow(Salv.WR)
head(Salv.WR)
unique(Salv.WR$Month)
unique(Salv.WR$Year)


#Salv.dat3.sum <- Salv.dat3 %>%
#  group_by(Year, Run) %>%
#  summarise(Total.Salvage = sum(Salvage.Total)) %>%
#  as.data.frame(Salv.dat3.sum)
#Salv.dat3.sum <- Salv.dat3.sum[order(Salv.dat3.sum$Run),]
#Salv.dat3.sum

#hist(Salv.dat3$Salvage.Total, breaks = 100)

#Salv.dat3.NoLF <- subset(Salv.dat3, Run !="LateFall")

#Salv.dat4 <- Salv.dat3.NoLF %>%
#  group_by(Year, Month) %>%
#  summarise(Salvage.Total = sum(Salvage.Total), Total.exports.z = unique(Total.exports.z), OMR.z = unique(OMR.z), Discharge.z = unique(Discharge.z)) %>%
#  as.data.frame(Salv.dat4)

#head(Salv.dat4)
#hist(Salv.dat4$Salvage.Total, breaks = 50)

##################################################################
##################################################################


#################### Winter Run ############################
head(Salv.WR)
min(Salv.WR$Year)
max(Salv.WR$Year)

Winter.salv <- subset(Salv.WR, Run == "Winter")
head(Winter.salv)
unique(Winter.salv$Month)

nrow(Winter.salv)

#jpeg(filename = "Salvage_Frequency.jpg", heigh=5, width = 6, units="in", bg="white", res=500)

par(mfrow=c(1,1), oma=c(3, 3.5, 0.5, 0.5), mar=c(0.5, 0.75, 0.5, 0.25), bty="l")

hist(Winter.salv$Salvage.Total, breaks = 100, las=1, xlab="Winter run salvage(n)", ylab="Frequency", main="")
mtext(1, text="Monthly Salvage (n)", line = 2.5)
mtext(2, text="Frequency", line = 3)

#dev.off()

Winter.salv <- transform(Winter.salv, P.A = ifelse(Salvage.Total > 0, 1, 0))

#1-(sum(Winter.salv$P.A)/nrow(Winter.salv))
# Note only 1.9% zeros so negbin should be fine #

############# Merging Escapement Data with Winter-run salvage #################

########### Escapement data #############

Escape.dat <- read.csv("Escapement.csv")
head(Escape.dat)

Winter.salv <- merge(Winter.salv, Escape.dat, by=c("Year", "Month"), all.x=TRUE)
Winter.salv <- transform(Winter.salv, Escape.z = scale(Total.Escapement))
head(Winter.salv)


############################ Clean Candidate Model Set #############################

# 1 Parameter Models
Mod.NegBin0 <- glm.nb(Salvage.Total~NULL, data=Winter.salv)
Mod.NegBin1 <- glm.nb(Salvage.Total~Month, data=Winter.salv)
Mod.NegBin2 <- glm.nb(Salvage.Total~Discharge.z, data=Winter.salv)
Mod.NegBin3 <- glm.nb(Salvage.Total~OMR.z, data=Winter.salv)
Mod.NegBin4 <- glm.nb(Salvage.Total~Total.exports.z, data=Winter.salv)
#Mod.NegBin5 <- glm.nb(Salvage.Total~Escape.z, data=Winter.salv)

# 2 Parameter Models
Mod.NegBin6 <- glm.nb(Salvage.Total~Month + Discharge.z, data=Winter.salv)
Mod.NegBin7 <- glm.nb(Salvage.Total~Month + OMR.z, data=Winter.salv)
Mod.NegBin8 <- glm.nb(Salvage.Total~Month + Total.exports.z, data=Winter.salv)
Mod.NegBin9 <- glm.nb(Salvage.Total~Month + Escape.z, data=Winter.salv)
Mod.NegBin10 <- glm.nb(Salvage.Total~OMR.z + Discharge.z, data=Winter.salv)
Mod.NegBin11 <- glm.nb(Salvage.Total~OMR.z + Total.exports.z, data=Winter.salv)
Mod.NegBin12 <- glm.nb(Salvage.Total~OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin13 <- glm.nb(Salvage.Total~Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin14 <- glm.nb(Salvage.Total~Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin15 <- glm.nb(Salvage.Total~Discharge.z + Escape.z, data=Winter.salv)

# 3 Parameter Models
Mod.NegBin16 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z, data=Winter.salv)
Mod.NegBin17 <- glm.nb(Salvage.Total~Month + Discharge.z + Total.exports.z, data=Winter.salv)
Mod.NegBin18 <- glm.nb(Salvage.Total~Month + Discharge.z + OMR.z, data=Winter.salv)
Mod.NegBin19 <- glm.nb(Salvage.Total~Month + Escape.z + OMR.z, data=Winter.salv)
Mod.NegBin20 <- glm.nb(Salvage.Total~Month + Escape.z + Total.exports.z, data=Winter.salv)
Mod.NegBin21 <- glm.nb(Salvage.Total~Month + Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin22 <- glm.nb(Salvage.Total~Total.exports.z + OMR.z + Discharge.z, data=Winter.salv)
Mod.NegBin23 <- glm.nb(Salvage.Total~Total.exports.z + OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin24 <- glm.nb(Salvage.Total~Discharge.z + OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin25 <- glm.nb(Salvage.Total~Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)

# 4 Parameter Models
Mod.NegBin26 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin27 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin28 <- glm.nb(Salvage.Total~Month + OMR.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin29 <- glm.nb(Salvage.Total~Month + Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin30 <- glm.nb(Salvage.Total~OMR.z + Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)

# 5 Parameter Models
Mod.NegBin31 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)

# 2 Parameter Models + 1 interaction
Mod.NegBin32 <- glm.nb(Salvage.Total~Month*Discharge.z, data=Winter.salv)
Mod.NegBin33 <- glm.nb(Salvage.Total~Month*OMR.z, data=Winter.salv)
Mod.NegBin34 <- glm.nb(Salvage.Total~Month*Total.exports.z, data=Winter.salv)
Mod.NegBin35 <- glm.nb(Salvage.Total~Month*Escape.z, data=Winter.salv)
Mod.NegBin36 <- glm.nb(Salvage.Total~OMR.z*Discharge.z, data=Winter.salv)
Mod.NegBin37 <- glm.nb(Salvage.Total~OMR.z*Total.exports.z, data=Winter.salv)
Mod.NegBin38 <- glm.nb(Salvage.Total~OMR.z*Escape.z, data=Winter.salv)
Mod.NegBin39 <- glm.nb(Salvage.Total~Total.exports.z*Discharge.z, data=Winter.salv)
Mod.NegBin40 <- glm.nb(Salvage.Total~Total.exports.z*Escape.z, data=Winter.salv)
Mod.NegBin41 <- glm.nb(Salvage.Total~Discharge.z*Escape.z, data=Winter.salv)

# 3 Parameter Models + 1 interaction
Mod.NegBin42 <- glm.nb(Salvage.Total~Month*OMR.z + Total.exports.z, data=Winter.salv)
Mod.NegBin43 <- glm.nb(Salvage.Total~Month*Total.exports.z + OMR.z, data=Winter.salv)
Mod.NegBin44 <- glm.nb(Salvage.Total~Month + OMR.z*Total.exports.z, data=Winter.salv)
Mod.NegBin45 <- glm.nb(Salvage.Total~Month*Discharge.z + Total.exports.z, data=Winter.salv)
Mod.NegBin46 <- glm.nb(Salvage.Total~Month*Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin47 <- glm.nb(Salvage.Total~Month + Discharge.z*Total.exports.z, data=Winter.salv)
Mod.NegBin48 <- glm.nb(Salvage.Total~Month*Discharge.z + OMR.z, data=Winter.salv)
Mod.NegBin49 <- glm.nb(Salvage.Total~Month*OMR.z + Discharge.z, data=Winter.salv)
#Mod.NegBin50 <- glm.nb(Salvage.Total~Month + Discharge.z*OMR.z, data=Winter.salv)
Mod.NegBin51 <- glm.nb(Salvage.Total~Month*Escape.z + OMR.z, data=Winter.salv)
Mod.NegBin52 <- glm.nb(Salvage.Total~Month*OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin53 <- glm.nb(Salvage.Total~Month + Escape.z*OMR.z, data=Winter.salv)
Mod.NegBin54 <- glm.nb(Salvage.Total~Month*Escape.z + Total.exports.z, data=Winter.salv)
Mod.NegBin55 <- glm.nb(Salvage.Total~Month*Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin56 <- glm.nb(Salvage.Total~Month + Escape.z*Total.exports.z, data=Winter.salv)
Mod.NegBin57 <- glm.nb(Salvage.Total~Month*Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin58 <- glm.nb(Salvage.Total~Month*Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin59 <- glm.nb(Salvage.Total~Month + Escape.z*Discharge.z, data=Winter.salv)
Mod.NegBin60 <- glm.nb(Salvage.Total~Total.exports.z*OMR.z + Discharge.z, data=Winter.salv)
Mod.NegBin61 <- glm.nb(Salvage.Total~Total.exports.z*Discharge.z + OMR.z, data=Winter.salv)
#Mod.NegBin62 <- glm.nb(Salvage.Total~Total.exports.z + OMR.z*Discharge.z, data=Winter.salv)
Mod.NegBin63 <- glm.nb(Salvage.Total~Total.exports.z*OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin64 <- glm.nb(Salvage.Total~Total.exports.z*Escape.z + OMR.z, data=Winter.salv)
Mod.NegBin65 <- glm.nb(Salvage.Total~Total.exports.z + OMR.z*Escape.z, data=Winter.salv)
Mod.NegBin66 <- glm.nb(Salvage.Total~Discharge.z*OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin67 <- glm.nb(Salvage.Total~Discharge.z*Escape.z + OMR.z, data=Winter.salv)
Mod.NegBin68 <- glm.nb(Salvage.Total~Discharge.z + OMR.z*Escape.z, data=Winter.salv)
Mod.NegBin69 <- glm.nb(Salvage.Total~Total.exports.z*Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin70 <- glm.nb(Salvage.Total~Total.exports.z*Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin71 <- glm.nb(Salvage.Total~Total.exports.z + Discharge.z*Escape.z, data=Winter.salv)

# 4 Parameter Models + 1 Interaction
Mod.NegBin72 <- glm.nb(Salvage.Total~Month*OMR.z + Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin73 <- glm.nb(Salvage.Total~Month*Total.exports.z + OMR.z + Discharge.z, data=Winter.salv)
Mod.NegBin74 <- glm.nb(Salvage.Total~Month*Discharge.z + OMR.z + Total.exports.z, data=Winter.salv)
Mod.NegBin75 <- glm.nb(Salvage.Total~Month + OMR.z*Total.exports.z + Discharge.z, data=Winter.salv)
#Mod.NegBin76 <- glm.nb(Salvage.Total~Month + OMR.z*Discharge.z + Total.exports.z, data=Winter.salv)
Mod.NegBin77 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z*Discharge.z, data=Winter.salv)
Mod.NegBin78 <- glm.nb(Salvage.Total~Month*OMR.z + Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin79 <- glm.nb(Salvage.Total~Month*Total.exports.z + OMR.z + Escape.z, data=Winter.salv)
#Mod.NegBin80 <- glm.nb(Salvage.Total~Month*Escape.z + OMR.z + Total.exports.z, data=Winter.salv)
Mod.NegBin81 <- glm.nb(Salvage.Total~Month + OMR.z*Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin82 <- glm.nb(Salvage.Total~Month + OMR.z*Escape.z + Total.exports.z, data=Winter.salv)
#Mod.NegBin83 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z*Escape.z, data=Winter.salv)
Mod.NegBin84 <- glm.nb(Salvage.Total~Month*OMR.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin85 <- glm.nb(Salvage.Total~Month*Discharge.z + OMR.z + Escape.z, data=Winter.salv)
Mod.NegBin86 <- glm.nb(Salvage.Total~Month*Escape.z + OMR.z + Discharge.z, data=Winter.salv)
#Mod.NegBin87 <- glm.nb(Salvage.Total~Month + OMR.z*Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin88 <- glm.nb(Salvage.Total~Month + OMR.z*Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin89 <- glm.nb(Salvage.Total~Month + OMR.z + Discharge.z*Escape.z, data=Winter.salv)
Mod.NegBin90 <- glm.nb(Salvage.Total~Month*Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin91 <- glm.nb(Salvage.Total~Month*Discharge.z + Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin92 <- glm.nb(Salvage.Total~Month*Escape.z + Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin93 <- glm.nb(Salvage.Total~Month + Total.exports.z*Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin94 <- glm.nb(Salvage.Total~Month + Total.exports.z*Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin95 <- glm.nb(Salvage.Total~Month + Total.exports.z + Discharge.z*Escape.z, data=Winter.salv)
Mod.NegBin96 <- glm.nb(Salvage.Total~OMR.z*Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin97 <- glm.nb(Salvage.Total~OMR.z*Discharge.z + Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin98 <- glm.nb(Salvage.Total~OMR.z*Escape.z + Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin99 <- glm.nb(Salvage.Total~OMR.z + Total.exports.z*Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin100 <- glm.nb(Salvage.Total~OMR.z + Total.exports.z*Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin101 <- glm.nb(Salvage.Total~OMR.z + Total.exports.z + Discharge.z*Escape.z, data=Winter.salv)

# 5 Parameter Models + 1 Interaction
Mod.NegBin102 <- glm.nb(Salvage.Total~Month*OMR.z + Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin103 <- glm.nb(Salvage.Total~Month*Total.exports.z + OMR.z + Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin104 <- glm.nb(Salvage.Total~Month*Discharge.z + OMR.z + Total.exports.z + Escape.z, data=Winter.salv)
#Mod.NegBin105 <- glm.nb(Salvage.Total~Month*Escape.z + OMR.z + Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin106 <- glm.nb(Salvage.Total~Month + OMR.z*Total.exports.z + Discharge.z + Escape.z, data=Winter.salv)
#Mod.NegBin107 <- glm.nb(Salvage.Total~Month + OMR.z*Discharge.z + Total.exports.z + Escape.z, data=Winter.salv)
Mod.NegBin108 <- glm.nb(Salvage.Total~Month + OMR.z*Escape.z + Total.exports.z + Discharge.z, data=Winter.salv)
Mod.NegBin109 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z*Discharge.z + Escape.z, data=Winter.salv)
Mod.NegBin110 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z*Escape.z + Discharge.z, data=Winter.salv)
Mod.NegBin111 <- glm.nb(Salvage.Total~Month + OMR.z + Total.exports.z + Discharge.z*Escape.z, data=Winter.salv)


#Model.list <- data.frame(Base.model = rep("Mod.NegBin", 69), Number = c(seq(1, 49, 1), seq(51, 61, 1), seq(63, 71, 1)))
#Model.list <- transform(Model.list, Model.number = paste(Base.model, Number, sep=""))

Mod.results  <- as.data.frame(AICc(Mod.NegBin0, Mod.NegBin1, Mod.NegBin2, Mod.NegBin3, Mod.NegBin4, Mod.NegBin5, Mod.NegBin6, Mod.NegBin7, Mod.NegBin8, Mod.NegBin9, Mod.NegBin10, Mod.NegBin11, Mod.NegBin12, Mod.NegBin13, Mod.NegBin14, Mod.NegBin15, Mod.NegBin16, Mod.NegBin17, Mod.NegBin18, Mod.NegBin19, Mod.NegBin20, Mod.NegBin21, Mod.NegBin22, Mod.NegBin23, Mod.NegBin24, Mod.NegBin25, Mod.NegBin26, Mod.NegBin27, Mod.NegBin28, Mod.NegBin29, Mod.NegBin30, Mod.NegBin31, Mod.NegBin32, Mod.NegBin33, Mod.NegBin34, Mod.NegBin35, Mod.NegBin36, Mod.NegBin37, Mod.NegBin38, Mod.NegBin39, Mod.NegBin40, Mod.NegBin41, Mod.NegBin42, Mod.NegBin43, Mod.NegBin44, Mod.NegBin45, Mod.NegBin46, Mod.NegBin47, Mod.NegBin48, Mod.NegBin49, Mod.NegBin51, Mod.NegBin52, Mod.NegBin53, Mod.NegBin54, Mod.NegBin55, Mod.NegBin56, Mod.NegBin57, Mod.NegBin58, Mod.NegBin59, Mod.NegBin60, Mod.NegBin61, Mod.NegBin63, Mod.NegBin64, Mod.NegBin65,Mod.NegBin66, Mod.NegBin67, Mod.NegBin68, Mod.NegBin69, Mod.NegBin70, Mod.NegBin71, Mod.NegBin72, Mod.NegBin73, Mod.NegBin74, Mod.NegBin75, Mod.NegBin77, Mod.NegBin78, Mod.NegBin79, Mod.NegBin81, Mod.NegBin82, Mod.NegBin84, Mod.NegBin85, Mod.NegBin86, Mod.NegBin88, Mod.NegBin89, Mod.NegBin90, Mod.NegBin91, Mod.NegBin92, Mod.NegBin93, Mod.NegBin94, Mod.NegBin95, Mod.NegBin96, Mod.NegBin97, Mod.NegBin98, Mod.NegBin99, Mod.NegBin100, Mod.NegBin101, Mod.NegBin102, Mod.NegBin103, Mod.NegBin104, Mod.NegBin106, Mod.NegBin108, Mod.NegBin109, Mod.NegBin110, Mod.NegBin111))

######## OMR -5000 ####################
Mod.results  <- as.data.frame(AICc(Mod.NegBin0, Mod.NegBin1, Mod.NegBin2, Mod.NegBin3, Mod.NegBin4, Mod.NegBin5, Mod.NegBin6, Mod.NegBin7, Mod.NegBin8, Mod.NegBin9, Mod.NegBin10, Mod.NegBin11, Mod.NegBin12, Mod.NegBin13, Mod.NegBin14, Mod.NegBin15, Mod.NegBin16, Mod.NegBin17, Mod.NegBin18, Mod.NegBin19, Mod.NegBin20, Mod.NegBin21, Mod.NegBin22, Mod.NegBin23, Mod.NegBin24, Mod.NegBin25, Mod.NegBin26, Mod.NegBin27, Mod.NegBin28, Mod.NegBin29, Mod.NegBin30, Mod.NegBin31, Mod.NegBin32, Mod.NegBin33, Mod.NegBin34, Mod.NegBin35, Mod.NegBin36, Mod.NegBin37, Mod.NegBin38, Mod.NegBin39, Mod.NegBin40, Mod.NegBin41, Mod.NegBin43, Mod.NegBin44, Mod.NegBin45, Mod.NegBin46, Mod.NegBin47, Mod.NegBin48, Mod.NegBin49, Mod.NegBin51, Mod.NegBin52, Mod.NegBin53, Mod.NegBin54, Mod.NegBin55, Mod.NegBin56, Mod.NegBin57, Mod.NegBin58, Mod.NegBin59, Mod.NegBin60, Mod.NegBin61, Mod.NegBin63, Mod.NegBin64, Mod.NegBin65,Mod.NegBin66, Mod.NegBin67, Mod.NegBin68, Mod.NegBin69, Mod.NegBin70, Mod.NegBin71, Mod.NegBin72, Mod.NegBin73, Mod.NegBin74, Mod.NegBin75, Mod.NegBin77, Mod.NegBin79, Mod.NegBin81, Mod.NegBin82, Mod.NegBin84, Mod.NegBin85, Mod.NegBin86, Mod.NegBin88, Mod.NegBin89, Mod.NegBin90, Mod.NegBin91, Mod.NegBin92, Mod.NegBin93, Mod.NegBin94, Mod.NegBin95, Mod.NegBin96, Mod.NegBin97, Mod.NegBin98, Mod.NegBin99, Mod.NegBin100, Mod.NegBin101, Mod.NegBin103, Mod.NegBin104, Mod.NegBin106, Mod.NegBin108, Mod.NegBin109, Mod.NegBin110, Mod.NegBin111))

######################################

######## OMR -3500 ####################
Mod.results  <- as.data.frame(AICc(Mod.NegBin0, Mod.NegBin1, Mod.NegBin2, Mod.NegBin3, Mod.NegBin4, Mod.NegBin5, Mod.NegBin6, Mod.NegBin7, Mod.NegBin8, Mod.NegBin9, Mod.NegBin10, Mod.NegBin11, Mod.NegBin12, Mod.NegBin13, Mod.NegBin14, Mod.NegBin15, Mod.NegBin16, Mod.NegBin17, Mod.NegBin18, Mod.NegBin19, Mod.NegBin20, Mod.NegBin21, Mod.NegBin22, Mod.NegBin23, Mod.NegBin24, Mod.NegBin25, Mod.NegBin26, Mod.NegBin27, Mod.NegBin28, Mod.NegBin29, Mod.NegBin30, Mod.NegBin31, Mod.NegBin32, Mod.NegBin33, Mod.NegBin34, Mod.NegBin35, Mod.NegBin36, Mod.NegBin37, Mod.NegBin38, Mod.NegBin39, Mod.NegBin40, Mod.NegBin41, Mod.NegBin42, Mod.NegBin43, Mod.NegBin44, Mod.NegBin45, Mod.NegBin46, Mod.NegBin47, Mod.NegBin48, Mod.NegBin49, Mod.NegBin52, Mod.NegBin53, Mod.NegBin54, Mod.NegBin55, Mod.NegBin56, Mod.NegBin57, Mod.NegBin58, Mod.NegBin59, Mod.NegBin60, Mod.NegBin61, Mod.NegBin63, Mod.NegBin64, Mod.NegBin65,Mod.NegBin66, Mod.NegBin67, Mod.NegBin68, Mod.NegBin69, Mod.NegBin70, Mod.NegBin71, Mod.NegBin72, Mod.NegBin73, Mod.NegBin74, Mod.NegBin75, Mod.NegBin77, Mod.NegBin78, Mod.NegBin79, Mod.NegBin81, Mod.NegBin82, Mod.NegBin84, Mod.NegBin85, Mod.NegBin88, Mod.NegBin89, Mod.NegBin90, Mod.NegBin91, Mod.NegBin93, Mod.NegBin94, Mod.NegBin95, Mod.NegBin96, Mod.NegBin97, Mod.NegBin98, Mod.NegBin99, Mod.NegBin100, Mod.NegBin101, Mod.NegBin102, Mod.NegBin103, Mod.NegBin104, Mod.NegBin106, Mod.NegBin108, Mod.NegBin109, Mod.NegBin110, Mod.NegBin111))

21.824458 - 5.157158

######################################

head(Mod.results)

Mod.results$Model.ID <- row.names(Mod.results)
Mod.results <- transform(Mod.results, Delta.AIC = AICc - min(Mod.results$AICc))
Mod.results <- Mod.results[order(Mod.results$Delta.AIC),]

Mod.results[1:11,]


# Top Model - Note 78 is < 1 AICc behind #

Mod.NegBin42 <- glm.nb(Salvage.Total~Month*OMR.z + Total.exports.z, data=Winter.salv)
summary(Mod.NegBin42)

####### Pseudo R^2 ###########
library("rsq")

rsq(Mod.NegBin42, adj=TRUE)
rsq.n(Mod.NegBin42)
?rsq.n

############# fit with Poisson distribution ###############

#Mod.pois21 <- glm(Salvage.Total~OMR.z + Total.exports.z*Month, family="poisson", data=Winter.salv)
#summary(Mod.pois21)
#hist(resid(Mod.pois21), breaks=10)

############# NegBin model coefficients #############

est <- cbind(Estimate = coef(Mod.NegBin42), confint(Mod.NegBin42))
exp(est)

#Interpretation: The "baseline" average salvage count is 91. The other exponentiated coefficients are interpreted multiplicatively. One unit increase (measured in Standard Diviations) in OMR.z increases the average salvage count by 0.36 times. Total exports (i.e., plus one SD) increases the average salvage count by 1.78 times. However, there are strong effects of month x OMR.z.

############# Model residulas ###############

Winter.salv$resid <- resid(Mod.NegBin42)

#boxplot(Winter.salv$resid)
hist(Winter.salv$resid, breaks = 20)

############## Looping prediction of Mod.NegBin21 ####################

Pred.dat <- data.frame(Observed.dat = rep(NA, nrow(Winter.salv)), Pred.NegBin = rep(NA, nrow(Winter.salv)))
head(Pred.dat)

for(i in 1:nrow(Winter.salv)){
  Hat.dat <- Winter.salv[-i, ]
  CV.dat <- Winter.salv[i, ]
  Temp.mod1 <- glm.nb(Salvage.Total~Month*OMR.z + Total.exports.z, data=Hat.dat)
  Pred.dat[i, "Observed.dat"] <- CV.dat[, "Salvage.Total"]
  Pred.dat[i, "Pred.NegBin"] <- predict(Temp.mod1, newdata = CV.dat, type="response")
}

head(Pred.dat)

plot(Pred.dat$Observed.dat, Pred.dat$Pred.NegBin)

loocv.lm <- lm(log10(Pred.NegBin+1)~log10(Observed.dat+1), data=Pred.dat)
summary(loocv.lm)

#jpeg(filename = "LOOCV_NegBin_Prediction.jpg", heigh=5, width = 6, units="in", bg="white", res=500)

par(mfrow=c(1,1), oma=c(3, 3.5, 0.5, 0.5), mar=c(0.5, 0.75, 0.5, 0.25), bty="l")

plot(log10(Pred.dat$Observed.dat), log10(Pred.dat$Pred.NegBin), pch=21, bg="darkgrey", cex=1.2, xlim=c(0, 4.5), ylim=c(0, 4.5), ylab="", xlab="", axes=FALSE, frame=TRUE)
abline(loocv.lm, col="red")
lines(c(0, 5), c(0, 5), lty=2)
axis(1, at=c(0, 1, 2, 3, 4))
axis(2, at=c(0, 1, 2, 3, 4), las=1)
mtext(1, text="Observed Monthy Salvage (log(n + 1))", line=2.5)
mtext(2, text="Predicted Monthy Salvage (log(n + 1))", line=3)
text(0, 0, labels="1:1 line", cex=0.65, pos=4)

#dev.off()

###################################################################
################# Prediction data frame ###########################
###################################################################
range(Winter.salv$Total.exports.z)
range(Winter.salv$OMR.z)

#length(seq(-1.2, 5.2, 0.1))
#length(rep(seq(-1.2, 5.2, 0.1), 4))

pred.dat <- data.frame(OMR.z  = rep(rep(seq(-1.2, 5.2, 0.1), 5), 2), Month = rep(c(rep("December", 65), rep("January", 65), rep("February", 65), rep("March", 65), rep("April", 65)), 2), Total.exports.z = c(rep(-1.5, 325), rep(1.5, 325)))

#pred.goods <- predict(Mod.NegBin42, newdata = pred.dat, type="link", se.fit=TRUE)
### OMR -5000
pred.goods <- predict(Mod.NegBin44, newdata = pred.dat, type="link", se.fit=TRUE)

pred.dat$fit <- pred.goods$fit
pred.dat$se.fit <- pred.goods$se.fit

pred.dat$count.fit <- exp(pred.dat$fit)
pred.dat$CIU <- exp(pred.dat$fit + 1.96*pred.dat$se.fit)
pred.dat$CIL <- exp(pred.dat$fit - 1.96*pred.dat$se.fit)

head(pred.dat)

library("RColorBrewer")
display.brewer.pal(11, "Spectral")
brewer.pal(11, "Spectral")

cols <- data.frame(Month = c("December", "January", "February", "March", "April"), col = c("#3288BD", "#5E4FA2", "#66C2A5", "#ABDDA4","#FEE08B"), poly.col=c("#3288BD35", "#5E4FA235", "#66C2A535", "#ABDDA435", "#FEE08B35"), pt.col=c("#3288BD65", "#5E4FA265", "#66C2A565", "#ABDDA465", "#FEE08B65"))
pred.dat <- merge(pred.dat, cols, by="Month")
pred.dat <- transform(pred.dat, col = as.character(col), poly.col = as.character(poly.col))

head(pred.dat)


#mean(Winter.salv$OMR)
#sd(Winter.salv$OMR)

#plot(Winter.salv$OMR.z, (Winter.salv$OMR-mean(Salv.dat2$OMR, na.rm=TRUE))/sd(Salv.dat2$OMR, na.rm=TRUE), xlim=c(-1, 5), ylim=c(-1, 5))
#lines(c(-1, 5), c(-1, 5), col="red")

##### Correcting OMR/exports back to natural scale ######

pred.dat <- transform(pred.dat, OMR = OMR.z*sd(Salv.dat2$OMR, na.rm=TRUE)+mean(Salv.dat2$OMR, na.rm=TRUE), Total.Exports = Total.exports.z*sd(Salv.dat2$Total.Exports, na.rm=TRUE)+mean(Salv.dat2$Total.Exports, na.rm=TRUE))

max(pred.dat$count.fit)

head(pred.dat)
unique(pred.dat$Total.Exports)

pred.dec.LE <- subset(pred.dat, Month == "December" & Total.Exports < 2000)
pred.dec.HE <- subset(pred.dat, Month == "December" & Total.Exports > 5000)

pred.jan.LE <- subset(pred.dat, Month == "January" & Total.Exports < 2000)
pred.jan.HE <- subset(pred.dat, Month == "January" & Total.Exports > 5000)

pred.feb.LE <- subset(pred.dat, Month == "February" & Total.Exports < 2000)
pred.feb.HE <- subset(pred.dat, Month == "February" & Total.Exports > 5000)

pred.mar.LE <- subset(pred.dat, Month == "March" & Total.Exports < 2000)
pred.mar.HE <- subset(pred.dat, Month == "March" & Total.Exports > 5000)

pred.apr.LE <- subset(pred.dat, Month == "April" & Total.Exports < 2000)
pred.apr.HE <- subset(pred.dat, Month == "April" & Total.Exports > 5000)

##################################################################################
##################################################################################
head(Winter.salv)

Winter.salv.plt <- merge(Winter.salv, cols, by="Month", all.x=TRUE)
Winter.salv.plt <- transform(Winter.salv.plt, pt.col=as.character(pt.col))
Winter.salv.plt$cex <- 0.5+(max(Winter.salv.plt$Total.Exports)-Winter.salv.plt$Total.Exports)/diff(range(Winter.salv.plt$Total.Exports))

head(Winter.salv.plt)

str(pred.dec.LE)

##################################################################################

#jpeg(filename = "Salvage_Model_LowExports.jpg", heigh=5, width = 6, units="in", bg="white", res=500)

par(mfrow=c(1,1), oma=c(3, 3.5, 0.5, 0), mar=c(0.5, 0.75, 0.5, 0.25), bty="l")

plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=TRUE)
polygon(c(pred.dec.LE$OMR, rev(pred.dec.LE$OMR)), c(pred.dec.LE$CIU, rev(pred.dec.LE$CIL)), border=FALSE, col=pred.dec.LE$poly.col)
par(new=TRUE)
plot(pred.dec.LE$OMR, pred.dec.LE$count.fit, type="l", lwd=2, col=pred.dec.LE$col, xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.jan.LE$OMR[order(pred.jan.LE$OMR)], rev(pred.jan.LE$OMR[order(pred.jan.LE$OMR)])), c(pred.jan.LE$CIU[order(pred.jan.LE$OMR)], rev(pred.jan.LE$CIL[order(pred.jan.LE$OMR)])), border=FALSE, col=pred.jan.LE$poly.col)
par(new=TRUE)
plot(pred.jan.LE$OMR[order(pred.jan.LE$OMR)], pred.jan.LE$count.fit[order(pred.jan.LE$OMR)], type="l", lwd=2, col=pred.jan.LE$col, xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.feb.LE$OMR[order(pred.feb.LE$OMR)], rev(pred.feb.LE$OMR[order(pred.feb.LE$OMR)])), c(pred.feb.LE$CIU[order(pred.feb.LE$OMR)], rev(pred.feb.LE$CIL[order(pred.feb.LE$OMR)])), border=FALSE, col=pred.feb.LE$poly.col)
par(new=TRUE)
plot(pred.feb.LE$OMR[order(pred.feb.LE$OMR)], pred.feb.LE$count.fit[order(pred.feb.LE$OMR)], type="l", lwd=2, col=pred.feb.LE$col, xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.mar.LE$OMR[order(pred.mar.LE$OMR)], rev(pred.mar.LE$OMR[order(pred.mar.LE$OMR)])), c(pred.mar.LE$CIU[order(pred.mar.LE$OMR)], rev(pred.mar.LE$CIL[order(pred.mar.LE$OMR)])), border=FALSE, col=pred.mar.LE$poly.col)
par(new=TRUE)
plot(pred.mar.LE$OMR[order(pred.mar.LE$OMR)], pred.mar.LE$count.fit[order(pred.mar.LE$OMR)], type="l", lwd=2, col=pred.mar.LE$col, xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.apr.LE$OMR[order(pred.apr.LE$OMR)], rev(pred.apr.LE$OMR[order(pred.apr.LE$OMR)])), c(pred.apr.LE$CIU[order(pred.apr.LE$OMR)], rev(pred.apr.LE$CIL[order(pred.apr.LE$OMR)])), border=FALSE, col=pred.apr.LE$poly.col)
par(new=TRUE)
plot(pred.apr.LE$OMR[order(pred.apr.LE$OMR)], pred.apr.LE$count.fit[order(pred.apr.LE$OMR)], type="l", lwd=2, col=pred.apr.LE$col, xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)
axis(1, at=c(-9000, -4500, 0, 4500, 9000))
axis(2, at=c(0, 500, 1000, 1500, 2000, 2500), las=1)
par(new=TRUE)
plot(Winter.salv.plt$OMR, Winter.salv.plt$Salvage.Total, pch=21, col="#00000085", bg=Winter.salv.plt$pt.col, cex=Winter.salv.plt$cex, xlim=c(-9000, 10000), ylim=c(0, 2500), xlab="", ylab="", axes=FALSE, frame=FALSE)
mtext(1, text="OMR (cfs)", line=2.25)
mtext(2, text = "Expected Monthly Winter-Run Salvage (n)", line=3.25)
text(c(6000, 6000, 6000, 6000, 6000), c(2500, 2375, 2250, 2125, 2000), labels=c("Dec", "Jan", "Feb", "Mar", "Apr"), pos=4, cex=0.8)
lines(c(7500, 9000), c(2500, 2500), col="#3288BD", lwd=2.5)
lines(c(7500, 9000), c(2375, 2375), col="#5E4FA2", lwd=2.5)
lines(c(7500, 9000), c(2250, 2250), col="#66C2A5", lwd=2.5)
lines(c(7500, 9000), c(2125, 2125), col="#ABDDA4", lwd=2.5)
lines(c(7500, 9000), c(2000, 2000), col="#FEE08B", lwd=2.5)

#dev.off()

############### High Exports ###############

jpeg(filename = "Salvage_Model_Exports5000_OMR3500.jpg", heigh=5, width = 6, units="in", bg="white", res=500)

par(mfrow=c(1,1), oma=c(3, 3.5, 0.5, 0), mar=c(0.5, 0.75, 0.5, 0.25), bty="l")

plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=TRUE)
polygon(c(pred.dec.HE$OMR[order(pred.dec.HE$OMR)], rev(pred.dec.HE$OMR[order(pred.dec.HE$OMR)])), c(pred.dec.HE$CIU[order(pred.dec.HE$OMR)], rev(pred.dec.HE$CIL[order(pred.dec.HE$OMR)])), border=FALSE, col=pred.dec.HE$poly.col)
par(new=TRUE)
plot(pred.dec.HE$OMR[order(pred.dec.HE$OMR)], pred.dec.HE$count.fit[order(pred.dec.HE$OMR)], type="l", lwd=2, col=pred.dec.HE$col, xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.jan.HE$OMR[order(pred.jan.HE$OMR)], rev(pred.jan.HE$OMR[order(pred.jan.HE$OMR)])), c(pred.jan.HE$CIU[order(pred.jan.HE$OMR)], rev(pred.jan.HE$CIL[order(pred.jan.HE$OMR)])), border=FALSE, col=pred.jan.HE$poly.col)
par(new=TRUE)
plot(pred.jan.HE$OMR[order(pred.jan.HE$OMR)], pred.jan.HE$count.fit[order(pred.jan.HE$OMR)], type="l", lwd=2, col=pred.jan.HE$col, xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.feb.HE$OMR[order(pred.feb.HE$OMR)], rev(pred.feb.HE$OMR[order(pred.feb.HE$OMR)])), c(pred.feb.HE$CIU[order(pred.feb.HE$OMR)], rev(pred.feb.HE$CIL[order(pred.feb.HE$OMR)])), border=FALSE, col=pred.feb.HE$poly.col)
par(new=TRUE)
plot(pred.feb.HE$OMR[order(pred.feb.HE$OMR)], pred.feb.HE$count.fit[order(pred.feb.HE$OMR)], type="l", lwd=2, col=pred.feb.HE$col, xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.mar.HE$OMR[order(pred.mar.HE$OMR)], rev(pred.mar.HE$OMR[order(pred.mar.HE$OMR)])), c(pred.mar.HE$CIU[order(pred.mar.HE$OMR)], rev(pred.mar.HE$CIL[order(pred.mar.HE$OMR)])), border=FALSE, col=pred.mar.HE$poly.col)
par(new=TRUE)
plot(pred.mar.HE$OMR[order(pred.mar.HE$OMR)], pred.mar.HE$count.fit[order(pred.mar.HE$OMR)], type="l", lwd=2, col=pred.mar.HE$col, xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)

par(new=TRUE)
plot(1, type="n", xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)
polygon(c(pred.apr.HE$OMR[order(pred.apr.HE$OMR)], rev(pred.apr.HE$OMR[order(pred.apr.HE$OMR)])), c(pred.apr.HE$CIU[order(pred.apr.HE$OMR)], rev(pred.apr.HE$CIL[order(pred.apr.HE$OMR)])), border=FALSE, col=pred.apr.HE$poly.col)
par(new=TRUE)
plot(pred.apr.HE$OMR[order(pred.apr.HE$OMR)], pred.apr.HE$count.fit[order(pred.apr.HE$OMR)], type="l", lwd=2, col=pred.apr.HE$col, xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)
axis(1, at=c(-9000, -4500, 0, 4500, 9000))
axis(2, at=c(0, 1000, 2000, 3000, 4000, 5000, 6000), las=1)

par(new=TRUE)
plot(Winter.salv.plt$OMR, Winter.salv.plt$Salvage.Total, pch=21, col="#00000085", bg=Winter.salv.plt$pt.col, cex=Winter.salv.plt$cex, xlim=c(-9000, 10000), ylim=c(0, 6000), xlab="", ylab="", axes=FALSE, frame=FALSE)
mtext(1, text="OMR (cfs)", line=2.25)
mtext(2, text = "Expected Monthly Winter-Run Salvage (n)", line=3.25)
text(c(6000, 6000, 6000, 6000), c(6000, 5700, 5400, 5100, 4800), labels=c("Dec", "Jan", "Feb", "Mar", "Apr"), pos=4, cex=0.8)
lines(c(7500, 9000), c(6000, 6000), col="#3288BD", lwd=2.5)
lines(c(7500, 9000), c(5700, 5700), col="#5E4FA2", lwd=2.5)
lines(c(7500, 9000), c(5400, 5400), col="#66C2A5", lwd=2.5)
lines(c(7500, 9000), c(5100, 5100), col="#ABDDA4", lwd=2.5)
lines(c(7500, 9000), c(4800, 4800), col="#FEE08B", lwd=2.5)

dev.off()

#cols <- data.frame(Month = c("December", "January", "February", "March"), col = c("#3288BD", "#5E4FA2", "#66C2A5", "#ABDDA4"), poly.col=c("#3288BD35", "#5E4FA235", "#66C2A535", "#ABDDA435"), pt.col=c("#3288BD65", "#5E4FA265", "#66C2A565", "#ABDDA465"))

##################################################################


############### Fitting Zero Inflated Model to All Salmon Data ####################

#Mod.ZeroInf0 <- zeroinfl(Salvage.Total~Month | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf1 <- zeroinfl(Salvage.Total~Total.exports.z | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf2 <- zeroinfl(Salvage.Total~OMR.z | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf3 <- zeroinfl(Salvage.Total~Total.exports.z+OMR.z | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf4 <- zeroinfl(Salvage.Total~Total.exports.z*OMR.z | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf5 <- zeroinfl(Salvage.Total~Discharge.z | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf6 <- zeroinfl(Salvage.Total~Discharge.z+Total.exports.z | 1, dist="negbin", data=Salv.dat4)
#Mod.ZeroInf7 <- zeroinfl(Salvage.Total~Discharge.z+OMR.z | 1, dist="negbin", data=Salv.dat4)


#AIC(Mod.ZeroInf0, Mod.ZeroInf1, Mod.ZeroInf2, Mod.ZeroInf3, Mod.ZeroInf4, Mod.ZeroInf5, Mod.ZeroInf6, Mod.ZeroInf7)

############### Fitting NegBin Model to All Salmon Data ####################

#Mod.NegBin0 <- glm.nb(Salvage.Total~NULL, data=Salv.dat4)
#Mod.NegBin00 <- glm.nb(Salvage.Total~Month, data=Salv.dat4)
#Mod.NegBin1 <- glm.nb(Salvage.Total~Total.exports.z, data=Salv.dat4)
#Mod.NegBin2 <- glm.nb(Salvage.Total~OMR.z, data=Salv.dat4)
#Mod.NegBin3 <- glm.nb(Salvage.Total~Total.exports.z+OMR.z, data=Salv.dat4)
#Mod.NegBin4 <- glm.nb(Salvage.Total~Total.exports.z*OMR.z, data=Salv.dat4)
#Mod.NegBin5 <- glm.nb(Salvage.Total~Discharge.z, data=Salv.dat4)
#Mod.NegBin6 <- glm.nb(Salvage.Total~Discharge.z+Total.exports.z, data=Salv.dat4)
#Mod.NegBin7 <- glm.nb(Salvage.Total~Discharge.z+OMR.z, data=Salv.dat4)

#AICc(Mod.NegBin0, Mod.NegBin00, Mod.NegBin1, Mod.NegBin2, Mod.NegBin3, Mod.NegBin4, Mod.NegBin5, Mod.NegBin6, Mod.NegBin7)

#summary(Mod.NegBin1)
#summary(Mod.NegBin2)
#summary(Mod.NegBin5)

############## Looping prediction of Mod.Zero ####################

#Pred.dat <- data.frame(Observed.dat = rep(NA, nrow(Salv.dat4)), Pred.OMR.dat = rep(NA, nrow(Salv.dat4)), Pred.Exports.dat = rep(NA, nrow(Salv.dat4)), Pred.Flow.dat = rep(NA, nrow(Salv.dat4)))

#head(Pred.dat)

#for(i in 1:nrow(Salv.dat4)){
#  Hat.dat <- Salv.dat4[-i, ]
#  CV.dat <- Salv.dat4[i, ]
#  Temp.mod1 <- glm.nb(Salvage.Total~Total.exports.z, data=Hat.dat)
#  Temp.mod2 <- glm.nb(Salvage.Total~OMR.z, data=Hat.dat)
#  Temp.mod3 <- glm.nb(Salvage.Total~Discharge.z, data=Hat.dat)
#  Pred.dat[i, "Observed.dat"] <- CV.dat[, "Salvage.Total"]
#  Pred.dat[i, "Pred.OMR.dat"] <- predict(Temp.mod2, newdata = CV.dat, type="response")
#  Pred.dat[i, "Pred.Exports.dat"] <- predict(Temp.mod1, newdata = CV.dat, type="response")
#  Pred.dat[i, "Pred.Flow.dat"] <- predict(Temp.mod3, newdata = CV.dat, type="response")
  #browser()
#}

#head(Pred.dat)

#plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.OMR.dat))
#plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.Exports.dat))
#plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.Flow.dat))


#plot(Pred.dat$Observed.dat, Pred.dat$Pred.OMR.dat)
#plot(Pred.dat$Observed.dat, Pred.dat$Pred.Exports.dat)
#plot(Pred.dat$Observed.dat, Pred.dat$Pred.Flow.dat)


##################################################################

#################### Spring Run ############################

#Spring.salv <- subset(Salv.dat3, Run == "Spring")
#head(Spring.salv)

#unique(Spring.salv$Month)

#nrow(Spring.salv2)
#hist(Spring.salv2$Salvage.Total, breaks = 100)

############### Fitting Zero Inflated Model to Spring Salmon Data ####################

#Mod.ZeroInf0 <- zeroinfl(Salvage.Total~NULL | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf00 <- zeroinfl(Salvage.Total~Month | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf1 <- zeroinfl(Salvage.Total~Total.exports.z | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf2 <- zeroinfl(Salvage.Total~OMR.z | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf3 <- zeroinfl(Salvage.Total~Total.exports.z+OMR.z | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf4 <- zeroinfl(Salvage.Total~Total.exports.z*OMR.z | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf5 <- zeroinfl(Salvage.Total~Discharge.z | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf6 <- zeroinfl(Salvage.Total~Discharge.z+Total.exports.z | 1, dist="negbin", data=Spring.salv)
#Mod.ZeroInf7 <- zeroinfl(Salvage.Total~Discharge.z+OMR.z | 1, dist="negbin", data=Spring.salv)


#AICc(Mod.ZeroInf0, Mod.ZeroInf00, Mod.ZeroInf1, Mod.ZeroInf2, Mod.ZeroInf3, Mod.ZeroInf4, Mod.ZeroInf5, Mod.ZeroInf6, Mod.ZeroInf7)

############### Fitting NegBin Model to Spring Salmon Data ####################

#Mod.NegBin0 <- glm.nb(Salvage.Total~NULL, data=Spring.salv)
#Mod.NegBin00 <- glm.nb(Salvage.Total~Month, data=Spring.salv)
#Mod.NegBin1 <- glm.nb(Salvage.Total~Total.exports.z, data=Spring.salv)
#Mod.NegBin2 <- glm.nb(Salvage.Total~OMR.z, data=Spring.salv)
#Mod.NegBin3 <- glm.nb(Salvage.Total~Total.exports.z+OMR.z, data=Spring.salv)
#Mod.NegBin4 <- glm.nb(Salvage.Total~Total.exports.z*OMR.z, data=Spring.salv)
#Mod.NegBin5 <- glm.nb(Salvage.Total~Discharge.z, data=Spring.salv)
#Mod.NegBin6 <- glm.nb(Salvage.Total~Discharge.z+Total.exports.z, data=Spring.salv)
#Mod.NegBin7 <- glm.nb(Salvage.Total~Discharge.z+OMR.z, data=Spring.salv)

#AICc(Mod.NegBin0, Mod.NegBin00, Mod.NegBin1, Mod.NegBin2, Mod.NegBin3, Mod.NegBin4, Mod.NegBin5, Mod.NegBin6, Mod.NegBin7)

#summary(Mod.NegBin1)
#summary(Mod.NegBin2)
#summary(Mod.NegBin5)

############## Looping prediction of Mod.Zero ####################

#Pred.dat <- data.frame(Observed.dat = rep(NA, nrow(Spring.salv2)), Pred.OMR.dat = rep(NA, nrow(Spring.salv2)), Pred.Exports.dat = rep(NA, nrow(Spring.salv2)))

#head(Pred.dat)

#for(i in 1:nrow(Spring.salv2)){
#  Hat.dat <- Spring.salv2[-i, ]
#  CV.dat <- Spring.salv2[i, ]
#  Temp.mod1 <- zeroinfl(Salvage.Total~Total.exports.z | 1, dist="negbin", data=Hat.dat)
#  Temp.mod2 <- zeroinfl(Salvage.Total~OMR.z | 1, dist="negbin", data=Hat.dat)
#  Pred.dat[i, "Observed.dat"] <- CV.dat[, "Salvage.Total"]
#  Pred.dat[i, "Pred.OMR.dat"] <- predict(Temp.mod2, newdata = CV.dat, type="response")
#  Pred.dat[i, "Pred.Exports.dat"] <- predict(Temp.mod1, newdata = CV.dat, type="response")
  #browser()
#}

#head(Pred.dat)

#plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.OMR.dat))
#plot(log(Pred.dat$Observed.dat), log(Pred.dat$Pred.Exports.dat))


#################### Fall Run ############################

Fall.salv <- subset(Salv.dat3, Run == "Fall")

head(Fall.salv)
nrow(Fall.salv)
hist(Fall.salv$Salvage.Total, breaks = 100)


############### Fitting Zero Inflated Model to Winter Salmon Data ####################

Mod.ZeroInf0 <- zeroinfl(Salvage.Total~NULL | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf1 <- zeroinfl(Salvage.Total~Total.exports.z | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf2 <- zeroinfl(Salvage.Total~OMR.z | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf3 <- zeroinfl(Salvage.Total~Total.exports.z+OMR.z | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf4 <- zeroinfl(Salvage.Total~Total.exports.z*OMR.z | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf5 <- zeroinfl(Salvage.Total~Discharge.z | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf6 <- zeroinfl(Salvage.Total~Discharge.z+Total.exports.z | 1, dist="negbin", data=Fall.salv)
Mod.ZeroInf7 <- zeroinfl(Salvage.Total~Discharge.z+OMR.z | 1, dist="negbin", data=Fall.salv)


AICc(Mod.ZeroInf0, Mod.ZeroInf1, Mod.ZeroInf2, Mod.ZeroInf3, Mod.ZeroInf4, Mod.ZeroInf5, Mod.ZeroInf6, Mod.ZeroInf7)

############### Fitting NegBin Model to Winter Salmon Data ####################

Mod.NegBin0 <- glm.nb(Salvage.Total~NULL, data=Fall.salv)
Mod.NegBin1 <- glm.nb(Salvage.Total~Total.exports.z, data=Fall.salv)
Mod.NegBin2 <- glm.nb(Salvage.Total~OMR.z, data=Fall.salv)
Mod.NegBin3 <- glm.nb(Salvage.Total~Total.exports.z+OMR.z, data=Fall.salv)
Mod.NegBin4 <- glm.nb(Salvage.Total~Total.exports.z*OMR.z, data=Fall.salv)
Mod.NegBin5 <- glm.nb(Salvage.Total~Discharge.z, data=Fall.salv)
Mod.NegBin6 <- glm.nb(Salvage.Total~Discharge.z+Total.exports.z, data=Fall.salv)
Mod.NegBin7 <- glm.nb(Salvage.Total~Discharge.z+OMR.z, data=Fall.salv)

AICc(Mod.NegBin0, Mod.NegBin1, Mod.NegBin2, Mod.NegBin3, Mod.NegBin4, Mod.NegBin5, Mod.NegBin6, Mod.NegBin7)

summary(Mod.NegBin5)







