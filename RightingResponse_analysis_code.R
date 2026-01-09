### Righting trials analysis for both 2019 and 2020 hatchling trials

rm(list = ls())

### TODO: set the working directory for the location of the data files
setwd('<<insert working directory path here>>')


library(lmerTest)
library(emmeans)
library(DHARMa)

###################################################
#########LOAD data files##########
# Saved as TAB delimited text files

###1 row per individual
hatchlingdata1row <- read.table(file = "hatchling_data.txt",
 header = TRUE, sep = "\t")

hatchlingdata1row <- within(hatchlingdata1row,{
  Year <- as.factor(Year)
  Clutch <- as.factor(Clutch)
  HatchlingID <- as.factor(HatchlingID)
  HatchlingNumber <- as.factor(HatchlingNumber)
  WinterHous <- as.factor(WinterHous) 
})

###3 rows per individual, with righting times
RRdata3rows <- read.table(file = "hatchlingAndRightingResponse_data.txt",
  header = TRUE, sep = "\t")

RRdata3rows <- within(RRdata3rows,{
  Year <- as.factor(Year)
  Clutch <- as.factor(Clutch)
  HatchlingID <- as.factor(HatchlingID)
  HatchlingNumber <- as.factor(HatchlingNumber)
  WinterHous <- as.factor(WinterHous) 
  Round <- as.factor(Round)
})

###incubation lengths
inclengths <- read.table(file = "incubationLengths_data.txt",
  header = TRUE, sep = "\t")

inclengths <- within(inclengths,{
  Year <- as.factor(Year)
  Clutch <- as.factor(Clutch)
  HatchlingID <- as.factor(HatchlingID)
  HatchlingNumber <- as.factor(HatchlingNumber)
})

###egg survival
eggsurvival <- read.table(file = "eggSurvival_data.txt",
  header = TRUE, sep = "\t")

eggsurvival <- within(eggsurvival,{
  Year <- as.factor(Year)
  Clutch <- as.factor(Clutch)
  EggID <- as.factor(EggID)
})

####################################################
###PCA on post-hatching measurements
###dataset with one row per individual


pca_posthatch <- prcomp( ~ SCL1 + CW1 + PL1 + PW1 + BD1 + NV1 + VT1 + Mass1,
                         data = hatchlingdata1row,
                         center = TRUE, scale. = TRUE)

summary(pca_posthatch)
plot(pca_posthatch)

barplot(summary(pca_posthatch)$importance["Proportion of Variance", ], main = "",
        ylab = "Variance (proportion)", ylim = c(0, 1))

pca_posthatch



###PCA on pre-trial measurements

pca_pretrial <- prcomp( ~ SCL2 + CW2 + PL2 + PW2 + BD2 + NV2 + VT2 + Mass2,
                        data = hatchlingdata1row,
                        center = TRUE, scale. = TRUE)

summary(pca_pretrial)
plot(pca_pretrial)

barplot(summary(pca_pretrial)$importance["Proportion of Variance", ], main = "",
        ylab = "Variance (proportion)", ylim = c(0, 1))

pca_pretrial



###add pc scores to dataframe

##make 2 new columns in dataframe that PC scores from post hatch PCA can be added to, then add scores
hatchlingdata1row[ ,c("PC1ph","PC2ph")]<-pca_posthatch$x[ ,c("PC1","PC2")]

##make 2 new columns in dataframe that PC scores from pre trial PCA can be added to, then add scores
hatchlingdata1row[ ,c("PC1pt","PC2pt")]<-pca_pretrial$x[ ,c("PC1","PC2")]

##make 2 new columns in dataframe that PC scores from pre trial PCA can be added
RRdata3rows[ ,c("PC1pt","PC2pt")] <-
hatchlingdata1row[match(RRdata3rows$HatchlingID,
                        hatchlingdata1row$HatchlingID), c("PC1pt", "PC2pt")]


###PCA on "growth": change in each phenotype from fall to spring (spring-fall)
###dataset with one row per individual
###10/28/2024
### growth calculated as:
#### (spring measurement) - (fall measurement) / (days between measurements)  
hatchlingdata1row[, c("SCLgrD", "CWgrD", "PLgrD", "PWgrD", "MassgrD")] <-
  (hatchlingdata1row[, c("SCL2", "CW2", "PL2", "PW2", "Mass2")] -
  hatchlingdata1row[, c("SCL1", "CW1", "PL1", "PW1", "Mass1")]) / hatchlingdata1row$DaysBwMeas
  

###PCA including SCL, CW, PL, PW, Mass, 
###standardized by # days between measurements (not including body depth or tail)
pca_growthD <- prcomp( ~ SCLgrD + CWgrD + PLgrD + PWgrD + MassgrD,
                       data = hatchlingdata1row,
                       center = TRUE, scale. = TRUE)

summary(pca_growthD)
plot(pca_growthD)

barplot(summary(pca_growthD)$importance["Proportion of Variance", ], main = "",
        ylab = "Variance (proportion)", ylim = c(0, 1))

pca_growthD

##make new column in dataframe that growth PC scores can be added to, then add scores
hatchlingdata1row[ ,c("PC1growthD")]<-pca_growthD$x[ ,c("PC1")]




##Use incubation temp as a continuous variable (IncTemp, centered in R)



###table 1

###EGG SURVIVAL

###center incubation temperature variable
eggsurvival$centeredIncTemp=scale(eggsurvival$IncTemp,center=TRUE,scale=FALSE)

eggsurv=glmer(EggSurvival~centeredIncTemp+Egg.Mass+(1|Clutch),data=eggsurvival,family=binomial)
summary(eggsurv)

###INCUBATION LENGTH

###center Incubation Temperature variable 
inclengths$centeredIncTemp=scale(inclengths$IncTemp,center=TRUE,scale=FALSE)

incdur=glmer(IncLength~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=inclengths,family=poisson(link ="log"))
summary(incdur)



###effect of inc temp on egg incubation and post hatch measurements (TABLE 2)

###center Incubation Temperature variable 
hatchlingdata1row$centeredIncTemp=scale(hatchlingdata1row$IncTemp,center=TRUE,scale=FALSE)

scl1=lmer(SCL1~centeredIncTemp+EggMass+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl1)

cw1=lmer(CW1~centeredIncTemp+EggMass+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(cw1)

pl1=lmer(PL1~centeredIncTemp+EggMass+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl1)

pw1=lmer(PW1~centeredIncTemp+EggMass+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pw1)

bd1=lmer(BD1~centeredIncTemp+EggMass+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(bd1)

nv1=lmer(NV1~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(nv1)

vt1=lmer(VT1~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(vt1)

mass1=lmer(Mass1~centeredIncTemp+EggMass+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(mass1)



###effect of inc temp, ovw location, and time in winter housing on pre trial measurements (TABLE 3)
###includes days in winter housing covariate (d)
###all include previous measurement as a fixed effect, except for growth PC

scl2d=lmer(SCL2~centeredIncTemp+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

cw2d=lmer(CW2~centeredIncTemp+WinterHous+CW1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(cw2d)

pl2d=lmer(PL2~centeredIncTemp+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

pw2d=lmer(PW2~centeredIncTemp+WinterHous+PW1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pw2d)

bd2d=lmer(BD2~centeredIncTemp+WinterHous+BD1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(bd2d)

mass2d=lmer(Mass2~centeredIncTemp+WinterHous+Mass1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(mass2d)

growthDd=lmer(PC1growthD~centeredIncTemp+WinterHous+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(growthDd)



#####################################################################
###analyses using 3 row per individual righting response dataset

RRdata3rows$centeredIncTemp3rows=RRdata3rows$IncTemp - 28#,center=TRUE,scale=FALSE)


###stepwise models(including year)  


######response = total RR time######
rrmodel1=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+PC1pt+PC2pt+as.factor(Year)+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel1)

###remove Days in winter housing
rrmodel2=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+PC1pt+PC2pt+as.factor(Year)+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel2)

###remove PC1
rrmodel3=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+PC2pt+as.factor(Year)+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel3)

###remove Year
rrmodel4=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+PC2pt+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel4)

###remove Age at Trials
rrmodel5=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PC2pt+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel5)

###remove DMIL
rrmodel6=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PC2pt+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel6)

###remove PC2
rrmodel7=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodel7)


###########################


######response = latency #####
latency1=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+PC1pt+PC2pt+ScaledDevMaxIncLength+as.factor(Year)+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latency1)

###remove Year
latency2=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+PC1pt+PC2pt+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latency2)

###remove Age
latency3=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+DaysWintHous+PC1pt+PC2pt+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latency3)

###remove Days in winter housing
latency4=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PC1pt+PC2pt+ScaledDevMaxIncLength+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latency4)

###remove DMIL
latency5=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PC1pt+PC2pt+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latency5)

###remove PC1
latency6=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PC2pt+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latency6)


######################



################response = Active righting response time### 
activerrtime1=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+PC1pt+PC2pt+ScaledDevMaxIncLength+as.factor(Year)+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activerrtime1)


##nothing left to remove after stepwise - nothing significant


################models with hatchling measurements instead of PC's###########


###TOTAL RR TIME###

rrmodelphen1=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen1)

###remove VT2
rrmodelphen2=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen2)

###remove DaysWintHous
rrmodelphen3=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen3)

###remove BD2
rrmodelphen4=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen4)

###remove Year
rrmodelphen5=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen5)

###remove Age at trials
rrmodelphen6=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen6)

###remove DMIL
rrmodelphen7=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+CW2+PL2+PW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen7)

###remove PL2
rrmodelphen8=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+CW2+PW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen8)

###remove Mass2
rrmodelphen9=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+CW2+PW2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen9)

###remove PW2
rrmodelphen10=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+CW2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen10)

###remove SCL2
rrmodelphen11=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+CW2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen11)

###remove CW2
rrmodelphen12=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen12)

###remove NV2
rrmodelphen13=lmer(log10(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(rrmodelphen13)

###LATENCY TO RIGHT###
latencyphen1=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen1)

###remove DMIL
latencyphen2=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+as.factor(Year)+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen2)

###remove Days wint hous (don't remove inc temp)
latencyphen3=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen3)

###remove BD2
latencyphen4=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+SCL2+CW2+PL2+PW2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen4)

###remove year
latencyphen5=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+SCL2+CW2+PL2+PW2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen5)

###remove PL2
latencyphen6=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+SCL2+CW2+PW2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen6)

###remove VT2
latencyphen7=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+SCL2+CW2+PW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen7)

###remove PW2
latencyphen8=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+SCL2+CW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen8)

###remove Age at trials
latencyphen9=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+CW2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen9)

###remove CW2
latencyphen10=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+NV2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(latencyphen10)

### ACTIVE RIGHTING RESPONSE TIME ###
activephen1=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+Mass2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen1)

###remove Mass2
activephen2=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+DaysWintHous+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen2)

###remove Days wint hous
activephen3=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+VT2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen3)

###remove VT2
activephen4=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+ScaledDevMaxIncLength+SCL2+CW2+PL2+PW2+BD2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen4)

###remove DMIL
activephen5=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+SCL2+CW2+PL2+PW2+BD2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen5)

###remove BD2
activephen6=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+SCL2+CW2+PL2+PW2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen6)

###remove CW2
activephen7=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+AgeAtTrials+as.factor(Year)+SCL2+PL2+PW2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen7)

###remove Age at trials
activephen8=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+as.factor(Year)+SCL2+PL2+PW2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen8)

###remove PW2
activephen9=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+as.factor(Year)+SCL2+PL2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen9)

###remove Year
activephen10=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+PL2+NV2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen10)

###remove NV2
activephen11=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+SCL2+PL2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen11)

###remove SCL2
activephen12=lmer(log10(ActiveRRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PL2+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)
summary(activephen12)













######################### Figures  #############################################
######################### Figures  #############################################
######################### Figures  #############################################
######################### Figures  #############################################

tempvals = seq(from = -4, to = 4, length.out = 100)

#tempdataf = expand.grid(WinterHous = unique(hatchlingdata1row$WinterHous),
#centeredIncTemp = tempvals)        ##this line not used below

###use the dataframes created below for model prediction lines in figure code

###CARAPACE LENGTH
emmipSCL2 = emmip(scl2d, WinterHous ~ centeredIncTemp, CIs = TRUE, cov.reduce = range, at = list(centeredIncTemp = tempvals), plotit=FALSE)

plot(yvar ~ xvar, data = emmipSCL2)

###CARAPACE WIDTH
emmipCW2 = emmip(cw2d, WinterHous ~ centeredIncTemp, CIs = TRUE, cov.reduce = range, at = list(centeredIncTemp = tempvals), plotit=FALSE)

plot(yvar ~ xvar, data = emmipCW2)

###PLASTRON LENGTH
emmipPL2 = emmip(pl2d, WinterHous ~ centeredIncTemp, CIs = TRUE, cov.reduce = range, at = list(centeredIncTemp = tempvals), plotit=FALSE)

plot(yvar ~ xvar, data = emmipPL2)

###PLASTRON WIDTH
emmipPW2 = emmip(pw2d, WinterHous ~ centeredIncTemp, CIs = TRUE, cov.reduce = range, at = list(centeredIncTemp = tempvals), plotit=FALSE)

plot(yvar ~ xvar, data = emmipPW2)

###BODY DEPTH
emmipBD2 = emmip(bd2d, WinterHous ~ centeredIncTemp, CIs = TRUE, cov.reduce = range, at = list(centeredIncTemp = tempvals), plotit=FALSE)

plot(yvar ~ xvar, data = emmipBD2)

###MASS
emmipMass2 = emmip(mass2d, WinterHous ~ centeredIncTemp, CIs = TRUE, cov.reduce = range, at = list(centeredIncTemp = tempvals), plotit=FALSE)

plot(yvar ~ xvar, data = emmipMass2)


#################################################################

# Now how to backtransform with emmeans if the response was log transformed

###############################################

######TOTAL RR TIME###########

# first with log transformation IN the model
rrmodel7=lmer(log(RRTimeMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)


# do emmeans without any change to the predictions for log transformation
RRtotal <- emmip(rrmodel7, WinterHous ~ centeredIncTemp3rows, CIs = TRUE,
                 cov.reduce = range,
                 at = list(centeredIncTemp3rows = tempvals),
                 plotit = FALSE)

## Note the message about results being on the log (not response) scale
head(RRtotal)                    

# now tell emmeans to backtransform when predicting emmeans
# Note from the help page of `ref_grid()`
## There is a subtle difference between specifying ‘⁠type = "response"⁠’ and ‘⁠regrid = "response"⁠’. While the summary statistics for the grid itself are the same, subsequent use in emmeans will yield different results if there is a response transformation or link function. With ‘⁠type = "response"⁠’, EMMs are computed by averaging together predictions on the linear-predictor scale and then back-transforming to the response scale; while with ‘⁠regrid = "response"⁠’, the predictions are already on the response scale so that the EMMs will be the arithmetic means of those response-scale predictions. To add further to the possibilities, geometric means of the response-scale predictions are obtainable via ‘⁠regrid = "log", type = "response"⁠’. See also the help page for regrid. 

## backtransform predictions from model scale THEN calculate EMMs on response scale
RRtotal_regrid <- emmip(rrmodel7, WinterHous ~ centeredIncTemp3rows, CIs = TRUE,
                        regrid = "response",
                        cov.reduce = range,
                        at = list(centeredIncTemp3rows = tempvals),
                        plotit = FALSE)

## Note the messages about the scale on which these predictions are made
head(RRtotal_regrid)


plot(yvar ~ xvar, data = RRtotal_regrid)  ###why are yvar values only from 3-7?

#######LATENCY TO RIGHT#########
# first with log transformation IN the model
latency6=lmer(log10(LatencyMin)~centeredIncTemp3rows+WinterHous+centeredIncTemp3rows*WinterHous+PC2pt+(1|HatchlingID)+(1|Clutch),data=RRdata3rows,na.action=na.omit)

# do emmeans without any change to the predictions for log transformation
Latency <- emmip(latency6, WinterHous ~ centeredIncTemp3rows, CIs = TRUE,
                 cov.reduce = range,
                 at = list(centeredIncTemp3rows = tempvals),
                 plotit = FALSE)

## Note the message about results being on the log (not response) scale
head(Latency)                    

# now tell emmeans to backtransform when predicting emmeans
# Note from the help page of `ref_grid()`
## There is a subtle difference between specifying ‘⁠type = "response"⁠’ and ‘⁠regrid = "response"⁠’. While the summary statistics for the grid itself are the same, subsequent use in emmeans will yield different results if there is a response transformation or link function. With ‘⁠type = "response"⁠’, EMMs are computed by averaging together predictions on the linear-predictor scale and then back-transforming to the response scale; while with ‘⁠regrid = "response"⁠’, the predictions are already on the response scale so that the EMMs will be the arithmetic means of those response-scale predictions. To add further to the possibilities, geometric means of the response-scale predictions are obtainable via ‘⁠regrid = "log", type = "response"⁠’. See also the help page for regrid. 

## backtransform predictions from model scale THEN calculate EMMs on response scale
Latency_regrid <- emmip(latency6, WinterHous ~ centeredIncTemp3rows, CIs = TRUE,
                        regrid = "response",
                        cov.reduce = range,
                        at = list(centeredIncTemp3rows = tempvals),
                        plotit = FALSE)

## Note the messages about the scale on which these predictions are made
head(Latency_regrid)


plot(yvar ~ xvar, data = Latency_regrid)






#############################
# Manuscript figures ########
#############################

###things to define first for all plots
xlab_in <- "Incubation Temperature (\u00B0C)"
xaxis <- seq(from = -6, to = 6, by = 2)  #use centered values here, re-label further down
# what should we set x-axis limits at
range(emmipSCL2$xvar, na.rm = TRUE)
xlim_in <- c(-4.5, 4.5)

degCexpr <- "(\u00B0C)"  #<-- degrees Celsius expression to paste in
TerrPtSymb <- 21
TerrPtCols <- c(bg = "#e86823", brd = "grey20")  #<-- "ORANGE"
TerrPtCx <- 0.8
AqSymb <- 22
AqPtCols <- c(bg = "#03244d", brd = "grey40")  #<-- "BLUE"
AqCx <- 0.8
jitfac <- 1.2  #<-- jitter factor
reglinewd <- 3.2  #<-- line width of all regression lines
ptLwd <- 0.6  #<-- line width of point border




#XXX for saving figures - use:
pdf(file = "Fig1_6morphologicalTraits_byIncubationAndOverwinter.pdf",
    width = 9, height = 12)

# use `par()` to set up some features of the entire figure
## `mfrow` designates the number of rows x columns to create in the figure
par(mfrow = c(3, 2),  #<-- (rows, columns)
    mar = c(5, 6.2, 4.5, 0.5),  #<-- space around each panel (bottom, left, top, right)
    mgp = c(3, 1, 0), #<-- adjustment of axis title, labels, and line
    cex.axis = 1.2, cex.lab = 1.3)  #<-- scaling for axis labels and axis title

################################
# CARAPACE LENGTH 
################################
plot(yvar ~ xvar, data = emmipSCL2, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(26, 49),
     xlab = xlab_in,
     ylab = "Carapace length (mm)")

# points first (to put in background)
## Terrestrial first
#points(SCL2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row, ##use this line to jitter
points(SCL2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(SCL2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row,
points(SCL2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = emmipSCL2, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = emmipSCL2, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=emmipSCL2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=emmipSCL2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=emmipSCL2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=emmipSCL2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])


# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(24, 52, 4))

mtext(text = expression((bold(A))),
      side = 3, line = 0.2, adj = -0.2, cex = 1.1)

################################
# LEGEND
################################
legend(x = "topleft",
       inset = c(-0.02, -0.2),  # tweak to shift above plotting area
       xpd = NA,                 # allow drawing outside plot region
       legend = c("Terrestrial", "Aquatic"),
       pch = c(TerrPtSymb, AqSymb),
       pt.bg = c(TerrPtCols["bg"], AqPtCols["bg"]),
       col = c(TerrPtCols["brd"], AqPtCols["brd"]),
       pt.cex = c(TerrPtCx, AqCx),
       pt.lwd = ptLwd,
       bty = "o",                # adds a box around the legend
       box.lwd = 1,              # thickness of legend border
       title = "Overwinter Treatment",
       horiz = TRUE,             # makes it horizontal (optional)
       cex = 1.1)                # text size



################################
# CARAPACE WIDTH
################################
plot(yvar ~ xvar, data = emmipCW2, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(27, 49),
     xlab = xlab_in,
     ylab = "Carapace width (mm)")

# points first (to put in background)
## Terrestrial first
#points(CW2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row,
points(CW2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(CW2 ~ jitter(centeredIncTemp, jitfac), data = hatchlingdata1row,
points(CW2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = emmipCW2, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = emmipCW2, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=emmipCW2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=emmipCW2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=emmipCW2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=emmipCW2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(24, 52, 4))

mtext(text = expression((bold(B))),
      side = 3, line = 0.2, adj = -0.2, cex = 1.1)

################################
# PLASTRON LENGTH
################################
plot(yvar ~ xvar, data = emmipPL2, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(25, 45),
     xlab = xlab_in,
     ylab = "Plastron length (mm)")

# points first (to put in background)
## Terrestrial first
#points(PL2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row,
points(PL2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(PL2 ~ jitter(centeredIncTemp, jitfac), data = hatchlingdata1row,
points(PL2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = emmipPL2, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = emmipPL2, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=emmipPL2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=emmipPL2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=emmipPL2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=emmipPL2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(24, 48, 4))

mtext(text = expression((bold(C))),
      side = 3, line = 0.2, adj = -0.2, cex = 1.1)

################################
# PLASTRON WIDTH
################################
plot(yvar ~ xvar, data = emmipPW2, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(21, 38),
     xlab = xlab_in,
     ylab = "Plastron width (mm)")

# points first (to put in background)
## Terrestrial first
#points(PW2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row,
points(PW2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(PW2 ~ jitter(centeredIncTemp, jitfac), data = hatchlingdata1row,
points(PW2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = emmipPW2, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = emmipPW2, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=emmipPW2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=emmipPW2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=emmipPW2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=emmipPW2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(18, 42, 4))

mtext(text = expression((bold(D))),
      side = 3, line = 0.2, adj = -0.2, cex = 1.1)

################################
# BODY DEPTH
################################
plot(yvar ~ xvar, data = emmipBD2, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(13, 23),
     xlab = xlab_in,
     ylab = "Body depth (mm)")

# points first (to put in background)
## Terrestrial first
#points(BD2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row,
points(BD2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(BD2 ~ jitter(centeredIncTemp, jitfac), data = hatchlingdata1row,
points(BD2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = emmipBD2, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = emmipBD2, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=emmipBD2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=emmipBD2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=emmipBD2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=emmipBD2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(10, 26, 2))

mtext(text = expression((bold(E))),
      side = 3, line = 0.2, adj = -0.2, cex = 1.1)

################################
# MASS
################################
plot(yvar ~ xvar, data = emmipMass2, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(5, 25),
     xlab = xlab_in,
     ylab = "Mass (g)")

# points first (to put in background) 
## Terrestrial first
#points(Mass2 ~ jitter(centeredIncTemp,jitfac), data = hatchlingdata1row,
points(Mass2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(Mass2 ~ jitter(centeredIncTemp, jitfac), data = hatchlingdata1row,
points(Mass2 ~ centeredIncTemp, data = hatchlingdata1row,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = emmipMass2, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = emmipMass2, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=emmipMass2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=emmipMass2, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=emmipMass2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=emmipMass2, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(2, 28, 4))

mtext(text = expression((bold(F))),
      side = 3, line = 0.2, adj = -0.2, cex = 1.1)



dev.off()









#XXX for saving figures - use:
pdf(file = "Fig3_RightingResponses_byIncubationAndOverwinter.pdf",
    width = 9, height = 4)

# use `par()` to set up some features of the entire figure
## `mfrow` designates the number of rows x columns to create in the figure
par(mfrow = c(1, 2),  #<-- (rows, columns)
    mar = c(5, 6.2, 4.5, 0.5),  #<-- space around each panel (bottom, left, top, right)
    mgp = c(3, 1, 0), #<-- adjustment of axis title, labels, and line
    cex.axis = 1.1, cex.lab = 1.2)  #<-- scaling for axis labels and axis title


################################
# TOTAL RIGHTING RESPONSE TIME
################################
plot(yvar ~ xvar, data = RRtotal_regrid, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(0, 20),
     xlab = xlab_in,
     ylab = "Righting response time (min)")

# points first (to put in background)
## Terrestrial first
#points(RRTimeMin ~ jitter(centeredIncTemp3rows,jitfac), data = RRdata3rows,
points(RRTimeMin ~ centeredIncTemp3rows, data = RRdata3rows,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(RRTimeMin ~ jitter(centeredIncTemp3rows, jitfac), data = RRdata3rows,
points(RRTimeMin ~ centeredIncTemp3rows, data = RRdata3rows,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = RRtotal_regrid, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = RRtotal_regrid, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=RRtotal_regrid, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=RRtotal_regrid, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=RRtotal_regrid, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=RRtotal_regrid, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(-2, 24, 4))

mtext(text = expression((bold(A))),
      side = 3, line = 0.9, adj = -0.2, cex = 1.1)

################################
# LEGEND
################################

legend(x = "topleft",
       inset = c(-0.02, -0.35),  # tweak to shift above plotting area
       xpd = NA,                 # allow drawing outside plot region
       legend = c("Terrestrial", "Aquatic"),
       pch = c(TerrPtSymb, AqSymb),
       pt.bg = c(TerrPtCols["bg"], AqPtCols["bg"]),
       col = c(TerrPtCols["brd"], AqPtCols["brd"]),
       pt.cex = c(TerrPtCx, AqCx),
       pt.lwd = ptLwd,
       bty = "o",                # adds a box around the legend
       box.lwd = 1,              # thickness of legend border
       title = "Overwinter Treatment",
       horiz = TRUE,             # makes it horizontal (optional)
       cex = 1.1)                # text size


###################

################################
# LATENCY TO RIGHT
################################
plot(yvar ~ xvar, data = Latency_regrid, type = "n", #<-- just set up
     axes = FALSE,  #<-- make our own fancy ones
     xlim = xlim_in, ylim = c(0, 20),
     xlab = xlab_in,
     ylab = "Latency to right (min)")

# points first (to put in background) 
## Terrestrial first
#points(LatencyMin ~ jitter(centeredIncTemp3rows,jitfac), data = RRdata3rows,
points(LatencyMin ~ centeredIncTemp3rows, data = RRdata3rows,
       subset = WinterHous == "Terrestrial",
       pch = TerrPtSymb, bg = TerrPtCols["bg"], col = TerrPtCols["brd"],
       cex = TerrPtCx, lwd = ptLwd)
## Aquatic second
#points(LatencyMin ~ jitter(centeredIncTemp3rows, jitfac), data = RRdata3rows,
points(LatencyMin ~ centeredIncTemp3rows, data = RRdata3rows,
       subset = WinterHous == "Aquatic",
       pch = AqSymb, bg = AqPtCols["bg"], col = AqPtCols["brd"],
       cex = AqCx, lwd = ptLwd)

# Lines from model
## Terrestrial first
lines(yvar ~ xvar, data = Latency_regrid, subset = tvar == "Terrestrial",
      lwd = reglinewd * 1.2, col = TerrPtCols["bg"])  #<-- make thicker so not covered
## Aquatic second
lines(yvar ~ xvar, data = Latency_regrid, subset = tvar == "Aquatic",
      lwd = reglinewd, col = AqPtCols["bg"])

##Confidence Intervals
#Terrestrial upper
lines(UCL~xvar,data=Latency_regrid, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Terrestrial lower
lines(LCL~xvar,data=Latency_regrid, subset = tvar == "Terrestrial",
      lwd = 1.0 * 1.2, lty="dashed", col = TerrPtCols["bg"])

#Aquatic upper
lines(UCL~xvar,data=Latency_regrid, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

#Aquatic lower
lines(LCL~xvar,data=Latency_regrid, subset = tvar == "Aquatic",
      lwd = 1.0, lty="dashed", col = AqPtCols["bg"])

# X-axis followed by Y-axis
axis(1, xaxis,labels=seq(from = 22, to = 34, by = 2))  
axis(2, seq(-2, 24, 4))

mtext(text = expression((bold(B))),
      side = 3, line = 0.9, adj = -0.2, cex = 1.1)



dev.off() 

#END OF FIGUREs#
#############################################################








#######################Dharma inspecting assumptions of linear models


simulationOutput <- simulateResiduals(fittedModel = eggsurv)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = scl1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = cw1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = pl1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = pw1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = bd1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = nv1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = vt1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = mass1)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = cw2d)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = pw2d)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = bd2d)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = mass2d)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = growthDd)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = rrmodel7)
plot(simulationOutput)

testQuantiles(simulationOutput)

###
simulationOutput <- simulateResiduals(fittedModel = latency6)
plot(simulationOutput)

testQuantiles(simulationOutput)



###DHARMa plot fixes 

##########INCUBATION DURATION MODEL###################

##run as glmm (didn't work)
incdur=glmer(IncLength~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=inclengths,na.action=na.omit,family=poisson)
summary(incdur)

simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)

###glmm and quadratic term (didn't work)
incdur=glmer(IncLength~centeredIncTemp+centeredIncTemp^2+ScaledDevMaxIncLength+(1|Clutch),data=inclengths,na.action=na.omit,family=poisson)
summary(incdur)

simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)

###try log transformation (didn't work)
incdur=lmer(log(IncLength)~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=inclengths,na.action=na.omit)
summary(incdur)

simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)

###try log transformation and quadratic term for inc temp (didn't work)
incdur=lmer(log(IncLength)~centeredIncTemp+centeredIncTemp^2+ScaledDevMaxIncLength+(1|Clutch),data=inclengths,na.action=na.omit)
summary(incdur)

simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)


###look at outliers (does this output row numbers? i'm assuming)
outliers(simulationOutput, lowerQuantile = 0, upperQuantile = 1,
         return = c("index", "logical"))

###remove outliers, lines 41, 42, 46, 47 from dataframe
rows_to_remove = c(41, 42, 46, 47)
inclengths_no_outliers=inclengths[-rows_to_remove,]

###try again with new dataframe with rows removed
incdur=lmer(log(IncLength)~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=inclengths_no_outliers,na.action=na.omit)
summary(incdur)

simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)

###try again with glm, no outliers, 
incdur=glmer(IncLength~centeredIncTemp+ScaledDevMaxIncLength+(1|Clutch),data=inclengths_no_outliers,na.action=na.omit,family=poisson)
summary(incdur)

simulationOutput <- simulateResiduals(fittedModel = incdur)
plot(simulationOutput)

##########SCL2 MODEL###################

###original
scl2d=lmer(SCL2~centeredIncTemp+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

###try adding DMIL variable (didn't help)
scl2d=lmer(SCL2~centeredIncTemp+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

###try log10 transformation (didn't help much)
scl2d=lmer(log10(SCL2)~centeredIncTemp+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

###try log transformation (didn't help much, same as log10)
scl2d=lmer(log(SCL2)~centeredIncTemp+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

###try adding log transf. and quadratic inc temp term
scl2d=lmer(log(SCL2)~centeredIncTemp+centeredIncTemp^2+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

###return outliers
scl2d=lmer(SCL2~centeredIncTemp+WinterHous+SCL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(scl2d)

simulationOutput <- simulateResiduals(fittedModel = scl2d)
plot(simulationOutput)

outliers(simulationOutput, lowerQuantile = 0, upperQuantile = 1,
         return = c("index", "logical"))

##########PL2 MODEL###################

###original
pl2d=lmer(PL2~centeredIncTemp+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)

###add DMIL variable (didn't help)
pl2d=lmer(PL2~centeredIncTemp+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)

###log10 transformation (helped some)
pl2d=lmer(log10(PL2)~centeredIncTemp+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)

###quadratic term (didn't help as much as log10 transf)
pl2d=lmer(PL2~centeredIncTemp+centeredIncTemp^2+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)


###quadratic term and log10 transf (not any better than log10transf)
pl2d=lmer(log(PL2)~centeredIncTemp+centeredIncTemp^2+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)


###return outliers of original
pl2d=lmer(PL2~centeredIncTemp+WinterHous+PL1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(pl2d)

simulationOutput <- simulateResiduals(fittedModel = pl2d)
plot(simulationOutput)

outliers(simulationOutput, lowerQuantile = 0, upperQuantile = 1,
         return = c("index", "logical"))

##########Mass2 MODEL###################

###original
mass2d=lmer(Mass2~centeredIncTemp+WinterHous+Mass1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(mass2d)

simulationOutput <- simulateResiduals(fittedModel = mass2d)
plot(simulationOutput)

###add DMIL (didn't help)
mass2d=lmer(Mass2~centeredIncTemp+WinterHous+Mass1+centeredIncTemp*WinterHous+DaysWintHous+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(mass2d)

###quadratic term (didn't help)
mass2d=lmer(Mass2~centeredIncTemp+WinterHous+Mass1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(mass2d)

simulationOutput <- simulateResiduals(fittedModel = mass2d)
plot(simulationOutput)

###log10 transformation (helped)
mass2d=lmer(log10(Mass2)~centeredIncTemp+WinterHous+Mass1+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(mass2d)

simulationOutput <- simulateResiduals(fittedModel = mass2d)
plot(simulationOutput)

###################growth model######################

###original
growthDd=lmer(PC1growthD~centeredIncTemp+WinterHous+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(growthDd)

simulationOutput <- simulateResiduals(fittedModel = growthDd)
plot(simulationOutput)

###DMIL - didn't help
growthDd=lmer(PC1growthD~centeredIncTemp+WinterHous+centeredIncTemp*WinterHous+DaysWintHous+ScaledDevMaxIncLength+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(growthDd)

simulationOutput <- simulateResiduals(fittedModel = growthDd)
plot(simulationOutput)


###quadratic - didn't help
growthDd=lmer(PC1growthD~centeredIncTemp+centeredIncTemp^2+WinterHous+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(growthDd)

simulationOutput <- simulateResiduals(fittedModel = growthDd)
plot(simulationOutput)

###log10/log - can't do log on PC scores because they are negative

###outliers
growthDd=lmer(PC1growthD~centeredIncTemp+WinterHous+centeredIncTemp*WinterHous+DaysWintHous+(1|Clutch),data=hatchlingdata1row,na.action=na.omit)
summary(growthDd)

simulationOutput <- simulateResiduals(fittedModel = growthDd)
plot(simulationOutput)

outliers(simulationOutput, lowerQuantile = 0, upperQuantile = 1,
         return = c("index", "logical"))













