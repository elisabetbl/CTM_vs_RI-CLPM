#'This script aims to compare how the relationship between cortical thickness and problem behavior
#'can be captured using a continuous time model with the R-package ctsem as compared to a 
#'(random intercept-)cross lagged panel model.
#'For a detailed description of the model and the different model set-ups, please read the original
#'paper and updates on the ctsem package manual (https://www.jstatsoft.org/article/view/v077i05/0).

#'Outline of the script:
#'1.     Continuous time model (CTM)
#'1.1    Model including time-dependent predictors only
#'1.1.1  Maximum likelihood, no priors specified
#'1.1.2  Maximum likelihood, including priors
#'1.2    Model including time-independent predictors
#'1.2.1  Maximum likelihood, no priors specified
#'1.2.2  Maximum likelihood, including priors
#'1.3    Time-independent predictors residualized before CTM
#'1.3.1  Maximum likelihood, no priors specified
#'1.3.2  Maximum likelihood, including priors
#'2.     Cross-lagged panel model (CLPM)
#'2.1    No random intercept (no distinction of within and between participant effects)
#'2.2    Random intercept CLPM

#load necessary libraries (if you have not installed those locally, you will need to do so before loading them)
library(dplyr)
library(tidyr)
library(ctsem)
library(xlsx)
library(lavaan)

#data
setwd("CTM_vs_RI-CLPM/data/")
df1 <- read.csv("CLPM_Simulation_Data_6y_9y_13y_Version_1.csv", sep = "|")
df2 <- read.csv("CLPM_Simulation_Data_6y_9y_13y_Version_2.csv", sep = "|")
df3 <- read.csv("CLPM_Simulation_Data_6y_9y_13y_Version_3.csv", sep = "|")
df4 <- read.csv("CLPM_Simulation_Data_6y_9y_13y_Version_4.csv", sep = "|")
df5 <- read.csv("CLPM_Simulation_Data_6y_9y_13y_Version_5.csv", sep = "|")

#Specify which dataset you want to run the analyses on
df <- df1
whichdf <- "df1"

####################################
##### 1. CONTINUOUS TIME MODEL #####
####################################

############################################################
##### 1.1 CTM INCLUDING TIME-DEPENDENT PREDICTORS ONLY #####
############################################################

#order data
df_order_nocovars <- select(df, ID, Beh_06, CT_06, Beh_09, CT_09, Beh_13, CT_13, Age_06, Age_09, Age_13)
colnames(df_order_nocovars) <- c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "T0", "T1", "T2")

#Change time variables
#Mean center T0
df_order_nocovars$MC_T0 <- df_order_nocovars$T0 - mean(df_order_nocovars$T0)
df_order_nocovars$time_T0 <- df_order_nocovars$MC_T0 - min(df_order_nocovars$MC_T0)
#Create time difference
df_order_nocovars$dT1 <- df_order_nocovars$T1 - df_order_nocovars$T0
df_order_nocovars$dT2 <- df_order_nocovars$T2 - df_order_nocovars$T1
#Calculate T1 & T2 
df_order_nocovars$time_T1 <- df_order_nocovars$time_T0 + df_order_nocovars$dT1
df_order_nocovars$time_T2 <- df_order_nocovars$time_T1 + df_order_nocovars$dT2

#Transform DP & CT To approach normality and then Z-transform
df_order_nocovars$Beh_T0 <- scale(sqrt(df_order_nocovars$Beh_T0))
df_order_nocovars$Beh_T1 <- scale(sqrt(df_order_nocovars$Beh_T1))
df_order_nocovars$Beh_T2 <- scale(sqrt(df_order_nocovars$Beh_T2))
df_order_nocovars$CT_T0 <- scale(df_order_nocovars$CT_T0)
df_order_nocovars$CT_T1 <- scale(df_order_nocovars$CT_T1)
df_order_nocovars$CT_T2 <- scale(df_order_nocovars$CT_T2)

#Select columns
df_order_nocovars_com <- select(df_order_nocovars, c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "time_T0", "time_T1", "time_T2"))

#create dataset appropriate for ctsem
df_long_nocovars <- pivot_longer(df_order_nocovars_com, -c(ID), names_to = c("variable", "timepoint"), values_to = c("value"), names_pattern = "(.*)_(.*)")
df_wide_nocovars <- pivot_wider(df_long_nocovars, names_from = c("variable"), values_from = c("value"))
df_fin_nocovars <- select(df_wide_nocovars, c(ID, Beh, CT, time))
df_fin_nocovars <- as.matrix(df_fin_nocovars)

#Random intercept continuous time model
model_nocovars <- ctModel(type = "stanct",
                          time = "time", id = "ID",
                          n.latent = 2, latentNames = c("Beh", "CT"),
                          n.manifest = 2, manifestNames = c("Beh", "CT"),
                          n.TDpred = 0,
                          n.TIpred = 0,
                          LAMBDA = diag(2))

###############################################
##### 1.1.1 MAXIMUM LIKELIHOOD, NO PRIORS #####
###############################################

fit_ML_nocovars<-ctStanFit(datalong = df_fin_nocovars, ctstanmodel = model_nocovars, nopriors = T, cores = "maxneeded")

#check results
#Mean dT1 = 3 years, mean 3T2 = 4 years, time intervals are chosen based on this
results_ML_nocovars_T1 <- summary(fit_ML_nocovars, timeinterval = 1)
results_ML_nocovars_T3 <- summary(fit_ML_nocovars, timeinterval = 3)
results_ML_nocovars_T7 <- summary(fit_ML_nocovars, timeinterval = 7)

coefs_ML_nocovars_T1 <- results_ML_nocovars_T1[[4]][c(33:36),c(1:7)]
coefs_ML_nocovars_T3 <- results_ML_nocovars_T3[[4]][c(33:36),c(1:7)]
coefs_ML_nocovars_T7 <- results_ML_nocovars_T7[[4]][c(33:36),c(1:7)]

setwd("CTM_vs_RI-CLPM/output/")

write.xlsx(coefs_ML_nocovars_T1, "results_ML_nocovars.xlsx", sheetName = paste0(whichdf,"_T1"), append = T)
write.xlsx(coefs_ML_nocovars_T3, "results_ML_nocovars.xlsx", sheetName = paste0(whichdf,"_T3"), append = T)
write.xlsx(coefs_ML_nocovars_T7, "results_ML_nocovars.xlsx", sheetName = paste0(whichdf,"_T7"), append = T)


#Plot results
jpeg(file=paste0(whichdf,"_plot_ML_nocovars.jpeg"))
ctStanDiscretePars(fit_ML_nocovars, plot = T, indices = "all", subjects = "all", times = seq(from = 0, to = 8, by = 1))
dev.off()



###############################################
##### 1.1.2 MAXIMUM LIKELIHOOD AND PRIORS #####
###############################################

fit_mid_nocovars <-ctStanFit(datalong = df_fin_nocovars, ctstanmodel = model_nocovars, nopriors = F, cores = "maxneeded")

#check results
#Mean dT1 = 3 years, mean 3T2 = 4 years, time intervals are chosen based on this
results_mid_nocovars_T1 <- summary(fit_mid_nocovars, timeinterval = 1)
results_mid_nocovars_T3 <- summary(fit_mid_nocovars, timeinterval = 3)
results_mid_nocovars_T7 <- summary(fit_mid_nocovars, timeinterval = 7)

coefs_mid_nocovars_T1 <- results_mid_nocovars_T1[[6]][c(33:36),c(1:7)]
coefs_mid_nocovars_T3 <- results_mid_nocovars_T3[[6]][c(33:36),c(1:7)]
coefs_mid_nocovars_T7 <- results_mid_nocovars_T7[[6]][c(33:36),c(1:7)]

write.xlsx(coefs_mid_nocovars_T1, "results_mid_nocovars.xlsx", sheetName = paste0(whichdf,"_T1"), append = T)
write.xlsx(coefs_mid_nocovars_T3, "results_mid_nocovars.xlsx", sheetName = paste0(whichdf,"_T3"), append = T)
write.xlsx(coefs_mid_nocovars_T7, "results_mid_nocovars.xlsx", sheetName = paste0(whichdf,"_T7"), append = T)


#Plot results
jpeg(file=paste0(whichdf,"_plot_mid_nocovars.jpeg"))
ctStanDiscretePars(fit_mid_nocovars, plot = T, indices = "all", subjects = "all", times = seq(from = 0, to = 8, by = 1))
dev.off()




#########################################################
##### 1.2 CTM INCLUDING TIME INDEPENDENT PREDICTORS #####
#########################################################

#Preparing the data so that ctsem recognizes it
#order and rename data
df_order <- select(df, ID, Beh_06, CT_06, Beh_09, CT_09, Beh_13, CT_13, Age_06, Age_09, Age_13, sex, education, ethnicity)
colnames(df_order) <- c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "T0", "T1", "T2", "sex", "education", "ethnicity")

#Change time variables
#Mean center T0
df_order$MC_T0 <- df_order$T0 - mean(df_order$T0)
df_order$time_T0 <- df_order$MC_T0 - min(df_order$MC_T0)
#Create time difference
df_order$dT1 <- df_order$T1 - df_order$T0
df_order$dT2 <- df_order$T2 - df_order$T1
#Calculate T1 & T2 
df_order$time_T1 <- df_order$time_T0 + df_order$dT1
df_order$time_T2 <- df_order$time_T1 + df_order$dT2

#Scale Time Dependent predictors in a way that 0 represents no effect
#Dummy coding of all categorical variables (sex, education & ethnicity)
df_order$girl <- ifelse(df_order$sex == 1, 1, 0)
df_order$ed_low <- ifelse(df_order$education == 0, 1, 0)
df_order$ed_mid <- ifelse(df_order$education == 1, 1, 0)
df_order$ed_high <- ifelse(df_order$education == 2, 1, 0)
df_order$eth_wes <- ifelse(df_order$ethnicity == 0, 1, 0)
df_order$eth_dutch <- ifelse(df_order$ethnicity == 1, 1, 0)
df_order$eth_nonwes <- ifelse(df_order$ethnicity == 2, 1, 0)

#Transform DP & CT To approach normality and then scale all time dependent predictors
df_order$Beh_T0 <- scale(sqrt(df_order$Beh_T0))
df_order$Beh_T1 <- scale(sqrt(df_order$Beh_T1))
df_order$Beh_T2 <- scale(sqrt(df_order$Beh_T2))
df_order$CT_T0 <- scale(df_order$CT_T0)
df_order$CT_T1 <- scale(df_order$CT_T1)
df_order$CT_T2 <- scale(df_order$CT_T2)

#Select columns
df_order_com <- select(df_order, c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "time_T0", "time_T1", "time_T2", "sex", "ed_low", "ed_mid", "ed_high", "eth_wes", "eth_dutch", "eth_nonwes"))

#create dataset appropriate for ctsem
df_long <- pivot_longer(df_order_com, -c(ID, sex, ed_low, ed_mid, ed_high, eth_wes, eth_dutch, eth_nonwes), names_to = c("variable", "timepoint"), values_to = c("value"), names_pattern = "(.*)_(.*)")
df_wide <- pivot_wider(df_long, names_from = c("variable"), values_from = c("value"))
df_fin <- select(df_wide, c(ID, Beh, CT, time, sex, ed_low, ed_mid, ed_high, eth_wes, eth_dutch, eth_nonwes))

#Random intercept continuous time model
model <- ctModel(type = "stanct",
                 time = "time", id = "ID",
                 n.latent = 2, latentNames = c("Beh", "CT"),
                 n.manifest = 2, manifestNames = c("Beh", "CT"),
                 n.TDpred = 0,
                 n.TIpred = 7, TIpredNames = c("sex", "ed_low", "ed_mid", "ed_high", "eth_wes", "eth_dutch", "eth_nonwes"),
                 LAMBDA = diag(2))

###############################################
##### 1.2.1 MAXIMUM LIKELIHOOD, NO PRIORS #####
###############################################

fit_ML<-ctStanFit(datalong = df_fin, ctstanmodel = model, nopriors = T, cores = "maxneeded")

#check results
#Mean dT1 = 3 years, mean 3T2 = 4 years, time intervals are chosen based on this
results_ML_T1 <- summary(fit_ML, timeinterval = 1)
results_ML_T3 <- summary(fit_ML, timeinterval = 3)
results_ML_T7 <- summary(fit_ML, timeinterval = 7)

coefs_ML_T1 <- results_ML_T1[[5]][c(33:36),c(1:7)]
coefs_ML_T3 <- results_ML_T3[[5]][c(33:36),c(1:7)]
coefs_ML_T7 <- results_ML_T7[[5]][c(33:36),c(1:7)]

write.xlsx(coefs_ML_T1, "results_ML.xlsx", sheetName = paste0(whichdf,"_T1"), append = T)
write.xlsx(coefs_ML_T3, "results_ML.xlsx", sheetName = paste0(whichdf,"_T3"), append = T)
write.xlsx(coefs_ML_T7, "results_ML.xlsx", sheetName = paste0(whichdf,"_T7"), append = T)


#Plot results
jpeg(file=paste0(whichdf,"_plot_ML.jpeg"))
ctStanDiscretePars(fit_ML, plot = T, indices = "all", subjects = "all", times = seq(from = 0, to = 8, by = 1))
dev.off()



###############################################
##### 1.1.2 MAXIMUM LIKELIHOOD AND PRIORS #####
###############################################

fit_mid <-ctStanFit(datalong = df_fin, ctstanmodel = model, nopriors = F, cores = "maxneeded")

#check results
#Mean dT1 = 3 years, mean dT2 = 4 years, time intervals are chosen based on this
results_mid_T1 <- summary(fit_mid, timeinterval = 1)
results_mid_T3 <- summary(fit_mid, timeinterval = 3)
results_mid_T7 <- summary(fit_mid, timeinterval = 7)

coefs_mid_T1 <- results_mid_T1[[7]][c(33:36),c(1:7)]
coefs_mid_T3 <- results_mid_T3[[7]][c(33:36),c(1:7)]
coefs_mid_T7 <- results_mid_T7[[7]][c(33:36),c(1:7)]

write.xlsx(coefs_mid_T1, "results_mid.xlsx", sheetName = paste0(whichdf,"_T1"), append = T)
write.xlsx(coefs_mid_T3, "results_mid.xlsx", sheetName = paste0(whichdf,"_T3"), append = T)
write.xlsx(coefs_mid_T7, "results_mid.xlsx", sheetName = paste0(whichdf,"_T7"), append = T)


#Plot results
jpeg(file=paste0(whichdf,"_plot_mid.jpeg"))
ctStanDiscretePars(fit_mid, plot = T, indices = "all", subjects = "all", times = seq(from = 0, to = 8, by = 1))
dev.off()





######################################################################
##### 1.3 TIME-INDEPENDENT PREDICTORS REGRESSED OUT PRIOR TO CTM #####
######################################################################

#order data
df_order <- select(df, ID, Beh_06, CT_06, Beh_09, CT_09, Beh_13, CT_13, Age_06, Age_09, Age_13, sex, education, ethnicity)
colnames(df_order) <- c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "T0", "T1", "T2", "sex", "education", "ethnicity")

#Change time variables
#Mean center T0
df_order$MC_T0 <- df_order$T0 - mean(df_order$T0)
df_order$time_T0 <- df_order$MC_T0 - min(df_order$MC_T0)
#Create time difference
df_order$dT1 <- df_order$T1 - df_order$T0
df_order$dT2 <- df_order$T2 - df_order$T1
#Calculate T1 & T2 
df_order$time_T1 <- df_order$time_T0 + df_order$dT1
df_order$time_T2 <- df_order$time_T1 + df_order$dT2

#Transform DP & CT To approach normality and then Z-transform
df_order$Beh_T0 <- scale(sqrt(df_order$Beh_T0))
df_order$Beh_T1 <- scale(sqrt(df_order$Beh_T1))
df_order$Beh_T2 <- scale(sqrt(df_order$Beh_T2))
df_order$CT_T0 <- scale(df_order$CT_T0)
df_order$CT_T1 <- scale(df_order$CT_T1)
df_order$CT_T2 <- scale(df_order$CT_T2)

#Select vars needed
df_order_com <- select(df_order, c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "time_T0", "time_T1", "time_T2", "sex", "education", "ethnicity"))


vars <- c("Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2")

#NOTE!! THIS OVERWRITES ORIGINAL VALUES!!
for (x in vars){
  f <- paste0(x, " ~ sex + education + ethnicity")
  df_order_com[x] <- residuals(lm(as.formula(f), data = df_order_com, na.action = na.exclude))
}

# Set up if you don't want to overwrite original values
#for (x in vars){
#  t <- df_order_com[x]
#  f <- paste0(x, " ~ sex + education + ethnicity")
#  colname <- paste0(x,"_resid")
#  df_order_com[,colname] <- residuals(lm(as.formula(f), data = df_order_com, na.action = na.exclude))
#}

#now proceed as if there are no covars/TIpreds
#Select vars needed
df_order_com_resid <- select(df_order_com, c("ID", "Beh_T0", "CT_T0", "Beh_T1", "CT_T1", "Beh_T2", "CT_T2", "time_T0", "time_T1", "time_T2"))

#create dataset appropriate for ctsem
df_long_resid <- pivot_longer(df_order_com_resid, -c(ID), names_to = c("variable", "timepoint"), values_to = c("value"), names_pattern = "(.*)_(.*)")
df_wide_resid <- pivot_wider(df_long_resid, names_from = c("variable"), values_from = c("value"))
df_fin_resid <- select(df_wide_resid, c(ID, Beh, CT, time))
df_fin_resid <- as.matrix(df_fin_resid)

#Random intercept continuous time model
model_resid <- ctModel(type = "stanct",
                       time = "time", id = "ID",
                       n.latent = 2, latentNames = c("Beh", "CT"),
                       n.manifest = 2, manifestNames = c("Beh", "CT"),
                       n.TDpred = 0,
                       n.TIpred = 0,
                       LAMBDA = diag(2))

###############################################
##### 1.3.1 MAXIMUM LIKELIHOOD, NO PRIORS #####
###############################################

fit_ML_resid<-ctStanFit(datalong = df_fin_resid, ctstanmodel = model_resid, nopriors = T, cores = "maxneeded")

#check results
#Mean dT1 = 3 years, mean 3T2 = 4 years, time intervals are chosen based on this
results_ML_resid_T1 <- summary(fit_ML_resid, timeinterval = 1)
results_ML_resid_T3 <- summary(fit_ML_resid, timeinterval = 3)
results_ML_resid_T7 <- summary(fit_ML_resid, timeinterval = 7)

coefs_ML_resid_T1 <- results_ML_resid_T1[[4]][c(33:36),c(1:7)]
coefs_ML_resid_T3 <- results_ML_resid_T3[[4]][c(33:36),c(1:7)]
coefs_ML_resid_T7 <- results_ML_resid_T7[[4]][c(33:36),c(1:7)]

write.xlsx(coefs_ML_resid_T1, "results_ML_resid.xlsx", sheetName = paste0(whichdf,"_T1"), append = T)
write.xlsx(coefs_ML_resid_T3, "results_ML_resid.xlsx", sheetName = paste0(whichdf,"_T3"), append = T)
write.xlsx(coefs_ML_resid_T7, "results_ML_resid.xlsx", sheetName = paste0(whichdf,"_T7"), append = T)


#Plot results
jpeg(file=paste0(whichdf,"_plot_ML_resid.jpeg"))
ctStanDiscretePars(fit_ML_resid, plot = T, indices = "all", subjects = "all", times = seq(from = 0, to = 8, by = 1))
dev.off()



###############################################
##### 1.3.1 MAXIMUM LIKELIHOOD AND PRIORS #####
###############################################

fit_mid_resid <-ctStanFit(datalong = df_fin_resid, ctstanmodel = model_resid, nopriors = F, cores = "maxneeded")

#check results
#Mean dT1 = 3 years, mean 3T2 = 4 years, time intervals are chosen based on this
results_mid_resid_T1 <- summary(fit_mid_resid, timeinterval = 1)
results_mid_resid_T3 <- summary(fit_mid_resid, timeinterval = 3)
results_mid_resid_T7 <- summary(fit_mid_resid, timeinterval = 7)

coefs_mid_resid_T1 <- results_mid_resid_T1[[6]][c(33:36),c(1:7)]
coefs_mid_resid_T3 <- results_mid_resid_T3[[6]][c(33:36),c(1:7)]
coefs_mid_resid_T7 <- results_mid_resid_T7[[6]][c(33:36),c(1:7)]

write.xlsx(coefs_mid_resid_T1, "results_mid_resid.xlsx", sheetName = paste0(whichdf,"_T1"), append = T)
write.xlsx(coefs_mid_resid_T3, "results_mid_resid.xlsx", sheetName = paste0(whichdf,"_T3"), append = T)
write.xlsx(coefs_mid_resid_T7, "results_mid_resid.xlsx", sheetName = paste0(whichdf,"_T7"), append = T)


#Plot results
jpeg(file=paste0(whichdf,"plot_mid_resid.jpeg"))
ctStanDiscretePars(fit_mid_resid, plot = T, indices = "all", subjects = "all", times = seq(from = 0, to = 8, by = 1))
dev.off()





##############################################
##### 2. CROSS-LAGGED PANEL MODEL (CLPM) #####
##############################################

###################################
##### 2.1 NO RANDOM INTERCEPT #####
###################################

#Model set-up performed like Ellen Hamaker
#https://github.com/ellenhamaker/RI-CLPM/blob/master/lavaan.Rmd

CLPM <- '

# Estimate the lagged effects between the observed variables.
Beh_09 + CT_09 ~ Beh_06 + CT_06
Beh_13 + CT_13 ~ Beh_09 + CT_09

# Estimate the covariance between the observed variables at the first wave. 
Beh_06 ~~ CT_06 # Covariance

# Estimate the covariances between the residuals of the observed variables.
Beh_09 ~~ CT_09
Beh_13 ~~ CT_13

# Estimate the (residual) variance of the observed variables.
Beh_06 ~~ Beh_06 # Variances
CT_06 ~~ CT_06 
Beh_09 ~~ Beh_09 # Residual variances
CT_09 ~~ CT_09 
Beh_13 ~~ Beh_13 
CT_13 ~~ CT_13

#Regression of observed variables on time independent predictors (sex, education, ethnicity)
Beh_06 + Beh_09 + Beh_13 ~ s1*sex + s2*education + s3*ethnicity
CT_06 + CT_09 + CT_13 ~ s4*sex + s5*education + s6*ethnicity

'
CLPM.fit <- lavaan(CLPM, data = df, missing = 'ML', meanstructure = T, int.ov.free = T) 
results_CLPM <- summary(CLPM.fit, standardized = T)

#Get values
#auto regressive coefficients
AR_Beh6_Beh9 <- results_CLPM[[1]][1,11]
AR_Beh9_Beh13 <- results_CLPM[[1]][5,11]
AR_CT6_CT9 <- results_CLPM[[1]][4,11]
AR_CT9_CT13 <- results_CLPM[[1]][8,11]
#cross lagged coefficients
CR_Beh6_CT9 <- results_CLPM[[1]][3,11]
CR_Beh9_CT13 <- results_CLPM[[1]][7,11]
CR_CT6_Beh9 <- results_CLPM[[1]][2,11]
CR_CT9_Beh13 <- results_CLPM[[1]][6,11]
#pvals (same order)
pvals_clpm <- c(results_CLPM[[1]][1,9],
                results_CLPM[[1]][5,9],
                results_CLPM[[1]][4,9],
                results_CLPM[[1]][8,9],
                results_CLPM[[1]][3,9],
                results_CLPM[[1]][7,9],
                results_CLPM[[1]][2,9],
                results_CLPM[[1]][6,9])

#create dataframe and write output
results_CLPM_combined <- data.frame()

results_CLPM_combined[c(1:8),1] <- c(AR_Beh6_Beh9, AR_Beh9_Beh13, AR_CT6_CT9, AR_CT9_CT13,
                                     CR_Beh6_CT9, CR_Beh9_CT13, CR_CT6_Beh9, CR_CT9_Beh13)
results_CLPM_combined[c(1:8),2] <- pvals_clpm
rownames(results_CLPM_combined) <- c("AR_Beh6_Beh9", "AR_Beh9_Beh13", "AR_CT6_CT9", "AR_CT9_CT13",
                                     "CR_Beh6_CT9", "CR_Beh9_CT13", "CR_CT6_Beh9", "CR_CT9_Beh13")
colnames(results_CLPM_combined) <- c("standardized_coefficient", "pval")

#Here all dataframes will be combined in one excel file, each sheet contains output from a different dataset if you run all
#Second time you run this, append = F needs to be switched to append = T
write.xlsx(results_CLPM_combined, "results_CLPM_combined.xlsx", sheetName = whichdf, append = F) 



#####################################
##### 2.2 RANDOM INTERCEPT CLPM #####
#####################################

#Model set-up performed like Ellen Hamaker
#https://github.com/ellenhamaker/RI-CLPM/blob/master/lavaan.Rmd

#Model with time invariant predictors that are modeled to have the same effect at all time points
RICLPM.covars <- '

#Random intercepts
RI_Beh =~ 1*Beh_06 + 1*Beh_09 + 1*Beh_13
RI_CT =~ 1*CT_06 + 1*CT_09 + 1*CT_13

#Create within-person centered variables
w_Beh_06 =~ 1*Beh_06
w_Beh_09 =~ 1*Beh_09
w_Beh_13 =~ 1*Beh_13
w_CT_06 =~ 1*CT_06
w_CT_09 =~ 1*CT_09
w_CT_13 =~ 1*CT_13

#Regression of observed variables on time independent predictors (sex, education, ethnicity)
Beh_06 + Beh_09 + Beh_13 ~ s1*sex + s2*education + s3*ethnicity
CT_06 + CT_09 + CT_13 ~ s4*sex + s5*education + s6*ethnicity

#Estimate the lagged effects between the within-person centered variables
w_Beh_09 + w_CT_09 ~ w_Beh_06 + w_CT_06
w_Beh_13 + w_CT_13 ~ w_Beh_09 + w_CT_09

#Estimate the covariance between the within-person centered variables at the first wave
w_Beh_06 ~~ w_CT_06

#Estimate the covariances between the residuals of the within-person centered variables (the innovations)
w_Beh_09 ~~ w_CT_09
w_Beh_13 ~~ w_CT_13

#Estimate the variance and covariance of the random intercepts
RI_Beh ~~ RI_Beh
RI_CT ~~ RI_CT
RI_Beh ~~ RI_CT

#Estimate the (residual) variance of the within-person centered variables
w_Beh_06 ~~ w_Beh_06
w_CT_06 ~~ w_CT_06
w_Beh_09 ~~ w_Beh_09
w_CT_09 ~~ w_CT_09
w_Beh_13 ~~ w_Beh_13
w_CT_13 ~~ w_CT_13
'

RICLPM.covars.fit <- lavaan(RICLPM.covars, data = df, missing = 'ML', meanstructure = T, int.ov.free = T) 
results_RICLPM <- summary(RICLPM.covars.fit, standardized = T)

#Get values
#auto regressive coefficients
AR_Beh6_Beh9 <- results_RICLPM[[1]][31,11]
AR_Beh9_Beh13 <- results_RICLPM[[1]][35,11]
AR_CT6_CT9 <- results_RICLPM[[1]][34,11]
AR_CT9_CT13 <- results_RICLPM[[1]][38,11]
#cross lagged coefficients
CR_Beh6_CT9 <- results_RICLPM[[1]][33,11]
CR_Beh9_CT13 <- results_RICLPM[[1]][37,11]
CR_CT6_Beh9 <- results_RICLPM[[1]][32,11]
CR_CT9_Beh13 <- results_RICLPM[[1]][36,11]
#pvals (same order)
pvals_RICLPM <- c(results_RICLPM[[1]][31,9],
                  results_RICLPM[[1]][35,9],
                  results_RICLPM[[1]][34,9],
                  results_RICLPM[[1]][38,9],
                  results_RICLPM[[1]][33,9],
                  results_RICLPM[[1]][37,9],
                  results_RICLPM[[1]][32,9],
                  results_RICLPM[[1]][36,9])

#create dataframe and write output
results_RICLPM_combined <- data.frame()

results_RICLPM_combined[c(1:8),1] <- c(AR_Beh6_Beh9, AR_Beh9_Beh13, AR_CT6_CT9, AR_CT9_CT13,
                                       CR_Beh6_CT9, CR_Beh9_CT13, CR_CT6_Beh9, CR_CT9_Beh13)
results_RICLPM_combined[c(1:8),2] <- pvals_RICLPM
rownames(results_RICLPM_combined) <- c("AR_Beh6_Beh9", "AR_Beh9_Beh13", "AR_CT6_CT9", "AR_CT9_CT13",
                                       "CR_Beh6_CT9", "CR_Beh9_CT13", "CR_CT6_Beh9", "CR_CT9_Beh13")
colnames(results_RICLPM_combined) <- c("standardized_coefficient", "pval")

#Here all dataframes will be combined in one excel file, each sheet contains output from a different dataset if you run all
#Second time you run this, append = F needs to be switched to append = T
write.xlsx(results_RICLPM_combined, "results_RICLPM_combined.xlsx", sheetName = whichdf, append = F) 
