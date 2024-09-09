                          #-----------------------------------------------------------------------------#
                          # Author: Eugene Baraka                                                       #
                          # Purpose: Flexible Cox modelling with Time-dependent & non-linear effects    #
                          # Date: September 7, 2024                                                     #
                          #-----------------------------------------------------------------------------#

# Load libraries
rm(list = ls()) # clean environment
library(data.table)
library(tidyverse)
library(table1)
library(GGally)
library(survival)
# library(rms)
setwd("~/Desktop/Github/survival-analysis/")

## ------------------------------------------------------------------------------------------------------------------------------
# Load source library for fitting flexible cox PH models
# This was developed by https://github.com/mebeauchamp
source("https://raw.githubusercontent.com/mebeauchamp/CoxFlex/main/CoxFlex%20-%2020220330%20-%20to%20share.R")
# load("./CoxFlex/dat.RData")


## ------------------------------------------------------------------------------------------------------------------------------
# Load data
# Data obtained from SurvSet, a python open-source time-to-event data repository (check get_data.ipynb for download steps)
# Check the SurvSet repo here: https://github.com/ErikinBC/SurvSet
prostate <- read.csv("data/prostatecancer.csv") %>% setDT() # convert dataframe to data.table
setorder(prostate, pid, time) # order by patient id and time to event/censoring 


## ------------------------------------------------------------------------------------------------------------------------------
# Data cleaning
## rename columns (num_ implies a numerical feature, while fac_ implies a categorical feature)
names(prostate) <- gsub(pattern = "num_", replacement = "", x = names(prostate)) 
names(prostate) <- gsub(pattern = "fac_", replacement = "", x = names(prostate))


#! Info: the dataset has four treatment categories: the three above + placebo
prostate$event <- as.factor(prostate$event)

prostate[, c("ap", "sdate", "ekg", "sg", "pf", "dbp") := NULL] # delete unnecessary columns

# rename other columns for easy identification
setnames(prostate, c("time", "wt", "rx", "hg", "sz", "hx", "bm"), 
         c("fwup", "weight", "treatment", "serum.hb", 
           "tumor.size", "cvd.history", "bone.metastasis"))
# prostate$event <- factor(prostate$event, levels = c(0, 1), labels = c("Did not die", "Died"))


## ------------------------------------------------------------------------------------------------------------------------------
# data processing
is.numeric(prostate$pid) 
is.data.frame(prostate)
# drop rows with missing values
# complex data imputation methods could be used here
prostate <- na.omit(prostate)

## rename treatment categories 
prostate[treatment == "0.2 mg estrogen", treatment := "estrogen0.2"]
prostate[treatment == "1.0 mg estrogen", treatment := "estrogen1.0"]
prostate[treatment == "5.0 mg estrogen", treatment := "estrogen5.0"]
prostate$treatment <- factor(prostate$treatment, levels = c("placebo", "estrogen0.2", "estrogen1.0", "estrogen5.0"))

# create treatment dummy variables to have 3 levels (i.e., n - 1)
trt.dummies <- model.matrix(~ treatment, data = prostate)
colnames(trt.dummies) <- sub('treatment', "", colnames(trt.dummies))
prostate_new <- cbind(prostate[, treatment := NULL], trt.dummies[, 2:4])

## cvd history
prostate_new[, cvd.history := ifelse(cvd.history == "N", 0, 1)]
## bone metastasis
prostate_new[, bone.metastasis := ifelse(bone.metastasis == "N", 0, 1)]
## cancer stage
prostate_new[, stage := ifelse(stage == 3, 0, 1)] # stage 3 == 0, stage 4 == 1

# event (death)
prostate_new$event <- as.numeric(prostate_new$event)
prostate_new[, event := ifelse(event == 1, 0, 1)]


## ------------------------------------------------------------------------------------------------------------------------------
modelling_dt <- prostate_new[, pid := pid + 1]
modelling_dt <- modelling_dt[, fwup := fwup + 1]
modelling_dt <- rename(modelling_dt, time = fwup)
modelling_dt <- rename(modelling_dt, id = pid)


## ------------------------------------------------------------------------------------------------------------------------------
# Model 1: Simple CoxPH model (COXPH Univariate Model)
# only treatment was included in the survival model, with no adjustment for the covariates;
model1 <- coxph(Surv(time, event) ~ estrogen0.2 + estrogen1.0 + estrogen5.0, data = modelling_dt, id = id)
summary(model1)


## ------------------------------------------------------------------------------------------------------------------------------
# Model 2: COXPH Multivariable Model
# modelling proportional hazard and linear effects of all the covariates
model2 <- coxph(Surv(time, event) ~ estrogen0.2 + estrogen1.0 + estrogen5.0 + age + weight + sbp + serum.hb + tumor.size + stage + cvd.history + bone.metastasis,  data = modelling_dt, id = id)
summary(model2)


## ------------------------------------------------------------------------------------------------------------------------------
# Model 3: Full flexible multivariable model
# consists of modelling the time dependent effects of exposure (treatment),
# time dependent effects of all covariates and 
# non-linear effects of continuous covariates
model3 <- CoxFlex(data = modelling_dt, Type = c("time", "event"), 
                  variables = c("estrogen0.2", "estrogen1.0", "estrogen5.0", 
                                "age", "weight", "sbp", "serum.hb", "tumor.size", 
                                "stage", "cvd.history", "bone.metastasis"), 
                  TD = c(1, 1, 1,
                         1, 1, 1, 1, 1,
                         1, 1, 1), 
                  NL = c(0, 0, 0,
                         1, 1, 1, 1, 1, 
                         0, 0, 0), 
                  m=1, p=2, knots=-999)


## ------------------------------------------------------------------------------------------------------------------------------
# Model 4: Flexible multivariable model with backward selection
# models TD/NL of all the covariates with stepwise elimination of non-significant effects (p > 0.05)
model4 <- backward_selection2(data=modelling_dt,
        Type = c("time", "event"),
        variables = c("estrogen0.2", "estrogen1.0", "estrogen5.0", "age", "weight", 
                      "sbp", "serum.hb", "tumor.size", "stage","cvd.history", "bone.metastasis"),
        continuous= c(0, 0, 0, 1, 1, 
                      1, 1, 1, 0, 0, 0),
        TD = c(0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0),
        NL = c(0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0),
        m=1, p=2, alpha_back=0.05, knots=-999)
# variables = c("estrogen0.2", "estrogen1.0", "estrogen5.0", 
#                                 "age", "weight", "sbp", "serum.hb", "tumor.size", 
#                                 "stage", "cvd.history", "bone.metastasis")


## ------------------------------------------------------------------------------------------------------------------------------
# Create table for model 1 summary for report
tabl <- data.frame(matrix(ncol = 2, nrow = 12, 
                           dimnames = list(c("Estrogen0.2mg", "Estrogen1.0mg", "Estrogen5.0mg", "Age", 
                                             "Weight", "SBP", "Serum", "Tumor size", "Stage", "CVD history", "Bone metastasis", "AIC"), 
                                           c("Model1", "Model2"))))
mod1_summary <- paste(round(exp(coef(model1)), 3), " [", round(exp(confint(model1))[1:3], 3), ", ", round(exp(confint(model1))[4:6], 3), "]", sep = "")
mod1_summary <- c(mod1_summary, "-", "-", "-", "-", "-", "-", "-", "-", round(AIC(model1), 2))
  

mod2_summary <- paste(round(exp(coef(model2)), 3), " [", round(exp(confint(model2))[1:11], 3), ", ", round(exp(confint(model2))[12:22], 3), "]", sep = "")
mod2_summary <- c(mod2_summary, round(AIC(model2), 2))

tabl$Model1 <- mod1_summary
tabl$Model2 <- mod2_summary

library(xtable)
latex_tabl2 <- xtable(tabl, caption = "Comparison between Model 1 and Model 2")
print(latex_tabl2, include.rownames = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
# Create table 2 for report
tab2 <- data.frame(matrix(ncol = 2, nrow = 11, 
                           dimnames = list(c("Estrogen0.1mg", "Estrogen1.0mg", "Estrogen5.0mg", "Age", 
                                             "Weight", "SBP", "Serum", "Tumor size", "Stage", "CVD history", "Bone metastasis"), 
                                           c("Model3", "Model4"))))
mod3_summs <- c("TD", "TD", "TD", "TD/NL", "TD/NL", "TD/NL", "TD/NL", "TD/NL", "TD", "TD", "TD")
mod4_sums <- c("-", "TD", "TD", "PH/LL", "NL", "-", "NL", "PH/LL", "-", "PH/LL", "PH/LL")

tab2$Model3 <- mod3_summs
tab2$Model4 <- mod4_sums

latex_tabl3 <- xtable(tab2, caption = "Comparison of forced TD/NL effects in Model 3 with selected effects in backward elimination in Model 4")
print(latex_tabl3, include.rownames = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
# Create table 3 for report
rows1 <- c("TD-Estrogen0.1mg", "TD-Estrogen1.0mg", "TD-Estrogen5.0mg", "NL-Age", "TD-Age", 
  "NL-Weight", "TD-Weight", "NL-SBP", "TD-SBP", "NL-Serum", 
  "TD-Serum", "NL-TumorSize", "TD-TumorSize", "TD-Stage4", "TD-CVDHistory", 
  "TD-BoneMetastasis")  
tabl_3 <- data.frame(matrix(ncol = 1, nrow = 16, 
                         dimnames = list(rows1, c("p-value"))))

tabl_3$p.value <- model3$pvalue
latex_tabl_3 <- xtable(tabl_3, caption = "P-values of variables in Model 3")
print(latex_tabl_3, include.rownames = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
# Figures for flexible model 3
par(mfrow=c(2,2))
plot.FlexSurv(model.FlexSurv = model3, variable="estrogen0.2", TD=1, NL=0, col="red", xlab="Time", ylab="log(HR)",
main="TD effect of 0.2mg Estrogen", type="l")
plot.FlexSurv(model.FlexSurv = model3, variable="estrogen1.0", TD=1, NL=0, col="red", xlab="Time", ylab="log(HR)",
main="TD effect of 1.0mg Estrogen", type="l")
plot.FlexSurv(model.FlexSurv = model3, variable="estrogen5.0", TD=1, NL=0, col="red", xlab="Time", ylab="log(HR)",
main="TD effect of 5.0mg Estrogen", type="l")
plot.FlexSurv(model.FlexSurv = model3, variable="weight", TD=0, NL=1, col="blue", ref.value.NL = min(modelling_dt$weight), xlab="Weight Index", ylab="log(HR)",
main="NL effect of Weight", type="l")


## ------------------------------------------------------------------------------------------------------------------------------
# Figures for flexible model 4
rows2 <- c("PH/LL-Age", "NL-Weight", "NL-Serum","PH/LL-TumorSize", "PH/LL-estrogen1.0", "TD-estrogen5.0", "PH/LL-CVDHistory", 
"PH/LL-BoneMetastasis")  
tabl_4 <- data.frame(matrix(ncol = 2, nrow = 8, 
                         dimnames = list(rows2, c("Estimate", "p-value"))))
tabl_4$Estimate <- round(exp(model4$final_model$coefficients),3)
tabl_4$p.value <- model4$final_model$pvalue

latex_tabl_4 <- xtable(tabl_4, caption = "Model 4 Results with alpha = 0.05 p-values")
print(latex_tabl_4, include.rownames = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))
plot.FlexSurv(model.FlexSurv = model4$final_model, variable="estrogen5.0", TD=1, NL=0,
         col="blue", xlab="Time", ylab="log(HR)",
         main="TD effect of 5.0mg Estrogen", type="l")

plot.FlexSurv(model.FlexSurv = model4$final_model, variable="serum.hb", TD=0, NL=1, ref.value.NL = min(modelling_dt$serum.hb),
         col="red", xlab="Serum", ylab="log(HR)",
         main="NL effect of Serum Hemoglobin", type="l")

plot.FlexSurv(model.FlexSurv = model4$final_model, variable="weight", TD=0, NL=1, ref.value.NL = min(modelling_dt$weight),
         col="red", xlab="Weight", ylab="log(HR)",
         main="NL effect of Weight", type="l")


## ------------------------------------------------------------------------------------------------------------------------------
# Models comparison using AIC
AIC1 <- AIC(model1)
AIC2 <- AIC(model2)
AIC3 <-  -2 * model3$Partial_Log_Likelihood + 2 * model3$Number_of_parameters
AIC4 <-  -2 * model4$final_model$Partial_Log_Likelihood + 2 * model4$final_model$Number_of_parameters

tabl_AIC <- data.frame(matrix(ncol = 1, nrow = 4, 
                         dimnames = list(c("Model 1", "Model 2", "Model 3", "Model 4"), c("AIC"))))
tabl_AIC$AIC <- c(AIC1, AIC2, AIC3, AIC4)

latex_tabl_AIC <- xtable(tabl_AIC, caption = "Comparison of the 4 Models")
print(latex_tabl_AIC, include.rownames = TRUE)




## ------------------------------------------------------------------------------------------------------------------------------
# Create table 1 for report
label(prostate$age) <- "Age"
units(prostate$age) <- "years"
label(prostate$weight) <- "Weight Index"
label(prostate$fwup) <- "Follow-up"
units(prostate$fwup) <- "months"
label(prostate$treatment) <- "Treatment"
label(prostate$sbp) <- "SBP"
# label(prostate$dbp) <- "DBP"
label(prostate$serum.hb) <- "Serum Hemoglobin"
label(prostate$tumor.size) <- "Tumor Size"
label(prostate$cvd.history) <- "Has CVD History"
label(prostate$performance.rating) <- "Performance Rating"
label(prostate$bone.metastasis) <- "Bone Metastatis"
table1(~ age + weight + fwup + treatment + sbp | event, data = prostate, caption = "Distribution of events across chosen predictors")
table1(~ serum.hb + tumor.size + performance.rating + cvd.history + bone.metastasis | event, data = prostate, caption = "Distribution of events across chosen predictors..cnt'd")


## ------------------------------------------------------------------------------------------------------------------------------

## SENSITIVITY ANALYSES-----------------------------------------------------
# exploratory data analysis of variables
# ggpairs(prostate[, c("event", "fwup", "age", "weight", "sbp", "dbp", "serum.hb", "tumor.size", "stage")])


## ------------------------------------------------------------------------------------------------------------------------------
fit1 <- coxph(Surv(fwup, event) ~ treatment, data = prostate, id = pid)
summary(fit1)


## ------------------------------------------------------------------------------------------------------------------------------
prostate$fwup_centered <- prostate$fwup - median(prostate$fwup)
fit2 <- coxph(Surv(fwup, event) ~ treatment*fwup_centered, data = prostate, id = pid)
summary(fit2)


## ------------------------------------------------------------------------------------------------------------------------------
fit3 <- coxph(Surv(fwup, event) ~ treatment + age  + weight + sbp + serum.hb + tumor.size + stage + performance.rating + cvd.history + bone.metastasis, data = prostate, id = pid)
summary(fit3)


## ------------------------------------------------------------------------------------------------------------------------------
fit3 <- coxph(Surv(fwup, event) ~ treatment + rcs(age, 3)  + weight + sbp + serum.hb + tumor.size + stage + performance.rating + cvd.history + bone.metastasis, data = prostate, id = pid)
summary(fit3)

