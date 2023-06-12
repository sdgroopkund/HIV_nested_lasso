library(xlsx)
library(stringr)
library(survival)
library(glmnet)
library(dplyr)
library(data.table)
library(ggplot2)
setwd('C:/Users/sdasgup2/Desktop/Work/Rachel/MSD_Data')

## DATA READING, PRE-PROCESSING STEPS

msd.dat <- read.csv('MSD.csv')
dates.dat <- read.csv('Dates.csv')
demo.dat <- read.csv('Demo.csv')
demofull.dat <- merge(demo.dat, dates.dat, by = intersect(names(demo.dat), names(dates.dat)))
msddemo.dat <- merge(demofull.dat, msd.dat, by='ptid')

# Transformed data to log10
combined10 <- msddemo.dat
combined10[, 34:44] <- log(combined10[,34:44], 10)

# Outcome
combined10$outcome = 1*(combined10$Case_or_Control=='Case')

combined.diff.abs.dat <- data.frame(summarize(group_by(combined10, ptid),
                                              TNF_a_enr = ifelse(sum(Time_Point == 'ENR') == 1, TNF_a[Time_Point == 'ENR'], NA),
                                              MIP_1a_enr = ifelse(sum(Time_Point == 'ENR') == 1, MIP_1a[Time_Point == 'ENR'], NA),
                                              IP_10_enr = ifelse(sum(Time_Point == 'ENR') == 1, IP_10[Time_Point == 'ENR'], NA),
                                              IL_7_enr = ifelse(sum(Time_Point == 'ENR') == 1, IL_7[Time_Point == 'ENR'], NA),
                                              IL_6_enr = ifelse(sum(Time_Point == 'ENR') == 1, IL_6[Time_Point == 'ENR'], NA),
                                              IL_2_enr = ifelse(sum(Time_Point == 'ENR') == 1, IL_2[Time_Point == 'ENR'], NA),
                                              IL_12p70_enr = ifelse(sum(Time_Point == 'ENR') == 1, IL_12p70[Time_Point == 'ENR'], NA),
                                              IFN_y_enr = ifelse(sum(Time_Point == 'ENR') == 1, IFN_y[Time_Point == 'ENR'], NA),
                                              IL_10_enr = ifelse(sum(Time_Point == 'ENR') == 1, IL_10[Time_Point == 'ENR'], NA),
                                              
                                              TNF_a_x1 =  ifelse(sum(Time_Point == 'X-1') == 1, TNF_a[Time_Point == 'X-1'], NA),
                                              MIP_1a_x1 = ifelse(sum(Time_Point == 'X-1') == 1, MIP_1a[Time_Point == 'X-1'], NA),
                                              IP_10_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IP_10[Time_Point == 'X-1'], NA),
                                              IL_7_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IL_7[Time_Point == 'X-1'], NA),
                                              IL_6_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IL_6[Time_Point == 'X-1'], NA),
                                              IL_2_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IL_2[Time_Point == 'X-1'], NA),
                                              IL_12p70_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IL_12p70[Time_Point == 'X-1'], NA),
                                              IFN_y_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IFN_y[Time_Point == 'X-1'], NA),
                                              IL_10_x1 = ifelse(sum(Time_Point == 'X-1') == 1, IL_10[Time_Point == 'X-1'], NA),
                                              
                                              TNF_a_change = ifelse(n()>1, TNF_a[Time_Point=='X-1'] - TNF_a[Time_Point=='ENR'], NA),
                                              MIP_1a_change = ifelse(n()>1, MIP_1a[Time_Point=='X-1'] - MIP_1a[Time_Point=='ENR'], NA),
                                              IP_10_change = ifelse(n()>1, IP_10[Time_Point=='X-1'] - IP_10[Time_Point=='ENR'], NA),
                                              IL_7_change = ifelse(n()>1, IL_7[Time_Point=='X-1'] - IL_7[Time_Point=='ENR'], NA),
                                              IL_6_change = ifelse(n()>1, IL_6[Time_Point=='X-1'] - IL_6[Time_Point=='ENR'], NA),
                                              IL_2_change = ifelse(n()>1, IL_2[Time_Point=='X-1'] - IL_2[Time_Point=='ENR'], NA),
                                              IL_12p70_change = ifelse(n()>1, IL_12p70[Time_Point=='X-1'] - IL_12p70[Time_Point=='ENR'], NA),
                                              IFN_y_change = ifelse(n()>1, IFN_y[Time_Point=='X-1'] - IFN_y[Time_Point=='ENR'], NA),
                                              IL_10_change = ifelse(n()>1, IL_10[Time_Point=='X-1'] - IL_10[Time_Point=='ENR'], NA),
                                              
                                              PatType = if (sum(Time_Point == 'ENR') == 1) PatType[Time_Point == 'ENR'] else PatType[Time_Point == 'X-1'],
                                              enrolldt = if(sum(Time_Point == 'ENR') == 1) enrolldt[Time_Point == 'ENR'] else enrolldt[Time_Point == 'X-1'],
                                              priordt = if(sum(Time_Point == 'ENR') == 1) priordt[Time_Point == 'ENR'] else  priordt[Time_Point == 'X-1'],
                                              postdt = if(sum(Time_Point == 'ENR') == 1) postdt[Time_Point == 'ENR'] else postdt[Time_Point == 'X-1'],
                                              Gender_Identity_Step1 = if(sum(Time_Point == 'ENR') == 1) Gender_Identity_Step1[Time_Point == 'ENR'] else Gender_Identity_Step1[Time_Point == 'X-1'],
                                              grpid = if(sum(Time_Point == 'ENR') == 1) grpid[Time_Point == 'ENR'] else grpid[Time_Point == 'X-1'],
                                              x1.after.EDDI = if(sum(Time_Point == 'ENR') == 1) x1.after.EDDI[Time_Point == 'ENR'] else x1.after.EDDI[Time_Point == 'X-1'], 
                                              x1.after.EPDDI = if(sum(Time_Point == 'ENR') == 1) x1.after.EPDDI[Time_Point == 'ENR'] else x1.after.EPDDI[Time_Point == 'X-1']
))

basedat <- read.xlsx('Step2_Cases_and_Controls_CLAI_and_Alcohol_080620_V2.xlsx', sheetIndex = 2)
clai <- read.xlsx('Step2_Cases_and_Controls_CLAI_and_Alcohol_080620_V2.xlsx', sheetIndex = 1)
combined.diff.abs.base <- merge(basedat, combined.diff.abs.dat, by ='ptid')
combined.diff.lassocov <- merge(combined.diff.abs.base, clai, by.x = c('ptid', 'priordt', 'type'), 
                              by.y = c('ptid', 'VisitDate', 'type'), all.x=T)


## CREATE VARIABLES FOR LASSO

## gather all relevant variables
XX <- combined.diff.lassocov[, c(21:29, 34, 37, 40:44, 47, 54:55, 57:58)]
XX <- na.omit(XX)

covmat.all <- XX[, - c(10, 11)] ## all covariates
covmat.cyto <- XX[, 1:9] ## only cytokines (immune response)
covmat.cov <- XX[, 12:21] ## only sociodemographic/behavioral

x_train.all <- model.matrix( ~ ., covmat.all)
x_train.cyto <- model.matrix( ~ ., covmat.cyto)
x_train.cov <- model.matrix( ~ ., covmat.cov)

y <- XX$outcome
groupid <- XX$grpid


## FULL LASSO MODEL

cv.lasso.allcov <- cv.glmnet(x = x_train, y = y, intercept = FALSE, family = "binomial")
lasso.allcov.fit <- glmnet(x_train, y, lambda = cv.lasso.allcov$lambda.min,
                           intercept = F, family="binomial")

## The following variables were selected by LASSO, and now we subset the covariate data to select only these variables
covmat.allcov.lasso <- covmat.allcov[, c("age", "Post_Secondary_Education", "Sex_Work",
                                         "AUD", "CLAI.enr", "CLAI_count_Except_Mains", "role",
                                         "IP_10_change", "IL_7_change", "IL_12p70_change", "IL_10_change",
                                         "IL_2_change")]

lasso.allcov.dat <- data.frame(outcome = y, covmat.allcov.lasso, grpid = groupid)
lcm.allcov.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                         data = lasso.allcov.dat, method='exact')


## COVARIATE LASSO MODEL

cv.lasso.cov <- cv.glmnet(x = x_train.cov, y = y, intercept = FALSE, family =  "binomial")
lasso.cov.fit <- glmnet(x_train.cov, y, lambda = cv.lasso.cov$lambda.min,
                        intercept = F, family="binomial")

## The following variables were selected by LASSO, and now we subset the covariate data to select only these variables
covmat.cov.lasso <- covmat.cov[, c("age", "Post_Secondary_Education", "Sex_Work", "Transgender", "AUD", 
                                   "Alcohol_YN", "CLAI.enr", "CLAI_count_Except_Mains", "role")]


lasso.cov.dat <- data.frame(outcome = y, covmat.cov.lasso, grpid = groupid)
lcm.cov.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                         data = lasso.cov.dat, method='exact')


## CYTOKINE LASSO MODEL

cv.lasso.cyto <- cv.glmnet(x = x_train.cyto, y = y, intercept = FALSE, family =  "binomial")
lasso.cyto.fit <- glmnet(x_train.cyto, y, lambda = cv.lasso.cyto$lambda.min,
                        intercept = F, family="binomial")

## The following variables were selected by LASSO, and now we subset the covariate data to select only these variables
covmat.cyto.lasso <- covmat.cov[, c("TNF_a_change", "MIP_1a_change", "IL_12p70_change", "IP_10_change")]


lasso.cyto.dat <- data.frame(outcome = y, covmat.cyto.lasso, grpid = groupid)
lcm.cyto.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                      data = lasso.cyto.dat, method='exact')


## GATHERING RESULTS

cc.all = exp(confint(lcm.all.fit))
cc.cov = exp(confint(lcm.cov.fit))
cc.cyto = exp(confint(lcm.cyto.fit))

zz.all = data.frame(round(cbind(summary(lcm.all.fit)$coefficients[, 1:3], 
                                summary(lcm.all.fit)$coefficients[, 5], cc.all),3))
names(zz.all) <- c('Estimate', 'exp(Estimate)', 'se(Estimate)', 'Pvalue', 'Lower.CI', 'Upper.CI')

zz.cov = data.frame(round(cbind(summary(lcm.cov.fit)$coefficients[, 1:3], 
                                summary(lcm.cov.fit)$coefficients[, 5], cc.cov),3))
names(zz.cov) <- c('Estimate', 'exp(Estimate)', 'se(Estimate)', 'Pvalue', 'Lower.CI', 'Upper.CI')

zz.cyto = data.frame(round(cbind(summary(lcm.cyto.fit)$coefficients[, 1:3], 
                                summary(lcm.cyto.fit)$coefficients[, 5], cc.cyto),3))
names(zz.cyto) <- c('Estimate', 'exp(Estimate)', 'se(Estimate)', 'Pvalue', 'Lower.CI', 'Upper.CI')

zz.cov$Pctdec.Est <- NA
zz.cyto$Pctdec.Est <- NA
zz.all$Pctdec.Est <- NA

list1 <- intersect(rownames(zz.cov), rownames(zz.all))
zz.all[list1, 7] <- (zz.cov[list1, 1] - zz.all[list1, 1])/abs(zz.cov[list1, 1])*100

zz.all$model <- 'LASSO_marker.change_confounders'
zz.cov$model <- 'LASSO_confounders'
zz.cyto$model <- 'LASSO_marker.change'

lasso.models <- rbind(zz.cov, zz.cyto, zz.all)

write.xlsx(lasso.models, file = 'Results_LASSO.xlsx', sheetName = 'Sheet 1', row.names = F)

