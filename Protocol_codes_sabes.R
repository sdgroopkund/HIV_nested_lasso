## 1. Call relevant libraries into workspace

library(survival)
library(glmnet)
library(dplyr)
library(xlsx)

## 2. Set Working Directory

## This step (line) could be skipped if the data file is stored in the local working directory.
## Otherwise one has to manually edit setwd() prior to running these codes.
## Include within quotes path to the folder where the data file is stored.
## Each forward slash in the path command needs to be replaced with a backslash as shown here.

setwd('C:\Work\MSD data')

## 3. Read Data

## The data that we have provided in Mendeley (https://data.mendeley.com/datasets/8k97yxxghm) 
## constitute a substantial part of the case control samples that were selected from Sabes. 
## However, note that the process of case-control selection for our analysis (two controls 
## were matched on time under observation and one was matched on calendar time) were conducted
## in separately. The choice was very particular to our dataset and would not be useful
## for the protocol in hand. For someone to conduct matching to select case-control samples 
## from a larger cohort, please refer to our sample codes 'sample_codes_clogit.R' in Github.


lasso.dat <- read.csv('SabesCaseControlData.csv')
head(lasso.dat) ## View first few lines of the data

## 4. Gather all relevant variables

covmat.cyto <- lasso.dat[, 14:22] ## only cytokines (immune response)
covmat.cov <- lasso.dat[, 4:13] ## only sociodemographic/behavioral
covmat.all <- cbind(covmat.cyto, covmat.cov) ## all covariates

x_train.all <- model.matrix( ~ ., covmat.all)
x_train.cyto <- model.matrix( ~ ., covmat.cyto)
x_train.cov <- model.matrix( ~ ., covmat.cov)

y <- lasso.dat$outcome
groupid <- lasso.dat$grpid


## 5. Run the LASSO models

## COVARIATE LASSO MODEL

cv.lasso.cov <- cv.glmnet(x = x_train.cov, y = y, intercept = FALSE, family =  "binomial")
lasso.cov.fit <- glmnet(x_train.cov, y, lambda = cv.lasso.cov$lambda.min,
                        intercept = F, family="binomial")

## CYTOKINE LASSO MODEL

cv.lasso.cyto <- cv.glmnet(x = x_train.cyto, y = y, intercept = FALSE, family =  "binomial")
lasso.cyto.fit <- glmnet(x_train.cyto, y, lambda = cv.lasso.cyto$lambda.min,
                         intercept = F, family="binomial")

## FULL LASSO MODEL

cv.lasso.all <- cv.glmnet(x = x_train.all, y = y, intercept = FALSE, family = "binomial")
lasso.all.fit <- glmnet(x_train.all, y, lambda = cv.lasso.all$lambda.min,
                        intercept = F, family="binomial")

## 6. Collecting selected variables 
## The following blocks of code collects the selected variables from GLMNET (Lasso) and create new datasets
## to run conditional logistic regression

## COVARIATE LASSO MODEL

ind.cov = rownames(lasso.cov.fit$beta)[as.numeric(lasso.cov.fit$beta) != 0]
names.cov.lasso = NULL
for (i in 1:ncol(covmat.cov)){
  gg = grep(names(covmat.cov)[i], ind.cov)
  if (length(gg) > 0) {
    names.cov.lasso = c(names.cov.lasso, names(covmat.cov)[i])
  } 
}
covmat.cov.lasso = covmat.cov %>% select(all_of(names.cov.lasso)) %>% data.frame()

## CYTOKINE LASSO MODEL

ind.cyto = rownames(lasso.cyto.fit$beta)[as.numeric(lasso.cyto.fit$beta) != 0]
names.cyto.lasso = NULL
for (i in 1:ncol(covmat.cyto)){
  gg = grep(names(covmat.cyto)[i], ind.cyto)
  if (length(gg) > 0) {
    names.cyto.lasso = c(names.cyto.lasso, names(covmat.cyto)[i])
  } 
}
covmat.cyto.lasso = covmat.cyto %>% select(all_of(names.cyto.lasso)) %>% data.frame()

## FULL LASSO MODEL

ind.all = rownames(lasso.all.fit$beta)[as.numeric(lasso.all.fit$beta) != 0]
names.all.lasso = NULL
for (i in 1:ncol(covmat.all)){
  gg = grep(names(covmat.all)[i], ind.all)
  if (length(gg) > 0) {
    names.all.lasso = c(names.all.lasso, names(covmat.all)[i])
  } 
}
covmat.all.lasso = covmat.all %>% select(all_of(names.all.lasso)) %>% data.frame()


## 7. Run conditional logistic regression on the selected covariates

## COVARIATE LASSO MODEL

lasso.cov.dat <- data.frame(outcome = y, covmat.cov.lasso, grpid = groupid)
lcm.cov.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                      data = lasso.cov.dat, method='exact')

## CYTOKINE LASSO MODEL

lasso.cyto.dat <- data.frame(outcome = y, covmat.cyto.lasso, grpid = groupid)
lcm.cyto.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                       data = lasso.cyto.dat, method='exact')

## FULL LASSO MODEL

lasso.all.dat <- data.frame(outcome = y, covmat.all.lasso, grpid = groupid)
lcm.all.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                      data = lasso.all.dat, method='exact')


## 8. Gather results and write out to an excel sheet

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

lasso.models <- bind_rows(zz.cov, zz.cyto, zz.all)

write.xlsx(lasso.models, file = 'Results_LASSO.xlsx', sheetName = 'Sheet 1')
