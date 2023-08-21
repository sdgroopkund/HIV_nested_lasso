library(dplyr) 
library(MatchIt) 
library(clogitL1)
library(survival)

# If any of the above packages are not installed in R, please use the following command to install them
# For example, to install the package clogitL1, one needs to run the following commented line
# install.packages('clogitL1')

# Set seed for analysis for reproducibility and declare functions required for the analysis.

set.seed(5132)
expit <- function(x) 1 / (1 + exp( -x))

# Generate trial data (Two-stage sampling - Stage 1).
# In stage 1, N participants are enrolled into the trial. Let's say we have a set of p predictors, and that predictors 1 and 2 are significantly associated with the outcome. 
# Note that we may not have these predictors measured for every participant in stage 1, however, we simulate these values for everyone in stage 1 to generate the outcome for this toy example, which will be observed for all participants in stage 1.

N = 500
p = 20
pred.set = data.frame(matrix(rnorm(N * p, 0, 1), ncol = p))
prob = expit(- 2 + 3 * pred.set[, 1] + 4 * pred.set[, 2])
outcome = rbinom(N, 1, prob)

# Get hold of covariates that we want to employ for matching, which may be different from the set of predictors that we are primarily interested in. 
# Here we simulate some imaginary covariates that we will use for matching. Note that the matching variables should be available for all participants in stage 1 (barring missingness). 

match.var = data.frame(W1 = rnorm(N, 0, 1), W2 = rnorm(N, 0, 1), W3 = rbinom(N, 1, 0.5))
data.stage1 = data.frame(outcome = outcome, match.var)

# Generate case-control data (Two-stage sampling - Stage 2).
# In this step, we create the case-control dataset. We will sample all cases from stage 1, and additionally, for each case, we will sample 3 matched controls from the trial cohort. So number of strata is the same as the number of cases in the trial for this example.
# The matchIt library performs pairing, subset selection, and subclassification with the aim of creating treatment and control groups balanced on included covariates. The option distance = 'glm' estimates distance using propensity scores (with logistic regression), method = 'nearest' implements nearest neighbor matching (on these propensity scores), and ratio denotes the number of controls to sample for each case. The option replace = F means that in this example, the matched controls are selected from the trial without replacement. By specifying replace = T we can select controls with replacement, in which case, the same control can be matched to more than one case.

m.out <- matchit(outcome ~ W1 + W2 + W3, data = data.stage1, method = "nearest", distance = "glm", ratio = 3, replace = F)
m.data <- match.data(m.out)

# Get phase 2 covariates (predictors).
# Typically in 2-phase sampling design, we assume that the set of p predictors are measured only in those who are in the case-control sample. 

mpred.set <- pred.set %>% slice(as.numeric(row.names(m.data)))
m.data <- cbind(m.data, mpred.set)

# For a real-data scenario, measurement for these covariates, available only for stage 2 participants, needs to be collected, read in R, and merged with the rest of the case control data (consisting of outcome and matching variables).
# These coding steps are to implement steps 17 and 18 of the protocol (see 'step by step method details' above).Run LASSO with conditional logistic regression model using clogitL1 package. 
# For the full model, we not only select the set of predictors, but also the matching variables.

y <- m.data$outcome
strata <- m.data$subclass
mCov = m.data %>% select(X1 : eval(expression(paste0('X', p))), W1, W2, W3)

# Standardize covariates

mCov.scaled <- scale(mCov, center=TRUE, scale=TRUE)
mCov.scaled <- data.frame(mCov.scaled)
X.train <- model.matrix( ~ ., mCov.scaled)

# Run cross validation and select variables picked up by LASSO for the value of ?? that produces the best CV performance. 

clObj = clogitL1(x = X.train, y = y, strata)
clcvObj = cv.clogitL1(clObj)
which.beta.pos = which(clcvObj$beta[clcvObj$lambda == clcvObj$minCV_lambda,] > 0)
names.train = colnames(X.train)[which.beta.pos]
mCov.lassoset <- mCov %>% select(which(names(mCov) %in% names.train))

# Run conditional logistic regression on the selected features using survival package in R. 

m.data.lassoset <- data.frame(outcome = y, mCov.lassoset, grpid = strata)
clogit.fit <- clogit(outcome  ~ . -grpid  + strata(grpid), 
                      data = m.data.lassoset, method = 'exact')
summary(clogit.fit)
