rm(list = ls())

library(parallel)
library(MASS)
library(sl3)
library(data.table)
library(origami)
library(assertthat)
library(tidyverse)
library(SuperLearner)
set.seed(4883845)
setwd("~/Documents/K99R00/R00/mediation_monotonicity")
source("truthfun.R")
source("fit_mechanisms_updated.R")
source("cvfuns.R")


n<-10000
data <- simdata(n)
  
  A <- data[, "A"]
  M <- data[, substr(names(data), 1, 1) == "M"]
  Z <- data[, substr(names(data), 1, 1) == "Z"]
  Y <- data[, "Y"]
  W <- data[, substr(names(data), 1, 1) == "W"]

  obs_weights = rep(1, length(Y))

 data <- data.table::as.data.table(cbind(Y, M, Z, A, W, obs_weights))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  m_names <- paste("M", seq_len(dim(data.table::as.data.table(M))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", m_names, "Z", "A", w_names, "obs_weights"))

mean_lrnr<-Lrnr_mean$new()
fglm_lrnr<-Lrnr_glm_fast$new(family=binomial())
fglm_lrnr_uv<-Lrnr_glm_fast$new(family=gaussian())
mars_lrnr<-Lrnr_earth$new()
glmnet_lrnr<-Lrnr_glmnet$new(family=binomial())
xgboost_lrnr<-Lrnr_xgboost$new()

slbin<-Lrnr_sl$new(learners=list(fglm_lrnr,mean_lrnr, mars_lrnr, glmnet_lrnr), metalearner=Lrnr_nnls$new())
sluv<-Lrnr_sl$new(learners=list(fglm_lrnr_uv, mean_lrnr, mars_lrnr, glmnet_lrnr), metalearner=Lrnr_nnls$new())

slbin<-Lrnr_sl$new(learners=list(fglm_lrnr,mean_lrnr), metalearner=Lrnr_nnls$new())
sluv<-Lrnr_sl$new(learners=list(fglm_lrnr_uv, mean_lrnr), metalearner=Lrnr_nnls$new())

learners<-slbin


est10<-est_onestep(data=data, contrast=c(1,0), g_learners=learners, e_learners=learners,q_learners=learners, r_learners=learners,b_learners=learners, w_names=w_names, m_names=m_names, ext_weights = NULL,
                        cv_folds = 5)
est11<-est_onestep(data=data, contrast=c(1,1), g_learners=learners, e_learners=learners,q_learners=learners, r_learners=learners,b_learners=learners, w_names=w_names, m_names=m_names, ext_weights = NULL,
                        cv_folds = 5)
est00<-est_onestep(data=data, contrast=c(0,0), g_learners=learners, e_learners=learners,q_learners=learners, r_learners=learners,b_learners=learners, w_names=w_names, m_names=m_names, ext_weights = NULL,
                        cv_folds = 5)

 
eifnie<-est11$eif-est10$eif

eifnde<-est10$eif-est00$eif

nie<-est11$theta - est10$theta
nde<-est10$theta - est00$theta

nie + c(-1,1)*1.96*sqrt(var(eifnie)/length(eifnie))
nde + c(-1,1)*1.96*sqrt(var(eifnde)/length(eifnde))
