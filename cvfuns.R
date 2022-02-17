cv_eif <- function(fold,data_in,contrast,g_learners,e_learners,q_learners,r_learners,b_learners,w_names,m_names) {
 # make training and validation data
 train_data <- origami::training(data_in)
 valid_data <- origami::validation(data_in)

#for debugging 
#train_data<-data
#valid_data<-data
#contrast<-c(1,0)
aprime <- contrast[1]
astar  <- contrast[2]

 # 1) fit regression for propensity score regression
 g_out <- fit_treat_mech(
  train_data = train_data,
  valid_data = valid_data,
  contrast = contrast,
  learners = g_learners,
  w_names = w_names,
  m_names = m_names,
  type = "g"
 )

g_prime <- g_out$treat_est_valid$treat_pred_A_prime
g_star <- g_out$treat_est_valid$treat_pred_A_star

 # 2) fit clever regression for treatment, conditional on mediators and z
 e_out <- fit_treat_mech(
  train_data = train_data,
  valid_data = valid_data,
  contrast = contrast,
  learners = e_learners,
  w_names = w_names,
  m_names = m_names,
  type = "e"
 )

e_prime_Z_1 <- e_out$treat_est_valid$treat_pred_A_prime_Z1
e_prime_Z_0 <- e_out$treat_est_valid$treat_pred_A_prime_Z0
e_star_Z_1 <- e_out$treat_est_valid$treat_pred_A_star_Z1
e_star_Z_0 <- e_out$treat_est_valid$treat_pred_A_star_Z0

 # 3) fit outcome regression
 mu_out <- fit_out_mech(
  train_data = train_data,
  valid_data = valid_data,
  contrast = contrast,
  learners = b_learners,
  m_names = m_names,
  w_names = w_names
 )

mu_pred <-mu_out$b_est_valid$b_pred_A_natural
mu_pred_aprime_z1 <- mu_out$b_est_valid$b_pred_A_prime_Z_1
mu_pred_aprime_z0 <- mu_out$b_est_valid$b_pred_A_prime_Z_0
m_z11<-mu_out$b_est_valid$m_z11
m_z10<-mu_out$b_est_valid$m_z10
m_z00<-mu_out$b_est_valid$m_z00

 # 4) fit mediator-outcome confounder regression, excluding mediator(s)
 q_out <- fit_moc_mech(
  train_data = train_data,
  valid_data = valid_data,
  contrast = contrast,
  learners = q_learners,
  m_names = m_names,
  w_names = w_names,
  type = "q"
 )
q_pred <- q_out$moc_est_valid_Z_one$moc_pred_A_natural
q_star <- q_out$moc_est_valid_Z_one$moc_pred_A_star
q_prime <- q_out$moc_est_valid_Z_one$moc_pred_A_prime

 r_out <- fit_moc_mech(
  train_data = train_data,
  valid_data = valid_data,
  contrast = contrast,
  learners = r_learners,
  m_names = m_names,
  w_names = w_names,
  type = "r"
 )

r_star <- r_out$moc_est_valid_Z_one$moc_pred_A_star
r_prime <- r_out$moc_est_valid_Z_one$moc_pred_A_prime

 eif_y_11 <- as.numeric(valid_data$A == aprime & valid_data$Z==1) / g_star * e_star_Z_1 / e_prime_Z_1 * (valid_data$Y - mu_pred)

 eif_y_10 <- as.numeric(valid_data$A == aprime & valid_data$Z==1) / (g_star * (1-q_star)) * e_star_Z_1 / e_prime_Z_1 * (1-r_star) / r_star * (q_prime - q_star) * (valid_data$Y - mu_pred)

 eif_y_00 <- as.numeric(valid_data$A == aprime & valid_data$Z==0) / g_prime * e_star_Z_0 / e_prime_Z_0 * (1-q_prime) / (1-q_star) * (valid_data$Y - mu_pred)

  eif_m_11 <- as.numeric(valid_data$A == astar & valid_data$Z == 1) / g_star * (mu_pred_aprime_z1 - m_z11)
 
  eif_m_10 <- as.numeric(valid_data$A == astar & valid_data$Z==0) / (g_star * (1-q_star)) * (q_prime - q_star) * (mu_pred_aprime_z1 - m_z10)

  eif_m_00 <- as.numeric(valid_data$A == astar & valid_data$Z==0) / (g_star * (1-q_star)) * (1-q_prime) * (mu_pred_aprime_z0 - m_z00)

  eif_z_11 <- as.numeric(valid_data$A == astar) / g_star * m_z11 * (valid_data$Z - q_pred)

  eif_z_10 <- (as.numeric(valid_data$A == aprime) / g_prime - as.numeric(valid_data$A == astar) / g_star) * m_z10 * (valid_data$Z - q_pred)

  eif_z_00 <- -1* as.numeric(valid_data$A == aprime) / g_prime * m_z00 * (valid_data$Z - q_pred)

  eif_w_11 <- m_z11 * q_star

  eif_w_10 <- m_z10 * (q_prime - q_star)

  eif_w_00 <- m_z00 * (1-q_prime)

   # un-centered efficient influence function
 eif_11 <- eif_y_11 + eif_m_11 + eif_z_11 + eif_w_11
 eif_10 <- eif_y_10 + eif_m_10 + eif_z_10 + eif_w_10
 eif_00 <- eif_y_00 + eif_m_00 + eif_z_00 + eif_w_00

 eif <- eif_11 + eif_10 + eif_00

 # output list
 out <- list(data.table::data.table(
  # components necessary for fluctuation step of TMLE
  g_prime = g_prime, g_star = g_star,

  # efficient influence function and fold IDs
  D_star = eif, fold = origami::fold_index()
 ))
 return(out)
}

est_onestep <- function(data,
            contrast,
            g_learners,
            e_learners,
            q_learners,
            r_learners,
            b_learners,
            w_names,
            m_names,
            ext_weights = NULL,
            cv_folds = 5) {

 # make sure that more than one fold is specified
 assertthat::assert_that(cv_folds > 1)

 # create cross-validation folds
 folds <- origami::make_folds(data,
  fold_fun = origami::folds_vfold,
  V = cv_folds
 )

 # estimate the EIF on a per-fold basis
 cv_eif_results <- origami::cross_validate(
  cv_fun = cv_eif,
  folds = folds,
  data_in = data,
  contrast = contrast,
  g_learners = g_learners,
  e_learners = e_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  b_learners = b_learners,
  w_names = w_names,
  m_names = m_names,
  use_future = FALSE,
  .combine = FALSE
 )

 # get estimated efficient influence function
 cv_eif_est <- do.call(c, lapply(cv_eif_results[[1]], `[[`, "D_star"))
 obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
 cv_eif_est <- cv_eif_est[order(obs_valid_idx)]

 # compute one-step estimate and variance from efficient influence function

  os_est <- mean(cv_eif_est)
  eif_est_out <- cv_eif_est

 os_var <- stats::var(eif_est_out) / length(eif_est_out)

 # output
 os_est_out <- list(
  theta = os_est,
  var = os_var,
  eif = eif_est_out,
  type = "onestep"
 )
 return(os_est_out)
}