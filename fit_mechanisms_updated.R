#utils::globalVariables(c("..w_names", "A", "Z"))

#contrast will be c(aprime, astar)
fit_treat_mech <- function(train_data,
                           valid_data = NULL,
                           contrast,
                           learners,
                           m_names,
                           w_names,
                           type = c("g", "e")) {
  if (type == "g") {
    cov_names <- w_names
  } else if (type == "e") {
    cov_names <- c(m_names, "Z", w_names)
  }

  aprime <- contrast[1]
  astar  <- contrast[2]
  ## construct task for treatment mechanism fit
  treat_task <- sl3::sl3_Task$new(
    data = train_data,
    weights = "obs_weights",
    covariates = cov_names,
    outcome = "A",
    outcome_type = "binomial"
  )

  ## fit and predict treatment mechanism
  treat_fit <- learners$train(treat_task)
  treat_pred <- treat_fit$predict()

    train_data_z1 <-  train_data_z0 <- data.table::copy(train_data)
    train_data_z1$Z<-1
    train_data_z0$Z<-0

    treat_task_z1 <- sl3::sl3_Task$new(
      data = train_data_z1,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "A",
      outcome_type = "binomial"
    )
    treat_z1_pred <- treat_fit$predict(treat_task_z1)

    treat_task_z0 <- sl3::sl3_Task$new(
      data = train_data_z0,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "A",
      outcome_type = "binomial"
    )
    treat_z0_pred <- treat_fit$predict(treat_task_z0)

  ## use full data for prediction if no validation data provided
  if (is.null(valid_data)) {
    treat_pred_A_prime <- aprime* treat_pred +
      (1 - astar) * (1 - treat_pred)
    treat_pred_A_star <- astar * treat_pred +
      (1 - astar) * (1 - treat_pred)
    #e_star_Z_1
    treat_pred_A_prime_Z1 <- aprime* treat_z1_pred +
      (1 - astar) * (1 - treat_z1_pred)
    treat_pred_A_star_Z1 <- astar * treat_z1_pred +
      (1 - astar) * (1 - treat_z1_pred)
    treat_pred_A_prime_Z0 <- aprime* treat_z0_pred +
      (1 - astar) * (1 - treat_z0_pred)
    treat_pred_A_star_Z0 <- astar * treat_z0_pred +
      (1 - astar) * (1 - treat_z0_pred)
    

    ## bounding to numerical precision and for positivity considerations
    out_treat_mat <- cbind(
      treat_pred_A_prime,
      treat_pred_A_star,
      treat_pred_A_prime_Z1,
      treat_pred_A_prime_Z0,
      treat_pred_A_star_Z1,
      treat_pred_A_star_Z0
    )
    out_treat_est <- apply(out_treat_mat, 2, function(x) {
      x_precise <- bound_precision(x)
      x_bounded <- bound_propensity(x_precise)
      return(x_bounded)
    })
    out_treat_est <- data.table::as.data.table(out_treat_est)
    data.table::setnames(out_treat_est, c(
      "treat_pred_A_prime",
      "treat_pred_A_star",
      "treat_pred_A_prime_Z1",
      "treat_pred_A_prime_Z0",
      "treat_pred_A_star_Z1",
      "treat_pred_A_star_Z0"
    ))

    ## output
    out <- list(
      treat_est = out_treat_est,
      treat_fit = treat_fit
    )
  } else {
    out_treat_est <- lapply(
      list(train_data, valid_data),
      function(data) {
        ## create task to generate contrast-specific predictions
        treat_task <- sl3::sl3_Task$new(
          data = data,
          weights = "obs_weights",
          covariates = cov_names,
          outcome = "A",
          outcome_type = "binomial"
        )

        ## predictions for training data
        treat_pred <- treat_fit$predict(treat_task)

        treat_pred_A_prime <- aprime * treat_pred +
          (1 - aprime) * (1 - treat_pred)
        treat_pred_A_star <- astar * treat_pred +
          (1 - astar) * (1 - treat_pred)


        data_z1 <-  data_z0 <- data.table::copy(data)
        data_z1$Z<-1
        data_z0$Z<-0

        treat_task_z1 <- sl3::sl3_Task$new(
          data = data_z1,
          weights = "obs_weights",
          covariates = c("A", cov_names),
          outcome = "A",
          outcome_type = "binomial"
        )
        treat_z1_pred <- treat_fit$predict(treat_task_z1)

        treat_task_z0 <- sl3::sl3_Task$new(
          data = data_z0,
          weights = "obs_weights",
          covariates = c("A", cov_names),
          outcome = "A",
          outcome_type = "binomial"
        )
        treat_z0_pred <- treat_fit$predict(treat_task_z0)

        treat_pred_A_prime_Z1 <- aprime* treat_z1_pred + (1 - astar) * (1 - treat_z1_pred)
        treat_pred_A_star_Z1 <- astar * treat_z1_pred + (1 - astar) * (1 - treat_z1_pred)
        treat_pred_A_prime_Z0 <- aprime* treat_z0_pred + (1 - astar) * (1 - treat_z0_pred)
        treat_pred_A_star_Z0 <- astar * treat_z0_pred + (1 - astar) * (1 - treat_z0_pred)
     
        ## bounding to numerical precision and for positivity considerations
        out_treat_mat <- cbind(
          treat_pred_A_prime,
          treat_pred_A_star,
          treat_pred_A_prime_Z1,
          treat_pred_A_prime_Z0,
          treat_pred_A_star_Z1,
          treat_pred_A_star_Z0
        )

        out_treat_est <- apply(out_treat_mat, 2, function(x) {
          x_precise <- bound_precision(x)
          x_bounded <- bound_propensity(x_precise)
          return(x_bounded)
        })
        out_treat_est <- data.table::as.data.table(out_treat_est)
        data.table::setnames(out_treat_est, c(
          "treat_pred_A_prime",
          "treat_pred_A_star",
          "treat_pred_A_prime_Z1",
          "treat_pred_A_prime_Z0",
          "treat_pred_A_star_Z1",
          "treat_pred_A_star_Z0"
        ))
      }
    )

    ## output
    out <- list(
      treat_est_train = out_treat_est[[1]],
      treat_est_valid = out_treat_est[[2]],
      treat_fit = treat_fit
    )
  }
  return(out)
}

fit_out_mech <- function(train_data,
                         valid_data = NULL,
                         contrast,
                         learners,
                         m_names,
                         w_names) {
  aprime <- contrast[1]
  astar  <- contrast[2]

  ##  construct task for propensity score fit
  b_natural_task <- sl3::sl3_Task$new(
    data = train_data,
    weights = "obs_weights",
    covariates = c(m_names, "Z", "A", w_names),
    outcome = "Y"
  )

  ## fit and predict
  b_natural_fit <- learners$train(b_natural_task)
  b_natural_pred <- b_natural_fit$predict()

  ## use full data for counterfactual prediction if no validation data given
  if (is.null(valid_data)) {
    ## set intervention to first contrast a_prime := contrast[1]
    train_data_aprime_z1 <- train_data_aprime_z0 <-train_data_astar_z1 <- train_data_astar_z0 <- data.table::copy(train_data)
    train_data_aprime_z1$A <- train_data_aprime_z0$A <- aprime
    train_data_astar_z1$A <- train_data_astar_z0$A <- astar
    train_data_aprime_z1$Z <- 1
    train_data_aprime_z0$Z <- 0
    train_data_astar_z1$Z <- 1
    train_data_astar_z0$Z <- 0

    ## predictions on observed data (i.e., under observed treatment status)
    b_natural_pred <- b_natural_fit$predict()

    ## create task for post-intervention outcome regression
    b_prime_z1_task <- sl3::sl3_Task$new(
      data = train_data_aprime_z1 ,
      weights = "obs_weights",
      covariates = c(m_names, "Z", "A", w_names),
      outcome = "Y"
    )

    ## predict from trained model on counterfactual data
    b_pred_A_prime_Z_1 <- b_natural_fit$predict(b_prime_z1_task)

    b_prime_z0_task <- sl3::sl3_Task$new(
      data = train_data_aprime_z0 ,
      weights = "obs_weights",
      covariates = c(m_names, "Z", "A", w_names),
      outcome = "Y"
    )

    ## predict from trained model on counterfactual data
    b_pred_A_prime_Z_0 <- b_natural_fit$predict(b_prime_z0_task)

    train_data$pseudoyz1<-b_pred_A_prime_Z_1

    muint_fit_task_z1 <- sl3::sl3_Task$new(
    data = train_data,
    weights = "obs_weights",
    covariates = c("Z", "A", w_names),
    outcome = "pseudoyz1",
    outcome_type="continuous"
    )

    muint_fit_z1 <- learners$train(muint_fit_task_z1)

    train_data_astar_z1$pseudoyz1<-train_data_astar_z0$pseudoyz1<-b_pred_A_prime_Z_1

    muint_star_z11_task <- sl3::sl3_Task$new(
      data = train_data_astar_z1 ,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "pseudoyz1"
    )
    
    m_z11 <- muint_fit_z1$predict(muint_star_z11_task)
    
    muint_star_z10_task <- sl3::sl3_Task$new(
      data = train_data_astar_z0 ,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "pseudoyz1"
    )
    
    m_z10 <- muint_fit_z1$predict(muint_star_z10_task)

    train_data$pseudoyz0<-train_data_astar_z0$pseudoyz0<-b_pred_A_prime_Z_0

    muint_fit_task_z0 <- sl3::sl3_Task$new(
    data = train_data,
    weights = "obs_weights",
    covariates = c("Z", "A", w_names),
    outcome = "pseudoyz0",
    outcome_type="continuous"
    )

    muint_fit_z0 <- learners$train(muint_fit_task_z0)

    muint_star_z00_task <- sl3::sl3_Task$new(
      data = train_data_astar_z0 ,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "pseudoyz0"
    )
    
    m_z00 <- muint_fit_z0$predict(muint_star_z00_task)

    ## output
    out_b_est <- data.table::as.data.table(cbind(
      b_natural_pred,
      b_pred_A_prime_Z_1,
      b_pred_A_prime_Z_0,
      m_z11,
      m_z10,
      m_z00
    ))
    data.table::setnames(out_b_est, c(
      "b_pred_A_natural",
      "b_pred_A_prime_Z_1",
      "b_pred_A_prime_Z_0",
      "m_z11",
      "m_z10",
      "m_z00"
    ))

    ## output
    out <- list(
      b_est = out_b_est,
      b_fit = b_natural_fit
    )
  } else {
    ## copy both training and validation data, once for each contrast
    train_data_intervene <- data.table::copy(train_data)
    valid_data_intervene <- data.table::copy(valid_data)

    ## predictions on observed data (i.e., under observed treatment status)
    b_natural_pred_train <- b_natural_fit$predict()
    b_natural_task_valid <- sl3::sl3_Task$new(
      data = valid_data,
      weights = "obs_weights",
      covariates = c(m_names, "Z", "A", w_names),
      outcome = "Y"
    )
    b_natural_pred_valid <- b_natural_fit$predict(b_natural_task_valid)

    ## set intervention to first contrast a' := contrast[1]
    out_b_est <- lapply(
      list(train_data_intervene, valid_data_intervene),
      function(data_intervene) {
        ## set intervention to first contrast a' := contrast[1]
        data_intervene_aprime_z1 <- data_intervene_aprime_z0 <-data_intervene_astar_z1 <- data_intervene_astar_z0 <- data.table::copy(data_intervene)
        data_intervene_aprime_z1$A <- data_intervene_aprime_z0$A <- aprime
        data_intervene_astar_z1$A <- data_intervene_astar_z0$A <- astar
        data_intervene_aprime_z1$Z <- 1
        data_intervene_aprime_z0$Z <- 0
        data_intervene_astar_z1$Z <- 1
        data_intervene_astar_z0$Z <- 0

        b_intervened_prime_z1_task <- sl3::sl3_Task$new(
          data = data_intervene_aprime_z1,
          weights = "obs_weights",
          covariates = c(m_names, "Z", "A", w_names),
          outcome = "Y"
        )

        ## predict from trained model on counterfactual data
        b_pred_A_prime_Z_1 <-
          b_natural_fit$predict(b_intervened_prime_z1_task)

        b_intervened_prime_z0_task <- sl3::sl3_Task$new(
          data = data_intervene_aprime_z0,
          weights = "obs_weights",
          covariates = c(m_names, "Z", "A", w_names),
          outcome = "Y"
        )

        ## predict from trained model on counterfactual data
        b_pred_A_prime_Z_0 <-
          b_natural_fit$predict(b_intervened_prime_z0_task)
        
        data_intervene$pseudoyz1<-data_intervene_astar_z1$pseudoyz1<-data_intervene_astar_z0$pseudoyz1<-b_pred_A_prime_Z_1

    muint_fit_task_z1 <- sl3::sl3_Task$new(
    data = data_intervene,
    weights = "obs_weights",
    covariates = c("Z", "A", w_names),
    outcome = "pseudoyz1",
    outcome_type="continuous"
    )

    muint_fit_z1 <- learners$train(muint_fit_task_z1)

    muint_star_z11_task <- sl3::sl3_Task$new(
      data = data_intervene_astar_z1 ,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "pseudoyz1"
    )
    
    m_z11 <- muint_fit_z1$predict(muint_star_z11_task)
    
    muint_star_z10_task <- sl3::sl3_Task$new(
      data = data_intervene_astar_z0 ,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "pseudoyz1"
    )
    
    m_z10 <- muint_fit_z1$predict(muint_star_z10_task)


    data_intervene$pseudoyz0<-data_intervene_astar_z0$pseudoyz0<-b_pred_A_prime_Z_0

    muint_fit_task_z0 <- sl3::sl3_Task$new(
    data = data_intervene,
    weights = "obs_weights",
    covariates = c("Z", "A", w_names),
    outcome = "pseudoyz0",
    outcome_type="continuous"
    )

    muint_fit_z0 <- learners$train(muint_fit_task_z0)

    muint_star_z00_task <- sl3::sl3_Task$new(
      data = data_intervene_astar_z0 ,
      weights = "obs_weights",
      covariates = c("Z", "A", w_names),
      outcome = "pseudoyz0"
    )
    
    m_z00 <- muint_fit_z0$predict(muint_star_z00_task)

        ## output
        out_b_est <- data.table::as.data.table(cbind(
          b_pred_A_prime_Z_1,
          b_pred_A_prime_Z_0,
          m_z11,
          m_z10,
          m_z00
        ))
        return(out_b_est)
      }
    )

    ## add natural treatment estimates to post-intervention predictions
    out_b_est[[1]] <- cbind(b_natural_pred_train, out_b_est[[1]])
    out_b_est[[2]] <- cbind(b_natural_pred_valid, out_b_est[[2]])
    lapply(out_b_est, function(x) {
      data.table::setnames(x, c(
        "b_pred_A_natural",
        "b_pred_A_prime_Z_1",
        "b_pred_A_prime_Z_0",
        "m_z11",
        "m_z10",
        "m_z00"
      ))
    })

    ## output
    out <- list(
      b_est_train = out_b_est[[1]],
      b_est_valid = out_b_est[[2]],
      b_fit = b_natural_fit
    )
  }
  return(out)
}

fit_moc_mech <- function(train_data,
                         valid_data = NULL,
                         contrast,
                         learners,
                         m_names,
                         w_names,
                         type = c("q", "r")) {
  
  aprime <- contrast[1]
  astar  <- contrast[2]
  ## construct task for nuisance parameter fit
  if (type == "q") {
    cov_names <- w_names
  } else if (type == "r") {
    cov_names <- c(m_names, w_names)
  }

  moc_task <- sl3::sl3_Task$new(
    data = train_data,
    weights = "obs_weights",
    covariates = c("A", cov_names),
    outcome = "Z",
    outcome_type = "binomial"
  )

  ## fit model on observed data
  moc_fit <- learners$train(moc_task)

  ## use full data for counterfactual prediction if no validation data given
  if (is.null(valid_data)) {
    ## set intervention to first contrast a_prime := contrast[1]
    train_data_aprime <- train_data_astar <- data.table::copy(train_data)
    #train_data_aprime_z1 <- train_data_aprime_z0 <- train_data_astar_z1 <- train_data_astar_z1 <- data.table::copy(train_data)

    train_data_aprime$A <- aprime
    train_data_astar$A <- astar

    #train_data_aprime$A <- train_data_aprime_z1$A <- train_data_aprime_z0$A <- aprime
    #train_data_astar$A <- train_data_astar_z1$A <- train_data_astar_z0$A <- astar
    #train_data_aprime_z1$Z <- 1
    #train_data_aprime_z0$Z <- 0
    #train_data_astar_z1$Z <- 1
    #train_data_astar_z0$Z <- 0
    
    ## predictions on observed data (i.e., under observed treatment status)
    moc_pred_A_natural <- moc_fit$predict()

    ## create task for post-intervention outcome regression
    moc_prime_task <- sl3::sl3_Task$new(
      data = train_data_aprime,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "Z",
      outcome_type = "binomial"
    )

    ## predict from trained model on counterfactual data
    moc_pred_A_prime <- moc_fit$predict(moc_prime_task)

    moc_star_task <- sl3::sl3_Task$new(
      data = train_data_astar,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "Z",
      outcome_type = "binomial"
    )

    ## predict from trained model on counterfactual data
    moc_pred_A_star <- moc_fit$predict(moc_star_task)

    ## output
    out_moc_est <- data.table::as.data.table(cbind(
      moc_pred_A_natural,
      moc_pred_A_prime,
      moc_pred_A_star
    ))
    data.table::setnames(out_moc_est, c(
      "moc_pred_A_natural",
      "moc_pred_A_prime",
      "moc_pred_A_star"
    ))

    ## output
    out <- list(
      moc_est = out_moc_est,
      moc_fit = moc_fit
    )
  } else {
    ## copy both training and validation data, once for each contrast
    train_data_intervene <- data.table::copy(train_data)
    valid_data_intervene <- data.table::copy(valid_data)

    ## predictions on observed data (i.e., under observed treatment status)
    moc_pred_A_natural_train <- moc_fit$predict()

    ## create task for post-intervention outcome regression
    moc_task_valid <- sl3::sl3_Task$new(
      data = valid_data,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "Z",
      outcome_type = "binomial"
    )

    ## prediction on observed data, in validation set
    moc_pred_A_natural_valid <- moc_fit$predict(moc_task_valid)

    ## set intervention to first contrast a_prime := contrast[1]
    out_moc_est <- lapply(
      list(train_data_intervene, valid_data_intervene),
      function(data_intervene) {
        ## intervene to set treatment to first contrast (A prime)
        data_intervene_aprime <- data_intervene_astar <- data.table::copy(data_intervene)
  
    data_intervene_aprime$A <- aprime
    data_intervene_astar$A <- astar
    
    ## create task for post-intervention outcome regression
    moc_prime_task <- sl3::sl3_Task$new(
      data = data_intervene_aprime,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "Z",
      outcome_type = "binomial"
    )

    ## predict from trained model on counterfactual data
    moc_pred_A_prime <- moc_fit$predict(moc_prime_task)

    moc_star_task <- sl3::sl3_Task$new(
      data = data_intervene_astar,
      weights = "obs_weights",
      covariates = c("A", cov_names),
      outcome = "Z",
      outcome_type = "binomial"
    )

    ## predict from trained model on counterfactual data
    moc_pred_A_star <- moc_fit$predict(moc_star_task)

        ## output
        out_moc_est <-
          data.table::as.data.table(cbind(
            moc_pred_A_prime,
            moc_pred_A_star
          ))
      }
    )

    ## add natural treatment estimates to post-intervention predictions
    out_moc_est[[1]] <- cbind(moc_pred_A_natural_train, out_moc_est[[1]])
    out_moc_est[[2]] <- cbind(moc_pred_A_natural_valid, out_moc_est[[2]])
    lapply(out_moc_est, function(x) {
      data.table::setnames(x, c(
        "moc_pred_A_natural",
        "moc_pred_A_prime",
        "moc_pred_A_star"
      ))
    })

    ## output
    out <- list(
      moc_est_train_Z_one = out_moc_est[[1]],
      moc_est_valid_Z_one = out_moc_est[[2]],
      moc_est_train_Z_natural = out_moc_est[[1]] * train_data$Z +
        (1 - out_moc_est[[1]]) * (1 - train_data$Z),
      moc_est_valid_Z_natural = out_moc_est[[2]] * valid_data$Z +
        (1 - out_moc_est[[2]]) * (1 - valid_data$Z),
      moc_fit = moc_fit
    )
  }
  return(out)
}