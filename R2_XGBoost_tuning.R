########################################################
############## SET Pathway and data input ##############


######## Bayesian Optimization ################################## 
## global set 
set.seed(123)

optimize_function <- function(learning_rate, max_depth, subsample, 
                              colsample_bytree, colsample_bylevel, 
                              alpha, gamma, lambda, nrounds) {
  
  current_train_cindex <- numeric()
  current_valid_cindex <- numeric()
  
  for (i in 1:20) {
    for (fold in 1:5) {
      
      set.seed(123)

      train_data <- train[[i]][[fold]]
      train.x <- as.matrix(as.data.frame(lapply(train_data[,c("age","bmi","mense","FH","breast_surgery","LN_surgery",
                                                              "histology_group","grade_group","ER","PR","HER2_IHC","Ki67","T_stage","N_stage","Ki67_fp1","Ki67_fp2")], as.numeric)))
      
      train.y <- ifelse(train_data$iDFS_status == 1, train_data$iDFS, -train_data$iDFS)
      dtrain <- xgb.DMatrix(data = train.x, label = train.y)
      
      valid_data <- validation[[i]][[fold]]
      valid_x <- as.matrix(as.data.frame(lapply(valid_data[,c("age","bmi","mense","FH","breast_surgery","LN_surgery",
                                                              "histology_group","grade_group","ER","PR","HER2_IHC","Ki67","T_stage","N_stage","Ki67_fp1","Ki67_fp2")], as.numeric)))
      valid_y <- ifelse(valid_data$iDFS_status == 1, valid_data$iDFS, -valid_data$iDFS)
      dvalid <- xgb.DMatrix(data = valid_x, label = valid_y)

      param <- list(objective = 'survival:cox',
                    booster = "gbtree",
                    eval_metric = 'cox-nloglik',
                    learning_rate = learning_rate,
                    max_depth = as.integer(max_depth),
                    subsample = subsample,
                    colsample_bytree = colsample_bytree,
                    colsample_bylevel = colsample_bylevel,
                    alpha = alpha,
                    gamma = gamma,
                    lambda = lambda)
      
      # Determine best iteration round
      watchlist <- list(train = dtrain, eval = dvalid)
      model <- xgb.train(params = param, data = dtrain, nrounds = as.integer(nrounds), 
                         watchlist = watchlist, early_stopping_rounds = 10, verbose = 0)
      best_nrounds <- model$best_iteration
      

      model <- xgb.train(params = param, data = dtrain, nrounds = best_nrounds)

      train_risk_score <- predict(model, newdata = dtrain)
      valid_risk_score <- predict(model, newdata = dvalid)
      
      train_data$xgb.risk_score <- train_risk_score
      valid_data$xgb.risk_score <- valid_risk_score
      
      train[[i]][[fold]] <- train_data
      validation[[i]][[fold]] <- valid_data
    }
  }
}


## define tuning bounds 
bounds <- list(learning_rate = c(0.0001, 0.1),
               max_depth = c(1L, 6L),
               subsample = c(0.1, 0.5),
               colsample_bytree = c(0.1, 0.8),
               colsample_bylevel = c(0.1, 0.8),
               alpha = c(0, 20),
               gamma = c(0, 20),
               lambda = c(0, 20),
               nrounds = c(1L, 500L))

opt_results <- bayesOpt(FUN = optimize_function, bounds = bounds, 
                        initPoints = 20, iters.n = 50, acq = "ucb", kappa = 2.576)
xgb.param.results = opt_results$scoreSummary

best_params <- getBestPars(opt_results)
best_params
best_params <- list(learning_rate = 0.01657858,
                    max_depth = 5,
                    subsample = 0.4044451,
                    colsample_bytree = 0.5425658,
                    colsample_bylevel = 0.5066444,
                    alpha = 2.660112,
                    gamma = 1.993577,
                    lambda = 15.41508,
                    nrounds = 410)

