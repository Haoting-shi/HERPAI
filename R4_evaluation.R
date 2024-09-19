########################################################
############## SET Pathway and data input ##############

############ 1. C index ############
library(Hmisc)
# C index estimation function
c_index <- function(c_index_result){       
  
  c_indices <- c_index_result$concordance
  SEs <- c_index_result$se
  n <- as.numeric(count(c_index_result))
  
  # combined C index
  combined_C <- mean(c_indices)
  Vw_C = sum(SEs^2) / n
  Vb_C = sum( ((c_indices - mean(c_indices))^2) )/(n-1)
  combined_C_se = sqrt(Vw_C + Vb_C + Vb_C/n)

  
  CI_lower <- combined_C - 1.96 * combined_C_se
  CI_upper <- min((combined_C + 1.96 * combined_C_se), 1)
  
  c_index_summary <- data.frame(
    mean_c_index = round(combined_C, 3), 
    CI = paste0(round(CI_lower, 3), " to ", round(CI_upper, 3))
  )
  
  return(c_index_summary)
}

# C index output function
c_index_merge <- function(risk_score_var, outcome_time, outcome_status, 
                          train, validation, internal.test, external.test, 
                          subgroup_var = NULL, subgroup_condition = NULL) {

  # subgroup definition
  filter_data <- function(data, subgroup_var, subgroup_condition) {
    if (!is.null(subgroup_var) && !is.null(subgroup_condition)) {
      return(subset(data, eval(parse(text=subgroup_condition))))
    }
    return(data)
  }
  
  # c index calculation
  calculate_for_subgroup <- function(subgroup_condition) {
    
    train_cindex <- numeric()
    val_cindex <- numeric()
    internal.test_cindex <- numeric() 
    external.test_cindex <- numeric() 
    
    for (i in 1:20) {
      for (fold in 1:5) {
        
        # train
        train_fold <- filter_data(train[[i]][[fold]], subgroup_var, subgroup_condition)
        fit_train <- coxph(Surv(train_fold[[outcome_time]] , train_fold[[outcome_status]]) ~ train_fold[[risk_score_var]], method="efron")
        train_cindex <- rbind(train_cindex, data.frame(concordance = concordance(fit_train)$concordance, se = sqrt(concordance(fit_train)$var)))
        
        # validation
        valid_fold <- filter_data(validation[[i]][[fold]], subgroup_var, subgroup_condition)
        fit_valid <- coxph(Surv(valid_fold[[outcome_time]] , valid_fold[[outcome_status]]) ~ valid_fold[[risk_score_var]], method="efron")
        val_cindex <- rbind(val_cindex, data.frame(concordance = concordance(fit_valid)$concordance, se = sqrt(concordance(fit_valid)$var)))
      }
      
      # internal test
      internal.test_data <- filter_data(internal.test[[i]], subgroup_var, subgroup_condition)
      fit_internal.test <- coxph(Surv(internal.test_data[[outcome_time]] , internal.test_data[[outcome_status]]) ~ internal.test_data[[risk_score_var]], method="efron")
      internal.test_cindex <- rbind(internal.test_cindex, data.frame(concordance = concordance(fit_internal.test)$concordance, se = sqrt(concordance(fit_internal.test)$var)))
      
      # external test
      external.test_data <- filter_data(external.test[[i]], subgroup_var, subgroup_condition)
      fit_external.test <- coxph(Surv(external.test_data[[outcome_time]] , external.test_data[[outcome_status]]) ~ external.test_data[[risk_score_var]], method="efron")
      external.test_cindex <- rbind(external.test_cindex, data.frame(concordance = concordance(fit_external.test)$concordance, se = sqrt(concordance(fit_external.test)$var)))
      
    }
    
    train_summary <- c_index(train_cindex)
    val_summary <- c_index(val_cindex)
    internal.test_summary <- c_index(internal.test_cindex)
    external.test_summary <- c_index(external.test_cindex)
    
    result <- data.frame(
      train_cindex = train_summary$mean_c_index, train_95CI = train_summary$CI,
      val_cindex = val_summary$mean_c_index, val_95CI = val_summary$CI,
      internal.test_cindex = internal.test_summary$mean_c_index, internal.test_95CI = internal.test_summary$CI, 
      external.test_cindex = external.test_summary$mean_c_index, external.test_95CI = external.test_summary$CI
    )
    
    return(result)
  }
  
  if (is.null(subgroup_var) || is.null(subgroup_condition)) {

    results_all <- calculate_for_subgroup(NULL)
    results_all$subgroup_var <- "all comers"
    results_all$subgroup_value <- "all comers"
    results_all$model <- risk_score_var
    
    return(results_all)
    
  } else {
    
    results_subgroup <- calculate_for_subgroup(subgroup_condition)
    results_subgroup$subgroup_var <- subgroup_var
    results_subgroup$subgroup_value <- subgroup_condition
    results_subgroup$model <- risk_score_var   

    return(results_subgroup)
  }
}


############ 2. AUC ############
# AUC estimation function
AUC <- function(AUC_result){       
  
  AUCs <- AUC_result$AUC
  SEs <-  AUC_result$se
  n <- as.numeric(count(AUC_result))
  
  combined_AUC <- mean(AUCs)
  Vw_C = sum(SEs^2) / n
  Vb_C = sum( ((AUCs - mean(AUCs))^2) )/(n-1)
  combined_AUC_se = sqrt(Vw_C + Vb_C + Vb_C/n)
  
  CI_lower <- combined_AUC - 1.96 * combined_AUC_se
  CI_upper <- min((combined_AUC + 1.96 * combined_AUC_se), 1)
  
  AUC_summary <- data.frame(
    mean_AUC = combined_AUC,
    
    AUC = paste0(round(combined_AUC, 3), " (",
                 round(CI_lower, 3), " to ", 
                 round(CI_upper, 3), ")"))
  
  return(AUC_summary)
}

# AUC output function
AUC_merge <- function(risk_score_var, outcome_time, outcome_status, 
                      validation, internal.test, external.test, 
                      subgroup_var = NULL, subgroup_condition = NULL) {
  
  filter_data <- function(data, subgroup_var, subgroup_condition) {
    if (!is.null(subgroup_var) && !is.null(subgroup_condition)) {
      return(subset(data, eval(parse(text=subgroup_condition))))
    }
    return(data)
  }
  
  
  
  calculate_for_subgroup <- function(subgroup_condition) {
    
    library(timeROC)
    
    val_auc1 <- numeric()
    val_auc3 <- numeric()
    val_auc5 <- numeric()
    
    internal.test_auc1 <- numeric()
    internal.test_auc3 <- numeric()
    internal.test_auc5 <- numeric()
    
    external.test_auc1 <- numeric()
    external.test_auc3 <- numeric() 
    external.test_auc5 <- numeric() 
    
    for (i in 1:20) {
      for (fold in 1:5) {
        
        # validation
        valid_fold <- filter_data(validation[[i]][[fold]], subgroup_var, subgroup_condition)
        valid.roc <- timeROC(T = valid_fold[[outcome_time]], delta = valid_fold[[outcome_status]], marker = valid_fold[[risk_score_var]],
                             cause = 1, times = c(1,3,5), iid = TRUE)
        
        val_auc1 <- rbind(val_auc1, data.frame(AUC = valid.roc$AUC["t=1"], se = valid.roc$inference$vect_sd_1["t=1"]))
        val_auc3 <- rbind(val_auc3, data.frame(AUC = valid.roc$AUC["t=3"], se = valid.roc$inference$vect_sd_1["t=3"]))
        val_auc5 <- rbind(val_auc5, data.frame(AUC = valid.roc$AUC["t=5"], se = valid.roc$inference$vect_sd_1["t=5"]))
      }
      
      # internal test
      internal.test_data <- filter_data(internal.test[[i]], subgroup_var, subgroup_condition)
      internal.test.roc <- timeROC(T = internal.test_data[[outcome_time]], delta = internal.test_data[[outcome_status]], marker = internal.test_data[[risk_score_var]],
                                   cause = 1, times = c(1,3,5), iid = TRUE)
      
      internal.test_auc1 <- rbind(internal.test_auc1, data.frame(AUC = internal.test.roc$AUC["t=1"], se = internal.test.roc$inference$vect_sd_1["t=1"]))
      internal.test_auc3 <- rbind(internal.test_auc3, data.frame(AUC = internal.test.roc$AUC["t=3"], se = internal.test.roc$inference$vect_sd_1["t=3"]))
      internal.test_auc5 <- rbind(internal.test_auc5, data.frame(AUC = internal.test.roc$AUC["t=5"], se = internal.test.roc$inference$vect_sd_1["t=5"]))
      
      # external test
      external.test_data <- filter_data(external.test[[i]], subgroup_var, subgroup_condition)
      external.test.roc <- timeROC(T = external.test_data[[outcome_time]], delta = external.test_data[[outcome_status]], marker = external.test_data[[risk_score_var]],
                                   cause = 1, times = c(1,3,5), iid = TRUE)
      
      external.test_auc1 <- rbind(external.test_auc1, data.frame(AUC = external.test.roc$AUC["t=1"], se = external.test.roc$inference$vect_sd_1["t=1"]))
      external.test_auc3 <- rbind(external.test_auc3, data.frame(AUC = external.test.roc$AUC["t=3"], se = external.test.roc$inference$vect_sd_1["t=3"]))
      external.test_auc5 <- rbind(external.test_auc5, data.frame(AUC = external.test.roc$AUC["t=5"], se = external.test.roc$inference$vect_sd_1["t=5"]))
      
    }
    
    val_auc1_summary <- AUC(val_auc1)
    val_auc3_summary <- AUC(val_auc3)
    val_auc5_summary <- AUC(val_auc5)
    
    internal.test_auc1_summary <- AUC(internal.test_auc1)
    internal.test_auc3_summary <- AUC(internal.test_auc3)
    internal.test_auc5_summary <- AUC(internal.test_auc5)
    
    external.test_auc1_summary <- AUC(external.test_auc1)
    external.test_auc3_summary <- AUC(external.test_auc3)
    external.test_auc5_summary <- AUC(external.test_auc5)
    
    
    result <- data.frame(
      
      val_auc1 = val_auc1_summary$AUC,
      val_auc3 = val_auc3_summary$AUC, 
      val_auc5 = val_auc5_summary$AUC, 
      
      internal.test_auc1 = internal.test_auc1_summary$AUC, 
      internal.test_auc3 = internal.test_auc3_summary$AUC, 
      internal.test_auc5 = internal.test_auc5_summary$AUC, 
      
      external.test_auc1 = external.test_auc1_summary$AUC, 
      external.test_auc3 = external.test_auc3_summary$AUC,
      external.test_auc5 = external.test_auc5_summary$AUC 
    )
    
    return(result)
  }
  
  if (is.null(subgroup_var) || is.null(subgroup_condition)) {
    
    results_all <- calculate_for_subgroup(NULL)
    results_all$subgroup_var <- "all comers"
    results_all$subgroup_value <- "all comers"
    results_all$model <- risk_score_var
    
    return(results_all)
    
  } else {
    
    results_subgroup <- calculate_for_subgroup(subgroup_condition)
    results_subgroup$subgroup_var <- subgroup_var
    results_subgroup$subgroup_value <- subgroup_condition
    results_subgroup$model <- risk_score_var   
    
    return(results_subgroup)
  }
}




