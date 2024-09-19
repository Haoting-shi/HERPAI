########################################################
############## SET Pathway and data input ##############

######## 1. Fractional polynomial terms identification ################
## Identify in complete data, for age at diagnosis, BMI at diagnosis, ER, PR, and Ki67 
library(mfp)
fit_fp <- mfp(Surv(iDFS, iDFS_status) ~ fp(age) + fp(bmi) + fp(ER) + fp(PR) + fp(Ki67),
              family = cox, 
              data = complete_dt, select=0.05, verbose=TRUE)
summary(fit_fp)
rm(complete_dt, fit_fp)

# Ki67 FP1 = ((Ki67 + 0.5)/10)^-0.5
# Ki67 FP2 = ((Ki67 + 0.5)/10)^-0.5 * log(((Ki67 + 0.5)/10)) 
# Rest keep linear

######## 2. Select variables in training cohort ###############

cox_model_select <- list()
for (i in seq_along(train)) {
  
  cox_model_select[[i]] <- list()
  
  for (fold in seq_along(train[[i]])) {
    train_set = train[[i]][[fold]]
    
    train_set$ER = train_set$ER/10
    train_set$PR = train_set$PR/10
    train_set$Ki67 = train_set$Ki67/10
    
    cox_model <- coxph(Surv(iDFS, iDFS_status) ~ age + bmi + factor(mense) + factor(FH) + factor(breast_surgery) + factor(LN_surgery) + 
                         factor(histology_group) + as.numeric(grade_group) + ER + PR + factor(HER2_IHC) + Ki67 + as.numeric(T_stage) + Ki67_fp1 + Ki67_fp2 + 
                         as.numeric(N_stage) + 
                         age * factor(FH),  ## Cox Full Model for Variable Selection
                       data = train_set)
    cox_model_select[[i]][[fold]] <- summary(cox_model)
  }
}

# Summarize HR using Rubin's Rules 
coefficients_select <- data.frame()
variables <- c("age", "bmi", "factor(mense)1", "factor(FH)1", "factor(breast_surgery)1","factor(LN_surgery)1",
               "factor(histology_group)1","factor(histology_group)2", 
               "as.numeric(grade_group)", "ER", "PR", "factor(HER2_IHC)1", "Ki67", "Ki67_fp1", "Ki67_fp2",
               "as.numeric(T_stage)", "as.numeric(N_stage)",
               "age:factor(FH)1") 

for (var in variables) {
  coef_estimates <- c()
  coef_se <- c()
  coef_p_values <- c()
  
  for (i in seq_along(cox_model_select)) {
    for (fold in seq_along(cox_model_select[[i]])) {
      coef_estimates <- c(coef_estimates, cox_model_select[[i]][[fold]]$coefficients[var, "coef"])
      coef_se <- c(coef_se, cox_model_select[[i]][[fold]]$coefficients[var, "se(coef)"])
      coef_p_values <- c(coef_p_values, cox_model_select[[i]][[fold]]$coefficients[var, "Pr(>|z|)"])
    }
  }
  
  mean_estimate <- mean(coef_estimates)
  Vw <- sum(coef_se^2)/100
  Vb <- sum((coef_estimates - mean_estimate)^2)/99
  total_SE <- sqrt(Vw + Vb + Vb/100)
  
  coefficients_select <- rbind(coefficients_select, data.frame(
    Variable = var,
    HR = exp(mean_estimate),
    Lower_95CI = exp(mean_estimate - 1.96 * total_SE),
    Upper_95CI = exp(mean_estimate + 1.96 * total_SE),
    P_Value = round(2*(1 - pnorm(abs(mean_estimate/total_SE))), 4)
  ))
}

coefficients_select



################## 3. Refit in training cohort to generate final Cox model #####################

cox_model_final <- list()
for (i in seq_along(train)) {
  
  cox_model_final[[i]] <- list()
  
  for (fold in seq_along(train[[i]])) {

    train_set = train[[i]][[fold]]

    train_set$PR = train_set$PR/10
    
    cox_model <- coxph(Surv(iDFS, iDFS_status) ~ age + as.numeric(grade_group) + PR + Ki67_fp1 + Ki67_fp2 + as.numeric(T_stage) + 
                         as.numeric(N_stage),  ## Cox Full Model for Variable Selection
                       data = train_set)
    cox_model_final[[i]][[fold]] <- summary(cox_model)
  }
}

# Summarize HR using Rubin's Rules 
coefficients_final <- data.frame()
variables2 <- c("age", "as.numeric(grade_group)", "PR", "Ki67_fp1", "Ki67_fp2",
                "as.numeric(T_stage)", "as.numeric(N_stage)") 

for (var in variables2) {
  coef_estimates <- c()
  coef_se <- c()
  coef_p_values <- c()
  
  for (i in seq_along(cox_model_final)) {
    for (fold in seq_along(cox_model_final[[i]])) {
      coef_estimates <- c(coef_estimates, cox_model_final[[i]][[fold]]$coefficients[var, "coef"])
      coef_se <- c(coef_se, cox_model_final[[i]][[fold]]$coefficients[var, "se(coef)"])
      coef_p_values <- c(coef_p_values, cox_model_final[[i]][[fold]]$coefficients[var, "Pr(>|z|)"])
    }
  }
  
  mean_estimate <- mean(coef_estimates)
  Vw <- sum(coef_se^2)/100
  Vb <- sum((coef_estimates - mean_estimate)^2)/99
  total_SE <- sqrt(Vw + Vb + Vb/100)
  
  coefficients_final <- rbind(coefficients_final, data.frame(
    Variable = var,
    HR = exp(mean_estimate),
    Lower_95CI = exp(mean_estimate - 1.96 * (total_SE)),
    Upper_95CI = exp(mean_estimate + 1.96 * (total_SE)),
    P_Value = round(2*(1 - pnorm(abs(mean_estimate/total_SE))), 4)
  ))
}

coefficients_final


################## 4. Refit in validation set and internal and external set ########################

# validation set
for (i in seq_along(validation)) {
  
  for (fold in seq_along(validation[[i]])) {
    validation_set <- validation[[i]][[fold]] 
    validation_set$PR = validation_set$PR/10
    
    validation_set$grade_group <- as.numeric(validation_set$grade_group) 
    validation_set$N_stage <- as.numeric(validation_set$N_stage)  
    validation_set$T_stage <- as.numeric(validation_set$T_stage)   
    
    # Estimate the risk score, refit
    validation_set$risk_score.cox <- exp(log(1.3836280) * validation_set$grade_group + 
                                         log(0.2693276) * validation_set$Ki67_fp1 + 
                                         log(0.6514327) * validation_set$Ki67_fp2 +
                                         log(0.9426962) * validation_set$PR +
                                         log(1.2568085) * validation_set$T_stage +  
                                         log(0.9948714) * validation_set$age + 
                                         log(1.0052645) * validation_set$N_stage)
  }
}

# internal test set
for (i in seq_along(internal.test)) {
  
    internal.test.set <- internal.test[[i]] 
    internal.test.set$PR = internal.test.set$PR/10
    
    internal.test.set$grade_group <- as.numeric(internal.test.set$grade_group)  
    internal.test.set$N_stage <- as.numeric(internal.test.set$N_stage)  
    internal.test.set$T_stage <- as.numeric(internal.test.set$T_stage)  
    
    internal.test.set$risk_score.cox <- exp(log(1.3836280) * internal.test.set$grade_group + 
                                            log(0.2693276) * internal.test.set$Ki67_fp1 + 
                                            log(0.6514327) * internal.test.set$Ki67_fp2 +
                                            log(0.9426962) * internal.test.set$PR +
                                            log(1.2568085) * internal.test.set$T_stage +  
                                            log(0.9948714) * internal.test.set$age + 
                                            log(1.0052645) * internal.test.set$N_stage)
  }


for (i in seq_along(external.test)) {
  
  external.test.set <- external.test[[i]] 
  external.test.set$PR = external.test.set$PR/10
  
  external.test.set$grade_group <- as.numeric(external.test.set$grade_group)  
  external.test.set$N_stage <- as.numeric(external.test.set$N_stage)  
  external.test.set$T_stage <- as.numeric(external.test.set$T_stage)  
  
  external.test.set$risk_score.cox <- exp(log(1.3836280) * external.test.set$grade_group + 
                                          log(0.2693276) * external.test.set$Ki67_fp1 + 
                                          log(0.6514327) * external.test.set$Ki67_fp2 +
                                          log(0.9426962) * external.test.set$PR +
                                          log(1.2568085) * external.test.set$T_stage +  
                                          log(0.9948714) * external.test.set$age + 
                                          log(1.0052645) * external.test.set$N_stage)
}


################## 5. Justification for using continuous variable (stage and grade) ################## 
# Comparing AIC in training set

cox_model_train_AIC <- list()
new_cox_model_train_AIC <- list()

for (i in seq_along(train)) {
  cox_model_train_AIC[[i]] <- list()
  new_cox_model_train_AIC[[i]] <- list()
  
  for (fold in seq_along(train[[i]])) {
    train_set = train[[i]][[fold]]
    
    train_set$ER = train_set$ER/10
    train_set$PR = train_set$PR/10
    train_set$Ki67 = train_set$Ki67/10
    
    train_set$grade_group1 = factor(train_set$grade_group, ordered = T)
    train_set$T_stage1 = factor(train_set$T_stage, ordered = T)
    train_set$N_stage1 = factor(train_set$N_stage, ordered = T) 
    
    
    cox_model <- coxph(Surv(iDFS, iDFS_status) ~ age + bmi + mense + FH + breast_surgery + LN_surgery + 
                         histology_group + as.numeric(grade_group) + ER + PR + HER2_IHC + Ki67 + as.numeric(T_stage) + Ki67_fp1 + Ki67_fp2 + 
                         as.numeric(N_stage) + age*FH,  
                       data = train_set)
    cox_model_train_AIC[[i]][[fold]] <- AIC(cox_model)
    

    new_cox_model <- coxph(Surv(iDFS, iDFS_status) ~ age + bmi + mense + FH + breast_surgery + LN_surgery + 
                             histology_group + grade_group1 + ER + PR + HER2_IHC + Ki67 + T_stage1 + Ki67_fp1 + Ki67_fp2 + 
                             N_stage1 + age*FH,  
                           data = train_set)
    new_cox_model_train_AIC[[i]][[fold]] <- AIC(new_cox_model)
  }
}

mean(unlist(cox_model_train_AIC))
mean(unlist(new_cox_model_train_AIC))



















