########################################################
############## SET Pathway and data input ##############


############ 1. Cut off definition: 5-year risk > 10 % ############################### 

cutoff_valid <- bind_rows(validation) # use stacked validation data
DL_cox <- coxph(Surv(iDFS, iDFS_status) ~ DL.risk_score, data = cutoff_valid)
DL_surv_prob <- survfit(DL_cox, newdata = cutoff_valid)
cutoff_valid$DL_surv_prob_year5 <- t(as.data.frame(summary(DL_surv_prob, times = 5)$surv))

cutoff_90 <- cutoff_valid %>% filter(DL_surv_prob_year5 > 0.9)
cutoff_90 <- min(cutoff_90$DL.risk_score)

# Cut-off: 0.2587269


############ 2. Risk group analysis ##################################
## Function definition ======================

# All comers 

# A. Get HR estimate
# For validation set
HR_estimate1 <- function(test_data, cutoff) {           
  
  survival <- data.frame(i = integer(), fold = integer(), time = numeric(), group = character(), rate = numeric(), rate_se = numeric())
  group_counts <- data.frame(i = integer(), fold = integer(), group = character(), count = integer())
  hr_results <- data.frame(i = integer(), fold = integer(), est = numeric(), est_se = numeric(), p = numeric())
  
  for (i in 1:20) {
    
    for (fold in 1:5) {
      
      # fit cox model
      data <- test_data[[i]][[fold]]; data$group <- ifelse(data$DL.risk_score >= cutoff, 0, 1)
      if (length(unique(data$group)) < 2) {
        next  
      }
      
      fit <- coxph(Surv(iDFS, iDFS_status) ~ group, data = data)
      cox_summary <- summary(fit)
      
      # get HR results
      est <- cox_summary$coefficients[1, "coef"]
      est_se <- cox_summary$coefficients[1, "se(coef)"]
      p <- cox_summary$coefficients[1, "Pr(>|z|)"]
      
      hr_results <- rbind(hr_results, data.frame(i = i, fold = fold, est = est, est_se = est_se, p = p))
      
      # get 1-, 3-, 5-year survival rate
      surv_fit <- survfit(Surv(iDFS, iDFS_status) ~ group, data = data)
      for (time in c(1, 3, 5)) {
        tryCatch({
          surv_summary <- summary(surv_fit, times = time, extend = TRUE)
          for (group in unique(data$group)) {
            surv_rate <- surv_summary$surv[surv_summary$strata == paste("group=", group, sep="")]
            surv_se <- surv_summary$std.err[surv_summary$strata == paste("group=", group, sep="")]
            survival <- rbind(survival, data.frame(i = i, fold = fold, time = time, group = group, rate = surv_rate, rate_se = surv_se))
          }
        }, error = function(e) {

          for (group in unique(data$group)) {
            survival <- rbind(survival, data.frame(i = i, fold = fold, time = time, group = group, rate = NA, rate_se = NA))
          }
        })
      }
      
      # get group counts and proportions
      group_count <- data %>% group_by(group) %>% summarise(count = n()) %>% mutate(i = i, fold = fold)
      total_count <- sum(group_count$count)
      group_count <- group_count %>% mutate(proportion = count / total_count)
      group_counts <- rbind(group_counts, group_count)
    }
  }
  
  return(list(hr_results = hr_results, survival = survival, group_counts = group_counts))
  
}

# For test set
HR_estimate2 <- function(test_data, cutoff) {         # for internal and external test data
  
  survival <- data.frame(i = integer(), fold = integer(), time = numeric(), group = character(), rate = numeric(), rate_se = numeric())
  group_counts <- data.frame(i = integer(), fold = integer(), group = character(), count = integer())
  hr_results <- data.frame(i = integer(), fold = integer(), hr = numeric(), hr_se = numeric(), hr_lower = numeric(), hr_upper = numeric(), p = numeric())
  
  for (i in 1:20) {
    
    data <- test_data[[i]]
    
    data$group <- ifelse(data$DL.risk_score >= cutoff, 0, 1)
    if (length(unique(data$group)) < 2) {
      next  
    }
    
    fit <- coxph(Surv(iDFS, iDFS_status) ~ group, data = data)
    cox_summary <- summary(fit)
    
    # Get HR results
    est <- cox_summary$coefficients[1, "coef"]
    est_se <- cox_summary$coefficients[1, "se(coef)"]
    p <- cox_summary$coefficients[1, "Pr(>|z|)"]
    
    hr_results <- rbind(hr_results, data.frame(i = i, fold = fold, est = est, est_se = est_se, p = p))
    
    # get 1-, 3-, 5-year survival rate
    surv_fit <- survfit(Surv(iDFS, iDFS_status) ~ group, data = data)
    for (time in c(1, 3, 5)) {
      surv_summary <- summary(surv_fit, times = time)
      for (group in unique(data$group)) {
        surv_rate <- surv_summary$surv[surv_summary$strata == paste("group=", group, sep="")]
        surv_se <- surv_summary$std.err[surv_summary$strata == paste("group=", group, sep="")]
        survival <- rbind(survival, data.frame(i = i, fold = fold, time = time, group = group, rate = surv_rate, rate_se = surv_se))
      }
    }
    
    # get group counts and proportions 
    group_count <- data %>% group_by(group) %>% summarise(count = n()) %>% mutate(i = i, fold = fold)
    total_count <- sum(group_count$count)
    group_count <- group_count %>% mutate(proportion = count / total_count)
    group_counts <- rbind(group_counts, group_count)
    
  }
  
  return(list(hr_results = hr_results, survival = survival, group_counts = group_counts))

}

# B. Get Group statistics
group_statistics <- function(HR_result, 
                             survival_rate_low1, survival_rate_high1, 
                             survival_rate_low3, survival_rate_high3, 
                             survival_rate_low5, survival_rate_high5, 
                             group_count_low, group_count_high) {
  
  # HR
  HR <- HR_result$est; HR_se <- HR_result$est_se; n_HR <- as.numeric(count(HR_result))
  
  # 1-yr rate
  rate_low1  <- survival_rate_low1$rate;   rate_se_low1 <- survival_rate_low1$rate_se;   n_rate_low1 <- as.numeric(count(survival_rate_low1))
  rate_high1 <- survival_rate_high1$rate;  rate_se_high1 <- survival_rate_high1$rate_se; n_rate_high1 <- as.numeric(count(survival_rate_high1))
  
  # 3-yr rate
  rate_low3  <- survival_rate_low3$rate;   rate_se_low3 <- survival_rate_low3$rate_se;   n_rate_low3 <- as.numeric(count(survival_rate_low3))
  rate_high3 <- survival_rate_high3$rate;  rate_se_high3 <- survival_rate_high3$rate_se; n_rate_high3 <- as.numeric(count(survival_rate_high3))
  
  # 5-yr rate
  rate_low5  <- survival_rate_low5$rate;   rate_se_low5 <- survival_rate_low5$rate_se;   n_rate_low5 <- as.numeric(count(survival_rate_low5))
  rate_high5 <- survival_rate_high5$rate;  rate_se_high5 <- survival_rate_high5$rate_se; n_rate_high5 <- as.numeric(count(survival_rate_high5))
  
  # proportion
  group_proportion_low <- group_count_low$proportion
  group_proportion_high <- group_count_high$proportion
  
  # Point estimate
  combined_HR <- exp(mean(HR))
  combined_rate_low1 <- mean(rate_low1);  combined_rate_high1 <- mean(rate_high1) 
  combined_rate_low3 <- mean(rate_low3);  combined_rate_high3 <- mean(rate_high3) 
  combined_rate_low5 <- mean(rate_low5);  combined_rate_high5 <- mean(rate_high5) 
  combined_group_proportion_low <- mean(group_proportion_low); combined_group_proportion_high <- mean(group_proportion_high)
  
  # Variance = average SE
   #HR
    Vw_HR = sum(HR_se^2) / n_HR
    Vb_HR = sum( ((HR - mean(HR))^2) )/(n_HR-1)
      combined_HR_SE = sqrt(Vw_HR + Vb_HR + Vb_HR/n_HR)
      z = mean(HR)/combined_HR_SE
   #Rate 1
    Vw_rate_low1 = sum(rate_se_low1^2) / n_rate_low1
    Vb_rate_low1 = sum((rate_low1 - combined_rate_low1)^2)/(n_rate_low1-1)
      combined_rate_low1_SE = sqrt(Vw_rate_low1 + Vb_rate_low1 + Vb_rate_low1/n_rate_low1)
    Vw_rate_high1 = sum(rate_se_high1^2) / n_rate_high1
    Vb_rate_high1 = sum((rate_high1 - combined_rate_high1)^2)/(n_rate_high1-1)
      combined_rate_high1_SE = sqrt(Vw_rate_high1 + Vb_rate_high1 + Vb_rate_high1/n_rate_high1)      
    # Rate 3
    Vw_rate_low3 = sum(rate_se_low3^2) / n_rate_low3
    Vb_rate_low3 = sum((rate_low3 - combined_rate_low3)^2)/(n_rate_low3-1)
      combined_rate_low3_SE = sqrt(Vw_rate_low3 + Vb_rate_low3 + Vb_rate_low3/n_rate_low3)
    Vw_rate_high3 = sum(rate_se_high3^2) / n_rate_high3
    Vb_rate_high3 = sum((rate_high3 - combined_rate_high3)^2)/(n_rate_high3-1)
      combined_rate_high3_SE = sqrt(Vw_rate_high3 + Vb_rate_high3 + Vb_rate_high3/n_rate_high3)       
    # Rate 3
    Vw_rate_low5 = sum(rate_se_low5^2) / n_rate_low5
    Vb_rate_low5 = sum((rate_low5 - combined_rate_low5)^2)/(n_rate_low5-1)
      combined_rate_low5_SE = sqrt(Vw_rate_low5 + Vb_rate_low5 + Vb_rate_low5/n_rate_low5)
    Vw_rate_high5 = sum(rate_se_high5^2) / n_rate_high5
    Vb_rate_high5 = sum((rate_high5 - combined_rate_high5)^2)/(n_rate_high5-1)
      combined_rate_high5_SE = sqrt(Vw_rate_high5 + Vb_rate_high5 + Vb_rate_high5/n_rate_high5)           
  
  # Output CI 
  # HR
  HR_CI_lower <- exp(mean(HR) - 1.96 * combined_HR_SE)
  HR_CI_upper <- exp(mean(HR) + 1.96 * combined_HR_SE)
  
  # 1-yr rate
  rate_low1_CI_lower  <- combined_rate_low1 - 1.96 * combined_rate_low1_SE
  rate_low1_CI_upper  <- min(combined_rate_low1 + 1.96 * combined_rate_low1_SE, 1) 
  rate_high1_CI_lower  <- combined_rate_high1 - 1.96 * combined_rate_high1_SE
  rate_high1_CI_upper  <- min(combined_rate_high1 + 1.96 * combined_rate_high1_SE, 1) 
  
  # 3-yr rate
  rate_low3_CI_lower   <- combined_rate_low3 - 1.96 * combined_rate_low3_SE
  rate_low3_CI_upper   <- min(combined_rate_low3 + 1.96 * combined_rate_low3_SE, 1) 
  rate_high3_CI_lower  <- combined_rate_high3 - 1.96 * combined_rate_high3_SE
  rate_high3_CI_upper  <- min(combined_rate_high3 + 1.96 * combined_rate_high3_SE, 1) 
  
  # 3-yr rate
  rate_low5_CI_lower   <- combined_rate_low5 - 1.96 * combined_rate_low5_SE
  rate_low5_CI_upper   <- min(combined_rate_low5 + 1.96 * combined_rate_low5_SE, 1) 
  rate_high5_CI_lower  <- combined_rate_high5 - 1.96 * combined_rate_high5_SE
  rate_high5_CI_upper  <- min(combined_rate_high5 + 1.96 * combined_rate_high5_SE, 1) 
  
  
  # Output summary table
  Group_summary <- data.frame(
    # group proportion
    percentage_high_risk = round(combined_group_proportion_high * 100, 2),
    
    # 1-yr survival rate
    survival_rate_low_1yr  = paste0(round(combined_rate_low1 * 100,2), 
                                    " (", round(rate_low1_CI_lower * 100, 2)," to ", round(rate_low1_CI_upper * 100, 2), ")"),
    
    survival_rate_high_1yr = paste0(round(combined_rate_high1 * 100,2), 
                                    " (", round(rate_high1_CI_lower * 100, 2)," to ", round(rate_high1_CI_upper * 100, 2), ")"),
    
    
    # 3-yr survival rate
    survival_rate_low_3yr  = paste0(round(combined_rate_low3 * 100, 2), 
                                    " (", round(rate_low3_CI_lower * 100, 2)," to ", round(rate_low3_CI_upper * 100, 2), ")"),
    
    survival_rate_high_3yr = paste0(round(combined_rate_high3 * 100,2), 
                                    " (", round(rate_high3_CI_lower * 100, 2)," to ", round(rate_high3_CI_upper * 100, 2), ")"),
    
    
    # 5-yr survival rate
    survival_rate_low_5yr  = paste0(round(combined_rate_low5 * 100, 2), 
                                    " (", round(rate_low5_CI_lower * 100, 2)," to ", round(rate_low5_CI_upper * 100, 2), ")"),
    
    survival_rate_high_5yr = paste0(round(combined_rate_high5 * 100, 2), 
                                    " (", round(rate_high5_CI_lower * 100, 2)," to ", round(rate_high5_CI_upper * 100, 2), ")"),
    
    
    # Hazard ratio
    HR = paste0(round(combined_HR,3), " (", round(HR_CI_lower, 3), " to ", round(HR_CI_upper, 3), ")"),
    HR_p.value = round(2 * (1-pnorm(abs(z))), 4)
  )
  
  return(Group_summary)
}

# C. Output result table: HR, survival rate, and group count
group_result = function(data1, data2, data3){
  # Validation
  HR_validation <- HR_estimate1(data1, cutoff_90)$hr_results
  survival_validation <- HR_estimate1(data1, cutoff_90)$survival
  survival_validation_low_1yr  <- subset(survival_validation, group == 0 & time == 1); survival_validation_high_1yr <- subset(survival_validation, group == 1 & time == 1)
  survival_validation_low_3yr  <- subset(survival_validation, group == 0 & time == 3); survival_validation_high_3yr <- subset(survival_validation, group == 1 & time == 3)
  survival_validation_low_5yr  <- subset(survival_validation, group == 0 & time == 5); survival_validation_high_5yr <- subset(survival_validation, group == 1 & time == 5)
  
  group_count_validation <- HR_estimate1(data1, cutoff_90)$group_counts   
  group_count_validation_low  <- subset(group_count_validation, group == 0); group_count_validation_high <- subset(group_count_validation, group == 1)
  
  # Internal test 
  HR_internal <- HR_estimate2(data2, cutoff_90)$hr_results
  survival_internal <- HR_estimate2(data2, cutoff_90)$survival
  survival_internal_low_1yr  <- subset(survival_internal, group == 0 & time == 1); survival_internal_high_1yr <- subset(survival_internal, group == 1 & time == 1)
  survival_internal_low_3yr  <- subset(survival_internal, group == 0 & time == 3); survival_internal_high_3yr <- subset(survival_internal, group == 1 & time == 3)
  survival_internal_low_5yr  <- subset(survival_internal, group == 0 & time == 5); survival_internal_high_5yr <- subset(survival_internal, group == 1 & time == 5)
  
  group_count_internal <- HR_estimate2(data2, cutoff_90)$group_counts   
  group_count_internal_low  <- subset(group_count_internal, group == 0); group_count_internal_high <- subset(group_count_internal, group == 1)
  
  
  # External test 
  HR_external <- HR_estimate2(data3, cutoff_90)$hr_results
  survival_external <- HR_estimate2(data3, cutoff_90)$survival
  survival_external_low_1yr  <- subset(survival_external, group == 0 & time == 1); survival_external_high_1yr <- subset(survival_external, group == 1 & time == 1)
  survival_external_low_3yr  <- subset(survival_external, group == 0 & time == 3); survival_external_high_3yr <- subset(survival_external, group == 1 & time == 3)
  survival_external_low_5yr  <- subset(survival_external, group == 0 & time == 5); survival_external_high_5yr <- subset(survival_external, group == 1 & time == 5)
  
  group_count_external <- HR_estimate2(data3, cutoff_90)$group_counts   
  group_count_external_low  <- subset(group_count_external, group == 0); group_count_external_high <- subset(group_count_external, group == 1)
  
  Group_result_summary <- cbind(
    as.data.frame(t(group_statistics(HR_validation,
                                     survival_validation_low_1yr,survival_validation_high_1yr,
                                     survival_validation_low_3yr,survival_validation_high_3yr,
                                     survival_validation_low_5yr,survival_validation_high_5yr,
                                     group_count_validation_low,group_count_validation_high))),
    
    as.data.frame(t(group_statistics(HR_internal,
                                     survival_internal_low_1yr,survival_internal_high_1yr,
                                     survival_internal_low_3yr,survival_internal_high_3yr,
                                     survival_internal_low_5yr,survival_internal_high_5yr,
                                     group_count_internal_low,group_count_internal_high))),
    
    as.data.frame(t(x = group_statistics(HR_external,
                                         survival_external_low_1yr,survival_external_high_1yr,
                                         survival_external_low_3yr,survival_external_high_3yr,
                                         survival_external_low_5yr,survival_external_high_5yr,
                                         group_count_external_low,group_count_external_high))))
  
  colnames(Group_result_summary) <- c("Validation","Internal test","External test")
  
  return(Group_result_summary)
}


########################################################################
################ Subgroup: Get HR, P, and P interaction ################
########################################################################

# Validation cohort
subgroup_estimate1 <- function(test_data, cutoff, factor_var) {           
  
  hr_results <- data.frame(i = integer(), fold = integer(), subgroup_level = character(), 
                           est = numeric(), est_se = numeric(), 
                           est_interaction = numeric(), se_interaction = numeric())
 
  for (i in 1:20) {
    for (fold in 1:5) {
      data <- test_data[[i]][[fold]]

      data[[factor_var]] <- relevel(factor(data[[factor_var]]), ref = "0")
      data$group <- factor(ifelse(data$DL.risk_score >= cutoff, 0, 1))

      formula <- as.formula(paste("Surv(iDFS, iDFS_status) ~ group *", factor_var))
      fit <- coxph(formula, data = data)
      cox_summary <- summary(fit)
      vcov_matrix <- vcov(fit)
      
      est <- cox_summary$coefficients["group1", "coef"]
      est_se <- sqrt(vcov_matrix["group1", "group1"]) 
      
      levels <- unique(data[[factor_var]])
      for (level in levels) {
        
        est_interaction <- NA
        se_interaction <- NA
        
        if (level != "0") { 
          
          interaction_term <- paste0("group1:", factor_var, level)
          est_interaction <- cox_summary$coefficients[interaction_term, "coef"]
          
          est_level <- est + est_interaction
          est_se_level <- sqrt(vcov_matrix["group1", "group1"] + 
                                 vcov_matrix[interaction_term, interaction_term] + 
                                 2 * vcov_matrix["group1", interaction_term])
          
          se_interaction <- cox_summary$coefficients[interaction_term, "se(coef)"]
          
        } else {  
          est_level <- est
          est_se_level <- est_se
        }
        
        hr_results <- rbind(hr_results, 
                            data.frame(i = i, 
                                       subgroup_level = level, 
                                       est = est_level, est_se = est_se_level, 
                                       est_interaction = est_interaction, se_interaction = se_interaction))
      }
    }
    return(hr_results)
  }
}

subgroup_summary1 <- function(test_data, cutoff, factor_var) {
  
  data <- subgroup_estimate1(test_data, cutoff, factor_var)
  
  unique_subgroups <- unique(data$subgroup_level)
  
  result <- data.frame(
    subgroup_level = unique_subgroups,
    HR = numeric(length(unique_subgroups)),
    HR_lower = numeric(length(unique_subgroups)),
    HR_upper = numeric(length(unique_subgroups)),
    HR_summary = character(length(unique_subgroups)),
    HR_p_value = numeric(length(unique_subgroups)),
    p.interaction = numeric(length(unique_subgroups)),
    stringsAsFactors = FALSE
  )
  

  for (subgroup in unique(data$subgroup_level)) {

    subgroup_data <- data[data$subgroup_level == subgroup, ]
    
    HR <- subgroup_data$est
    HR_se <- subgroup_data$est_se
    n_HR <- nrow(subgroup_data)  
    
    Vw_HR <- sum(HR_se^2) / n_HR
    Vb_HR <- sum((HR - mean(HR))^2) / (n_HR - 1)
    combined_HR_SE <- sqrt(Vw_HR + Vb_HR + Vb_HR / n_HR)
    
    combined_HR <- exp(mean(HR))
    HR_lower <- exp(mean(HR) - 1.96 * combined_HR_SE)
    HR_upper <- exp(mean(HR) + 1.96 * combined_HR_SE)
    
    HR_z_value <- mean(HR) / combined_HR_SE
    

    interaction_effect <- subgroup_data$est_interaction
    interaction_se <- subgroup_data$se_interaction
    Vw_interaction <- sum(interaction_se^2) / n_HR
    Vb_interaction <- sum((interaction_effect - mean(interaction_effect))^2) / (n_HR - 1)
    combined_interaction_SE <- sqrt(Vw_interaction + Vb_interaction + Vb_interaction / n_HR)
    combined_interaction <- mean(interaction_effect)
    interaction_z_value <- mean(interaction_effect) / combined_interaction_SE
    
    result[result$subgroup_level == subgroup, ] <- list(
      subgroup_level = subgroup,
      HR = combined_HR, 
      HR_lower = HR_lower,
      HR_upper = HR_upper,
      HR_summary = paste0(round(combined_HR,3), " (", 
                          round(HR_lower,3), " to ", 
                          round(HR_upper,3),")"),
      HR_p_value = round(2 * (1 - pnorm(abs(HR_z_value))), 4),
      p.interaction = round(2 * (1 - pnorm(abs(interaction_z_value))), 4)
    )
  }
  
  return(result)
}


# Test cohort
subgroup_estimate2 <- function(test_data, cutoff, factor_var) {           
  
  hr_results <- data.frame(i = integer(), subgroup_level = character(), 
                           est = numeric(), est_se = numeric(), 
                           est_interaction = numeric(), se_interaction = numeric())
  
  for (i in 1:20) {
    data <- test_data[[i]]
    
    data[[factor_var]] <- relevel(factor(data[[factor_var]]), ref = "0")
    data$group <- factor(ifelse(data$DL.risk_score >= cutoff, 0, 1))
    
    formula <- as.formula(paste("Surv(iDFS, iDFS_status) ~ group *", factor_var))
    fit <- coxph(formula, data = data)
    cox_summary <- summary(fit)
    vcov_matrix <- vcov(fit)
    
    est <- cox_summary$coefficients["group1", "coef"]
    est_se <- sqrt(vcov_matrix["group1", "group1"]) 
    
    levels <- unique(data[[factor_var]])
    for (level in levels) {
      
      est_interaction <- NA
      se_interaction <- NA
      
      if (level != "0") {  
        
        interaction_term <- paste0("group1:", factor_var, level)
        est_interaction <- cox_summary$coefficients[interaction_term, "coef"]
        
        est_level <- est + est_interaction
        est_se_level <- sqrt(vcov_matrix["group1", "group1"] + 
                               vcov_matrix[interaction_term, interaction_term] + 
                               2 * vcov_matrix["group1", interaction_term])
        
        se_interaction <- cox_summary$coefficients[interaction_term, "se(coef)"]
        
      } else { 
        est_level <- est
        est_se_level <- est_se
      }
      
      hr_results <- rbind(hr_results, 
                          data.frame(i = i, 
                                     subgroup_level = level, 
                                     est = est_level, est_se = est_se_level, 
                                     est_interaction = est_interaction, se_interaction = se_interaction))
    }
  }
  return(hr_results)
}

subgroup_summary2 <- function(test_data, cutoff, factor_var) {
  
  data <- subgroup_estimate2(test_data, cutoff, factor_var)
  
  unique_subgroups <- unique(data$subgroup_level)
  
  result <- data.frame(
    subgroup_level = unique_subgroups,
    HR = numeric(length(unique_subgroups)),
    HR_lower = numeric(length(unique_subgroups)),
    HR_upper = numeric(length(unique_subgroups)),
    HR_summary = character(length(unique_subgroups)),
    HR_p_value = numeric(length(unique_subgroups)),
    p.interaction = numeric(length(unique_subgroups)),
    stringsAsFactors = FALSE
  )
  
  
  for (subgroup in unique(data$subgroup_level)) {
    
    subgroup_data <- data[data$subgroup_level == subgroup, ]
    
    HR <- subgroup_data$est
    HR_se <- subgroup_data$est_se
    n_HR <- nrow(subgroup_data)  
    
    Vw_HR <- sum(HR_se^2) / n_HR
    Vb_HR <- sum((HR - mean(HR))^2) / (n_HR - 1)
    combined_HR_SE <- sqrt(Vw_HR + Vb_HR + Vb_HR / n_HR)
    
    combined_HR <- exp(mean(HR))
    HR_lower <- exp(mean(HR) - 1.96 * combined_HR_SE)
    HR_upper <- exp(mean(HR) + 1.96 * combined_HR_SE)
    
    HR_z_value <- mean(HR) / combined_HR_SE
    
    
    interaction_effect <- subgroup_data$est_interaction
    interaction_se <- subgroup_data$se_interaction
    Vw_interaction <- sum(interaction_se^2) / n_HR
    Vb_interaction <- sum((interaction_effect - mean(interaction_effect))^2) / (n_HR - 1)
    combined_interaction_SE <- sqrt(Vw_interaction + Vb_interaction + Vb_interaction / n_HR)
    combined_interaction <- mean(interaction_effect)
    interaction_z_value <- mean(interaction_effect) / combined_interaction_SE
    
    result[result$subgroup_level == subgroup, ] <- list(
      subgroup_level = subgroup,
      HR = combined_HR, 
      HR_lower = HR_lower,
      HR_upper = HR_upper,
      HR_summary = paste0(round(combined_HR,3), " (", 
                          round(HR_lower,3), " to ", 
                          round(HR_upper,3),")"),
      HR_p_value = round(2 * (1 - pnorm(abs(HR_z_value))), 4),
      p.interaction = round(2 * (1 - pnorm(abs(interaction_z_value))), 4)
    )
  }
  
  return(result)
}


