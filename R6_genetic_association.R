########################################################
############## SET Pathway and data input ##############

###########################################################################
############# I. association between gene data and risk score #############
###########################################################################

reference_genes <- c("ACTIN_CT", "GAPDH_CT", "GUS_CT", "RPLPO_CT", "TFRC_CT")

gene_var <- c("GRB7_CT","HER2_CT","ER_CT", "PR_CT", "Bcl2_CT", "CEGP1_CT",
              "CCNB1_CT", "Ki67_CT", "MYBL2_CT", "STK15_CT", "SURV_CT", "CTSL2_CT",
              "STMY3_CT", "CD68_CT", "GSTM1_CT", "BAG1_CT", "RS",
              "AR_CT", "CCND1_CT", "CDK6_CT","E2F1_CT",  "GluT1_CT",  "HK2_CT", "S6K2_CT")

#############################################
## Regression 1: Crude results  #############
#############################################

association_crude_results <- list()

for (i in seq_along(data)) {
  
  fold_results <- list()
  
  for (fold in seq_along(data[[i]])) {
    
    dt <- data[[i]][[fold]]

    genes <- intersect(names(dt), gene_var) 
    ref_genes <- intersect(names(dt), reference_genes)
    dt$ref_mean_CT <- rowMeans(dt[ref_genes], na.rm = TRUE)
 
    coefficients <- numeric(length(genes))
    standard_errors <- numeric(length(genes))
    linear.p.value <- numeric(length(genes))
    
    for (j in seq_along(genes)) {
      gene <- genes[j]
      dt$normalized_gene_CT <- dt[[gene]] - dt$ref_mean_CT
      
      fit <- lm(DL.risk_score ~ scale(dt$normalized_gene_CT), data = dt)
      summary_fit <- summary(fit)
      coefficients[j] <- summary_fit$coefficients[2, 1] 
      standard_errors[j] <- summary_fit$coefficients[2, 2] 
      linear.p.value[j] <- summary_fit$coefficients[2, 4]  
    }
    
    fold_df <- data.frame(Gene = genes, 
                          Coefficient = coefficients, Standard_Error = standard_errors, 
                          linear.p.value = linear.p.value)
    fold_results[[fold]] <- fold_df
  }
  association_crude_results[[i]] <- fold_results
}



association_crude_summary <-  do.call(rbind, lapply(association_crude_results, function(fold_list) do.call(rbind, fold_list)))

association_crude_summary <- association_crude_summary %>%
  group_by(Gene) %>%
  summarise(est = mean(Coefficient),
            Vw = sum(Standard_Error^2) / 100,
            Vb = sum( ((Coefficient - mean(Coefficient))^2) )/(100-1),
            combined_SE = sqrt(Vw + Vb + Vb/100),
            
            est_lower = est - 1.96*combined_SE,
            est_upper = est + 1.96*combined_SE) %>%
  
  mutate(z.sta = est/combined_SE) %>%
  mutate(p.value = round(2 * (1-pnorm(abs(z.sta))), 4))


association_crude_summary <- association_crude_summary %>%  
  mutate(linear.coefficient = round(est, 4)) %>%
  mutate(linear.lower = round(est - 1.96*combined_SE, 4)) %>%
  mutate(linear.upper = round(est + 1.96*combined_SE, 4)) %>%
  mutate(linear.association = paste0(linear.coefficient, "(", linear.lower, " to ", linear.upper, ")")) %>%
  mutate(linear.p.value = round(p.value, 3)) %>%
  select(Gene, linear.coefficient, linear.lower, linear.upper, linear.association, linear.p.value)

###################################################################
## Regression 2: Adjusted for age and family history ##############
###################################################################

association_adj_results1 <- list()

for (i in seq_along(data)) {

  fold_results <- list()
  
  for (fold in seq_along(data[[i]])) {

    dt <- data[[i]][[fold]]
    
    genes <- intersect(names(dt), gene_var) 
    ref_genes <- intersect(names(dt), reference_genes)
    dt$ref_mean_CT <- rowMeans(dt[ref_genes], na.rm = TRUE)
    
    coefficients <- numeric(length(genes))
    standard_errors <- numeric(length(genes))
    linear.p.value <- numeric(length(genes))
    
    for (j in seq_along(genes)) {
      gene <- genes[j]
      dt$normalized_gene_CT <- dt[[gene]] - dt$ref_mean_CT
      
      fit <- lm(DL.risk_score ~ scale(dt$normalized_gene_CT) + age + FH, data = dt)
      summary_fit <- summary(fit)
      coefficients[j] <- summary_fit$coefficients[2, 1] 
      standard_errors[j] <- summary_fit$coefficients[2, 2] 
      linear.p.value[j] <- summary_fit$coefficients[2, 4]  
    }
    
    fold_df <- data.frame(Gene = genes, 
                          Coefficient = coefficients, Standard_Error = standard_errors, 
                          linear.p.value = linear.p.value)
    fold_results[[fold]] <- fold_df
  }
  association_adj_results1[[i]] <- fold_results
}



association_adj_summary1 <-  do.call(rbind, lapply(association_adj_results1, function(fold_list) do.call(rbind, fold_list)))

association_adj_summary1 <- association_adj_summary1 %>%
  group_by(Gene) %>%
  summarise(est = mean(Coefficient),
            Vw = sum(Standard_Error^2) / 100,
            Vb = sum( ((Coefficient - mean(Coefficient))^2) )/(100-1),
            combined_SE = sqrt(Vw + Vb + Vb/100),
            
            est_lower = est - 1.96*combined_SE,
            est_upper = est + 1.96*combined_SE) %>%
  
  mutate(z.sta = est/combined_SE) %>%
  mutate(p.value = round(2 * (1-pnorm(abs(z.sta))), 4))
         

association_adj_summary1 <- association_adj_summary1 %>%  
      mutate(linear.coefficient = round(est, 4)) %>%
      mutate(linear.lower = round(est - 1.96*combined_SE, 4)) %>%
      mutate(linear.upper = round(est + 1.96*combined_SE, 4)) %>%
      mutate(linear.association = paste0(linear.coefficient, "(", linear.lower, " to ", linear.upper, ")")) %>%
      mutate(linear.p.value = round(p.value, 3)) %>%
    select(Gene, linear.coefficient, linear.lower, linear.upper, linear.association, linear.p.value)



######## Regression 3: Adjusted for age, menopausal status, and family history ##############

association_adj_results2 <- list()

for (i in seq_along(data)) {
  
  fold_results <- list()
  
  for (fold in seq_along(data[[i]])) {
    
    dt <- data[[i]][[fold]]
    
    genes <- intersect(names(dt), gene_var) 
    ref_genes <- intersect(names(dt), reference_genes)
    dt$ref_mean_CT <- rowMeans(dt[ref_genes], na.rm = TRUE)
    
    coefficients <- numeric(length(genes))
    standard_errors <- numeric(length(genes))
    linear.p.value <- numeric(length(genes))
    
    for (j in seq_along(genes)) {
      gene <- genes[j]
      dt$normalized_gene_CT <- dt[[gene]] - dt$ref_mean_CT
      
      fit <- lm(DL.risk_score ~ scale(dt$normalized_gene_CT) + age + mense + FH, data = dt)
      summary_fit <- summary(fit)
      coefficients[j] <- summary_fit$coefficients[2, 1] 
      standard_errors[j] <- summary_fit$coefficients[2, 2] 
      linear.p.value[j] <- summary_fit$coefficients[2, 4]  
    }
    
    fold_df <- data.frame(Gene = genes, 
                          Coefficient = coefficients, Standard_Error = standard_errors, 
                          linear.p.value = linear.p.value)
    fold_results[[fold]] <- fold_df
  }
  association_adj_results2[[i]] <- fold_results
}



association_adj_summary2 <-  do.call(rbind, lapply(association_adj_results2, function(fold_list) do.call(rbind, fold_list)))

association_adj_summary2 <- association_adj_summary2 %>%
  group_by(Gene) %>%
  summarise(est = mean(Coefficient),
            Vw = sum(Standard_Error^2) / 100,
            Vb = sum( ((Coefficient - mean(Coefficient))^2) )/(100-1),
            combined_SE = sqrt(Vw + Vb + Vb/100),
            
            est_lower = est - 1.96*combined_SE,
            est_upper = est + 1.96*combined_SE) %>%
  
  mutate(z.sta = est/combined_SE) %>%
  mutate(p.value = round(2 * (1-pnorm(abs(z.sta))), 4))


association_adj_summary2 <- association_adj_summary2 %>%  
  mutate(linear.coefficient = round(est, 4)) %>%
  mutate(linear.lower = round(est - 1.96*combined_SE, 4)) %>%
  mutate(linear.upper = round(est + 1.96*combined_SE, 4)) %>%
  mutate(linear.association = paste0(linear.coefficient, "(", linear.lower, " to ", linear.upper, ")")) %>%
  mutate(linear.p.value = round(p.value, 3)) %>%
  select(Gene, linear.coefficient, linear.lower, linear.upper, linear.association, linear.p.value)


