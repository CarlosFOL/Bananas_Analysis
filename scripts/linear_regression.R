library(car)
library(MASS)
library(nortest)
library(leaps)
library(lmtest)


###################### Search the best model ###################################
get_AIC <- function(models, dataset){
  ########################################
  # Compute the AIC for each model       #
  ########################################
  AIC_val <- c()
  for (i in 1:length(models)){
    model <- lm(models[i], data = dataset)
    AIC_val <- c(AIC_val, AIC(model))
  }
  return (AIC_val)
}

get_expl_var <- function(which_obj, response){
  #############################################
  # Get the explanatory variables that belong #
  # to the model.                             #
  #############################################
  #k <- ncol(which_obj)-1 # Not include \beta0
  expl_vars <- c()
  for (i in 1:nrow(which_obj)){
    vars <- names(which(which_obj[i,] == TRUE)[-1])
    value <- ""
    for (j in 1:length(vars)){
      value <- paste0(value, vars[j], sep="+")
    }
    # To remove the last +
    value <- paste0(response, "~", substring(value, 1, nchar(value)-1))
    expl_vars <- c(expl_vars, value)
  }
  return(expl_vars)
}


search_best_model <- function(dataset, response_idx){
  ###########################################################
  # Calculate the best model by using different number of   #
  # explanatory variables                                   #
  #
  # Parameters
  # ----------
  # dataset: data.frame
  #    Dataset = Matrix with the features and response variable
  # response_idx: numeric
  #    Matrix's index of the response variable.
  ##########################################################
  # Build the feature matrix, but it does not take into account the categorical
  # variables.
  num_col <- sapply(dataset, function(col) class(col) == "numeric" || class(col) == "integer")
  X <- dataset[, num_col][, -response_idx] # Feature matrix
  # Response variable:
  Y <- dataset[, response_idx]; response <- colnames(dataset)[response_idx]
  k <- ncol(X) # Numero de explicativas
  # Exhaustive search.
  exh_model <- regsubsets(x = X, 
                          y = Y, method = "exhaustive")
  # Summary de la busqueda exhaustiva
  summ_exh <- summary(exh_model)
  # Compute the important statistic to measure the goodness and complexity of
  # a model.
  stats <- attributes(summ_exh)$names[2:6]
  # Stats Table
  table <- c() 
  for (i in 1:length(stats)){
    row <- unlist(summ_exh[stats[i]])
    table <- rbind(table, row)
  }
  table <- data.frame(t(matrix(table, nrow = 5)), row.names = 1:(ncol(X)-1))
  colnames(table) <- c("R^2", "SS(R)", "R^2_adj", "Cp", "BIC")
  # Add the features variables for each model
  table$Modelo = get_expl_var(summ_exh$which, response)
  # Add AIC column
  table$AIC = get_AIC(table$Modelo, dataset)
  table <- table[, c(2, 1, 3, 7, 5, 4, 6)]
  return(table)
}

###################### INFLUENCE ANALYSIS ####################################
get_influential_points <- function(values, threshold, statistic){
  indexes <- which(values >= threshold) 
  cat("Observations with high", statistic, ":", indexes, "\n")
}

influential_analysis <- function(model, alpha){
  # Highest leverage
  leverages <- hatvalues(model)
  get_influential_points(leverages, 2 * mean(leverages), "leverage")
  # Highest DFFIT
  n <- length(model$fittet.values); k <- nrow(model$coefficients) - 1
  dffits_val <- abs(dffits(model))
  get_influential_points(dffits_val, 2/sqrt( (k+1)/n ), "DFFIT")
  # Highest DFBETA
  th_dfbeta <- 2/sqrt(n)
  dfbetas_infl <- apply(dfbetas(model), 
                       MARGIN = 1, 
                       FUN = function(row) any(abs(row) > th_dfbeta) )
  cat("Observations with high DFBETA:", which(dfbetas_infl == T), "\n")
  # Highest Cook's distance
  dist_cook <- cooks.distance(model)
  get_influential_points(dist_cook, qf(alpha, k+1, n-k-1), "Cook's Distance")
}


################### DIAGNOSIS DEL MODELO ######################################

diag_model <- function(model){
  ##############################################################
  # Check the regression hypotheses of our model.              #
  ##############################################################
  # Use the Student's residuals
  std_residuals <- stdres(model)
  # Linearity 
  x11()
  plot(model, which=1)
  # Normal Distribution (H_0: e ~ N)
  p_value_lillie <- lillie.test(std_residuals)$p.value
  p_value_shapiro <- shapiro.test(std_residuals)$p.value
  cat("Normal hipothesis\nLillie:", p_value_lillie, "\nShapiro", p_value_shapiro, "\n\n")
  # Homocedasticity (H_0: Homocedasticity)
  p_value_bp <- bptest(model)$p.value
  cat("Homocedasticity hipothesis\nBreusch-Pagan:", p_value_bp, "\n\n")
  # Randomness (H_0: Residuals are independent)
  p_value_box <- Box.test(std_residuals, type = "Ljung-Box", lag = 20)$p.value
  cat("Randomness hypothesis\nLjung-Box:", p_value_box)
}



