library(car)
library(MASS)
library(nortest)
library(leaps)
library(lmtest)


###################### Search the best model ###################################
get_AIC <- function(models, dataset){
  ########################################
  # Calcular el AIC para cada modelo.    #
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
  # Obtener las variables explicativas que deben
  # ser incluidas en el modelo.
  #############################################
  k <- ncol(which_obj)-1 # Excluyendo el intercepto
  # Obtiene las explicativas que se necesitan para el modelo
  expl_vars <- c()
  for (i in 1:k){
    # Obtener los nombres de las variables explicativas sin tener en cuenta \beta0
    vars <- names(which(which_obj[i,] == TRUE)[-1])
    value <- ""
    for (j in 1:length(vars)){
      value <- paste0(value, vars[j], sep="+")
    }
    # Para remove el ultimo "+"
    value <- paste0(response, "~", substring(value, 1, nchar(value)-1))
    expl_vars <- c(expl_vars, value)
  }
  return(expl_vars)
}


search_best_model <- function(dataset, response_idx){
  ###########################################################
  # Calcula el mejor modelo que se puede formar a partir
  # de un numero dado de variables explicativas.
  #
  # Parameters
  # ----------
  # dataset: data.frame
  #    Conjunto de datos
  # response_idx: numeric
  #    Indice de la columna que indica la variable respuesta
  ##########################################################
  # Se define la matriz de features y el vector response
  X <- dataset[, -response_idx]
  Y <- dataset[, response_idx]; response <- colnames(dataset)[response_idx]
  k <- ncol(X) # Numero de explicativas
  # Busqueda exhaustiva del modelo.
  exh_model <- regsubsets(x = dataset[, -response_idx], 
                          y = dataset[, response_idx],
                          nbest = 1)
  # Summary de la busqueda exhaustiva
  summ_exh <- summary(exh_model)
  # Se extraen aquellos estadisticos de interes
  stats <- attributes(summ_exh)$names[2:6]
  # Creacion de tabla estadistica
  table <- c() 
  for (i in 1:length(stats)){
    row <- unlist(summ_exh[stats[i]])
    table <- rbind(table, row)
  }
  table <- data.frame(t(matrix(table, nrow = 5)),
                      row.names = 1: (ncol(dataset) - 1))
  colnames(table) <- c("R^2", "SS(R)", "R^2_adj", "Cp", "BIC")
  # Adicion de las explicativas pertinentes
  table$Modelo = get_expl_var(summ_exh$which, response)
  # Calculo del AIC
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
  # Mayor leverage
  leverages <- hatvalues(model)
  get_influential_points(leverages, 2 * mean(leverages), "leverage")
  # Mayor DFFIT
  n <- length(model$fittet.values); k <- nrow(model$coefficients) - 1
  dffits_val <- abs(dffits(model))
  get_influential_points(dffits_val, 2/sqrt( (k+1)/n ), "DFFIT")
  # Mayor DFBETA
  th_dfbeta <- 2/sqrt(n)
  dfbetas_infl <- apply(dfbetas(model), 
                       MARGIN = 1, 
                       FUN = function(row) any(abs(row) > th_dfbeta) )
  cat("Observations with high DFBETA:", which(dfbetas_infl == T), "\n")
  # Mayor Distancia de Cook
  dist_cook <- cooks.distance(model)
  get_influential_points(dist_cook, qf(alpha, k+1, n-k-1), "Cook's Distance")
}


################### DIAGNOSIS DEL MODELO ######################################

diag_model <- function(model){
  ##############################################################
  # Verificar que se cumplan las hipotesis estructurales del   #
  # modelo de regresion.                                       #
  ##############################################################
  # Se emplean los residuos studentizados
  std_residuals <- stdres(model)
  # Linealidad
  x11()
  plot(model, which=1)
  # Normalidad
  p_value_lillie <- lillie.test(std_residuals)$p.value
  p_value_shapiro <- shapiro.test(std_residuals)$p.value
  cat("HIPOTESIS DE NORMALIDAD:\nLillie:", p_value_lillie, "\nShapiro", p_value_shapiro, "\n")
  # Homocedasticidad (H_0: Homocedasticidad)
  p_value_bp <- bptest(model)$p.value
  cat("HIPOTESIS DE HOMOCEDASTICIDAD:\nBreusch-Pagan:", p_value_bp, "\n")
  # Aleatoriedad (H_0: Los residuos son independientes)
  p_value_box <- Box.test(std_residuals, type = "Ljung-Box", lag = 20)$p.value
  cat("HIPOTESIS DE ALEATORIEDAD:\nLjung-Box:", p_value_box)
}



