#'  lambdaCalculationClusteredOption.R    Function for Calculation of the penalty parameter - with a clustering option
#'  
#'  This function build over the function `hdm::lambdaCalculation` . Only the addition are commented.
#'  To compute the penalty level for clustered Lasso, cluster.presence must be set at TRUE but also homoscedastic to FALSE 
#'  and X.dependent.lambda to FALSE.
#'  
#'  
#'  @param penalty A list of parameters, with in addition cluster.presence.
#'  @param cluster_index A n-sized vector with identifiers of clusters.
#'  

lambdaCalculationClusteredOption<-function (
  penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, 
                 lambda.start = NULL, c = 1.1, gamma = 0.1, cluster.presence=FALSE),  
  y = NULL, x = NULL,
  cluster_index=NULL) {
  checkmate::checkChoice(penalty$X.dependent.lambda, c(TRUE, 
                                                       FALSE, NULL))
  checkmate::checkChoice(penalty$homoscedastic, c(TRUE, FALSE, 
                                                  "none"))
  if (!exists("homoscedastic", where = penalty)) 
    penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty)) 
    penalty$X.dependent.lambda = "FALSE"
  if (!exists("cluster.presence", where = penalty)) {
    penalty$cluster.presence = !is.null(cluster_index)
  }
  if (!exists("gamma", where = penalty) & penalty$homoscedastic != 
      "none") {
    penalty$gamma = 0.1
  }
  if (!exists("cluster", where = penalty) & penalty$homoscedastic != 
      "none") {
    penalty$gamma = 0.1
  }
  
  ## == HOMOSKEDACTIC CASE == ##
  if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == 
      FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 * 
                                                                    p))
    Ups0 <- sqrt(var(y))
    lambda <- rep(lambda0 * Ups0, p)
  }
  if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == 
      TRUE) {
    if (!exists("numSim", where = penalty)) {
      penalty$numSim = 5000
    }
    p <- dim(x)[2]
    n <- dim(x)[1]
    R <- penalty$numSim
    sim <- vector("numeric", length = R)
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each = p), ncol = p, byrow = TRUE)
      sim[l] <- n * max(2 * colMeans(x * g))
    }
    lambda0 <- penalty$c * quantile(sim, probs = 1 - penalty$gamma)
    Ups0 <- sqrt(var(y))
    lambda <- rep(lambda0 * Ups0, p)
  }
  
  
  ## == HETEROSKEDASTIC CASE == ##
  if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == 
      FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 * 
                                                                    p * 1))
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
    
    ## ======================  (addition) ============================= ##
    if(penalty$cluster.presence==TRUE){ 
      Ups0 <-NULL
      Ups0k<-0
      
      Kk<-ifelse(length(dim(y))==0,1,dim(y)[2]) # Number of endogeneous variables (in the IV Case), dependant variable size the linear model case (1).
      
      # x is n x p (number of penalized variables)
      # y is n x Kk (1 or number of endogeneous variables)
      for(k in 1:Kk){
        
        for(clus in unique(cluster_index)){
          
          # Residuals for the given cluster (sized n_c x Kk):
          if(length(dim(y))==0){y_c<-unlist(y[cluster_index==clus])}else{y_c<-unlist(y[cluster_index==clus,k])}
          
          # Residuals * explanatory (coefficient-by-coefficient) for the given cluster (sized n_c x p):
          xy_c<-x[cluster_index==clus,]*matrix(y_c,sum(cluster_index==clus),p) 
          
          # Intermediate for penalized-variable-specific penalty (size p):
          Ups0k <- Ups0k+ apply(xy_c,MARGIN=2,FUN=function(x) sum(x)^2)  
          
        }
        Ups0<-cbind(Ups0,sqrt(Ups0k))
      }
      
      Ups0<-1/sqrt(n)*Ups0
      
      lambda <- lambda0 * Ups0
    }
    ## ========================  (addition) ============================= ##
    
  }
  
  if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == 
      TRUE) {
    if (!exists("numSim", where = penalty)) {
      penalty$numSim = 5000
    }
    p <- dim(x)[2]
    n <- dim(x)[1]
    R <- penalty$numSim
    sim <- vector("numeric", length = R)
    lasso.x.y <- rlasso(y ~ x)
    eh <- lasso.x.y$residuals
    ehat <- matrix(rep(eh, each = p), ncol = p, byrow = TRUE)
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each = p), ncol = p, byrow = TRUE)
      sim[l] <- n * max(2 * colMeans(x * ehat * g))
    }
    lambda0 <- penalty$c * quantile(sim, probs = 1 - penalty$gamma)
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }
  if (!is.null(penalty$lambda.start)) {
    p <- dim(x)[2]
    if (length(penalty$lambda.start) == 1) {
      lambda.start <- rep(penalty$lambda.start, p)
    }
    lambda <- as.matrix(penalty$lambda.start)
  }
  if (penalty$homoscedastic == "none") {
    if (is.null(penalty$lambda.start) | !exists("lambda.start", 
                                                where = penalty)) 
      stop("For method \"none\" lambda.start must be provided")
    n <- dim(x)[1]
    lambda0 <- penalty$lambda.start
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }
  return(list(lambda0 = lambda0, lambda = lambda, Ups0 = Ups0, 
              method = penalty))
}
