#'  rlasso_unpenalized_controls.R    Function for rlasso estimation - with unpenalized coefficients and with a clustering option
#'  
#'  This function build over the function `hdm::rlasso` . Only the addition are commented.
#'  
#'  
#'  @param index_x_no_selection A p-sized vector storing the index where no selection is performed (after the rLasso estimation, they are added to the selection set for post Lasso estimation)
#'  @param clusters A n-sized vector with identifiers of clusters (to be passed as argument to `lambdaCalculationClusteredOption`).
#'  If not NULL, indicates to compute the penalty level with Clustered-Lasso
#'  

rlasso_unpenalized_controls<-function (x, y,index_x_no_selection,clusters, post = TRUE, intercept = TRUE, model = TRUE, 
                                       penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, 
                                                      lambda.start = NULL, c = 1.1, gamma = 0.1/log(n)), control = list(numIter = 15, 
                                                                                                                        tol = 10^-5, threshold = NULL), ...) 
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(colnames(x))) 
    colnames(x) <- paste("V", 1:p, sep = "")
  ind.names <- 1:p
  if (!exists("homoscedastic", where = penalty)) 
    penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty)) 
    penalty$X.dependent.lambda = "FALSE"
  if (!exists("gamma", where = penalty)) 
    penalty$gamma = 0.1/log(n)
  if (penalty$homoscedastic == "none" & !exists("lambda.start", 
                                                where = penalty)) 
    stop("lambda.start must be provided!")
  if (!exists("numIter", where = control)) {
    control$numIter = 15
  }
  if (!exists("tol", where = control)) {
    control$tol = 10^-5
  }
  if (post == FALSE & (!exists("c", where = penalty))) {
    penalty$c = 0.5
  }
  if (intercept) {
    meanx <- colMeans(x)
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- y - mu
  }
  else {
    meanx <- rep(0, p)
    mu <- 0
  }
  normx <- sqrt(apply(x, 2, var))
  ind <- rep(FALSE, p)
  XX <- crossprod(x)
  Xy <- crossprod(x, y)
  startingval <- init_values(x, y)$residuals
  
  ## === (addition) Starting penalty with Clustered Option ==== ##
  if(!is.null(clusters)){
    penalty$cluster.presence<-TRUE
    pen <- lambdaCalculationClusteredOption(penalty = penalty, y = startingval, 
                             x = x,cluster_index = clusters)
  }else{
    penalty$cluster.presence<-FALSE
    pen <- lambdaCalculationClusteredOption(penalty = penalty, y = startingval, 
                             x = x)
  }
  ##=============================================================== ##
  lambda <- pen$lambda
  Ups0 <- Ups1 <- pen$Ups0
  lambda0 <- pen$lambda0
  mm <- 1
  s0 <- sqrt(var(y))
  
  ## ============ Put the penalty at zero for chosen coefficients ============== ## 
  lambda[index_x_no_selection] <- 0
  ##=============================================================== ##
  
  while (mm <= control$numIter) {
    if (mm == 1 && post) {
      coefTemp <- LassoShooting.fit(x, y, lambda/2, XX = XX, 
                                    Xy = Xy)$coefficients
    }
    else {
      coefTemp <- LassoShooting.fit(x, y, lambda, XX = XX, 
                                    Xy = Xy)$coefficients
    }
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    
    ## =========== (Addition) Keep track of the selected variables ===== ##
    selected<-sum(ind1)
    #print(coefTemp[index_x_no_selection]) # Should be selected by construction as not penalized
    
    # a/ When 0 selected variables
    x1 <- as.matrix(x[, ind1, drop = FALSE])
    if (dim(x1)[2] == 0) {
      if (intercept) {
        intercept.value <- mean(y + mu)
        coef <- rep(0, p + 1)
        names(coef) <- c("intercept", colnames(x))
      }
      else {
        intercept.value <- mean(y)
        coef <- rep(0, p)
        names(coef) <- colnames(x)
      }
      est <- list(coefficients = coef, beta = rep(0, p), 
                  intercept = intercept.value, index = rep(FALSE, 
                                                           p), lambda = lambda, lambda0 = lambda0, loadings = Ups0, 
                  residuals = y - mean(y), sigma = var(y), iter = mm, 
                  call = match.call(), options = list(post = post, 
                                                      intercept = intercept, ind.scale = ind, control = control, 
                                                      mu = mu, meanx = meanx))
      if (model) {
        est$model <- x
      }
      else {
        est$model <- NULL
      }
      est$tss <- est$rss <- sum((y - mean(y))^2)
      est$dev <- y - mean(y)
      class(est) <- "rlasso"
      return(est)
    }
    
    # B/When at least one selected variables, iterate penalty computation based on lasso selection
    
    if (post) {
      reg <- lm(y ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y - x1 %*% coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- y - x1 %*% coefTemp[ind1]
    }
    s1 <- sqrt(var(e1))
    if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == 
        FALSE) {
      Ups1 <- s1 * normx
      lambda <- rep(pen$lambda0 * s1, p)
      lambda[index_x_no_selection] <- 0
    }
    if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == 
        TRUE) {
      Ups1 <- s1 * normx
      lambda <- rep(pen$lambda0 * s1, p)
      lambda[index_x_no_selection] <- 0
    }
    if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == 
        FALSE) {
      Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      lambda <- pen$lambda0 * Ups1
      lambda[index_x_no_selection] <- 0
      
      ## =========== (Addition) Next step penalty computation with clustered option ===== ##
      
      if(penalty$cluster.presence==TRUE){  
        Ups1 <-NULL
        Ups1k<-0
        for(k in 1:dim(e1)[2]){
          for(clus in unique(clusters)){
            # x is n x p (number of penalized variables)
            # e1 is n x kK
            xy_c<-x[clusters==clus,]*matrix(unlist(e1[clusters==clus,k]),sum(clusters==clus),p) #n_c x p
            Ups1k <- Ups1k+ apply(xy_c,MARGIN=2,FUN=function(x) sum(x)^2)
          }
          Ups1<-cbind(Ups1,sqrt(Ups1k))
        }
        Ups1<-1/sqrt(n)*Ups1
        
        lambda <- pen$lambda0 * Ups1
        lambda[index_x_no_selection] <- 0
      }
      ## ================================================================================== ##
      
    }
    if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == 
        TRUE) {
      lc <- lambdaCalculation(penalty, y = e1, x = x)
      Ups1 <- lc$Ups0
      lambda <- lc$lambda
      lambda[index_x_no_selection] <- 0
    }
    if (penalty$homoscedastic == "none") {
      if (is.null(penalty$lambda.start)) 
        stop("Argument lambda.start required!")
      Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      lambda <- pen$lambda0 * Ups1
      lambda[index_x_no_selection] <- 0
    }
    mm <- mm + 1
    if (abs(s0 - s1) < control$tol) {
      break
    }
    s0 <- s1
  }
  
  
  if (dim(x1)[2] == 0) {
    coefTemp = NULL
    ind1 <- rep(0, p)
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  coefTemp <- as.vector(as.vector(coefTemp))
  names(coefTemp) <- names(ind1) <- colnames(x)
  if (intercept) {
    if (is.null(mu)) 
      mu <- 0
    if (is.null(meanx)) 
      meanx <- rep(0, length(coefTemp))
    if (sum(ind) == 0) {
      intercept.value <- mu - sum(meanx * coefTemp)
    }
    else {
      intercept.value <- mu - sum(meanx * coefTemp)
    }
  }
  else {
    intercept.value <- NA
  }
  if (intercept) {
    beta <- c(intercept.value, coefTemp)
    names(beta)[1] <- "(Intercept)"
  }
  else {
    beta <- coefTemp
  }
  s1 <- sqrt(var(e1))
  est <- list(coefficients = beta, beta = coefTemp, intercept = intercept.value, 
              index = ind1, lambda = lambda, lambda0 = lambda0, loadings = Ups1, 
              residuals = as.vector(e1), sigma = s1, iter = mm, call = match.call(), 
              options = list(post = post, intercept = intercept, control = control, 
                             penalty = penalty, ind.scale = ind, mu = mu, meanx = meanx), 
              model = model,n_selected=selected)
  if (model) {
    x <- scale(x, -meanx, FALSE)
    est$model <- x
  }
  else {
    est$model <- NULL
  }
  est$tss <- sum((y - mean(y))^2)
  est$rss <- sum(est$residuals^2)
  est$dev <- y - mean(y)
  class(est) <- "rlasso"
  return(est)
}

init_values<-function (X, y, number = 5, intercept = TRUE) 
{
  suppressWarnings(corr <- abs(cor(y, X)))
  kx <- dim(X)[2]
  index <- order(corr, decreasing = T)[1:min(number, kx)]
  coefficients <- rep(0, kx)
  if (intercept == TRUE) {
    reg <- lm(y ~ X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)[-1]
  }
  else {
    reg <- lm(y ~ -1 + X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)
  }
  coefficients[is.na(coefficients)] <- 0
  res <- list(residuals = reg$residuals, coefficients = coefficients)
  return(res)
}
