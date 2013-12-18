polyCov <-
  function(x,y=NULL,thresholds=NULL,progress=TRUE,...)
  {
    if (!is.null(y))
    {
      polyCov(cbind(x,y), thresholds = thresholds)
    } else
    {
      # Load mvtnorm:
      library("mvtnorm")
      
      # Functions to optimize:
      fcov <- function(param,x,y,threshX,threshY,means,vars)
      {
        tab <- table(x,y)
        covmat <- matrix(1,2,2)
        covmat[1,2] <- covmat[2,1] <- param
        diag(covmat) <- vars
        tabExp <- matrix(1,nrow(tab),ncol(tab))
        for (i in 1:(length(threshX)-1))
        {
          for (j in 1:(length(threshY)-1))
          {
            tabExp[i,j] <- pmvnorm(c(threshX[i],threshY[j]),c(threshX[i+1],threshY[j+1]),means,sigma=covmat)
          }
        }
        return(sum(tab*log(tabExp)))
      }
      fmeanvar <- function(params,dat,thresh)
      {
        tab <- table(dat)
        exp <- numeric(0)
        for (i in 1:length(tab)) exp[i] <- pnorm(thresh[i+1],params[1],params[2]) - pnorm(thresh[i],params[1],params[2])
        return(sum(tab*log(exp)))
      }
      fmean <- function(params,dat,thresh)
      {
        tab <- table(dat)
        exp <- numeric(0)
        for (i in 1:length(tab)) exp[i] <- pnorm(thresh[i+1],params[1],1) - pnorm(thresh[i],params[1],1)
        return(sum(tab*log(exp)))
      }
      
      x <- as.data.frame(x)
      n <- ncol(x)
      
      if (is.null(thresholds))
      {
        newTresh <- TRUE
        thresholds <- lapply(lapply(x,table),function(y)c(-Inf,sapply(cumsum(y/sum(y)),qnorm)))
        means <- rep(0,n)
        vars <- rep(1,n)
      } else 
      {
        newTresh <- FALSE
        thresholds <- lapply(thresholds,function(x)c(-Inf,x,Inf))
      }
      # Number of Thresholds:
      Nthresh <- sapply(thresholds,function(x)length(x)-2)
      
      # Estimate means and vars:
      if (!newTresh)
      {
        means <- numeric(n)
        vars <- numeric(n)
        for (i in 1:n)
        {
          if (Nthresh[i]==1)
          {
            est <- suppressWarnings(optim(0,fmean,dat=x[,i],thresh=thresholds[[i]],control=list(fnscale=-1))$par )
            means[i] <- est
            vars[i] <- 1
          } else if (Nthresh[i] > 1)
          {
            est <- optim(c(0,1),fmeanvar,dat=x[,i],thresh=thresholds[[i]],method="L-BFGS-B",lower=c(-Inf,0.1),control=list(fnscale=-1))$par
            means[i] <- est[1]
            vars[i] <- est[2]^2
          } else stop ("Number of thresholds not correct")
        }
      }
      
      # Create Covariance matrix:
      covmat <- matrix(,n,n)
      diag(covmat) <- vars
      if (progress)
      {
        pb <- txtProgressBar(0,n*(n-1)/2,title="Computing polychoric covariances:",style=3)
        c <- 0
      }
      for (i in 1:n)
      {
        for (j in i:n)
        {
          if (j>i)
          {
            covmat[i,j] <- optimize(fcov,interval=c(-sqrt(vars[i]*vars[j]),sqrt(vars[i]*vars[j])),x=x[,i],y=x[,j],threshX=thresholds[[i]],threshY=thresholds[[j]],means=means[c(i,j)],vars=vars[c(i,j)],maximum=TRUE,...)$maximum
            if (progress)
            {
              c <- c + 1
              setTxtProgressBar(pb, c)
            }
          }
        }
      }
      if (progress)
      {
        close(pb)
      }
      thresholds <- lapply(thresholds,function(x)x[-c(1,length(x))])
      covmat[lower.tri(covmat,diag=FALSE)] <- t(covmat)[lower.tri(covmat,diag=FALSE)]
      output <- list(covmat=covmat,thresholds=thresholds,means=means,vars=vars)
      class(output) <- "polyCov"
      return(output)
    }
  }

