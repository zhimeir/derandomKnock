#' A stable variable selection procedure based on the knockoffs
#' 
#' The main function that implements the derandomized knockoffs procedure.
#'
#' @param X a n-by-p matrix of the covariates.
#' @param y the response vector of length in (can be continuous or binary).
#' @param type the type of error to control. Options include "pfer" and "kfwer".
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param tau a number betweem 0 and 1 indicating the selection frequency (default: 0.5).
#' @param seed an integer specifying the random seed used in the procedure.
#' @param v a positive numver indicating the PFER target (default: 1). Can be left NULL if using the kfwer error.
#' @param k a positive integer corresponding to k-FWER.
#' @param alpha a number between 0 and 1 indicating the target k-FWER level.
#' @param knockoff_method either "gaussian" or "hmm" (default: "gaussian").
#' @param mu a length-p mean vector of X if it follows a Gaussian distribution.
#' @param Sigma a p-by-p covariance matrix of X if it follows a Gaussian distribution.
#' @param pInit n array of length K, containing the marginal distribution of the states for the first variable, if X is sampled from an HMM.
#' @param Q an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K states of the Markov chain, if X is sampled from an HMM.
#' @param pEmit an array of size (p,M,K), containing the emission probabilities for each of the M possible emission states, from each of the K hidden states and the p variables, if X is sampled from an HMM.
#' @return S the selection set.
#' @return frequency the selection frequency of the selected variables. 
#'
#' @examples
#'  #Generate data
#'  n <- 100; p <- 50; s <- 10;
#'  rho <- 0.5;
#'  Sigma <- toeplitz(rho^(1:p-1))
#'  X <- matrix(rnorm(n*p),n,p)%*%chol(Sigma)
#'  beta <- rep(0,p)
#'  beta[1:s] <- 3.5/sqrt(n)
#'  y <- X%*%beta+rnorm(n)
#' 
#' # Control PFER at level v=1
#' res <- derandomKnock(X,y,type = "pfer",v=1, knockoff_method = "gaussian",
#'                    mu = rep(0,p),Sigma = Sigma)
#'
#' # Control 1-FWER at level alpha=0.1
#' res <- derandomKnock(X,y,type = "kfwer", k=1, alpha = 0.1, knockoff_method = "gaussian",
#'                    mu = rep(0,p),Sigma = Sigma)
#'
#' @export

derandomKnock <- function(X,y,type = "pfer",
                          M=30, tau = 0.5, seed = 24601,
                          v = 1, k=1, alpha= 0.05,
                          knockoff_method = "gaussian",
                          mu = NULL,Sigma =NULL,#parameter for gaussian knockoff
                          pInit = NULL, Q = NULL,pEmit = NULL #parameter for hmm knockoff 
                          ){

  ## Check the dependencies
  list_of_packages <- c("knockoff","glmnet","SNPknock","doMC")
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=list_of_packages,FUN=require,character.only=TRUE))

  ## Check the dimension
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(length(y)!=n) stop("The dimension of X and y does not match!")

  ## Check the type of error criterion
  if(type %in% c("pfer","kfwer") == 0) stop("The error criterion is not supported!")

  ## Check the parameters
  if(M<=0) stop("The number of knockoff runs should be positive!")
  if(tau>1 | tau<0) stop("The selection probability should be between 0 and 1!")

  ## Check the knockoff-related input
  if(knockoff_method %in% c("gaussian","hmm") == 0) stop("The type of knockoffs is not supported!")
  
  ## Determine the type of the response
  if(length(unique(y))==2){
    response <- "binary"
  }else{
    response <- "continuous"
  }

  ## Run the derandomized knockoffs procedure
  ## Controlling the PFER 
  if(type == "pfer"){
    if(v<0) stop("The PFER target v should be positive!")
    ## Continuous response
    if(response == "continuous"){
      knockoff_stat  <-  stat.glmnet_coefdiff
      res <- pfer_filter(X,y,v0 = v,M,tau,knockoff_method,
                         knockoff_stat,seed,mu,Sigma,pInit,Q,pEmit)
    }
    ## Binary response
    if(response == "binary"){
      knockoff_stat  <-  stat.lasso_coefdiff_bin
      res <- pfer_filter(X,y,v0 = v,M,tau,knockoff_method,
                         knockoff_stat,seed,mu,Sigma,pInit,Q,pEmit)
    }
  }

  ## Controlling the k-FWER
  if(type == "kfwer"){
    if(round(k)!=k | k<=0) stop("k should be a positive natural number!")
    if(alpha>1 | alpha<0) stop("alpha should be between 0 and 1!")
    
    ## Continuous response
    if(response == "continuous"){
      knockoff_stat  <-  stat.glmnet_coefdiff
      res <- fwer_filter(X,y,k,alpha,M,tau,knockoff_method,
                         knockoff_stat,seed,mu,Sigma,pInit,Q,pEmit)
    }
    ## Binary response
    if(response == "binary"){
      knockoff_stat  <-  stat.lasso_coefdiff_bin
      res <- fwer_filter(X,y,k,alpha,M,tau,knockoff_method,
                         knockoff_stat,seed,mu,Sigma,pInit,Q,pEmit)
    }
  }
  return(list(S = res$S, frequency = res$pi[res$S]))
}

