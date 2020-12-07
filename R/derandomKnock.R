#' A stable variable selection procedure based on the knockoffs
#' 
#' The main function that implements the derandomized knockoffs procedure.
#'
#' @param X a n-by-p matrix of the covariates.
#' @param y the response vector of length in (can be continuous or binary).
#' @param type the type of error to control. Options include "pfer" and "kfwer".
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param eta a number betweem 0 and 1 indicating the selection frequency (default: 0.5).
#' @param seed an integer specifying the random seed used in the procedure.
#' @param v a positive numver indicating the PFER target (default: 1). Can be left NULL if using the kfwer error.
#' @param k a positive integer corresponding to k-FWER.
#' @param alpha a number between 0 and 1 indicating the target k-FWER level.
#' @param knockoff_method either "gaussian" or "hmm" (default: "gaussian").
#' @param knockoff_stat Knockoff statistics to used. See \code{\link[knockoff]{knockoff_stat}}(default: "stat.glmnet_diff")
#' @param mu a length-p mean vector of X if it follows a Gaussian distribution.
#' @param Sigma a p-by-p covariance matrix of X if it follows a Gaussian distribution.
#' @param diags a length-p vector, containing the precomputed covariances between the original variables and the knockoffs.
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
#'  beta[1:s] <- 5/sqrt(n)
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

derandomKnock <- function(X,y,type,
                          M = NULL, eta = NULL, seed = 24601,
                          v = 1, k = 1, alpha = 0.05,
                          knockoff_method = "gaussian",
                          knockoff_stat = stat.glmnet_coefdiff,
                          mu = NULL,Sigma =NULL,diags = NULL,
                          pInit = NULL, Q = NULL,pEmit = NULL,
                          thres = 50, beta = 1
                          ){

  ## Check the dependencies
  list_of_packages <- c("knockoff","glmnet","SNPknock","doMC","CVXR")
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=list_of_packages,FUN=require,character.only=TRUE))

  ## Check the dimension
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(length(y)!=n) stop("The dimensions of X and y are not compatible!")

  ## Check the type of error criterion
  if(type %in% c("pfer","kfwer") == 0) stop("The error criterion is not supported!")

  ## Check the parameters
  if(!is.null(M)){
    if(M<=0) stop("The number of knockoff runs should be positive!")
  }
  if(!is.null(eta)){
    if(eta>1 | eta<0) stop("The selection probability should be between 0 and 1!")
  }

  ## Check the knockoff-related input
  if(knockoff_method %in% c("gaussian","hmm") == 0) stop("The type of knockoffs is not supported!")

  ## Check the error criterion target 
  if(type == "pfer"){
    if(v<0) stop("The PFER target v should be positive!")
  }

  if(type == "kfwer"){
    if(round(k)!=k | k<=0) stop("k should be a positive natural number!")
    if(alpha>1 | alpha<0) stop("alpha should be between 0 and 1!")  
  }

  ## Determine the type of the response
  if(length(unique(y))==2){
    response <- "binary"
  }else{
    response <- "continuous"
  }

  ## Determining knockoff statistics (if not specified by the input)
    ## Continuous response
    if(response == "continuous"){
      if(is.null(knockoff_stat)){
        knockoff_stat  <-  stat.glmnet_coefdiff
      }
    }
    ## Binary response
    if(response == "binary"){
      if(is.null(knockoff_stat)){
        knockoff_stat  <-  stat.lasso_coefdiff_bin
      }
    }

  ## Determining parameters
  params <- get_params(type = type, v = v,
                       k = k, alpha = alpha, 
                       M = M, eta = eta,
                       thres = thres, beta = beta)  
  M <- params$M
  eta <- params$eta
  v0 <- params$v0

  ## Run the derandomized knockoffs procedure
  res <- base_filter(X = X,y = y,v0 = v0,
                     M = M, tau = eta, 
                     knockoff_method = knockoff_method,
                     knockoff_stat = knockoff_stat,
                     seed = seed, mu = mu, Sigma = Sigma,
                     pInit = pInit, Q =Q, pEmit = pEmit)

  return(list(S = res$S, frequency = res$pi[res$S]))
}


