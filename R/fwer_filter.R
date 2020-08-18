#' Stable variable selectino with k-FWER control
#'
#' The function to implement the derandomized knockoffs procedure with k-FWER control.
#'
#' @param X a n-by-p matrix of the covariates.
#' @param y the response vector of length in (can be continuous or binary).
#' @param k a positive integer corresponding to k-FWER.
#' @param alpha a number between 0 and 1 indicating the target k-FWER level.
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param tau a number betweem 0 and 1 indicating the threshold of the selection frequency (default: 0.5).
#' @param seed an integer specifying the random seed used in the procedure.
#' @param knockoff_method either "gaussian" or "hmm" (default: "gaussian").
#' @param knockoff_stat a feature importance statistic in the \code{\link{knockoff}} pacakge (e.g., \code{\link{stat.glmnet_coefdiff}}).
#' @param mu a length-p mean vector of X if it follows a Gaussian distribution.
#' @param Sigma a p-by-p covariance matrix of X if it follows a Gaussian distribution.
#' @param pInit n array of length K, containing the marginal distribution of the states for the first variable, if X is sampled from an HMM.
#' @param Q an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K states of the Markov chain, if X is sampled from an HMM.
#' @param pEmit an array of size (p,M,K), containing the emission probabilities for each of the M possible emission states, from each of the K hidden states and the p variables, if X is sampled from an HMM.
#' @return S the selection set.
#' @return pi the selection frequency of all selected variables.
#' @return tau the selection threshold
#' @return W a M-by-p matrix of the computed knockoff feature importance statistics.
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
#' # Control 1-FWER at level alpha=0.1
#' res <- derandomKnock(X,y,type = "kfwer", k=1, alpha = 0.1, knockoff_method = "gaussian",
#'                    mu = rep(0,p),Sigma = Sigma)
#'
#' @export
fwer_filter<- function(X,y,k = 1,alpha=0.1,M = 50, tau =0.5, 
                       knockoff_method = "gaussian",
                       knockoff_stat = stat.glmnet_coefdiff,
                       seed = 24601,
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

  ## Check the parameters
  if(M<=0) stop("The number of knockoff runs should be positive!")
  if(tau>1 | tau<0) stop("The selection probability should be between 0 and 1!")
  if(round(k)!=k | k<=0) stop("k should be a positive natural number!")
  if(alpha>1 | alpha<0) stop("alpha should be between 0 and 1!")

  ## Check the knockoff-related input
  if(knockoff_method %in% c("gaussian","hmm") == 0) stop("The type of knockoffs is not supported!")
  
  #Initialization
  set.seed(seed)
  pi = rep(0,p)
  if(knockoff_method == "gaussian"){
    diags = knockoff::create.solve_asdp(Sigma)
  }

  ## Computing the parameter of the base procedure
  v0 = getV(k,alpha,accuracy = 0.01,tau,xi = 2*tau, nu = 1, h1 = k,h2 = 1/2*k)
  v0 = max(v0,tau)

  ## Multiple knockoff runs
  for (m in 1:M){
    v = floor(v0)+rbinom(1,1,v0-floor(v0))
    if(knockoff_method == "gaussian"){
      Xk = create.gaussian(X,mu,Sigma,diag_s = diags)}
    if(knockoff_method == "hmm"){
      Xk = knockoffHMM(X, pInit, Q,pEmit,seed = seed+m)
    }
    W = knockoff_stat(X,Xk,y)
    order_w = order(abs(W),decreasing = TRUE)
    sorted_w = W[order_w]
    negid = which(sorted_w<0)
    if(v>0){
      if(sum(W<0)<v){
        ## If the total number of negative W is less than v, select all positice W_j's
        S = which(W>0)
        pi[S] = pi[S]+1
      }else{
        ## Stop when seeing v negative W_j's
        TT = negid[v]
        S = which(sorted_w[1:TT]>0)
        S = order_w[S]
        pi[S] = pi[S]+1
      }
    }
  }
  pi = pi/M
  S = which(pi>=tau)

return(list(S=S,pi=pi,tau = tau,W = W))
}

## Computing the parameter of the base procedure
getV <- function(k,alpha,accuracy = 0.05,tau = 0.5,
                 h1,h2,nu,xi){
  #number of grids
  length = k*tau/accuracy
  
  #calculating roots
  vexp = root_exp(k = k,alpha=alpha,tau = tau,length = length,accuracy = accuracy,nu=nu)
  v1 = (k+h1)*(tau+xi)*alpha
  v2 = root_2(k = k,alpha=alpha,tau = tau,length = length,accuracy = accuracy,h2 = h2)
  v = max(vexp,v1,v2)
  return(v)
}


root_exp <- function(k,alpha,accuracy = 0.05,length,tau=0.5,nu){
  vfun = rep(0,length)
  for (i in 1:length){
    v = accuracy*i
    vfun[i] = (2*v/tau/(v/tau+k))^v*((2*k)/(v/tau+k))^(k*tau)
    vfun[i] = 1/vfun[i]/(1+nu)
  }
  v = suppressWarnings(max(which(vfun<=alpha)))*accuracy
  return(v)
}

root_2 <- function(k,alpha,accuracy = 0.05,length,tau=0.5,h2 = h2){
  vfun = rep(0,length)
  for (i in 1:length){
    v = accuracy*i
    vfun[i] =(v^2+2*v)/tau^2/k/(k+2*h2)
  }
  v =suppressWarnings(max(which(vfun<=alpha)))*accuracy
  return(v)
}
