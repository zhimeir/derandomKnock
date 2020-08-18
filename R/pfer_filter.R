#' Stable variable selection with PFER control
#'
#' The function to implement the derandomized knockoffs procedure with PFER control. 
#'
#' @param X a n-by-p matrix of the covariates.
#' @param y the response vector of length in (can be continuous or binary).
#' @param v0 a positive numver indicating the PFER target (default: 1). Can be left NULL if using the kfwer error.
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param tau a number betweem 0 and 1 indicating the selection frequency (default: 0.5).
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
#' #Generate data
#' n <- 100; p <- 50; s <- 10;
#' rho <- 0.5;
#' Sigma <- toeplitz(rho^(1:p-1))
#' X <- matrix(rnorm(n*p),n,p)%*%chol(Sigma)
#' beta <- rep(0,p)
#' beta[1:s] <- 3.5/sqrt(n)
#' y <- X%*%beta+rnorm(n)
#' 
#' # Control PFER at level v=1
#' res <- pfer_filter(X,y,v0=1, knockoff_method = "gaussian",
#'                  knockoff_stat = stat.glmnet_coefdiff,
#'                  mu = rep(0,p),Sigma = Sigma)
#'
#'
#' @export

pfer_filter <- function(X,y, v0 = 1, M = 30, tau =0.5, knockoff_method = "gaussian",
                            knockoff_stat = stat.glmnet_coefdiff,seed = 24601,
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
  if(v0<0) stop("The PFER target v0 should be positive!")

  ## Check the knockoff-related input
  if(knockoff_method %in% c("gaussian","hmm") == 0) stop("The type of knockoffs is not supported!")

  ## Initialization
  set.seed(seed)
  pi = rep(0,p)
  if(knockoff_method == "gaussian"){
    diags = knockoff::create.solve_asdp(Sigma)
  }

  ## Multiple knockoff runs
  for (m in 1:M){
    v <- floor(v0)+rbinom(1,1,v0-floor(v0))
    if(knockoff_method == "gaussian"){
      Xk <- create.gaussian(X,mu,Sigma,diag_s = diags)}
    if(knockoff_method == "hmm"){
      Xk <- knockoffHMM(X, pInit, Q,pEmit,seed = seed+m)
    }
    W <- knockoff_stat(X,Xk,y)
    order_w <- order(abs(W),decreasing = TRUE)
    sorted_w <- W[order_w]
    negid <- which(sorted_w<0)
    if(v>0){
      if(sum(W<0)<v){
        S <- which(W>0)
        pi[S] <- pi[S]+1
      }else{
        TT <- negid[v]
        S <- which(sorted_w[1:TT]>0)
        S <- order_w[S]
        pi[S] <- pi[S]+1
        }
      }
  }
  pi <- pi/M
  S <- which(pi>=tau)

return(list(S=S,pi=pi,tau=tau,W=W))
}


