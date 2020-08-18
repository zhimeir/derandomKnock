#' A filter for PFER control
#'
#' The function to implement the derandomized knockoffs procedure with PFER control. 
#'
#' @param X a n-by-p matrix of the covariates.
#' @param y the response vector of length in (can be continuous or binary).
#' @param type the type of error to control. Options include "pfer" and "kfwer".
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param tau a number betweem 0 and 1 indicating the selection frequency (default: 0.5).
#' @param seed an integer specifying the random seed used in the procedure.
#' @param k a positive integer corresponding to k-FWER.
#' @param alpha a number between 0 and 1 indicating the target k-FWER level.
#' @param knockoff_method either "gaussian" or "hmm" (default: "gaussian").
#' @param knockoff_stat An function in \code{\link[knockoff]{statistics}}
#' @param mu a length-p mean vector of X if it follows a Gaussian distribution.
#' @param Sigma a p-by-p covariance matrix of X if it follows a Gaussian distribution.
#' @param pInit n array of length K, containing the marginal distribution of the states for the first variable, if X is sampled from an HMM.
#' @param Q an array of size (p-1,K,K), containing a list of p-1 transition matrices between the K states of the Markov chain, if X is sampled from an HMM.
#' @param pEmit an array of size (p,M,K), containing the emission probabilities for each of the M possible emission states, from each of the K hidden states and the p variables, if X is sampled from an HMM.
#'
#' @examples
#' #Generate data
#' n <- 100; p <- 50; s <- 10;
#' rho <- 0.5;
#' Sigma <- toeplitz(rho^(1:p-1))
#'
#' @export

pfer_filter <- function(X,y, v0 = 1, M = 50, tau =0.5, knockoff_method = "gaussian",
                            knockoff_stat = stat.glmnet_coefdiff,seed = 24601,
                            mu = NULL,Sigma =NULL,#parameter for gaussian knockoff
                            pInit = NULL, Q = NULL,pEmit = NULL #parameter for hmm knockoff
                            ){
#check input
#Initialization
set.seed(seed)
n = dim(X)[1]
p = dim(X)[2]
pi = rep(0,p)
lambda = rep(0,p)
if(knockoff_method == "gaussian"){
  diags = knockoff::create.solve_asdp(Sigma)
}

#Main part
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
        S = which(W>0)
        pi[S] = pi[S]+1
      }else{
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


