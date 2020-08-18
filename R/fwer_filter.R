fwer_filter<- function(X,y,k = 1,alpha=0.1,M = 50, tau =0.5, 
                       knockoff_method = "gaussian",
                       knockoff_stat = stat.glmnet_coefdiff,
                       mu = NULL,Sigma =NULL,#parameter for gaussian knockoff
                       pInit = NULL, Q = NULL,pEmit = NULL, #parameter for hmm knockoff
                       seed = 24601){
#Initialization
set.seed(seed)
n = dim(X)[1]
p = dim(X)[2]
pi = rep(0,p)
if(knockoff_method == "gaussian"){
  diags = knockoff::create.solve_asdp(Sigma)
  Xk = create.gaussian(X,mu,Sigma,diag_s = diags)
}
if(knockoff_method == "hmm"){
  Xk = knockoffHMM(X, pInit, Q,pEmit,seed = seed)
}

v0 = getV(k,alpha,accuracy = 0.01,tau,xi = 2*tau, nu = 1, h1 = k,h2 = 1/2*k)
v0 = max(v0,tau)

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

