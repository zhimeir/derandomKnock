#' Stable knockoff(fdx)
#' internal function to determine v
#' 
#' @export

#Auxiliary functions
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
  v =max(which(vfun<=alpha))*accuracy
  return(v)
}

root_2 <- function(k,alpha,accuracy = 0.05,length,tau=0.5,h2 = h2){
  vfun = rep(0,length)
  for (i in 1:length){
    v = accuracy*i
    vfun[i] =(v^2+2*v)/tau^2/k/(k+2*h2)
  }
  v =max(which(vfun<=alpha))*accuracy
  return(v)
}
