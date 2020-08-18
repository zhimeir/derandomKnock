#' Stable knockoff(fdx)
#' 
#'Internal function to generate the output
#'
#' @export

getRejectionSet <- function(W,v){
  order = order(abs(W),decreasing = TRUE)
  ordw = W[order]
  negid = which(ordw<0)
  if(v>0){
    if(sum(ordw<0) == 0){
      S = which(W>0)
      TT = length(S)
    }else{
      TT = negid[min(v,length(negid))]-1
      S = which(ordw[1:TT]>0)
      S = order[S]
    }
  }else{S = NULL;TT = 0}
  return(list(S=S))
}

