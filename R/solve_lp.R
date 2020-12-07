#' Determine the parameter for PFER/k-FWER control
#'
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param eta a number betweem 0 and 1 indicating the selection frequency (default: 0.5).
#' @param beta a number between 0 and 1 specifying the vanishing speed of the density.
#'
#' @export


get_lp_bnd <- function(M,eta,beta = 1){
  one <- rep(1,M+1)
  one_eta <- rep(0,M+1)
  supp <- (0:M)/M
  one_eta[supp>=eta] <- 1
  
  y <- Variable(M+1)
  obj <- Maximize(t(one_eta)%*% y)
  constraints <- list(
    y>=0,
    t(supp)%*%y==1,
    y[2:(M+1)]-y[1:M] <= (beta-1)*y[-(M+1)]
  )
  prob <- Problem(obj, constraints)
  res <- psolve(prob,solver = "ECOS")
  y <- res$getValue(y)
  opt_ratio <- one_eta%*%y
  return(opt_ratio)
}

