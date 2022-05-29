#' Determining the parameter for the base filter
#'
#' A function to determine the parameters for the derandomized knockoffs procedure 
#'
#' @param type the type of error to control. Options include "pfer" and "kfwer".
#' @param v a positive numver indicating the PFER target (default: 1). Can be left NULL if using the kfwer error.
#' @param k a positive integer corresponding to k-FWER.
#' @param alpha a number between 0 and 1 indicating the target k-FWER level.
#' @param M an integer specifying the number of knockoff copies computed (default: 30).
#' @param eta a number betweem 0 and 1 indicating the selection frequency (default: 0.5).
#'
#' @export

get_params <- function(type, v = NULL,
                       k = NULL, alpha = NULL,
                       M = NULL,eta = NULL,
                       thres = 50, beta = 1){
  tol <- 1e-6
  load("sysdata.rda")
  ## PFER control
  M_eta_mat <- M_eta_mat[M_list <= thres,]
  if(type == "pfer"){
    if(is.null(M)){
      if(is.null(eta)){
        M <- 31
        eta <- 0.5
        v0 <- v
      }else{
        v0 <- v
        eta <- ceiling(eta*100)/100
        eta.ind <- which(eta_list == eta)
        M.ind <- suppressWarnings(max(which(M_eta_mat[,eta.ind] <= 1 + tol)))
        if(M.ind == -Inf){stop("No suitable M found!")}
        M <- M_list[M.ind]
      }
    }else{
      if(is.null(eta)){
        v0 <- v
        M.ind <- which(M_list == M)
        eta.ind <- suppressWarnings(min(which(M_eta_mat[M.ind,]<=1 + tol)))
        if(eta.ind == Inf){stop("No suitable eta found!")}
        eta <- eta_list[eta.ind]
      }else{
        ratio <- get_lp_bnd(M,eta,beta)
        v0 <- v/ratio
      }
    }
  }
    
  ## k-FWER control
  if(type == "kfwer"){
    if(is.null(M)){
      if(is.null(eta)){
        if(k==1){
          v0 <- 1
          res <- search_M_eta(M_eta_mat,alpha,M_list,eta_list)
          M <- res$M
          eta <- res$eta
        }else{
          if(4*k*alpha<=1){
            v0 <- 1
            res <- search_M_eta(M_eta_mat,2*k*alpha,M_list,eta_list)
            M <- res$M
            eta <- res$eta
          }else{
            M <- 31
            eta <- 0.5
            v0 <- 2*k*alpha
          }
        }
      }else{
        stop("The value of M cannot be found. Please either input both M and eta or leave both of them NULL.")
      }
    }else{
      if(is.null(eta)){
        stop("The value of M cannot be found. Please either input both M and eta or leave both of them NULL.")
      }else{
        ratio <- get_lp_bnd(M,eta,beta)
        if(k == 1){
          v0 <- alpha/ratio
        }else{
          v0 <- 2*k*alpha/ratio
        }
        if(v0<0){stop("There is no suitable parameter.")}
      }
    }
  }

  return(list(M = M, eta = eta, v0 = v0))
}

search_M_eta <- function(M_eta_mat,ratio,M_list,eta_list){
  len_M <- dim(M_eta_mat)[1]
  len_eta <- dim(M_eta_mat)[2]
  for(i in 1:len_eta){
    ind <- suppressWarnings(max(which(M_eta_mat[,i]<=ratio)))
    if(ind!=-Inf){
      M <- M_list[ind]
      eta <- eta_list[i]
      break
    }    
  }
  if(ind == -Inf){stop("No suitable parameters found.")}
  return(list(M = M, eta = eta))
}

