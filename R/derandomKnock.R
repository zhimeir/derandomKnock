#' A stable variable procedure based on knockoffs
#' 
#' @export

derandomKnock <- function(){
  ## Check the dependencies
  list_of_packages <- c("knockff","glmnet","SNPknock")
  new.packages <- list_of_packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=list_of_packages,FUN=require,character.only=TRUE))


}

