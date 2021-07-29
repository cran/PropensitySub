#' Calculate standardized difference
#' 
#' @param vars variables of interest. standardized difference of each variable \cr
#' listed here will be calculated.
#' @param data0 A \code{data.frame} which include \code{vars} as columns from reference arm. 
#' All data are expected to be numerical. If a column is
#' not numerical, it will be turned to numerical by \code{\link{model.matrix}}.
#' @param weight0 weights for each record in reference arm.
#' @param data1 A \code{data.frame} which include \code{vars} as columns from comparison arm. 
#' All data are expected to be numerical. If a column is
#' not numerical, it will be turned to numerical by \code{\link{model.matrix}}.
#' @param weight1 weights for each record in comparison arm.
#' 
#' @return return a `numeric` vector for standardized difference of each variable
#' 
#' @examples 
#' library(dplyr)
#' data0 <- clinical %>% filter(ARM == "experimental")
#' data1 <- clinical %>% filter(ARM == "control")
#' calc_std_diff(
#'   vars = c("BECOG", "SEX"),
#'   data0 = data0,
#'   weight0 = rep(1, nrow(data0)),
#'   data1 = data1,
#'   weight1 = rep(1, nrow(data1))
#'  )
#' @importFrom stats model.matrix
#' @export
#' @note Calculation from Austin and Stuart (2015)
#' 


calc_std_diff <- function(vars, data0, weight0, data1, weight1) {
  stopifnot(all(vars%in%colnames(data0)))
  stopifnot(all(vars%in%colnames(data1)))
  stopifnot(length(weight0)==nrow(data0))
  stopifnot(length(weight1)==nrow(data1))
  data.cb <- rbind(data0, data1)
  data.cb <- data.cb[,vars]
  colnames(data.cb) <- paste0(colnames(data.cb),".")
  for (i in colnames(data.cb)){
    if (class(data.cb[[i]]) %in% c("factor", "character")){
      if (length(unique(data.cb[[i]])) > 2){
        lvls <- as.character(unique(data.cb[[i]]))
        stopifnot(!("tmp" %in% lvls))
        data.cb[[i]] <- factor(data.cb[[i]], levels = c("tmp", lvls))
      }
    }
  }
  
  # for da
  num.unique <- apply(data.cb, 2, function(i)length(unique(i)))
  data.cb.1lev <- data.cb[,which(num.unique==1), drop = FALSE ] # model.matrix cannot hanlde variable with only 1 unique value
  data.cb.morelev <- data.cb[,which(num.unique>1), drop = FALSE]
  
    
  data.cb.morelev <- model.matrix(~., data.cb.morelev)[ , -1, drop = FALSE]
  vars2 <- colnames(data.cb.morelev)
  data.cb <- cbind(data.cb.morelev, data.cb.1lev)
  data0.use <- data.cb[1:nrow(data0), , drop = FALSE]
  data1.use <- data.cb[(nrow(data0)+1): nrow(data.cb), , drop = FALSE]
  
  mean.std.out <- rep(NA, length(vars2) )
  names(mean.std.out) <- vars2
  for (i in vars2){
    v0 <- data0.use[,i, drop = FALSE]
    v1 <- data1.use[,i, drop = FALSE]
    Z.wmean0 <- sum(v0*weight0)/sum(weight0)
    Z.wmean1 <- sum(v1*weight1)/sum(weight1)
    Z.ws0 <- sum(weight0)/(sum(weight0)^2 -sum(weight0^2))*sum(weight0*(v0-Z.wmean0)^2)
    Z.ws1 <- sum(weight1)/(sum(weight1)^2 -sum(weight1^2))*sum(weight1*(v1-Z.wmean1)^2)
    Z.wmd <- (Z.wmean1-Z.wmean0)/sqrt((Z.ws0+Z.ws1)/2)
    # if all values in trt and ctrl are the same, diff is 0
    if(Z.wmean0 == Z.wmean1 & length(unique(c(v0,v1)))==1) Z.wmd <- 0
    # if all values in trt are one level and all values in ctrl are another values, diff is Inf
    if(!is.na(Z.ws0) & Z.ws0 == 0 & !is.na(Z.ws1) & Z.ws1 == 0 & length(unique(c(v0,v1)))>1) Z.wmd <- Inf
    nz.0 <- weight0[which(v0!=0)]
    nz.1 <- weight1[which(v1!=0)]
    if(all(nz.0==0) & all(nz.1==0)) Z.wmd <- 0     
    
    mean.std.out[i] <- Z.wmd
  }
  mean.std.out
}
