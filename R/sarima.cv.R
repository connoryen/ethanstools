#' time series cross validation on multiple SARIMA models
#' 
#' @param x data
#' 
#' @importFrom astsa sarima.for
#' 
#' @export
sarima.cv <- function(x, models = c(p=0,d=0,q=1), 
                      n.frames = 10, n.points = 5) {
  # CHECK INPUT:
  # ----------------------------------------------------------------------------
  stopifnot("`x` must be a numeric vector" = 
              (is.numeric(x)))
  stopifnot("`models` must be a list or numeric" = 
              (is.list(models) | is.numeric(models)))
  stopifnot("the number of testing points is too large" = 
              (length(x) > n.frames+n.points))
  for(i in c(1:length(models))){
    stopifnot("models must be of length 3 or 7" = 
                (length(models[[i]]) %in% c(3,7)))
  }
  # HELPER FUNCTION:
  # ----------------------------------------------------------------------------
  # astsa::sarima.for wrapper to accommodate vector order input. 
  # 
  # Pass a SARIMA(p,d,q)x(P,D,Q)[S] model specified by a vector with 3 
  # elements: c(p,d,q), or a vector with 7 elements: c(p,d,q, P,D,Q,S) to 
  # astsa::sarima.for.
  # 
  # @param x time series data
  # @param ord vector of SARIMA orders
  # 
  # @return vector of predictions
  #
  to.sarima.for <- function(x, ord, n.ahead) {
    if (length(ord) == 3){  # ARIMA(p,d,q)
      return(astsa::sarima.for(xdata = x, n.ahead = n.ahead, plot = FALSE,
                               p = ord[1], d = ord[2], q = ord[3])$pred)
    } else if (length(models[[i]]) == 7){  # SARIMA(p,d,q)x(P,D,Q)[S]
      return(astsa::sarima.for(xdata = x, n.ahead = n.ahead, plot = FALSE,
                               p = ord[1], d = ord[2], q = ord[3], 
                               P = ord[4], D = ord[5], Q = ord[6],  S = ord[7])$pred)
    }
  }
  # PARSE MODELS:
  # ----------------------------------------------------------------------------
  my.models <- list()
  for (i in c(1:length(models))) {
    models[[i]]
  }
  # CROSS VALIDATE:
  # ----------------------------------------------------------------------------
  RMSE <- rep(0, length(models))
  for (f in c(0:(n.frames-1))){
    train <- x[1:(length(x)-n.points-f)]
    test <- x[(length(x)-n.points+1-f):(length(x)-f)]
    for (m in c(1:length(models))) {
      RMSE[m] <- RMSE[m] + sum(abs(to.sarima.for(train, models[[m]], n.ahead=n.points)-test)**2)
    }
  }
  RMSE <- base::sqrt(RMSE/(n.frames*n.points))
  names(RMSE) <- paste0("m", c(1:length(models)))
  return(RMSE)
}