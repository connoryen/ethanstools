#' Compare the behavior of multiple arima models to a reference model's residuals
#' 
#' Fits ARIMA(p,d,q) models to the data and plots their ACF/PACF, 
#' Ljung-Box Statistic p-values, and spectral density. Also plots for comparison
#' the residuals of a reference ARIMA model. 
#' 
#' @param x data
#' @param model List of vectors of ARIMA(p,d,q) model orders. 
#' @param ref Order of reference ARIMA model whose residuals will be used for comparison
#' @param max.lag number of lags to compute in the ACF/PACF and 
#' Ljung-Box Statistic plots. 
#' @param plot plots diagnostic
#' 
#' @return For the fitted ARIMA models, list of theoretical ACF/PACF values, 
#' theoretical spectral density values, and Ljung-Box Statistic p-values.
sarima.compare <- function(x, models = c(p=0,d=0,q=1), ref = c(p=0,d=0,q=0),
                           max.lag = 15, plot = TRUE, 
                           model.colors = c("dodgerblue", "hotpink")) {
  # CHECK INPUT:
  # ----------------------------------------------------------------------------
  stopifnot("`x` must be a numeric vector" = 
              (is.numeric(x)))
  stopifnot("`models` must be a list or numeric" = 
              (is.list(models) | is.numeric(models)))
  stopifnot("Reference ARIMA order vector (`ref`) must be numeric" = 
              is.numeric(ref))
  stopifnot("Reference ARIMA order vector (`ref`) must have length of 3 or 7" = 
              ((length(ref) == 3) | (length(ref) == 7)))
  
  p_ref = ref[1]; d_ref = ref[2]; q_ref = ref[3]
  # ref's seasonal parameters
  if(length(ref) == 7){
    P_ref = ref[4]; D_ref = ref[5]; Q_ref = ref[6]; S_ref = ref[7]
  } else {
    P_ref = 0; D_ref = 0; Q_ref = 0; S_ref = -1
  }
  # PASS TO `sarima.compare1`:
  # ----------------------------------------------------------------------------
  # only 1 arima model -- use default sarima.compare1
  if (is.numeric(models)) {  # 1 ARIMA model input as a vector
    return(sarima.compare1(x=x, model=models, ref=ref, max.lag=max.lag, plot=plot))
  } else if (length(models) == 1){  # 1 ARIMA model input as a list
    return(sarima.compare1(x=x, model=models[[1]], ref=ref, max.lag=max.lag, plot=plot))
  }
  
  my.models <- list()
  for (i in c(1:length(models))) {
    my.models[[i]] <- sarima.compare1(x=x, model=models[[i]], ref = ref,
                                      max.lag=max.lag, plot = FALSE)
  }
  # REFERENCE MODEL:
  # ----------------------------------------------------------------------------
  m_ref <- astsa::sarima(x, p=p_ref, d=d_ref, q=q_ref, 
                         P=P_ref, D=D_ref, Q=Q_ref, S=S_ref,
                         details=FALSE, model=FALSE)
  # ljung-box p-vals
  lb.pvals.ref <- c()
  for (i in c(1:max.lag)){
    lb.pvals.ref <- 
      c(lb.pvals.ref, stats::Box.test(m_ref$fit$residuals, lag = i, type = "Ljung-Box")$p.value)
  }
  # spectral density
  spec.ref <- astsa::mvspec(m_ref$fit$residuals, plot = FALSE)
  # PLOT
  # ----------------------------------------------------------------------------
  if (plot) {
    # get colors:
    if(length(model.colors)==length(models)) {my.colors = model.colors} else 
      {my.colors = rainbow(length(models))}
    par(mfrow=c(2,2), mar=c(2.5,2,2,2)+0.1)
    # ACF Plot:
    stats::acf(m_ref$fit$residuals, lag.max = max.lag, ylim = c(-0.75,1), xaxt = "n",
               main = '', xlab = '', ylab = '')
    for (i in c(1:length(models))) {
      points(c(0:max.lag), my.models[[i]]$acf.theoretical, col = my.colors[i], lwd = 2)
    }
    axis(side = 1, at = seq(0,max.lag, by=2))
    mtext(substitute(bold("ACF")), side=3)
    mtext(substitute(italic("Lag")), side=1, line = 1.5, cex = 0.85)
    # PACF Plot:
    stats::pacf(m_ref$fit$residuals, lag.max = max.lag, ylim = c(-0.75,1), xaxt = "n",
                main = '', xlab = '', ylab = '')
    for (i in c(1:length(models))) {
      points(c(1:max.lag), my.models[[i]]$pacf.theoretical, col = my.colors[i], lwd = 2)
    }
    axis(side = 1, at = seq(2,max.lag, by=2))
    mtext(substitute(bold("PACF")), side=3)
    mtext(substitute(italic("Lag")), side=1, line = 1.5, cex = 0.85)
    # Ljung-Box:
    plot(c(3:max.lag), lb.pvals.ref[3:length(lb.pvals.ref)], 
         ylim = c(0,1), lwd = 2, main = '', xlab = '', ylab = '')
    for (i in c(1:length(models))) {
      points(c(3:max.lag), my.models[[i]]$lb.pvals.model[3:length(my.models[[i]]$lb.pvals.model)], 
             ylim = c(0,1), col = my.colors[i], lwd = 2)
    }
    abline(h = 0.05, lty = "dashed", col = "blue")
    mtext(substitute(bold("Ljung-Box Statistic p-value")), side=3)
    mtext(substitute(italic("Lag")), side=1, line = 1.5, cex = 0.85)
    # Spectral Density:
    plot(spec.ref$freq, spec.ref$spec, type = "l", yaxt = "n", 
         main = '', xlab = '', ylab = '', lwd = 2)
    for (i in c(1:length(models))) {
      par(new = TRUE)
      plot(my.models[[i]]$spec.theoretical$freq, my.models[[i]]$spec.theoretical$spec, 
           xaxt = "n", yaxt = "n", ylab = '', xlab = , 
           main = "", type = "l", lwd = 2, col = my.colors[i])
    }
    mtext(substitute(bold("Spectral Density")), side=3)
    mtext(substitute(italic("Frequency")), side=1, line = 1.5, cex = 0.85)
    # reset par
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  }
  # name `my.models`
  my.names <- c()
  for (i in c(1:length(models))) {
    my.names <- c(my.names,paste0("m",i))
  }
  names(my.models) <- my.names
  return(invisible(my.models))
}