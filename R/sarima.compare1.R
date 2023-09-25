#' Diagnostic plots for a SARIMA model
#' 
#' Plots ACF/PACF, Ljung-Box Statistic p-values, and spectral density of the
#' ARIMA(p,d,q) or SARIMA(p,d,q, P,D,Q,S) process fitted to the inputted data. 
#' 
#' @param x data
#' @param model Vectors of ARIMA(p,d,q) or SARIMA(p,d,q, P,D,Q,S) model order. 
#' @param ref reference model. 
#' @param max.lag number of lags to compute in the ACF/PACF and 
#' Ljung-Box Statistic plots. 
#' @param plot plots diagnostic
#' 
#' @return list of fitted ACF, PACF, Ljung-Box Statistic p-values, and spectral
#' density of the fitted ARIMA model. 
#' 
#' @import dplyr
#' @importFrom astsa sarima
#' @importFrom astsa mvspec
#' @importFrom stats Box.test
#' @importFrom stats ARMAacf
#' @importFrom stats acf
#' @importFrom stats pacf
#' 
#' @export
sarima.compare1 <- function(x, model = c(p=0,d=0,q=1), ref = c(p=0,d=0,q=0),  
                            max.lag = 15, plot = TRUE) {
  # CHECK INPUT:
  # ----------------------------------------------------------------------------
  stopifnot("ARIMA order vector (`model`) must be numeric" = 
              is.numeric(model))
  stopifnot("Reference ARIMA order vector (`ref`) must be numeric" = 
              is.numeric(ref))
  stopifnot("ARIMA order vector (`model`) must have length of 3 or 7" = 
              ((length(model) == 3) | (length(model) == 7)))
  stopifnot("Reference ARIMA order vector (`ref`) must have length of 3 or 7" = 
              ((length(ref) == 3) | (length(ref) == 7)))
  stopifnot("`max.lag` must be numeric" = 
              is.numeric(max.lag))
  
  p_ref = ref[1]; d_ref = ref[2]; q_ref = ref[3]
  p = model[1]; d = model[2]; q = model[3]
  
  stopifnot("model cannot be white noise" = (p+q>0))
  
  # ref's seasonal parameters
  if(length(ref) == 7){
    P_ref = ref[4]; D_ref = ref[5]; Q_ref = ref[6]; S_ref = ref[7]
  } else {
    P_ref = 0; D_ref = 0; Q_ref = 0; S_ref = -1
  }
  # model's seasonal parameters
  if(length(model) == 7){
    P = model[4]; D = model[5]; Q = model[6]; S = model[7]
  } else {
    P = 0; D = 0; Q = 0; S = -1
  }
  
  if((d_ref != d) | (D_ref != D) | (S_ref != S)){
    warning("difference or seasonality values are not equal between model and ref")
  }
  
  # HELPER FUNCTION:
  # ----------------------------------------------------------------------------
  # Multiply two polynomials
  # 
  # @param a vector of coefficients for first polynomial
  # @param b vector of coefficients for second polynomial
  # 
  # @return vector of coefficients of product of inputs
  mult_poly <- function(a, b){
    m <- outer(a, b)
    return(as.vector(tapply(m, row(m) + col(m), sum)))
  }
  # FITTED SARIMA:
  # ----------------------------------------------------------------------------
  # reference:
  m_ref <- astsa::sarima(x, p=p_ref, d=d_ref, q=q_ref, 
                         P=P_ref, D=D_ref, Q=Q_ref, S=S_ref,
                         details=FALSE, model=FALSE)
  # model:
  m <- astsa::sarima(x, p=p, d=d, q=q, P=P, D=D, Q=Q, S=S, 
                     details=FALSE, model=FALSE)
  # AR & MA components:
  # ----------------------------------------------------------------------------
  coefs <- data.table(coef = m$fit$coef, param = names(m$fit$coef))
  
  # get ar parameters
  ars <- coefs %>%
    dplyr::filter(grepl("^ar", param)) %>%
    dplyr::arrange(param) %>%
    dplyr::select(coef)
  ars <- as.vector(ars$coef)
  # get sar parameters
  sars <- coefs %>%
    dplyr::filter(grepl("^sar", param)) %>%
    dplyr::arrange(param) %>%
    dplyr::select(coef)
  sars <- as.vector(sars$coef)
  # get ma parameters
  mas <- coefs %>%
    dplyr::filter(grepl("^ma", param)) %>%
    dplyr::arrange(param) %>%
    dplyr::select(coef)
  mas <- as.vector(mas$coef)
  # get sma parameters
  smas <- coefs %>%
    dplyr::filter(grepl("^sma", param)) %>%
    dplyr::arrange(param) %>%
    dplyr::select(coef)
  smas <- as.vector(smas$coef)
  
  # multiply ar
  sar_poly <- c(1)
  for (i in sars) {sar_poly <- c(sar_poly, rep(0, (S-1)), i)}
  ar_expanded <- mult_poly(c(1, ars), sar_poly)
  if (length(ar_expanded) == 1){ ar_expanded <- c() } else {
    ar_expanded <- ar_expanded[2:length(ar_expanded)]  # drop 0th order coef 
  }
  # multiply ma
  sma_poly <- c(1)
  for (i in smas) {sma_poly <- c(sar_poly, rep(0, (S-1)), i)}
  ma_expanded <- mult_poly(c(1, mas), sma_poly)
  if (length(ma_expanded) == 1){ ma_expanded <- c() } else {
    ma_expanded <- ma_expanded[2:length(ma_expanded)]  # drop 0th order coef
  }
  # ACF & PACF:
  # ----------------------------------------------------------------------------
  if(!is.null(ar_expanded) & !is.null(ma_expanded)){  # p+P>0 and q+Q>0
    acf.theoretical <- stats::ARMAacf(ar = ar_expanded, ma = ma_expanded, lag = max.lag)
    pacf.theoretical <- stats::ARMAacf(ar = ar_expanded, ma = ma_expanded, 
                                       lag = max.lag, pacf = TRUE)
  } else if(!is.null(ar_expanded) & is.null(ma_expanded)){  # p+P>0 and q+Q=0
    acf.theoretical <- stats::ARMAacf(ar = ar_expanded, lag = max.lag)
    pacf.theoretical <- stats::ARMAacf(ar = ar_expanded, lag = max.lag, pacf = TRUE)
  } else if(is.null(ar_expanded) & !is.null(ma_expanded)){  # p+P=0 and q+Q>0 
    acf.theoretical <- stats::ARMAacf(ma = ma_expanded, lag = max.lag)
    pacf.theoretical <- stats::ARMAacf(ma = ma_expanded, lag = max.lag, pacf = TRUE)
  }
  # Ljung-Box Statistic:
  # ----------------------------------------------------------------------------
  lb.pvals.ref <- c()
  lb.pvals.model <- c()
  for (i in c(1:max.lag)){
    lb.pvals.ref <- 
      c(lb.pvals.ref, stats::Box.test(m_ref$fit$residuals, lag = i, type = "Ljung-Box")$p.value)
    lb.pvals.model <- 
      c(lb.pvals.model, stats::Box.test(m$fit$residuals, lag = i, type = "Ljung-Box")$p.value)
  }
  # Spectral Density:
  # ----------------------------------------------------------------------------
  spec.ref <- astsa::mvspec(m_ref$fit$residuals, plot = FALSE)
  # theoretical spectral density of noise process:
  if(!is.null(ar_expanded) & !is.null(ma_expanded)){  # p+P>0 and q+Q>0
    spec.theoretical <- arma.spec2(ar = ar_expanded, ma = ma_expanded)
  } else if(!is.null(ar_expanded) & is.null(ma_expanded)){  # p+P>0 and q+Q=0
    spec.theoretical <- arma.spec2(ar = ar_expanded)
  } else if(is.null(ar_expanded) & !is.null(ma_expanded)){  # p+P=0 and q+Q>0 
    spec.theoretical <- arma.spec2(ma = ma_expanded)
  }
  # Plot:
  # ----------------------------------------------------------------------------
  if (plot) {
    par(mfrow=c(2,2), mar=c(2.5,2,2,2)+0.1)
    # ACF Plot:
    stats::acf(m_ref$fit$residuals, lag.max = max.lag, ylim = c(-0.75,1), xaxt = "n",
               main = '', xlab = '', ylab = '')
    points(c(0:max.lag), acf.theoretical, col = "dodgerblue", lwd = 2)
    axis(side = 1, at = seq(0,max.lag, by=2))
    mtext(substitute(bold("ACF")), side=3)
    mtext(substitute(italic("Lag")), side=1, line = 1.5, cex = 0.85)
    # PACF Plot:
    stats::pacf(m_ref$fit$residuals, lag.max = max.lag, ylim = c(-0.75,1), xaxt = "n",
                main = '', xlab = '', ylab = '')
    points(c(1:max.lag), pacf.theoretical, col = "dodgerblue", lwd = 2)
    axis(side = 1, at = seq(2,max.lag, by=2))
    mtext(substitute(bold("PACF")), side=3)
    mtext(substitute(italic("Lag")), side=1, line = 1.5, cex = 0.85)
    # Ljung-Box:
    plot(c(3:max.lag), lb.pvals.ref[3:length(lb.pvals.ref)], 
         ylim = c(0,1), lwd = 2,
         main = '', xlab = '', ylab = '')
    points(c(3:max.lag), lb.pvals.model[3:length(lb.pvals.model)], 
           ylim = c(0,1), col = "dodgerblue", lwd = 2)
    abline(h = 0.05, lty = "dashed", col = "blue")
    mtext(substitute(bold("Ljung-Box Statistic p-value")), side=3)
    mtext(substitute(italic("Lag")), side=1, line = 1.5, cex = 0.85)
    # Spectral Density:
    plot(spec.ref$freq, spec.ref$spec, type = "l", yaxt = "n", 
         lwd = 2, main = '', xlab = '', ylab = '')
    par(new = TRUE)
    plot(spec.theoretical$freq, spec.theoretical$spec, 
         xaxt = "n", yaxt = "n", type = "l", lwd = 2, col = "dodgerblue",
         main = '', xlab = '', ylab = '')
    mtext(substitute(bold("Spectral Density")), side=3)
    mtext(substitute(italic("Frequency")), side=1, line = 1.5, cex = 0.85)
    # reset par
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  # RETURN:
  # ----------------------------------------------------------------------------
  return(invisible(list("acf.theoretical" = acf.theoretical,
                        "pacf.theoretical" = pacf.theoretical,
                        "lb.pvals.model" = lb.pvals.model,
                        "spec.theoretical" = spec.theoretical)))
}