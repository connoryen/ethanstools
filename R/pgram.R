#' Returns the periodogram of a signal
#' 
#' @description computes and plots the periodogram of a signal. The periodogram
#' is an estimate of the spectral density of a signal and is used to examine
#' the dominant frequencies present in the data. 
#' 
#' @param x Data
#' @param plot If TRUE (default), a graph is produced. If FALSE, no graph is 
#' produced.
#' @param ... Additional arguments passed to \code{plot}
#' 
#' @author Jared Fisher
#' 
#' @examples
#' pgram(rnorm(100))
#' pgram(sin(2*pi*0.1*0:100) + 3*cos(2*pi*0.33*0:100))
#' 
#' @export
pgram <- function(x, plot = TRUE, ...){
  m <- floor(length(x)/2)
  # we start at 2 because j=0 tells us nothing about the frequency space: 
  pgram <- abs(fft(x)[2:(m+1)])^2/length(x)
  # fourier frequencies:
  f_freqs <- c(1:floor(length(x)/2))/length(x)
  
  # plot
  if (plot == TRUE){
    plot(x=f_freqs, y=pgram, type = "h", col = "black", lwd = 2,
         xlab = "fourier frequency", ...)
    abline(h=0)
  }
  
  return(data.frame(f_freq = f_freqs, pgram = pgram))
}