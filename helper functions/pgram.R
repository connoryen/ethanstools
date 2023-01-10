# Modified from Jared Fisher:
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