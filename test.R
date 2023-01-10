source("helper functions/load.R")

x <- JohnsonJohnson  # get test data
x <- as.vector(x)  # convert out of time series form
plot.ts(x)

sarima.compare1(x, c(p=4,d=1,q=0, P=0,D=1,Q=0,S=4), 
                ref = c(p=0,d=1,q=0, P=0,D=1,Q=0,S=4))

models <- list()
models[['model1']] <- c(p=4,d=1,q=0, P=0,D=1,Q=1,S=4)
models[['model2']] <- c(p=1,d=1,q=1, P=0,D=1,Q=1,S=4)

sarima.compare(x, models, ref = c(p=0,d=1,q=0, P=0,D=1,Q=0,S=4))

