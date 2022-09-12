library(parameters)

####################
#Create pseudo data#
####################

x <- 1:50
set.seed(666)
y <- 1.5*x+rnorm(length(x), 0, 10)

##############################
#Create linear model function#
##############################

linear_model <- function(x, y){
  
  if(length(x) != length(y)){stop("x is not the same length as y")}
  
  n  <- length(x)
  b1 <- (sum((x-mean(x))*(y-mean(y)))*1/n) / (sum((x-mean(x))^2)*1/n)
  b0 <- mean(y)-b1*mean(x)
  
  resid <- y-(b0+b1*x)
  sigma <- sqrt(sum((resid-mean(resid))^2/n))
  
  return(list(parameters=c(b0, b1, sigma), residuals=resid))}

####################################################
#Set number of ABC simulations and extract observed#
####################################################

nsim            <- 400000
qtol            <- 0.0025
m               <- linear_model(x,y)
observations    <- c(m$parameters[1], m$parameters[2], m$parameters[3])

################
#Set the priors#
################

#Uniform priors bleeh never do this for the parameters of interest (coefficients
#of fixed effects[prior1 and 2]) and often #sufficient information is available
priors    <- list(prior1  = runif(nsim, -15, 15),
                  prior2  = runif(nsim,  -5, 5),
                  prior3  = runif(nsim, 0, m$parameters[3]*2)) 

##################
#Simulation model#
##################

model     <- function(par1, par2, par3){
  
    mu      <- par1 + par2 * x
    samples <- rnorm(length(x), mu, par3)
    
    mod     <- linear_model(x, samples)
    distance<- sqrt(sum((mod$parameters-observations)^2))
  
  return(c(mod$parameters[1], mod$parameters[2], mod$parameters[3], distance))}

###################
#Simulate our data#
###################

sim_df <- matrix(ncol=4, nrow = nsim)

for(s in 1:nsim){
  
  print(s)
  
  sim_df[s,] <- model(priors[[1]][s], priors[[2]][s], priors[[3]][s])}

#######################
#Rejection < qtol*nsim#
#######################

  priors2               <- t(do.call(rbind, priors))
  statistics            <- setNames(cbind.data.frame(priors2, sim_df), c("prior1", "prior2", "prior3", "sim1", "sim2", "sim3", "epsilon"))
  estimates             <- statistics[order(statistics$epsilon),]
  rownames(estimates)   <- NULL
  estimates             <- estimates[1:c(qtol*nsim),]
  
  plot(estimates$sim1, estimates$epsilon)
  
###################
#Linear correction#
###################
  
  posterior <- data.frame(b0=mean(estimates$prior1)+linear_model(estimates$sim1, estimates$prior1)$residuals,
                          b1=mean(estimates$prior2)+linear_model(estimates$sim2, estimates$prior2)$residuals,
                          sigma=mean(estimates$prior3)+linear_model(estimates$sim3, estimates$prior3)$residuals)
  
################################
#Plot posteriors and regression#
################################
  
par(mfrow=c(2,2), mar=c(5,5,2,2))

#posterior mu and sigma
hist(posterior$b0, breaks = 10, xlab = "Posterior b0", main = "")
hist(posterior$b1, breaks = 10, xlab = "Posterior b1", main = "")
hist(posterior$sigma, breaks = 10, xlab = "Posterior sigma", main = "")

plot(x, y, cex=2, pch=19)
for(i in 1:1000){
    lines(seq(-5, 55, 2), posterior$b0[i]+posterior$b1[i]*seq(-5, 55, 2), col=alpha(rgb(0,0,0), 0.05))
  }
abline(coef = c(mean(posterior$b0), mean(posterior$b1)), col="red", lwd=2)  

plot_label <- paste0("\nIntercept=", round(mean(posterior$b0, na.rm = T),2), " (",round(quantile(posterior$b0,0.025, na.rm = T),2), "; ",round(quantile(posterior$b0,0.975, na.rm = T),2),")",
                     "\nRegression coefficient=", round(mean(posterior$b1, na.rm = T),2), " (",round(quantile(posterior$b1,0.025, na.rm = T),2), "; ",round(quantile(posterior$b1,0.975, na.rm = T),2),")")

text(x=30, y=10, labels=plot_label)

#Meewh close enough
parameters::parameters(summary(lm(y~x)))
plot_label
