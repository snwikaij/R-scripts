####################
#Create pseudo data#
####################

x <- 1:50
set.seed(666)
y <- 1.5*x+rnorm(length(x), 0, 10)

##############################
#Create linear model function#
##############################

linear_model <-  function(x,y){
  
  if(length(x) != length(y)){stop("x is not the same length as y")}
  
  n  <- length(x)
  xs <- sum(x)
  ys <- sum(y)
  xy <- sum(x*y)
  x2 <- sum(x^2)
  
  b1    <- ((n*xy)-(xs*ys))/((n*x2)-(xs^2))
  b0    <- (ys-b1*xs)/n
  resid <- y-(b0+b1*x)
  sigma <- sqrt(sum((resid-mean(resid))^2/n))

  return(list(parameters=c(b0, b1, sigma), residuals=resid))}

####################################################
#Set number of ABC simulations and extract observed#
####################################################

nsim            <- 200000
qtol            <- 0.005
m               <- linear_model(x,y)
observations    <- c(m$parameters[1], m$parameters[2], m$parameters[3])

################
#Set the priors#
################

priors    <- list(prior1  = rnorm(nsim, 1.5, 2),
                  prior2  = rnorm(nsim, 0, 5),
                  prior3  = runif(nsim, 0, m$parameters[3]*2)) 

##################
#Simulation model#
##################

model     <- function(par1, par2, par3){
  
    mu      <- par1 + par2 * x
    samples <- rnorm(length(x), mu, par3)
    
    mod     <- linear_model(x, samples)
    distance<- sum(abs(mod$parameters-observations))
  
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

  priors                <- t(do.call(rbind, priors))
  statistics            <- cbind(priors, sim_df)
  estimates             <- statistics[order(statistics[,7]),][1:c(qtol*nsim),]
  rownames(estimates)   <- NULL
  
###################
#Linear correction#
###################
  
  posterior <- data.frame(b0=mean(estimates[,1])+linear_model(estimates[,1], estimates[,4])$residuals,
                          b1=mean(estimates[,2])+linear_model(estimates[,2], estimates[,5])$residuals,
                          sigma=mean(estimates[,3])+linear_model(estimates[,3], estimates[,6])$residuals)
  
################################
#Plot posteriors and regression#
################################
  
par(mfrow=c(2,2), mar=c(5,5,2,2))

#posterior mu and sigma
hist(posterior$b0, breaks = 15, xlab = "Posterior b0", main = "")
hist(posterior$b1, breaks = 15, xlab = "Posterior b1", main = "")
hist(posterior$sigma, breaks = 15, xlab = "Posterior sigma", main = "")

plot(x, y, cex=2, pch=19)
for(i in 1:1000){
    lines(x, posterior$b0[i]+posterior$b1[i]*x, col=alpha(rgb(0,0,0), 0.05))
  }
abline(coef = c(mean(posterior$b0), mean(posterior$b1)), col="red", lwd=2)  

plot_label <- paste0("\nIntercept=", round(mean(posterior$b0, na.rm = T),2), " (",round(quantile(posterior$b0,0.025, na.rm = T),2), "; ",round(quantile(posterior$b0,0.975, na.rm = T),2),")",
                     "\nRegression coefficient=", round(mean(posterior$b1, na.rm = T),2), " (",round(quantile(posterior$b1,0.025, na.rm = T),2), "; ",round(quantile(posterior$b1,0.975, na.rm = T),2),")")

text(x=30, y=10, labels=plot_label)
