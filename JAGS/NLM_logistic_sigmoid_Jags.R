##########################
#Generate artificial data#
##########################
x <- seq(4, 10.5, 0.1)
A <- 150
M <- 6.2
S <- 0.3

set.seed(123)
y <- A / (1+exp((M-x)/S))+rgamma(length(x), 5, 0.08)
plot(x, y, pch=19)

df <- data.frame(x=x, y=y)

##############
#Create model#
##############

Asym_jags <- function() {
  
  # Likelihood
  for (i in 1:N){
    y[i]    ~ dnorm(mu[i], tau1)
    
    mu[i]  <- Asymptote / (1+exp((Midpoint-x[i])/Scale))}

  # Priors
  sigma1     ~ dexp(.001)
  tau1       <- 1/(sigma1 * sigma1)
  
  Asymptote ~ dgamma(1, .001)
  Midpoint  ~ dnorm(0, 1/10^2)
  Scale     ~ dt(0, 1/10^2, 1)}

##################################
#Fit initiation params for chains#
##################################

data_for_Asym.mod <- list(N=nrow(df), x=df$x, y=df$y)
initstart         <- function(){list(Asymptote=runif(1), Midpoint=runif(1), Scale=runif(1), sigma1=runif(1))}
save_params       <- c("Asymptote", "Midpoint", "Scale")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_for_Asym.mod,
                 model.file = Asym_jags,
                 inits = initstart,
                 parameters.to.save = save_params,
                 n.iter = 50000,
                 n.burnin = 10000,
                 n.thin = 10,
                 n.chains = 2)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)

#Estimation of Asymptote
print(mu_A <- mean(Jags_mod$BUGSoutput$sims.list$Asymptote))

#Estimation of Midpoint
print(mu_M <- mean(Jags_mod$BUGSoutput$sims.list$Midpoint))

#Estimation of scale
print(mu_S <- mean(Jags_mod$BUGSoutput$sims.list$Scale))

plot(df$x, df$y, pch=19, xlab="pH", ylab="Bicarbonate [mg/L]")
lines(x=df$x, y=mu_A/ (1+exp((mu_M-x)/mu_S)), col="red", lwd=2)
dev.off()
