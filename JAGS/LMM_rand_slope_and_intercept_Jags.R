library(R2jags)

###############################
#Using some ``real life`` data#
###############################

df              <- read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Hydrobiologia_Kaijser_et_al._2022.csv"), header = T)
df              <- df[df$Location == "Mountain",]
df$EC           <- log(df$EC)
df              <- df[df$EC < 8 & df$EC > 4.2,]
df              <- df[df$pH < 9 & df$pH > 7,]
df              <- na.omit(df[c("EC", "pH", "mittlere_tiefevop_opt")])
df              <- df[!df$mittlere_tiefevop_opt %in% c(0,3),]
colnames(df)    <- c("x", "y", "r")

##########################################################
#Linear model with random effects for slope and intercept#
##########################################################

lmm_jags <- function() {
  
  # Likelihood
  for (i in 1:N){
    y[i]    ~ dnorm(mu[i], tau1)
  
    mu[i]  <- alpha[random[i]] + beta[random[i]] * x[i]}
  
  # Random effects 
  for (r in 1:n_rand) {
    alpha[r]    ~ dnorm(mu_alpha, tau_alpha)
    beta[r]     ~ dnorm(mu_beta, tau_beta)
  }
  
  # Priors
  sigma1     ~ dexp(.01)
  tau1       <- 1/(sigma1 * sigma1)
  
  mu_alpha   ~ dnorm(0, 1/5)
  mu_beta    ~ dnorm(0.5, 1/0.5)
  
  sigma2     ~ dexp(.01)
  tau_alpha  <- 1/(sigma2 * sigma2)
  sigma3     ~ dexp(.01)
  tau_beta   <- 1/(sigma3 * sigma3)
  
}

##################################
#Fit initiation params for chains#
##################################

data_for_lmm.mod <- list(N=nrow(df), n_rand=length(unique(df$r)), random=as.factor(df$r), x=df$x, y=df$y)
initstart        <- function(){list(mu_alpha=runif(1), mu_beta=runif(1), sigma1=runif(1), sigma2=runif(1), sigma3=runif(1))}
save_params      <- c("alpha", "beta")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_for_lmm.mod,
                 model.file = lmm_jags,
                 inits = initstart,
                 parameters.to.save = save_params,
                 n.iter = 50000,
                 n.burnin = 10000,
                 n.thin = 100,
                 n.chains = 5)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)

#Pooled estimation for alpha
mean(Jags_mod$BUGSoutput$sims.list$alpha)

#Pooled estimation for beta
mean(Jags_mod$BUGSoutput$sims.list$beta)

plot(df$x, df$y, pch=19, xlab="EC mS/cm (log transformed)", ylab="pH", col=df$r)
abline(a=mean(Jags_mod$BUGSoutput$sims.list$alpha), b=mean(Jags_mod$BUGSoutput$sims.list$beta), col="red", lwd=2)
dev.off()
