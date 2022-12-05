library(readr)
library(R2jags)

df  <- read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Hydrobiologia_Kaijser_et_al._2022.csv"), header = T)
dfp <- df[df$Species == "Fontinalis antipyretica",]
dfa <- df[df$Species != "Fontinalis antipyretica",]
dfa <- dfa[!duplicated(dfa$ID),]

dfp$Abundance <- 1
dfa$Abundance <- 0
dfp   <- dfp[c("Abundance", "EC")]
dfa   <- dfa[c("Abundance", "EC")]
df    <- na.omit(rbind(dfp, dfa))
df$EC <- log(df$EC)

plot(df$EC, df$Abundance)

##############
#Create model#
##############

logistic_jags <- function() {
  
  # Likelihood
  for (i in 1:N){
    y[i]            ~ dbern(mu[i])
    logit(mu[i])   <- alpha + beta * x[i]}
  
  # Priors
  alpha ~ dnorm(0, 1/10^2)
  beta  ~ dnorm(0, 1/1^2)}

##################################
#Fit initiation params for chains#
##################################

data_for_log.mod  <- list(N=nrow(df), x=df$EC, y=df$Abundance)
initstart         <- function(){list(alpha=runif(1), beta=runif(1))}
save_params       <- c("alpha", "beta")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_for_log.mod,
                 model.file = logistic_jags,
                 inits = initstart,
                 parameters.to.save = save_params,
                 n.iter = 20000,
                 n.burnin = 5000,
                 n.thin = 5,
                 n.chains = 2)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)

#Estimation of alpha
print(mu_alpha <- mean(Jags_mod$BUGSoutput$sims.list$alpha))

#Estimation of beta
print(mu_beta <- mean(Jags_mod$BUGSoutput$sims.list$beta))

x  <- seq(4,9,.1)
plot(df$EC, df$Abundance, pch=19)
lines(x=x, y=exp(mu_alpha+mu_beta*x)/(1+exp(mu_alpha+mu_beta*x)), col="red", lwd=2)
