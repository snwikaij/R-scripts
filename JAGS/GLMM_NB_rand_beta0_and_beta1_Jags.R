library(R2jags)

###########################
#Using some simulated data#
###########################

x1          <- seq(4,10,0.1)
x2          <- seq(0,8,0.1)
y1          <- exp(.1*x1+1)+rpois(length(x1),2)*(1:length(x1))/2
y2          <- exp(0.15*x2)+rpois(length(x2),2)*(1:length(x2))/2
df          <- data.frame(x=c(x1, x2), y=c(y1, y2), r=c(rep("A", length(x1)), rep("B", length(x2))))
colnames(df)<- c("x", "y", "r")

plot(df$x, df$y, pch=19, col=as.factor(df$r))
##########################################################
#Linear model with random effects for slope and intercept#
##########################################################

glmm_jags     <- function(){
  
  #Likelihood
  for (i in 1:N){
    
    y[i]           ~ dnegbin(p[i], r)
    
    p[i]           <- r/(r+lambda[i])
    
    log(lambda[i]) <- mu[i]
    
    mu[i]          <- beta0[random[i]] + beta1[random[i]] * x[i]}
  
  for(r in 1:n_rand){
    beta0[r]         ~ dnorm(b0, 1/(t0*t0))
    beta1[r]         ~ dnorm(b1, 1/(t1*t1))}
  
    b0               ~ dnorm(1, 1/1)
    t0               ~ dexp(.01)
  
    b1               ~ dnorm(0.25,1/0.1)
    t1               ~ dexp(.01)
    
    r                ~ dexp(0.1)
}

##################################
#Fit initiation params for chains#
##################################

data_for_glmm.mod <- list(N=nrow(df), n_rand=length(unique(df$r)), random=as.factor(df$r), x=df$x, y=as.integer(df$y))
initstart         <- function(){list(b0=runif(1), t0=runif(1), b1=runif(1), t1=runif(1), r=runif(1))}
save_params       <- c("beta0", "beta1", "r")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_for_glmm.mod,
                 model.file = glmm_jags,
                 inits = initstart,
                 parameters.to.save = save_params,
                 n.iter = 10000,
                 n.burnin = 1000,
                 n.thin = 5,
                 n.chains = 5)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)

#Pooled estimation for beta0
mean(Jags_mod$BUGSoutput$sims.list$beta0)

#Pooled estimation for beta
mean(Jags_mod$BUGSoutput$sims.list$beta1)

plot(df$x, df$y, pch=19, ylab="Counts", xlab="Env-var", col=as.factor(df$r))
lines(x=seq(0,10,0.5), y=exp(mean(Jags_mod$BUGSoutput$sims.list$beta0)+seq(0,10,0.5)*mean(Jags_mod$BUGSoutput$sims.list$beta1)), col="red", lwd=2)
dev.off()
