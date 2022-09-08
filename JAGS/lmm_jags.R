library(R2jags)
library(brms)

###############################
#Using some ``real life`` data#
###############################

df              <- read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Hydrobiologia_Kaijser_et_al._2022.csv"), header = T)
df              <- df[df$Location == "Mountain",]
df$EC           <- log(df$EC)
df              <- df[df$EC < 8 & df$EC > 4.2,]
df              <- df[df$pH < 9 & df$pH > 7,]
df              <- na.omit(df[c("pH", "EC", "mittlere_tiefevop_opt")])
colnames(df)[3] <-"Width"

###########################################
#Linear model with one level random effect#
###########################################

lmm_Jags <- function(){
  
  #Likelihood
    for(i in 1:N){
    y[i]    ~ dnorm(mu[i], tau1)
    mu[i]   <- b0 + R[W[i]] + b1 * x[i]}
  
  #Random effect
    for(r in 1:nW){
    R[r] ~ dnorm(0, tau2)}
  
  #Priors
    b0      ~ dnorm(0, 1)
    b1      ~ dnorm(0.5, 0.5)
    
    sigma1  ~ dexp(1)
    tau1    <- 1/(sigma1 * sigma1)
    
    sigma2  ~ dexp(3)
    tau2    <- 1/(sigma2 * sigma2)
}

##################################
#Fit initiation params for chains#
##################################

data_for_lmm.mod <- list(N=nrow(df), nW=length(unique(df$Width)), x=df$pH, y=df$EC, W=as.factor(df$Width))
initstart        <- function(){list(b0=runif(1), b1=runif(1), sigma1=runif(1), sigma2=runif(1))}
save_params      <- c("b0", "b1", "R")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_for_lmm.mod,
     model.file = lmm_Jags,
     inits = initstart,
     parameters.to.save = save_params,
     n.iter = 2000,
     n.burnin = 1000,
     n.thin = 1,
     n.chains = 10)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)
dev.off()

#######################################
#Plot expected values in a scatterplot#
#######################################

plot(df$pH, df$EC, col=df$Width, type="p", pch=19, cex=2, xlab="pH", ylab="EC")
for(i in 1:300){
  lines(df$pH, as.numeric(Jags_mod$BUGSoutput$sims.list$b0[i])+as.numeric(Jags_mod$BUGSoutput$sims.list$b1[i])*df$pH, alpha=0.1, col=alpha(rgb(0,0,0), 0.05))
}
lines(df$pH, as.numeric(Jags_mod$BUGSoutput$mean$b0)+as.numeric(Jags_mod$BUGSoutput$mean$b1)*df$pH, col="red", lwd=2)

plot_label <- paste0("\nIntercept=", round(mean(b0, na.rm = T),2), " (",round(quantile(b0,0.025, na.rm = T),2), "; ",round(quantile(b0,0.975, na.rm = T),2),")",
                     "\nRegression coefficient=", round(mean(b1, na.rm = T),2), " (",round(quantile(b1,0.025, na.rm = T),2), "; ",round(quantile(b1,0.975, na.rm = T),2),")")

text(x=7.4, y=7, labels=plot_label)