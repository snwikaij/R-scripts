library(R2jags)

#######################
#Create pseudo-dataset#
#######################

x <- 1:50

set.seed(666)
df <- data.frame(
  x = x,
  y = as.integer(exp(rnorm(50, 2, .5)+x*rnorm(50, -0.05, 0.01))))

##################################################
#Build Jags negative binomial model with log link#
##################################################

nb.mod_log     <- function(){
  
  #Likelihood
  for (i in 1:N){
    
    y[i]           ~ dnegbin(p[i], r)#likelihood of negative binomial distribution
    
    p[i]           <- r/(r+lambda[i])#dispersion/shape parameter of negative binomial distribution
    
    log(lambda[i]) <- mu[i]#lambda (mu similar to poisson mu) log link
    
    mu[i]          <- b0 + b1 * x[i]#linear model 
    
  }
  
  #Priors
  b0     ~ dnorm(2.25, 0.02)
  b1     ~ dnorm(-0.05,0.05)
  r      ~ dnorm(0, 1);T(0,1)
  
}

###################################################################
#Fit the data to simulate in a file, start and parameters to store#
###################################################################

data_for_nb.mod <- list("N" = nrow(df), "y" = df$y, "x" = df$x) 
initsstart <- function (){list (b0=runif(1), b1=runif(1), r=runif(1))}
parameters <- c("b0", "b1")

################
#Run Jags model#
################

jags_nb.mod <- jags(data  =  data_for_nb.mod, 
                    model.file = nb.mod_log,
                    inits = initsstart, 
                    parameters.to.save = parameters,
                    n.chains = 2, 
                    n.thin = 5, 
                    jags.seed = 666,
                    n.iter = 10000,
                    n.burnin = 1000)

################################
#Plot the chains and posteriors#
################################

results<- as.mcmc(jags_nb.mod)
plot(results)
dev.off()

###############
#Calculated R2#
###############

n_draws  <- jags_nb.mod$BUGSoutput$n.keep*jags_nb.mod$BUGSoutput$n.chains
r2       <- rep(NA, jags_nb.mod$model$iter())

b0       <- jags_nb.mod$BUGSoutput$sims.list$b0
b1       <- jags_nb.mod$BUGSoutput$sims.list$b1

x        <- jags_nb.mod$model$data()$x
y        <- jags_nb.mod$model$data()$y

for(i in 1:n_draws){
  
  r2[i] <- cor(y, exp(b0[i]+b1[i]*x))^2
  
}

#######################################
#Plot expected values in a scatterplot#
#######################################

plot(df$x, df$y, type="p", pch=19, cex=2, xlab="x", ylab="y")
for(i in 1:500){
  lines(df$x, exp(as.numeric(jags_nb.mod$BUGSoutput$sims.list$b0[i])+as.numeric(jags_nb.mod$BUGSoutput$sims.list$b1[i])*df$x), alpha=0.1, col=alpha(rgb(0,0,0), 0.05))
}
lines(df$x, exp(as.numeric(jags_nb.mod$BUGSoutput$mean$b0)+as.numeric(jags_nb.mod$BUGSoutput$mean$b1)*df$x), col="red", lwd=2)

plot_label <- paste0("R-squared=", round(mean(r2, na.rm = T),2), " (",round(quantile(r2,0.025, na.rm = T),2), "; ",round(quantile(r2,0.975, na.rm = T),2),")",
                     "\nIntercept=", round(mean(b0, na.rm = T),2), " (",round(quantile(b0,0.025, na.rm = T),2), "; ",round(quantile(b0,0.975, na.rm = T),2),")",
                     "\nRegression coefficient=", round(mean(b1, na.rm = T),2), " (",round(quantile(b1,0.025, na.rm = T),2), "; ",round(quantile(b1,0.975, na.rm = T),2),")")

text(x=35, y=10, labels=plot_label)
