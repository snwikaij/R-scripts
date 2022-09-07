library(R2jags)

###########################
#Use ``real life`` dataset#
###########################

df<- read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Hydrobiologia_Kaijser_et_al._2022.csv"), header = T)

y1 <- na.omit(df$PPO4[df$Species == "Lemna minor"])#blabla species is so indicative for higher PPO4 NOT!!!
y2 <- na.omit(df$PPO4[df$Species == "Callitriche brutia"])#species is soo sensitive for higher PPO4 NOT!!!

#################################################################
#Build Jags two mean comparison (bit clumpsy with the two loops)#
#################################################################

m2c     <- function(){
  
      #Likelihood 
        for (i in 1:N1) {
            y1[i]   ~ dnorm(mu1, tau1)}
  
        for (i in 1:N2) {
            y2[i]   ~ dnorm(mu2, tau2)}
  
       #Priors
            sigma1   ~ dgamma(6.26, 25)#Set as relative insensitive 
            tau1     <- 1/(sigma1 * sigma1)
            mu1      ~ dnorm(0.25, tau1);T(0,100000)#set above critical boundary of 100 (see Kaisjer et al. 2022)
            
            #Priors
            sigma2   ~ dgamma(1, 2)#Set as relative sensitive
            tau2     <- 1/(sigma2 * sigma2)
            mu2      ~ dnorm(0.10, tau2);T(0,100000)#set mu as 50 below the critical boundary (see Kaijser et al. 2022) 
            #and at just above this of Poikane et al. (2021)

}

###########################################
#Setup the data, iterations and parameters#
###########################################

data_m2c <- list("N1" = length(y1), "y1"=y1, "N2" = length(y2), "y2"=y2)
initsstart <- function (){list (mu1=runif(1),  sigma1=runif(1), mu2=runif(1),  sigma2=runif(1))}
parameters <- c("mu1", "sigma1", "mu2", "sigma2")

###########
#Run model#
###########

jags_m2c   <-  jags(data  =  data_m2c, 
                    model.file = m2c,
                    inits = initsstart, 
                    parameters.to.save = parameters,
                    n.chains = 4, 
                    n.thin = 5, 
                    jags.seed = 666,
                    n.iter = 5000,
                    n.burnin = 1000)

########################
#Plot trace and density#
########################

plot(as.mcmc(jags_m2c), density = F)
plot(as.mcmc(jags_m2c), trace = F)

################
#Calculate CLES#
################

mc.samp  <- jags_m2c$BUGSoutput$n.chains*jags_m2c$BUGSoutput$n.keep
shape1   <- jags_m2c$BUGSoutput$sims.list$mu1^2/jags_m2c$BUGSoutput$sims.list$sigma1^2
rate1    <- jags_m2c$BUGSoutput$sims.list$mu1/jags_m2c$BUGSoutput$sims.list$sigma1^2

shape2   <- jags_m2c$BUGSoutput$sims.list$mu2^2/jags_m2c$BUGSoutput$sims.list$sigma2^2
rate2    <- jags_m2c$BUGSoutput$sims.list$mu2/jags_m2c$BUGSoutput$sims.list$sigma2^2

 cl <- rep(NA, mc.samp)
for(i in 1:mc.samp){
 cl[i] <- sum(rgamma(1000, shape1[i], rate1[i])-rgamma(1000, shape2[i], rate2[i]) >= 0)/1000}

####################
#Plot CLES and data#
####################
 
main.lab <- paste0("The probability a random sample from \nLemna minor would fall higher \nalong the P-PO4 gradient than Callitriche brutia = \nP(L.minor > C. brutia) = ", round(mean(cl, na.rm = T),2), " (",round(quantile(cl,0.025, na.rm = T),2), "; ",round(quantile(cl,0.975, na.rm = T),2),")")
 
par(mfrow=c(1,2))
hist(cl, main = main.lab, breaks=30, xlab = "", cex.main=0.8)
abline(v=mean(cl), col="red", lwd="3")
abline(v=quantile(cl, c(.025, .975)), col="red", lwd=2, lty=3)

df.box.plot <- cbind.data.frame(Species=c(rep("L. minor", length(y1)), rep("C. brutia", length(y2))),
                                Values=c(y1, y2))

boxplot(log(df.box.plot$Values)~df.box.plot$Species, ylab="log transformed for visualization", xlab="Species")