library(R2jags)

#############################
#Using real and awkward data#
#############################
df          <- read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Aquatic_Botany_Kaijser_et_al._2019.csv"), header = T)

df1 <- df[df$Soort %in% c("Elodea canadensis", "Elodea nuttallii", "Ceratophyllum demersum", 
                          "Ceratophyllum submersum", "Azolla filiculoides", "Glyceria fluitans",
                          "Lemna gibba", "Lemna minor", "Lemna minuta", "Lemna trisulca", "Myriophyllum spicatum", "Nuphar lutea",
                          "Numphaeae alba", "Nymphaea alba ", "Potamogeton crispus", "Stuckenia pectinata", "Potamogeton pusillus", "Potamogeton trichoides",
                          "Ranunculus aquatilis ", "Ranunculus aquatilis", "Ruppia cirrhosa", "Ruppia maritima", "Zannichellia spp.",
                          "Callitriche spp."),]

df1$Soort <- plyr::mapvalues(df1$Soort, from=c("Ranunculus aquatilis ", "Nymphaea alba "), to=c("Ranunculus aquatilis", "Nymphaea alba"))
df2       <- df1[c("MPN.Jaar.Maand", "Soort", "Mediaan")]
colnames(df2) <- c("Code", "Species", "Cl")

df3 <- read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Hydrobiologia_Kaijser_et_al._2022.csv"), header = T)
df3 <- df3[c("ID", "Species", "Cl")]
colnames(df3)[1] <- "Code"

df2$set <- "NL"
df3$set <- "GE"

df4          <- rbind(df2, df3)
df5          <- aggregate(data=df4, Species~set*Code*Cl, length)
df6          <- df5
df6$Cl       <- as.numeric(df6$Cl/1000)
df6$set      <- as.factor(df6$set)
x            <- log(df6$Cl)
y            <- df6$Species
r            <- as.factor(df6$set)

plot(x, y, pch=19, col=as.factor(r))

df <- data.frame(x=x, y=y, r=r)

##############
#Create model#
##############

Gaussian_curve <- function() {
  
  # Likelihood
  for (i in 1:N){
    
    y[i]    ~ dpois(mu[i])
    
    #dnorm(mu[i], tau)
    
    mu[i]  <- H[r[i]]/(SD[r[i]]*sqrt(2*3.14))*exp(-0.5*((x[i]-M[r[i]])/SD[r[i]])^2)+U[r[i]]
    
  }
  
  for (r in 1:n_random){
    H[r]   ~ dgamma(H_alpha, H_beta)
    M[r]   ~ dnorm(M_mu, 1/M_sd^2)
    SD[r]  ~ dgamma(SD_alpha, SD_beta)
    U[r]   ~ dgamma(U_alpha, U_beta)}
  
    H_alpha   ~ dgamma(4, 1)
    H_beta    ~ dgamma(0.4, 1)
  
    M_mu      ~ dnorm(-3, 1/1)
    M_sd      ~ dgamma(1, 1)
  
    SD_alpha  ~ dexp(1)
    SD_beta   ~ dexp(1)
  
    U_alpha   ~ dexp(1)
    U_beta    ~ dexp(1)}

##################################
#Fit initiation params for chains#
##################################

data_Gauss_mod    <- list(N=nrow(df), x=df$x, y=df$y, r=df$r, n_random=length(unique(df$r)))
initstart         <- function(){list(H_alpha=runif(1), H_beta=runif(1),
                                     M_mu=runif(1), M_sd=runif(1),
                                     SD_alpha=runif(1), SD_beta=runif(1),
                                     U_alpha=runif(1), U_beta=runif(1))}
save_params       <- c("M", "H", "SD", "U")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_Gauss_mod,
                 model.file = Gaussian_curve,
                 inits = initstart,
                 parameters.to.save = save_params,
                 n.iter = 5000,
                 n.burnin = 1000,
                 n.thin = 5,
                 n.chains = 5)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)

#Estimation of Height
print(mu_H <- mean(Jags_mod$BUGSoutput$sims.list$H))

#Estimation of Mu
print(mu_M <- mean(Jags_mod$BUGSoutput$sims.list$M))

#Estimation of SD
print(mu_SD<- mean(Jags_mod$BUGSoutput$sims.list$SD))

#Estimation of U
print(mu_U<- mean(Jags_mod$BUGSoutput$sims.list$U))

par(mfrow=c(1,3))
plot(x=df$x, df$y, pch=19, col=df$r, xlim=c(-5.5,3), main="log", xlab="Chloride g/L [Log transformed]", ylab="Species richness")
lines(x=df$x, y= mu_H/(mu_SD*sqrt(2*3.14))*exp(-0.5*((df$x-mu_M)/mu_SD)^2)+mu_U, col="blue", lwd=2)
plot(x=exp(df$x), df$y, pch=19, col=df$r, main="normal", xlab="Chloride g/L", ylab="Species richness")
lines(x=exp(df$x), y= mu_H/(mu_SD*sqrt(2*3.14))*exp(-0.5*((df$x-mu_M)/mu_SD)^2)+mu_U, col="blue", lwd=2)
plot(mu_H/(mu_SD*sqrt(2*3.14))*exp(-0.5*((df$x-mu_M)/mu_SD)^2)+mu_U, df$y-mu_H/(mu_SD*sqrt(2*3.14))*exp(-0.5*((df$x-mu_M)/mu_SD)^2)+mu_U,
     main="Predicted~Fitted", xlab="Predicted", ylab = "Fitted", col=df$r, pch=19)
dev.off()

cor(df$y, mu_H/(mu_SD*sqrt(2*3.14))*exp(-0.5*((df$x-mu_M)/mu_SD)^2)+mu_U)^2



