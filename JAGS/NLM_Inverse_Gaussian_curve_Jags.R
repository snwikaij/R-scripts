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
x            <- df6$Cl
y            <- df6$Species
r            <- as.factor(df6$set)

plot(x, y, pch=19, col=as.factor(r))

df <- data.frame(x=x, y=y, r=r)

##############
#Create model#
##############

Inverse_Gaussian_curve <- function() {
  
  # Likelihood
  for (i in 1:N){
    
    y[i]    ~ dpois(mu[i])
    
    mu[i]  <- sqrt(lambda[r[i]]/(2*3.14*x[i]^3))*exp(-(lambda[r[i]]*(x[i]-M[r[i]])^2)/(2*M[r[i]]^2*x[i]))+U[r[i]]
    
  }
  
  for (r in 1:n_random){
    lambda[r]   ~ dgamma(lambda_alpha, lambda_beta)
    M[r]        ~ dgamma(M_alpha, M_beta)
    U[r]        ~ dexp(U_rate)}
  
  lambda_alpha  ~ dgamma(9, 1)
  lambda_beta   ~ dgamma(9, 0.3)
  
  M_alpha       ~ dgamma(1.5626, 0.0625)
  M_beta        ~ dgamma(2.5626, 0.0625)
  
  U_rate        ~ dexp(1)
  
  }

##################################
#Fit initiation params for chains#
##################################

data_Inverse_mod  <- list(N=nrow(df), x=df$x, y=df$y, r=df$r, n_random=length(unique(df$r)))
initstart             <- function(){list(lambda_alpha=runif(1), lambda_beta=runif(1),
                                         M_alpha=runif(1), M_beta=runif(1),
                                         U_rate=runif(1))}
save_params           <- c("lambda", "M", "U")

################
#Run Jags model#
################

Jags_mod <- jags(data = data_Inverse_mod,
                 model.file = Inverse_Gaussian_curve,
                 inits = initstart,
                 parameters.to.save = save_params,
                 n.iter = 10000,
                 n.burnin = 3000,
                 n.thin = 5,
                 n.chains = 5)

################################
#Plot the chains and posteriors#
################################

plot(as.mcmc(Jags_mod), density=F)
plot(as.mcmc(Jags_mod), trace=F)

#Estimation of lambda
print(mu_lambda <- mean(Jags_mod$BUGSoutput$sims.list$lambda))

#Estimation of mean
print(mu_M <- mean(Jags_mod$BUGSoutput$sims.list$M))

#Estimation of uppy
print(mu_U <- mean(Jags_mod$BUGSoutput$sims.list$U))

mean_list <- rowMeans(Jags_mod$BUGSoutput$sims.matrix[,c(1,2)])
lamb_list <- rowMeans(Jags_mod$BUGSoutput$sims.matrix[,c(6,7)])
uppy_list <- rowMeans(Jags_mod$BUGSoutput$sims.matrix[,c(3,4)])

yep <- cbind(lamb_list, mean_list, uppy_list)

hop <- yep[sample(1:nrow(yep), 100, T),]

par(mfrow=c(1,2))
plot(x=df$x, df$y, ylim=c(0,15), pch=19, col=df$r, xlab="Chloride g/L", ylab="Species richness")
for(i in 1:100){
  lines(df$x, sqrt(hop[i,1]/(2*3.14*df$x^3))*exp(-(hop[i,1]*(df$x-hop[i,2])^2)/(2*hop[i,2]^2*df$x))+hop[i,3])
}
lines(x=df$x, y= sqrt(mu_lambda/(2*3.14*df$x^3))*exp(-(mu_lambda*(df$x-mu_M)^2)/(2*mu_M^2*df$x))+mu_U, col="blue", lwd=2)

plot(x=log(df$x), df$y, ylim=c(0,15), pch=19, col=df$r, xlab="Chloride g/L [log transformed]", ylab="Species richness")
for(i in 1:100){
  lines(log(df$x), sqrt(hop[i,1]/(2*3.14*df$x^3))*exp(-(hop[i,1]*(df$x-hop[i,2])^2)/(2*hop[i,2]^2*df$x))+hop[i,3])
}
lines(x=log(df$x), y= sqrt(mu_lambda/(2*3.14*df$x^3))*exp(-(mu_lambda*(df$x-mu_M)^2)/(2*mu_M^2*df$x))+mu_U, col="blue", lwd=2)

dev.off()




