
  m <- c(40, 82, 47)
  v <- c(32, 11, 92)
  n <- c(404, 336, 1474)
  
  distribution <- "gamma"
  
  nsim                <- 250000
  quantile_acceptance <- 0.05
  mm                  <- 50
  vm                  <- 20
  r                   <- NULL
  
  #Group the observations
  observations <- data.frame(m=m, v=v, n=n)

  #Choose distribution for prior
  if(distribution == "normal"){
    prior1 <- rnorm(nsim, mm, vm)
  }else if(distribution == "gamma"){
    prior1 <- rgamma(nsim, mm^2/vm^2, mm/vm^2)}else{stop("No distribution selected")}
  
  #Set the priors 
  prior    <- list(prior1  = prior1,
                   prior2  = rexp(nsim, ifelse(is.null(r), 1/(sum(observations[,2]*observations[,3])/sum(observations[,3])), r))) 

  #Simulation model
  model     <- function(par1, par2){
    
    vec        <- rep(NA, sum(observations[,3]))
    
    for(i in 1:nrow(observations)){
    
    start          <- ifelse(i == 1, 0, sum(observations[c(1:(i-1)),3]))+1
    end            <- sum(observations[c(1:i),3])
      
    sim            <- rnorm(observations[i,3], par1, par2)
    vec[start:end] <- sim}
    
    return(c(mean(vec), sd(vec)))}
  
  sim_df <- array(dim = c(nsim, 2))
  
  #Simulate our data
  for(s in 1:nsim){
    
  print(s)
    
  sim_df[s,] <- model(prior[[1]][s], prior[[2]][s])}
  
  #Distance
  m_all      <- sum(observations[,1]*observations[,3])/sum(observations[,3])
  v_all      <- sum(observations[,2]*observations[,3])/sum(observations[,3])
  priors     <- do.call(rbind, prior)
  distance   <- rowSums(abs(sim_df-c(m_all, v_all)))
  
  sim_df     <- rbind(sim_df)
  priors     <- t(rbind(priors))
  statistics <- cbind(sim_df, priors)

  #Select closest < mean(sd) or quantile acceptance
  if(is.null(quantile_acceptance)){
  parameters             <- statistics[distance<v_all,]
  rownames(parameters)   <- NULL}else{
  parameters             <- statistics[order(distance),][1:(quantile_acceptance*nsim),]
  rownames(parameters)   <- NULL}
  
  #Linear adjustment
  estimates <- array(dim = c(nsim, 4))          
  estimates <- cbind(mean(parameters[,1])+resid(lm(parameters[,1]~parameters[,3])), 
                     mean(parameters[,2])+resid(lm(parameters[,2]~parameters[,4])))

  #Plot the linear correction and posterior
  par(mfrow=c(2,2), mar=c(5,5,2,2))
  
  #mu correction
  plot(parameters[,1], parameters[,3], main = "Posterior mu", xlab="Simulated mu", ylab = "Prior mu", col="grey80")
  abline(lm(parameters[,3]~parameters[,1]), col="red", lwd=2)
  points(parameters[,3], estimates[,1])
  abline(lm(estimates[,1]~parameters[,3]), col="blue", lwd=2)
  
  #sigma correction
  plot(parameters[,2], parameters[,4], main = "Posterior sigma", xlab="Simulated sigma", ylab = "Prior sigma", col="grey80")
  abline(lm(parameters[,4]~parameters[,2]), col="red", lwd=2)
  points(parameters[,4], estimates[,2])
  abline(lm(estimates[,2]~parameters[,4]), col="blue", lwd=2) 
  
  #posterior mu and sigma
  hist(estimates[,1], breaks = 30, xlab = "Posterior mu", main = "",
       xlim = c(quantile(estimates[,1], .005), quantile(estimates[,1], .995)))
  hist(estimates[,2], breaks = 20, xlab = "Posterior sigma", main = "",
       xlim = c(quantile(estimates[,2], .005), quantile(estimates[,2], .995)))

  qgamma(0.05, mean(estimates[,1])^2/mean(estimates[,2])^2, mean(estimates[,1])/mean(estimates[,2])^2)
  qgamma(0.95, mean(estimates[,1])^2/mean(estimates[,2])^2, mean(estimates[,1])/mean(estimates[,2])^2)