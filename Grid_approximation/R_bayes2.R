library(vegan)

data("dune")
data("dune.env")

x <- dune
g <- dune.env$Management

#####################################################################
#R calculation function from anosim vegan with credibility intervals#
#####################################################################

R_anosim_bayes <- function(data, group, prior_R=c(1,5), prior_NULL=c(1.5,1), 
                           ci=0.95, method = "euclidean", dist="base", 
                           nboot=5000, ndraws=10000, plot=T, show_progress = T,
                           transform=NULL){
  df       <- data
  df$group <- group
  
  param_B <- function(vec){
  
    vec   <- vec[vec>0 & vec<1]
    mu    <- mean(vec)
    v     <- var(vec)
    a     <- ((1-mu)/v-1/mu)*mu^2
    b     <- a*(1/mu-1)
  
  return(c(a, b))}
  
  R_function <- function(x, group, randomize=F){
    
      x         <- as.dist(x)
      if(randomize == F){grouping  <- as.factor(group)}else{grouping <- as.factor(sample(group))}
      x.rank    <- rank(x)
      div       <- length(x)/2
      n         <- length(group)
  
      in_row    <- as.vector(as.dist(row(matrix(nrow = n, ncol = n))))
      in_col    <- as.vector(as.dist(col(matrix(nrow = n, ncol = n))))
  
      within    <- grouping[in_row] == grouping[in_col]
  
      aver      <- tapply(x.rank, within, mean)
      
  return(-diff(aver)/div)}

  vec_ano <- matrix(ncol=2, nrow=nboot)
  
  set.seed(123)
  
  for(t in 1:nrow(vec_ano)){
    
    if(show_progress == T){print(t)}
    
    dfboot       <- df[sample(1:nrow(df), nrow(df), replace=T),]
    group        <- dfboot$group
    
    if(is.null(transform)){
      dfboot       <- dfboot[,-ncol(dfboot)]}else{
      dfboot       <- decostand(dfboot[,-ncol(dfboot)], method = transform)}
    
    if(dist == "base"){
      x <- dist(dfboot, method = method)}
    else if(dist == "vegan"){
      x <- vegdist(dfboot, method = method)}else{stop("No correct method selected for `dist`.")}
       
    vec_ano[t,] <- c(R_function(x, group = group, randomize = F), R_function(x, group = group, randomize = T))
    
    }
  
  param1 <- param_B(vec_ano[,1])
  param2 <- param_B(vec_ano[,2])
  
  grid         <- seq(0.0001, 0.9999, 0.001)
  
  grid.approx1 <- data.frame(grid = grid,
                            likelihood = dbeta(grid, param1[1], param1[2]),
                            prior = dbeta(grid, prior_R[1], prior_R[2]))
  grid.approx1$posterior <- grid.approx1$likelihood*grid.approx1$prior/max(grid.approx1$likelihood*grid.approx1$prior)
  
  grid.approx0 <- data.frame(grid = grid,
                             likelihood = dbeta(grid, 1, param2[2]),
                             prior = dbeta(grid, prior_NULL[1], prior_NULL[2]))
  grid.approx0$posterior <-  grid.approx0$likelihood*grid.approx0$prior/max(grid.approx0$likelihood*grid.approx0$prior)
  
  if(plot == T){
    
    par(mfrow=c(1, 3))

    plot(grid.approx1$grid, grid.approx1$likelihood, main="P(theta1|Data)~P(Data|theta1)*P(theta1)", type="l", col="green", xlab="", ylab="density")
    lines(grid.approx1$grid, grid.approx1$prior, type="l", col="red")
    lines(grid.approx1$grid, grid.approx1$posterior, type="l", col= "blue")
    legend(0.6, max(grid.approx1$likelihood)*0.9, lty = 1, col=c("green", "red", "blue"),
           legend = c("likelihood", "prior", "posterior"))
    
    plot(grid.approx0$grid, grid.approx0$likelihood, main="P(theta0|Data)~P(Data|theta0)*P(theta0)", type="l", col="green", xlab="", ylab="density")
    lines(grid.approx0$grid, grid.approx0$prior, type="l", col="red")
    lines(grid.approx0$grid, grid.approx0$posterior, type="l", col= "blue")
    legend(0.6, max(grid.approx0$likelihood)*0.9, lty = 1, col=c("green", "red", "blue"),
           legend = c("likelihood", "prior", "posterior"))
    
    plot(grid.approx1$grid, grid.approx1$posterior, main="P(theta1|Data)/P(theta0|Data)", type="l", col="blue", xlab="", ylab="density")
    lines(grid.approx0$grid, grid.approx0$posterior, type="l", col="red")
    legend(0.5, 0.9, lty = 1, col=c("blue", "red"),
           legend = c("P(theta1|Data)", "P(theta0|Data)"))
    
    }
  
  post_samp1   <- sample(grid.approx1$grid, ndraws,  T, prob = grid.approx1$posterior)
  post_samp0   <- sample(grid.approx0$grid, ndraws,  T, prob = grid.approx0$posterior)
  
  results      <- setNames(cbind.data.frame(mean(post_samp1), quantile(post_samp1, 1-(ci/2+0.5)), quantile(post_samp1, ci/2+0.5),
                                            mean(post_samp0), quantile(post_samp0, 1-(ci/2+0.5)), quantile(post_samp0, ci/2+0.5),
                                            (mean(post_samp1/post_samp0))), 
                           c("theta1", "ll1", "ul1", "theta0", "ll0", "ul0", "posterior_odds"))
  
  rownames(results) <- NULL
  return(round(results,2))}

##############
#Run function#
##############

res <- R_anosim_bayes(x, g, method = "bray", dist = "vegan", transform = "hellinger")
res
#check anosim
anosim(x = x, grouping = g, distance = "bray")
