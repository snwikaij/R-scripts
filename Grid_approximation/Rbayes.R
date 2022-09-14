library(vegan)

data("dune")
data("dune.env")

x <- dune
g <- dune.env$Management

#####################################################################
#R calculation function from anosim vegan with credibility intervals#
#####################################################################

Rbayes <- function(data, group, prior_beta=c(1,5), ci=0.95, method = "euclidean", dist="base", 
                   nboot=5000, ndraws=10000, plot=T, show_progress = T,
                   transform=NULL){
      df       <- data
      df$group <- group
  
vec_ano <- rep(NA, nboot)

for(t in 1:length(vec_ano)){
  
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
      
  x         <- as.dist(x)
  grouping  <- as.factor(group)
  x.rank    <- rank(x)
  div       <- length(x)/2
  n         <- nrow(df)
  
  in_row    <- as.vector(as.dist(row(matrix(nrow = n, ncol = n))))
  in_col    <- as.vector(as.dist(col(matrix(nrow = n, ncol = n))))
  
  within    <- grouping[in_row] == grouping[in_col]
  
  aver      <- tapply(x.rank, within, mean)
  vec_ano[t]<- -diff(aver)/div}

  mu    <- mean(vec_ano)
  v     <- var(vec_ano)
  a     <- ((1-mu)/v-1/mu)*mu^2
  b     <- a*(1/mu-1)
  param <- c(a, b)

  grid        <- seq(0.0001, 0.9999, 0.001)
  grid.approx <- data.frame(grid = grid,
                             likelihood = dbeta(grid, param[1], param[2]),
                             prior = dbeta(grid, prior_beta[1], prior_beta[2]))
  grid.approx$posterior  <-  grid.approx$likelihood*grid.approx$prior/max(grid.approx$likelihood*grid.approx$prior)

if(plot == T){
  plot(grid.approx$grid, grid.approx$likelihood, type="l", col="green", xlab="", ylab="density")
  lines(grid.approx$grid, grid.approx$prior, type="l", col="red")
  lines(grid.approx$grid, grid.approx$posterior, type="l", col= "blue")
  legend(0.8, max(grid.approx$likelihood)*0.9, lty = 1, col=c("green", "red", "blue"),
         legend = c("likelihood", "prior", "posterior"))}

post_samp   <- sample(grid.approx$grid, ndraws,  T, prob = grid.approx$posterior)

results     <- setNames(cbind.data.frame(mean(post_samp), 
                                         quantile(post_samp, 1-(ci/2+0.5)), 
                                         quantile(post_samp, ci/2+0.5)), c("mu", "ll", "ul"))

rownames(results) <- NULL
return(results)}

##############
#Run function#
##############

res <- Rbayes(x, g, prior_beta = c(1,5), method = "euclidean", dist="vegan")
res
#check anosim
anosim(x=x, grouping=g, distance = "euclidean")

