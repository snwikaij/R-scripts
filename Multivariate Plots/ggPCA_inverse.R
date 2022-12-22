library(vegan)
data("varechem")

#Set pH >= 3 as High pH and < 3 as Low pH
colour_dataframe <- data.frame(group=ifelse(varechem$pH>=3, "High pH", "Low pH"), colour=ifelse(varechem$pH>=3, "pink", "green"))

#Remove pH vector from varechem 
pca_df <- varechem[-ncol(varechem)]

ggPCA_inverse <- function(data, rank_inverse=T, scale_arrows=1, text_pos=0.075, colour_groups=NULL){
  
  if(any(apply(data, 2, is.numeric)) == F){stop("Data needs to be all numeric")}
  if(any(apply(data, 2, is.na)) == T){stop("Cannot have NAs in data")}
  
  if(is.null(colour_groups)){
    groups  <- as.factor(1)
    scm     <- scale_colour_manual(breaks=1, values="grey75")
    lps     <- "none"}else{
      if(nrow(colour_groups) != nrow(data)){stop("Lenght of the colours should match with the number of rows of the data")}
      groups  <- colour_groups[,1]
      scm     <- scale_colour_manual(breaks = colour_groups[,1], values=colour_groups[,2])
      lps     <- "bottom"}
  
  if(rank_inverse == T){
    m <- apply(data, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))}else{
      m <- data}
  
  pl      <- vegan::rda(m)
  p1      <- as.data.frame(scores(pl)$sites)
  p2      <- as.data.frame(scores(pl)$species)
  
  ev      <- eigenvals(pl)[1:2]/sum(eigenvals(pl))
  ev2     <- paste0("PC (", round(ev*100),"%)")
  
  ggplot(p1, aes(x=PC1, y=PC2, col=groups))+xlab(ev2[1])+ylab(ev2[2])+
    geom_point(pch=19, size=3)+scm+
    geom_segment(data=p2, aes(x = 0, y = 0, xend = PC1*scale_arrows, yend = PC2*scale_arrows), 
                 col="black", lwd=1.05,
                 arrow = arrow(length = unit(0.35, "cm")))+
    annotate("text", x=ifelse(p2$PC1 > 0, p2$PC1*scale_arrows+text_pos, p2$PC1*scale_arrows-text_pos),
             y=ifelse(p2$PC2  > 0, p2$PC2*scale_arrows+text_pos, p2$PC2*scale_arrows-text_pos),
             label=rownames(p2), col="black", size=5, fontface=2)+
    labs(col="")+
    geom_vline(xintercept = 0, lty=2)+
    geom_hline(yintercept = 0, lty=2)+
    theme_classic()+
    theme(legend.position = lps)}

ggPCA_inverse(pca_df, rank_inverse = T, colour_groups = colour_dataframe, text_pos = 0.1)
