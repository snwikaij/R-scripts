#load vegan and varespec dataset
library(vegan)
data(varespec)

#Run nmds
set.seed(1)
my_nmds <- metaMDS(varespec) 

#Store function as ggmds_fun
ggMDS_fun <- function(my_nmds, groups=NULL, colours=NULL, alpha_group=0.1, size_species=2,
                      size_point=2, pch_point=19){

if(is.null(groups)){warning("No groups are displayed only the row names. To change this give a character vector
                            with group names using the argument `group`.")}
if(!is.null(groups) && nrow(my_nmds$points) != length(groups)){stop("The sum of the groups needs 
                                                 to be the same length as the number 
                                                 of samples given in the NMDS")}
if(!is.null(colours) && length(unique(groups)) != length(unique(colours))){stop("Number of groups is not the same as number of unique colours")}
  
  positions_nmds <- list()
  
  positions_nmds[["species"]]              <- cbind.data.frame(my_nmds$species)
  colnames(positions_nmds[["species"]])    <- c("MDS1", "MDS2")
  
  if(!is.null(groups)){
    positions_nmds[["points"]]                    <- cbind.data.frame(my_nmds$points, groups)}else{
    positions_nmds[["points"]]                    <- setNames(cbind.data.frame(my_nmds$points, rownames(my_nmds$points)), c("MDS1", "MDS2", "groups"))}

sites_df   <- positions_nmds$points
species_df <- positions_nmds$species

mds_hull   <- split(sites_df, sites_df$groups)
mds_hull2  <- do.call(rbind, lapply(mds_hull, function(x) x[chull(x[c("MDS1", "MDS2")]),]))
rownames(mds_hull2) <- NULL

if(!is.null(colours)){coliboli <- colours}else{coliboli <- rainbow(length(unique(groups)))}

plot <- ggplot(sites_df, aes(MDS1, MDS2, col=groups))+
  annotate("text", x=species_df$MDS1,  y=species_df$MDS2, label=rownames(species_df), size=size_species)+
  geom_point(size=size_point, pch=pch_point)+
  scale_color_manual(values=coliboli, breaks=levels(unique(sites_df$groups)))+
  geom_polygon(data=mds_hull2, aes(x=MDS1, y=MDS2, group=groups), inherit.aes=F, alpha=alpha_group)+
  theme_classic()+labs(col="")+
  theme(axis.line = element_line(colour = "black", size = .4),
        axis.text = element_text(colour = "black"),
        legend.position = "bottom",
        legend.spacing = unit(0.01, "cm"),
        legend.text = element_text(size=6),
        legend.margin=margin(t=-15),
        panel.border = element_rect(colour = "black", fill=NA, size=.6),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title =  element_text(size=6))

return(plot)}

set.seed(1)
ggMDS_fun(my_nmds = my_nmds, groups = factor(sample(c("B", "C"), 24, T), levels = c("B", "C")), colours=c("black", "pink"))
