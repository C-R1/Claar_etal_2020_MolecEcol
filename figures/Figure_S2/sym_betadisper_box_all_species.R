# Load necessary data
load("analyses/ITS2/sym_betadisper.RData")

# Make figure
pdf("figures/Figure_S2/sym_betadisper_box_all_species.pdf",height=12, width=12,useDingbats = FALSE)
par(mar=c(2.5,4,1,1),mfrow=c(3,3))
par(cex.lab=1.25) # is for y-axis
par(cex.axis=1.25) # is for x-axis
boxplot(Montipora.bd.dist,col=distcols[c(1:3,5)])
boxplot(Pocillopora.bd.dist,col=distcols[c(1:3,5)])
boxplot(Porites.bd.dist,col=distcols[c(1:3,5)])
boxplot(Hydnophora.bd.dist,col=distcols[c(1:3,5)])
boxplot(Platygyra.bd.dist,col=distcols[c(2,3,5)])
boxplot(Favites.bd.dist,col=distcols[c(2,3,5)])
boxplot(FaviaM.bd.dist,col=distcols[c(2,3,5)])
dev.off()

jpeg("figures/Figure_S2/sym_betadisper_box_all_species.jpg",height=12, width=12,units="in",res=300)
par(mar=c(2.5,4,1,1),mfrow=c(3,3))
par(cex.lab=1.25) # is for y-axis
par(cex.axis=1.25) # is for x-axis
boxplot(Montipora.bd.dist,col=distcols[c(1:3,5)])
boxplot(Pocillopora.bd.dist,col=distcols[c(1:3,5)])
boxplot(Porites.bd.dist,col=distcols[c(1:3,5)])
boxplot(Hydnophora.bd.dist,col=distcols[c(1:3,5)])
boxplot(Platygyra.bd.dist,col=distcols[c(2,3,5)])
boxplot(Favites.bd.dist,col=distcols[c(2,3,5)])
boxplot(FaviaM.bd.dist,col=distcols[c(2,3,5)])
dev.off()
