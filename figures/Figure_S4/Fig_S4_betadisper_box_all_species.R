# Load necessary libraries
library(ggplot2)
library(phyloseq)
library(vegan)

#Load necessary data
load(file="analyses/16S/mic_betadisper.RData")

# Make figure
jpeg("figures/Figure_S4/Fig_S4_betadisper_box_all_species.jpg",height=8, width=12, unit="in",res=300)
par(mar=c(2.5,4,1,1),mfrow=c(2,3))
boxplot(phy.f.f.Mfol.bd.braydist,col=distcols[c(2:3,5)])
boxplot(phy.f.f.Plob.bd.braydist,col=distcols[c(2:3,5)])
boxplot(phy.f.f.Fpent.bd.braydist,col=distcols[c(3,5)])
boxplot(phy.f.f.Favmat.bd.braydist,col=distcols[c(3,5)])
boxplot(phy.f.f.Hydno.bd.braydist,col=distcols[c(2:3,5)])
boxplot(phy.f.f.Platy.bd.braydist,col=distcols[c(2:3,5)])
dev.off()

pdf("figures/Figure_S4/Fig_S4_betadisper_box_all_species.pdf",height=8, width=12,useDingbats = FALSE)
par(mar=c(2.5,6,1,1),mfrow=c(2,3))
par(cex.lab=2) # is for y-axis
par(cex.axis=2) # is for x-axis
boxplot(phy.f.f.Mfol.bd.braydist,col=distcols[c(2:3,5)])
boxplot(phy.f.f.Plob.bd.braydist,col=distcols[c(2:3,5)])
boxplot(phy.f.f.Fpent.bd.braydist,col=distcols[c(3,5)])
boxplot(phy.f.f.Favmat.bd.braydist,col=distcols[c(3,5)])
boxplot(phy.f.f.Hydno.bd.braydist,col=distcols[c(2:3,5)])
boxplot(phy.f.f.Platy.bd.braydist,col=distcols[c(2:3,5)])
dev.off()
