# Load in necessary RData files
load("analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData")
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped_byspecies.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Load necessary packages
library(vegan)
library(phyloseq)

# For plotting boxplots
# From https://stackoverflow.com/questions/29943251/displaying-values-from-a-character-vector-as-italic-labels-in-boxplot-in-r
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

# Remove low sample size coral species not included in analysis
phyASV.f.c.p <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species!="Favia speciosa")
phyASV.f.c.p <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species!="Favia sp")
phyASV.f.c.p <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species!="Acropora sp")
phyASV.f.c.p <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species!="Favites halicora")

# Fix names for plotting
sample_data(phyASV.f.c.p)$Coral_Species[sample_data(phyASV.f.c.p)$Coral_Species=="Platygyra rykyuensis"]<- "P. ryukyuensis"
sample_data(phyASV.f.c.p)$Coral_Species[sample_data(phyASV.f.c.p)$Coral_Species=="Favia matthaii"]<- "D. matthaii"

# Reorder factor levels
sample_data(phyASV.f.c.p)$Coral_Species <- factor(sample_data(phyASV.f.c.p)$Coral_Species,levels=c("Montipora aequituberculata","Pocillopora grandis","Porites lobata","Hydnophora microconos","P. ryukyuensis","Favites pentagona","D. matthaii"))

# Reorder factor levels
sample_data(phyASV.f.c.p)$Dist <- factor(sample_data(phyASV.f.c.p)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.FaviaM)$Dist <- factor(sample_data(phyASV.f.c.p.FaviaM)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.Favites)$Dist <- factor(sample_data(phyASV.f.c.p.Favites)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.Hydnophora)$Dist <- factor(sample_data(phyASV.f.c.p.Hydnophora)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.Montipora)$Dist <- factor(sample_data(phyASV.f.c.p.Montipora)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.Platygyra)$Dist <- factor(sample_data(phyASV.f.c.p.Platygyra)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.Pocillopora)$Dist <- factor(sample_data(phyASV.f.c.p.Pocillopora)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))
sample_data(phyASV.f.c.p.Porites)$Dist <- factor(sample_data(phyASV.f.c.p.Porites)$Dist, levels=c("VeryLow","Low","Medium","VeryHigh"))

# Calculate Unifrac distances for each compartment
coral.ufdist <- UniFrac(phyASV.f.c.p, weighted=T, 
                        normalized=F, parallel=F, fast=T)
Pocillopora.ufdist <- UniFrac(phyASV.f.c.p.Pocillopora, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Porites.ufdist <- UniFrac(phyASV.f.c.p.Porites, weighted=T, 
                       normalized=F, parallel=F, fast=T)
Montipora.ufdist <- UniFrac(phyASV.f.c.p.Montipora, weighted=T, 
                       normalized=F, parallel=F, fast=T)
FaviaM.ufdist <- UniFrac(phyASV.f.c.p.FaviaM, weighted=T, 
                            normalized=F, parallel=F, fast=T)
Hydnophora.ufdist <- UniFrac(phyASV.f.c.p.Hydnophora, weighted=T, 
                            normalized=F, parallel=F, fast=T)
Platygyra.ufdist <- UniFrac(phyASV.f.c.p.Platygyra, weighted=T, 
                            normalized=F, parallel=F, fast=T)
Favites.ufdist <- UniFrac(phyASV.f.c.p.Favites, weighted=T, 
                            normalized=F, parallel=F, fast=T)


# Run betadisper for each compartment
coral.bd.dist <- betadisper(d=coral.ufdist, 
                            group=sample_data(phyASV.f.c.p)$Dist,
                            type="centroid", bias.adjust=FALSE)
Pocillopora.bd.dist <- betadisper(d=Pocillopora.ufdist, 
                            group=sample_data(phyASV.f.c.p.Pocillopora)$Dist,
                            type="centroid", bias.adjust=FALSE)
Porites.bd.dist <- betadisper(d=Porites.ufdist, 
                                  group=sample_data(phyASV.f.c.p.Porites)$Dist,
                                  type="centroid", bias.adjust=FALSE)
Montipora.bd.dist <- betadisper(d=Montipora.ufdist, 
                                  group=sample_data(phyASV.f.c.p.Montipora)$Dist,
                                  type="centroid", bias.adjust=FALSE)
FaviaM.bd.dist <- betadisper(d=FaviaM.ufdist, 
                                  group=sample_data(phyASV.f.c.p.FaviaM)$Dist,
                                  type="centroid", bias.adjust=FALSE)
Hydnophora.bd.dist <- betadisper(d=Hydnophora.ufdist, 
                                  group=sample_data(phyASV.f.c.p.Hydnophora)$Dist,
                                  type="centroid", bias.adjust=FALSE)
Platygyra.bd.dist <- betadisper(d=Platygyra.ufdist, 
                                  group=sample_data(phyASV.f.c.p.Platygyra)$Dist,
                                  type="centroid", bias.adjust=FALSE)
Favites.bd.dist <- betadisper(d=Favites.ufdist, 
                                  group=sample_data(phyASV.f.c.p.Favites)$Dist,
                                  type="centroid", bias.adjust=FALSE)
coral.bd.spp <- betadisper(d=coral.ufdist, 
                            group=sample_data(phyASV.f.c.p)$Coral_Species,
                            type="centroid", bias.adjust=FALSE)


# Test betadisper results
betadisper.coral <- anova(coral.bd.dist)
betadisper.Pocillopora <- anova(Pocillopora.bd.dist) #
betadisper.Porites <- anova(Porites.bd.dist)
betadisper.Montipora <- anova(Montipora.bd.dist)
betadisper.FaviaM <- anova(FaviaM.bd.dist) #
betadisper.Hydnophora <- anova(Hydnophora.bd.dist) #
betadisper.Platygyra <- anova(Platygyra.bd.dist)
betadisper.Favites <- anova(Favites.bd.dist) #
betadisper.coral.spp <- anova(coral.bd.spp)

# Betadisper Tukey tests
TukeyHSD(coral.bd.spp)
TukeyHSD(coral.bd.dist)
TukeyHSD(Pocillopora.bd.dist)
TukeyHSD(FaviaM.bd.dist)
TukeyHSD(Hydnophora.bd.dist)
TukeyHSD(Favites.bd.dist)
TukeyHSD(Platygyra.bd.dist)

# Save 
save.image(file="analyses/ITS2/sym_betadisper.RData")
save(coral.bd.dist,coral.bd.spp,file="figures/sym_beta_plots.RData")

####
jpeg("figures/sym_betadisper_box_all.jpg",height=3.2, width=5.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(coral.bd.dist,col=distcols[1:5])
dev.off()

pdf("figures/sym_betadisper_box_all.pdf",height=3.2, width=4,useDingbats = FALSE)
par(mar=c(2.5,4,1,1))
boxplot(coral.bd.dist,col=distcols[1:5])
dev.off()


jpeg("figures/sym_betadisper_box_spp.jpg",height=5.5, width=4.5, unit="in",res=300)
par(mar=c(11,4,1,1))
boxplot(coral.bd.spp,col=speccols,las=2, names = make.italic(rownames(coral.bd.spp$centroids)),ylim=c(0,0.15))
mtext(side = 1,text = "a",adj = 0.1,padj=-1.5,col="darkgray")
mtext(side = 1,text = "b",adj = 0.22,padj=-1.5,col="darkgray")
mtext(side = 1,text = "a,b",adj = 0.37,padj=-1.5,col="darkgray")
mtext(side = 1,text = "b",adj = 0.5,padj=-1.5,col="darkgray")
mtext(side = 1,text = "d",adj = 0.64,padj=-1.5,col="darkgray")
mtext(side = 1,text = "c,d",adj = 0.78,padj=-1.5,col="darkgray")
mtext(side = 1,text = "c",adj = 0.91,padj=-1.5,col="darkgray")
dev.off()

pdf("figures/sym_betadisper_box_spp.pdf",height=5.5, width=4.5,useDingbats = FALSE)
par(mar=c(11,4,1,1))
boxplot(coral.bd.spp,col=speccols,las=2, names = make.italic(rownames(coral.bd.spp$centroids)),ylim=c(0,0.15))
mtext(side = 1,text = "a",adj = 0.1,padj=-1.5,col="darkgray")
mtext(side = 1,text = "b",adj = 0.22,padj=-1.5,col="darkgray")
mtext(side = 1,text = "a,b",adj = 0.37,padj=-1.5,col="darkgray")
mtext(side = 1,text = "b",adj = 0.5,padj=-1.5,col="darkgray")
mtext(side = 1,text = "d",adj = 0.64,padj=-1.5,col="darkgray")
mtext(side = 1,text = "c,d",adj = 0.78,padj=-1.5,col="darkgray")
mtext(side = 1,text = "c",adj = 0.91,padj=-1.5,col="darkgray")
dev.off()


