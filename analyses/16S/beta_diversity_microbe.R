# Load necessary libraries
library(ggplot2)
library(phyloseq)
library(vegan)

# Set random seed for reproducibility
set.seed(1010)

# Load necessary data
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# For plotting boxplots
# From https://stackoverflow.com/questions/29943251/displaying-values-from-a-character-vector-as-italic-labels-in-boxplot-in-r
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))


# Transform sample counts to relative abundances (McMurdie and Holmes 2014)
phy.f.f <- transform_sample_counts(phy.f.f,function(x) x / sum(x))

# Set up function to easly use phyloseq objects in vegan
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Use function to convert phyloseq to vegan-friendly object
phy.f.f.veg <- vegan_otu(phy.f.f)
# Extract metadata
phy.f.f.metadata <- as(sample_data(phy.f.f), "data.frame")

# Calculate Bray-Curtis distance
phy.f.f.braydist <- phyloseq::distance(phy.f.f, method="bray")
# Run betadispersion calculation
phy.f.f.bd.braydist <- betadisper(d=phy.f.f.braydist, 
                              group=sample_data(phy.f.f)$human_disturbance,
                              bias.adjust=FALSE)
# Test for significant model
anova(phy.f.f.bd.braydist) # sig
# Test for differences between levels
TukeyHSD(phy.f.f.bd.braydist) # high-low, high-medium, veryhigh-high
# Plot betadispersion results
plot(phy.f.f.bd.braydist)
boxplot(phy.f.f.bd.braydist)

##
table(sample_data(phy.f.f)$field_host_name) # check how many samples from each species of coral. Favia speciosa (n=1), Acropora sp (n=3), Favites halicora (n=3), Favia sp (n=4) Removing all species that have super low sample size, as they are not used in this manuscript
phy.f.f.spp <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name!="Favia speciosa")
phy.f.f.spp <- subset_samples(phy.f.f.spp,sample_data(phy.f.f.spp)$field_host_name!="Favia sp")
phy.f.f.spp <- subset_samples(phy.f.f.spp,sample_data(phy.f.f.spp)$field_host_name!="Acropora sp")
phy.f.f.spp <- subset_samples(phy.f.f.spp,sample_data(phy.f.f.spp)$field_host_name!="Favites halicora")


sample_data(phy.f.f.spp)$field_host_name[sample_data(phy.f.f.spp)$field_host_name=="Platygyra daedalea"] <- "Platygyra ryukyuensis"
sample_data(phy.f.f.spp)$field_host_name[sample_data(phy.f.f.spp)$field_host_name=="Favia matthaii"] <- "D. matthaii"

# Order species levels
sample_data(phy.f.f.spp)$field_host_name <- factor(sample_data(phy.f.f.spp)$field_host_name,levels=c("Montipora aequituberculata","Pocillopora grandis","Porites lobata","Hydnophora microconos","Platygyra ryukyuensis","Favites pentagona","D. matthaii"))

# Calculate distance matrix
phy.f.f.brayspp <- phyloseq::distance(phy.f.f.spp, method="bray")
# Run betadisper
phy.f.f.bd.brayspp <- betadisper(d=phy.f.f.brayspp, 
                                  group=sample_data(phy.f.f.spp)$field_host_name,
                                  bias.adjust=FALSE)
# Run anova on betadisper
anova(phy.f.f.bd.brayspp)
# Run posthoc on betadisper
TukeyHSD(phy.f.f.bd.brayspp)
# Plot betadisper results
plot(phy.f.f.bd.brayspp)
# Boxplot of betadisper results
boxplot(phy.f.f.bd.brayspp)

# Plot
plot(phy.f.f.bd.braydist, hull=F, label=F, 
     main=NULL,col=distcols,
     xlab="PCoA 1", ylab="PCoA 2", sub="")
ordihull(phy.f.f.bd.braydist, sample_data(phy.f.f)$human_disturbance,  
         draw = c("polygon"), col=distcols,
         alpha=0.2, lwd=0.05)
################

phy.f.f.Mfol <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name=="Montipora aequituberculata")
phy.f.f.Mfol.veg <- vegan_otu(phy.f.f.Mfol)
phy.f.f.Mfol.metadata <- as(sample_data(phy.f.f.Mfol), "data.frame")
phy.f.f.Mfol.braydist <- phyloseq::distance(phy.f.f.Mfol, method="bray")
phy.f.f.Mfol.bd.braydist <- betadisper(d=phy.f.f.Mfol.braydist, 
                                  group=sample_data(phy.f.f.Mfol)$human_disturbance,
                                  bias.adjust=FALSE)
anova(phy.f.f.Mfol.bd.braydist) # Not Sig
TukeyHSD(phy.f.f.Mfol.bd.braydist) # 
plot(phy.f.f.Mfol.bd.braydist,col=distcols)
boxplot(phy.f.f.Mfol.bd.braydist,col=distcols)

##
phy.f.f.Plob <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name=="Porites lobata")
phy.f.f.Plob.veg <- vegan_otu(phy.f.f.Plob)
phy.f.f.Plob.metadata <- as(sample_data(phy.f.f.Plob), "data.frame")
phy.f.f.Plob.braydist <- phyloseq::distance(phy.f.f.Plob, method="bray")
phy.f.f.Plob.bd.braydist <- betadisper(d=phy.f.f.Plob.braydist, 
                                       group=sample_data(phy.f.f.Plob)$human_disturbance,
                                       bias.adjust=FALSE)
anova(phy.f.f.Plob.bd.braydist) # Not Sig
TukeyHSD(phy.f.f.Plob.bd.braydist) 
plot(phy.f.f.Plob.bd.braydist,col=distcols)
boxplot(phy.f.f.Plob.bd.braydist,col=distcols)

# Too few samples for Acropora sp., Favites halicora, Favia sp., and P. grandis

phy.f.f.Fpent <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name=="Favites pentagona")
phy.f.f.Fpent.veg <- vegan_otu(phy.f.f.Fpent)
phy.f.f.Fpent.metadata <- as(sample_data(phy.f.f.Fpent), "data.frame")
sample_data(phy.f.f.Fpent)$human_disturbance
phy.f.f.Fpent <- subset_samples(phy.f.f.Fpent,sample_data(phy.f.f.Fpent)$human_disturbance!="Low") # Remove low disturbance sample since there's only one

phy.f.f.Fpent.braydist <- phyloseq::distance(phy.f.f.Fpent, method="bray")
phy.f.f.Fpent.bd.braydist <- betadisper(d=phy.f.f.Fpent.braydist, 
                                       group=sample_data(phy.f.f.Fpent)$human_disturbance,
                                       bias.adjust=FALSE)
anova(phy.f.f.Fpent.bd.braydist) # Not sig
plot(phy.f.f.Fpent.bd.braydist,col=distcols)
boxplot(phy.f.f.Fpent.bd.braydist,col=distcols)


phy.f.f.Favmat <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name=="Favia matthaii")
sample_data(phy.f.f.Favmat)$human_disturbance # Only one low and one high disturbance, so remove these two disturbance levels from analysis.
phy.f.f.Favmat <- subset_samples(phy.f.f.Favmat,sample_data(phy.f.f.Favmat)$human_disturbance=="Medium"|sample_data(phy.f.f.Favmat)$human_disturbance=="VeryHigh")
phy.f.f.Favmat.veg <- vegan_otu(phy.f.f.Favmat)
phy.f.f.Favmat.metadata <- as(sample_data(phy.f.f.Favmat), "data.frame")
sample_data(phy.f.f.Favmat)$human_disturbance # Good
phy.f.f.Favmat.braydist <- phyloseq::distance(phy.f.f.Favmat, method="bray")
phy.f.f.Favmat.bd.braydist <- betadisper(d=phy.f.f.Favmat.braydist, 
                                        group=sample_data(phy.f.f.Favmat)$human_disturbance,
                                        bias.adjust=FALSE)
anova(phy.f.f.Favmat.bd.braydist) # Not Sig
plot(phy.f.f.Favmat.bd.braydist,col=distcols)
boxplot(phy.f.f.Favmat.bd.braydist,col=distcols)

##
phy.f.f.Hydno <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name=="Hydnophora microconos")
phy.f.f.Hydno.veg <- vegan_otu(phy.f.f.Hydno)
phy.f.f.Hydno.metadata <- as(sample_data(phy.f.f.Hydno), "data.frame")
sample_data(phy.f.f.Hydno)$human_disturbance #Good
phy.f.f.Hydno.braydist <- phyloseq::distance(phy.f.f.Hydno, method="bray")
phy.f.f.Hydno.bd.braydist <- betadisper(d=phy.f.f.Hydno.braydist, 
                                         group=sample_data(phy.f.f.Hydno)$human_disturbance,
                                         bias.adjust=FALSE)
anova(phy.f.f.Hydno.bd.braydist) # Not sig
plot(phy.f.f.Hydno.bd.braydist,col=distcols)
boxplot(phy.f.f.Hydno.bd.braydist,col=distcols)

##
phy.f.f.Platy <- subset_samples(phy.f.f,sample_data(phy.f.f)$field_host_name=="Platygyra daedalea")
phy.f.f.Platy.veg <- vegan_otu(phy.f.f.Platy)
phy.f.f.Platy.metadata <- as(sample_data(phy.f.f.Platy), "data.frame")
sample_data(phy.f.f.Platy)$human_disturbance # Good
phy.f.f.Platy.braydist <- phyloseq::distance(phy.f.f.Platy, method="bray")
phy.f.f.Platy.bd.braydist <- betadisper(d=phy.f.f.Platy.braydist, 
                                        group=sample_data(phy.f.f.Platy)$human_disturbance,
                                        bias.adjust=FALSE)
anova(phy.f.f.Platy.bd.braydist) # Not sig
plot(phy.f.f.Platy.bd.braydist,col=distcols)
boxplot(phy.f.f.Platy.bd.braydist,col=distcols)

# Save 
save.image(file="analyses/16S/mic_betadisper.RData")
save(phy.f.f.bd.braydist,phy.f.f.bd.brayspp, file="figures/mic_beta_plots.RData")

####
# Make boxplots
jpeg("figures/betadisper_box_all.jpg",height=3.2, width=5.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.bd.braydist,col=distcols[2:5])
dev.off()

pdf("figures/betadisper_box_all.pdf",height=3.2, width=4)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.bd.braydist,col=distcols[2:5])
dev.off()

jpeg("figures/betadisper_box_all_spp.jpg",height=5.5, width=4.5, unit="in",res=300)
par(mar=c(11.5,4,1,1))
boxplot(phy.f.f.bd.brayspp,col=speccols,las=2, names = make.italic(rownames(phy.f.f.bd.brayspp$centroids)),ylim=c(0,0.8))
dev.off()

pdf("figures/betadisper_box_all_spp.pdf",height=5.5, width=4.5)
par(mar=c(11.5,4,1,1))
boxplot(phy.f.f.bd.brayspp,col=speccols,las=2, names = make.italic(rownames(phy.f.f.bd.brayspp$centroids)),ylim=c(0,0.8))
dev.off()


jpeg("figures/betadisper_box.jpg",height=11, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1),mfrow=c(2,1))
layout(matrix(c(1,2)), heights=c(1,1.25))
boxplot(phy.f.f.bd.braydist,col=distcols[2:5])
par(mar=c(11,4,1,1))
boxplot(phy.f.f.bd.brayspp,col=rainbow(7),las=2)
dev.off()

pdf("figures/betadisper_box.pdf",height=11, width=4.5)
par(mar=c(2.5,4,1,1),mfrow=c(2,1))
layout(matrix(c(1,2)), heights=c(1,1.25))
boxplot(phy.f.f.bd.braydist,col=distcols[2:5])
par(mar=c(11,4,1,1))
boxplot(phy.f.f.bd.brayspp,col=rainbow(7),las=2)
dev.off()

###

jpeg("figures/betadisper_box_Mfol.jpg",height=3.5, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.Mfol.bd.braydist,col=distcols[2:5])
dev.off()

jpeg("figures/betadisper_box_Plob.jpg",height=3.5, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.Plob.bd.braydist,col=distcols[2:5])
dev.off()

jpeg("figures/betadisper_box_FPenta.jpg",height=3.5, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.Fpent.bd.braydist,col=distcols[3:5])
dev.off()

jpeg("figures/betadisper_box_FaviaM.jpg",height=3.5, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.Favmat.bd.braydist,col=distcols[c(3,5)])
dev.off()

jpeg("figures/betadisper_box_Hydno.jpg",height=3.5, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.Hydno.bd.braydist,col=distcols[2:5])
dev.off()

jpeg("figures/betadisper_box_Platy.jpg",height=3.5, width=4.5, unit="in",res=300)
par(mar=c(2.5,4,1,1))
boxplot(phy.f.f.Platy.bd.braydist,col=distcols[2:5])
dev.off()