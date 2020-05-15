# Load necessary packages
library(vegan)
library(gridExtra)
library(metagMisc)
library(phyloseq)
library(ggplot2)
library(data.table)

# Load necessary data
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData")

# Set random seed for reproducibility
set.seed(1010)

# Print number of OTUS
phyASV.f.c # 1073

# This is observed raw ASVs, this is what is in the manuscript
phyASV.f.c.spp <- merge_samples(phyASV.f.c, "Coral_Species")
plot_richness(phyASV.f.c.spp,measures=c("Observed"))
rich.spp <- estimate_richness(phyASV.f.c,measures="Observed")
rich.spp$Coral_Species <- sample_data(phyASV.f.c)$Coral_Species
spp.obs.mean <- aggregate(rich.spp, by=list(rich.spp$Coral_Species), 
                          FUN=mean, na.rm=TRUE)
spp.obs.sd <- aggregate(rich.spp, by=list(rich.spp$Coral_Species), 
                        FUN=sd, na.rm=TRUE)
spp.obs <- data.frame("Coral_Species"=spp.obs.mean$Group.1, "spp_obs_mean"=spp.obs.mean$Observed,"spp_obs_sd"=spp.obs.sd$Observed)
spp.obs
# Calculate observed ASV richness for each site
sample_data(phyASV.f.c)$Site <- as.character(sample_data(phyASV.f.c)$Site)
phyASV.f.c.site <- merge_samples(phyASV.f.c, "Site")
plot_richness(phyASV.f.c.site,measures=c("Observed"))
rich.site <- estimate_richness(phyASV.f.c,measures="Observed")
rich.site$Site <- sample_data(phyASV.f.c)$Site
site.obs.mean <- aggregate(rich.site, by=list(rich.site$Site), 
                           FUN=mean, na.rm=TRUE)
site.obs.sd <- aggregate(rich.site, by=list(rich.site$Site), FUN=sd, na.rm=TRUE)
site.obs <- data.frame("Site"=site.obs.mean$Group.1, "site_obs_mean"=site.obs.mean$Observed,"site_obs_sd"=site.obs.sd$Observed)
site.obs
# Calculate observed ASV richness for each disturbance level
phyASV.f.c.dist <- merge_samples(phyASV.f.c, "Dist")
plot_richness(phyASV.f.c.dist,measures=c("Observed"))
rich.dist <- estimate_richness(phyASV.f.c,measures="Observed")
rich.dist$Dist <- sample_data(phyASV.f.c)$Dist
dist.obs.mean <- aggregate(rich.dist, by=list(rich.dist$Dist), 
                           FUN=mean, na.rm=TRUE)
dist.obs.sd <- aggregate(rich.dist, by=list(rich.dist$Dist), FUN=sd, na.rm=TRUE)
dist.obs <- data.frame("Dist"=dist.obs.mean$Group.1, "dist_obs_mean"=dist.obs.mean$Observed,"dist_obs_sd"=dist.obs.sd$Observed)
dist.obs

rich.all <- estimate_richness(phyASV.f.c,measures="Observed")
all.obs.mean <- mean(rich.all$Observed)
all.obs.sd <- sd(rich.all$Observed)
all.obs <- data.frame("mean"=all.obs.mean,"sd"=all.obs.sd)
all.obs
