# Load necessary packages
library(vegan)
library(gridExtra)
library(metagMisc)
library(phyloseq)
library(ggplot2)
library(data.table)

# Load necessary data
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData")

# Look at histogram of OTU counts
tdt = data.table(tax_table(phyASV.f.c),
                 TotalCounts = taxa_sums(phyASV.f.c),
                 OTU = taxa_names(phyASV.f.c))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts") + xlim(c(0,1000))

# Estimate richness - Chao1
rich <- estimate_richness(phyASV.f.c,measures = c("Observed","Chao1"))

# See how far off observed is from Chao1 estimate
sort((rich$Chao1 - rich$Observed)/rich$Observed)

# Plot Chao1 by site
plot_richness(phyASV.f.c, measures=c("Chao1"), x="Site")

# Merge by species, and plot Chao1
phyASV.f.c.spp <- merge_samples(phyASV.f.c, "Coral_Species")
plot_richness(phyASV.f.c.spp, measures=c("Chao1"))
# Merge by site, and plot Chao1
phyASV.f.c.site <- merge_samples(phyASV.f.c, "Site")
plot_richness(phyASV.f.c.site, measures=c("Chao1"))
# Merge by disturbance, and plot Chao1
sample_data(phyASV.f.c)$Dist <- as.character(sample_data(phyASV.f.c)$Dist)
phyASV.f.c.dist <- merge_samples(phyASV.f.c, "Dist")
plot_richness(phyASV.f.c.dist, measures=c("Chao1"))
# Merge all together and plot Chao1
phyASV.f.c.all <- merge_samples(phyASV.f.c,"SampleType")
plot_richness(phyASV.f.c.all, measures=c("Chao1"))
# Extract Chao1 estimate
estimate_richness(phyASV.f.c.all, measures=c("Chao1"))

# Append Chao1 to sample data table
rich.Chao1 <- cbind(sample_data(phyASV.f.c),rich)
# Plot Chao1 using ggplot
ggplot(rich.Chao1, mapping=aes(x=Dist, y=Chao1, color=Coral_Species)) +
  geom_point()
# Aggregate Chao by human disturbance - calculate mean
Chao.mean <- aggregate(rich.Chao1, by=list(rich.Chao1$Dist), 
                    FUN=mean, na.rm=TRUE)
# Aggregate Chao by human disturbance - calculate sd
Chao.sd <- aggregate(rich.Chao1$Chao1, by=list(rich.Chao1$Dist), 
                  FUN=sd, na.rm=TRUE)
Chao.mean$Group.1 <- as.character(Chao.mean$Group.1)
# Make data frame with disturbance level, mean, sd
Chao.Dist <- data.frame("Dist"=Chao.mean$Group.1, 
                        "Chao1_mean"=Chao.mean$Chao1, 
                        "Chao1_sd"=Chao.sd$x)
# Aggregate Chao by site - calculate mean
Chao.mean <- aggregate(rich.Chao1, by=list(rich.Chao1$Site), 
                       FUN=mean, na.rm=TRUE)
# Aggregate Chao by site - calculate sd
Chao.sd <- aggregate(rich.Chao1$Chao1, by=list(rich.Chao1$Site), 
                     FUN=sd, na.rm=TRUE)
Chao.mean$Group.1 <- as.character(Chao.mean$Group.1)
# Make data frame with site, mean, sd
Chao.Site <- data.frame("Site"=Chao.mean$Group.1, 
                        "Chao1_mean"=Chao.mean$Chao1, 
                        "Chao1_sd"=Chao.sd$x)

# Aggregate Chao by coral species - calculate mean
Chao.mean <- aggregate(rich.Chao1, by=list(rich.Chao1$Coral_Species), 
                       FUN=mean, na.rm=TRUE)
# Aggregate Chao by coral species - calculate sd
Chao.sd <- aggregate(rich.Chao1$Chao1, by=list(rich.Chao1$Coral_Species), 
                     FUN=sd, na.rm=TRUE)
Chao.mean$Group.1 <- as.character(Chao.mean$Group.1)
# Make data frame with disturbance level, mean, sd
Chao.Spp <- data.frame("Coral_Species"=Chao.mean$Group.1, 
                       "Chao1_mean"=Chao.mean$Chao1, 
                       "Chao1_sd"=Chao.sd$x)
