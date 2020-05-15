# Load necessary packages
library(vegan)
library(gridExtra)
library(metagMisc)
library(phyloseq)
library(ggplot2)
library(data.table)

# Load necessary data
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")

# Look at histogram of OTU counts
tdt = data.table(tax_table(phy.f.f),
                 TotalCounts = taxa_sums(phy.f.f),
                 OTU = taxa_names(phy.f.f))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts") + xlim(c(0,1000))

# Estimate richness - Chao1
rich <- estimate_richness(phy.f.f,measures = c("Observed","Chao1"))

# See how far off observed is from Chao1 estimate
sort((rich$Chao1 - rich$Observed)/rich$Observed)

# Plot Chao1 by site
plot_richness(phy.f.f, measures=c("Chao1"), x="reef_name")

# Merge by species, and plot Chao1
phy.f.f.spp <- merge_samples(phy.f.f, "field_host_name")
plot_richness(phy.f.f.spp, measures=c("Chao1"))
# Merge by site, and plot Chao1
phy.f.f.site <- merge_samples(phy.f.f, "reef_name")
plot_richness(phy.f.f.site, measures=c("Chao1"))
# Merge by disturbance, and plot Chao1
sample_data(phy.f.f)$human_disturbance <- as.character(sample_data(phy.f.f)$human_disturbance)
phy.f.f.dist <- merge_samples(phy.f.f, "human_disturbance")
plot_richness(phy.f.f.dist, measures=c("Chao1"))
# Merge all together and plot Chao1
phy.f.f.all <- merge_samples(phy.f.f,"center_project_name")
plot_richness(phy.f.f.all, measures=c("Chao1"))
# Extract Chao1 estimate
estimate_richness(phy.f.f.all, measures=c("Chao1"))

# Append Chao1 to sample data table
rich.Chao1 <- cbind(sample_data(phy.f.f),rich)
# Plot Chao1 using ggplot
ggplot(rich.Chao1, mapping=aes(x=human_disturbance, y=Chao1, color=field_host_name)) +
  geom_point()
# Aggregate Chao by human disturbance - calculate mean
Chao.mean <- aggregate(rich.Chao1, by=list(rich.Chao1$human_disturbance), 
                    FUN=mean, na.rm=TRUE)
# Aggregate Chao by human disturbance - calculate sd
Chao.sd <- aggregate(rich.Chao1$Chao1, by=list(rich.Chao1$human_disturbance), 
                  FUN=sd, na.rm=TRUE)
Chao.mean$Group.1 <- as.character(Chao.mean$Group.1)
# Make data frame with disturbance level, mean, sd
Chao.Dist <- data.frame("human_disturbance"=Chao.mean$Group.1, 
                        "Chao1_mean"=Chao.mean$Chao1, 
                        "Chao1_sd"=Chao.sd$x)
# Aggregate Chao by site - calculate mean
Chao.mean <- aggregate(rich.Chao1, by=list(rich.Chao1$reef_name), 
                       FUN=mean, na.rm=TRUE)
# Aggregate Chao by site - calculate sd
Chao.sd <- aggregate(rich.Chao1$Chao1, by=list(rich.Chao1$reef_name), 
                     FUN=sd, na.rm=TRUE)
Chao.mean$Group.1 <- as.character(Chao.mean$Group.1)
# Make data frame with site, mean, sd
Chao.Site <- data.frame("reef_name"=Chao.mean$Group.1, 
                        "Chao1_mean"=Chao.mean$Chao1, 
                        "Chao1_sd"=Chao.sd$x)

# Aggregate Chao by coral species - calculate mean
Chao.mean <- aggregate(rich.Chao1, by=list(rich.Chao1$field_host_name), 
                       FUN=mean, na.rm=TRUE)
# Aggregate Chao by coral species - calculate sd
Chao.sd <- aggregate(rich.Chao1$Chao1, by=list(rich.Chao1$field_host_name), 
                     FUN=sd, na.rm=TRUE)
Chao.mean$Group.1 <- as.character(Chao.mean$Group.1)
# Make data frame with disturbance level, mean, sd
Chao.Spp <- data.frame("field_host_name"=Chao.mean$Group.1, 
                       "Chao1_mean"=Chao.mean$Chao1, 
                       "Chao1_sd"=Chao.sd$x)
