# Load necessary libraries
library(ggplot2)
library(phyloseq)
library(vegan)

# Set random seed for reproducibility
set.seed(1010)

# Load necessary data
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Function for easily switching from phyloseq to vegan-friendly format
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Switch from phyloseq to vegan-friendly object
phy.f.f.p.veg <- vegan_otu(phy.f.f.p)
# Extract metadata
phy.f.f.p.metadata <- as(sample_data(phy.f.f.p), "data.frame")

# Run adonis test, with Bray-Curtis similarity, full model
ado <- adonis(phyloseq::distance(phy.f.f.p, method="bray") ~ human_disturbance/reef_name+field_host_name, data = phy.f.f.p.metadata,permutations = 9999)
ado$aov.tab

### Montipora aequitberculata
# Subset this species
phy.f.f.p.MAeq <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Montipora aequituberculata")
# Ordinate
ord3 <- ordinate(phy.f.f.p.MAeq,method="PCoA",formula = ~human_disturbance)
# Plot ordination
plot_ordination(physeq = phy.f.f.p.MAeq, ordination = ord3, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
# Extract metadata
phy.f.f.p.MAeq.metadata <- as(sample_data(phy.f.f.p.MAeq), "data.frame")
# Run adonis
ado.MAeq <- adonis(phyloseq::distance(phy.f.f.p.MAeq, method="bray") ~ human_disturbance/reef_name, data = phy.f.f.p.MAeq.metadata,permutations = 9999)
ado.MAeq$aov.tab

# Porites lobata
phy.f.f.p.Plob <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Porites lobata")
ord4 <- ordinate(phy.f.f.p.Plob,method="PCoA",formula = ~human_disturbance)
plot_ordination(physeq = phy.f.f.p.Plob, ordination = ord4, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
phy.f.f.p.Plob.metadata <- as(sample_data(phy.f.f.p.Plob), "data.frame")
ado.Plob <- adonis(phyloseq::distance(phy.f.f.p.Plob, method="bray") ~ human_disturbance/reef_name, data = phy.f.f.p.Plob.metadata,permutations = 9999)
ado.Plob$aov.tab

# Pocillopora grandis
phy.f.f.p.Pgrand <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Pocillopora grandis")
phy.f.f.p.Pgrand <- subset_samples(phy.f.f.p.Pgrand,sample_data(phy.f.f.p.Pgrand)$human_disturbance!="Low") # Remove low, since there's only one sample
ord4 <- ordinate(phy.f.f.p.Pgrand,method="PCoA",formula = ~human_disturbance)
plot_ordination(physeq = phy.f.f.p.Pgrand, ordination = ord4, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
phy.f.f.p.Pgrand.metadata <- as(sample_data(phy.f.f.p.Pgrand), "data.frame")

# Favites pentagona
phy.f.f.p.Favites <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Favites pentagona")
phy.f.f.p.Favites <- subset_samples(phy.f.f.p.Favites,sample_data(phy.f.f.p.Favites)$human_disturbance!="Low") # Remove low, since there's only one sample
ord4 <- ordinate(phy.f.f.p.Favites,method="PCoA",formula = ~human_disturbance)
plot_ordination(physeq = phy.f.f.p.Favites, ordination = ord4, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
phy.f.f.p.Favites.metadata <- as(sample_data(phy.f.f.p.Favites), "data.frame")
ado.Favites <- adonis(phyloseq::distance(phy.f.f.p.Favites, method="bray") ~ human_disturbance/reef_name, data = phy.f.f.p.Favites.metadata,permutations = 9999)
ado.Favites$aov.tab

# Favia matthaii
phy.f.f.p.FaviaM <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Favia matthaii")
phy.f.f.p.FaviaM <- subset_samples(phy.f.f.p.FaviaM,sample_data(phy.f.f.p.FaviaM)$human_disturbance!="Low") # Remove low, since there's only one sample
phy.f.f.p.FaviaM <- subset_samples(phy.f.f.p.FaviaM,sample_data(phy.f.f.p.FaviaM)$human_disturbance!="High") # Remove high, since there's only one sample
ord4 <- ordinate(phy.f.f.p.FaviaM,method="PCoA",formula = ~human_disturbance)
plot_ordination(physeq = phy.f.f.p.FaviaM, ordination = ord4, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
phy.f.f.p.FaviaM.metadata <- as(sample_data(phy.f.f.p.FaviaM), "data.frame")
ado.FaviaM <- adonis(phyloseq::distance(phy.f.f.p.FaviaM, method="bray") ~ human_disturbance/reef_name, data = phy.f.f.p.FaviaM.metadata,permutations = 9999)
ado.FaviaM$aov.tab

# Hydnophora microconos
phy.f.f.p.Hydno <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Hydnophora microconos")
ord4 <- ordinate(phy.f.f.p.Hydno,method="PCoA",formula = ~human_disturbance)
plot_ordination(physeq = phy.f.f.p.Hydno, ordination = ord4, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
phy.f.f.p.Hydno.metadata <- as(sample_data(phy.f.f.p.Hydno), "data.frame")
ado.Hydno <- adonis(phyloseq::distance(phy.f.f.p.Hydno, method="bray") ~ human_disturbance/reef_name, data = phy.f.f.p.Hydno.metadata,permutations = 9999)
ado.Hydno$aov.tab

# Platygyra daedalea
phy.f.f.p.Platy <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Platygyra daedalea")
ord4 <- ordinate(phy.f.f.p.Platy,method="PCoA",formula = ~human_disturbance)
plot_ordination(physeq = phy.f.f.p.Platy, ordination = ord4, color="human_disturbance")+
  theme(panel.background = element_blank())+
  scale_color_manual(values=distcols)+ 
  scale_fill_manual(values=distcols)+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=human_disturbance))
phy.f.f.p.Platy.metadata <- as(sample_data(phy.f.f.p.Platy), "data.frame")
ado.Platy <- adonis(phyloseq::distance(phy.f.f.p.Platy, method="bray") ~ human_disturbance/reef_name, data = phy.f.f.p.Platy.metadata,permutations = 9999)
ado.Platy$aov.tab

