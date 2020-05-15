# Load necessary packages
library(vegan)

# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData")
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped_byspecies.RData")

# Extract metadata
coral.sp.metadata <- as(sample_data(phyASV.f.c.p), "data.frame")
# Run adonis analysis - weighted unifrac
ado1 <- adonis(phyloseq::distance(phyASV.f.c.p, method="wunifrac") ~ Coral_Species + Dist/Site, data = coral.sp.metadata,permutations = 9999)
ado1$aov.tab
# Run adonis analysis - bray-curtis
ado2 <- adonis(phyloseq::distance(phyASV.f.c.p, method="bray") ~ Coral_Species + Dist/Site, data = coral.sp.metadata,permutations = 9999)
ado2$aov.tab

# Extract metadata
Pocillopora.metadata <- as(sample_data(phyASV.f.c.p.Pocillopora), "data.frame")
Montipora.metadata <- as(sample_data(phyASV.f.c.p.Montipora), "data.frame")
Porites.metadata <- as(sample_data(phyASV.f.c.p.Porites), "data.frame")
Favites.metadata <- as(sample_data(phyASV.f.c.p.Favites), "data.frame")
FaviaM.metadata <- as(sample_data(phyASV.f.c.p.FaviaM), "data.frame")
Hydnophora.metadata <- as(sample_data(phyASV.f.c.p.Hydnophora), "data.frame")
Platygyra.metadata <- as(sample_data(phyASV.f.c.p.Platygyra), "data.frame")

# Run adonis analysis for each species
ado1.Pocillopora <- adonis(phyloseq::distance(phyASV.f.c.p.Pocillopora, 
                    method="wunifrac") ~ Dist/Site, 
                    data = Pocillopora.metadata,permutations = 9999)
ado1.Pocillopora$aov.tab

ado1.Porites <- adonis(phyloseq::distance(phyASV.f.c.p.Porites, 
                                              method="wunifrac") ~ Dist/Site, 
                           data = Porites.metadata,permutations = 9999)
ado1.Porites$aov.tab

ado1.Montipora <- adonis(phyloseq::distance(phyASV.f.c.p.Montipora, 
                                              method="wunifrac") ~ Dist/Site, 
                           data = Montipora.metadata,permutations = 9999)
ado1.Montipora$aov.tab

ado1.Favites <- adonis(phyloseq::distance(phyASV.f.c.p.Favites, 
                                              method="wunifrac") ~ Dist/Site, 
                           data = Favites.metadata,permutations = 9999)
ado1.Favites$aov.tab

ado1.FaviaM <- adonis(phyloseq::distance(phyASV.f.c.p.FaviaM, 
                                              method="wunifrac") ~ Dist/Site, 
                           data = FaviaM.metadata,permutations = 9999)
ado1.FaviaM$aov.tab

ado1.Hydnophora <- adonis(phyloseq::distance(phyASV.f.c.p.Hydnophora, 
                                              method="wunifrac") ~ Dist/Site, 
                           data = Hydnophora.metadata,permutations = 9999)
ado1.Hydnophora$aov.tab

ado1.Platygyra <- adonis(phyloseq::distance(phyASV.f.c.p.Platygyra, 
                                              method="wunifrac") ~ Dist/Site, 
                           data = Platygyra.metadata,permutations = 9999)
ado1.Platygyra$aov.tab

# Run adonis analysis for each species
# ado1.Pocillopora.B <- adonis(phyloseq::distance(phyASV.f.c.p.Pocillopora, 
#                                               method="bray") ~ Dist/Site, 
#                            data = Pocillopora.metadata,permutations = 9999)
# ado1.Pocillopora.B$aov.tab
# 
# ado1.Porites.B <- adonis(phyloseq::distance(phyASV.f.c.p.Porites, 
#                                           method="bray") ~ Dist/Site, 
#                        data = Porites.metadata,permutations = 9999)
# ado1.Porites.B$aov.tab
# 
# ado1.Montipora.B <- adonis(phyloseq::distance(phyASV.f.c.p.Montipora, 
#                                             method="bray") ~ Dist/Site, 
#                          data = Montipora.metadata,permutations = 9999)
# ado1.Montipora.B$aov.tab
# 
# ado1.Favites.B <- adonis(phyloseq::distance(phyASV.f.c.p.Favites, 
#                                           method="bray") ~ Dist/Site, 
#                        data = Favites.metadata,permutations = 9999)
# ado1.Favites.B$aov.tab
# 
# ado1.FaviaM.B <- adonis(phyloseq::distance(phyASV.f.c.p.FaviaM, 
#                                          method="bray") ~ Dist/Site, 
#                       data = FaviaM.metadata,permutations = 9999)
# ado1.FaviaM.B$aov.tab
# 
# ado1.Hydnophora.B <- adonis(phyloseq::distance(phyASV.f.c.p.Hydnophora, 
#                                              method="bray") ~ Dist/Site, 
#                           data = Hydnophora.metadata,permutations = 9999)
# ado1.Hydnophora.B$aov.tab
# 
# ado1.Platygyra.B <- adonis(phyloseq::distance(phyASV.f.c.p.Platygyra, 
#                                             method="bray") ~ Dist/Site, 
#                          data = Platygyra.metadata,permutations = 9999)
# ado1.Platygyra.B$aov.tab

# Save output
save.image(file="analyses/ITS2/sym_adonis.RData")
