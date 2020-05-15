# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData")

# Calculate number of samples by species and site
phyASV.mic.Montipora.3 <- subset_samples(phyASV.mic.Montipora,
                                         site=="Site_3") #12
phyASV.mic.Montipora.8 <- subset_samples(phyASV.mic.Montipora,
                                         site=="Site_8") #6
phyASV.mic.Montipora.14 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_14") #8
phyASV.mic.Montipora.25 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_25") #5
phyASV.mic.Montipora.27 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_27") #10
phyASV.mic.Montipora.30 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_30") #5
phyASV.mic.Montipora.32 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_32") #0
phyASV.mic.Montipora.34 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_34") #5
phyASV.mic.Montipora.35 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_35") #6
phyASV.mic.Montipora.38 <- subset_samples(phyASV.mic.Montipora,
                                          site=="Site_38") #9

phyASV.mic.Pocillopora.3 <- subset_samples(phyASV.mic.Pocillopora,
                                           site=="Site_3") #1
phyASV.mic.Pocillopora.8 <- subset_samples(phyASV.mic.Pocillopora,
                                           site=="Site_8") #7
phyASV.mic.Pocillopora.14 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_14") #0
phyASV.mic.Pocillopora.25 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_25") #0
phyASV.mic.Pocillopora.27 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_27") #0
phyASV.mic.Pocillopora.30 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_30") #0
phyASV.mic.Pocillopora.32 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_32") #0
phyASV.mic.Pocillopora.34 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_34") #0
phyASV.mic.Pocillopora.35 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_35") #3
phyASV.mic.Pocillopora.38 <- subset_samples(phyASV.mic.Pocillopora,
                                            site=="Site_38") #4

phyASV.mic.Porites.3 <- subset_samples(phyASV.mic.Porites,site=="Site_3") #7
phyASV.mic.Porites.8 <- subset_samples(phyASV.mic.Porites,site=="Site_8") #3
phyASV.mic.Porites.14 <- subset_samples(phyASV.mic.Porites,site=="Site_14") #4
phyASV.mic.Porites.25 <- subset_samples(phyASV.mic.Porites,site=="Site_25") #3
phyASV.mic.Porites.27 <- subset_samples(phyASV.mic.Porites,site=="Site_27") #7
phyASV.mic.Porites.30 <- subset_samples(phyASV.mic.Porites,site=="Site_30") #4
phyASV.mic.Porites.32 <- subset_samples(phyASV.mic.Porites,site=="Site_32") #0
phyASV.mic.Porites.34 <- subset_samples(phyASV.mic.Porites,site=="Site_34") #3
phyASV.mic.Porites.35 <- subset_samples(phyASV.mic.Porites,site=="Site_35") #6
phyASV.mic.Porites.38 <- subset_samples(phyASV.mic.Porites,site=="Site_38") #7

phyASV.mic.Hydnophora.3 <- subset_samples(phyASV.mic.Hydnophora,
                                          site=="Site_3") #5
phyASV.mic.Hydnophora.8 <- subset_samples(phyASV.mic.Hydnophora,
                                          site=="Site_8") #4
phyASV.mic.Hydnophora.14 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_14") #4
phyASV.mic.Hydnophora.25 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_25") #5
phyASV.mic.Hydnophora.27 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_27") #3
phyASV.mic.Hydnophora.30 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_30") #7
phyASV.mic.Hydnophora.32 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_32") #1
phyASV.mic.Hydnophora.34 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_34") #5
phyASV.mic.Hydnophora.35 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_35") #4
phyASV.mic.Hydnophora.38 <- subset_samples(phyASV.mic.Hydnophora,
                                           site=="Site_38") #4

phyASV.mic.Platygyra.3 <- subset_samples(phyASV.mic.Platygyra,
                                         site=="Site_3") #5
phyASV.mic.Platygyra.8 <- subset_samples(phyASV.mic.Platygyra,
                                         site=="Site_8") #0
phyASV.mic.Platygyra.14 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_14") #4
phyASV.mic.Platygyra.25 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_25") #4
phyASV.mic.Platygyra.27 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_27") #8
phyASV.mic.Platygyra.30 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_30") #3
phyASV.mic.Platygyra.32 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_32") #2
phyASV.mic.Platygyra.34 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_34") #3
phyASV.mic.Platygyra.35 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_35") #1
phyASV.mic.Platygyra.38 <- subset_samples(phyASV.mic.Platygyra,
                                          site=="Site_38") #1

phyASV.mic.Favites.3 <- subset_samples(phyASV.mic.Favites,site=="Site_3") #1
phyASV.mic.Favites.8 <- subset_samples(phyASV.mic.Favites,site=="Site_8") #1
phyASV.mic.Favites.14 <- subset_samples(phyASV.mic.Favites,site=="Site_14") #2
phyASV.mic.Favites.25 <- subset_samples(phyASV.mic.Favites,site=="Site_25") #3
phyASV.mic.Favites.27 <- subset_samples(phyASV.mic.Favites,site=="Site_27") #7
phyASV.mic.Favites.30 <- subset_samples(phyASV.mic.Favites,site=="Site_30") #3
phyASV.mic.Favites.32 <- subset_samples(phyASV.mic.Favites,site=="Site_32") #0
phyASV.mic.Favites.34 <- subset_samples(phyASV.mic.Favites,site=="Site_34") #2
phyASV.mic.Favites.35 <- subset_samples(phyASV.mic.Favites,site=="Site_35") #1
phyASV.mic.Favites.38 <- subset_samples(phyASV.mic.Favites,site=="Site_38") #3

phyASV.mic.FaviaM.3 <- subset_samples(phyASV.mic.FaviaM,site=="Site_3") #1
phyASV.mic.FaviaM.8 <- subset_samples(phyASV.mic.FaviaM,site=="Site_8") #1
phyASV.mic.FaviaM.14 <- subset_samples(phyASV.mic.FaviaM,site=="Site_14") #1
phyASV.mic.FaviaM.25 <- subset_samples(phyASV.mic.FaviaM,site=="Site_25") #0
phyASV.mic.FaviaM.27 <- subset_samples(phyASV.mic.FaviaM,site=="Site_27") #2
phyASV.mic.FaviaM.30 <- subset_samples(phyASV.mic.FaviaM,site=="Site_30") #1
phyASV.mic.FaviaM.32 <- subset_samples(phyASV.mic.FaviaM,site=="Site_32") #0
phyASV.mic.FaviaM.34 <- subset_samples(phyASV.mic.FaviaM,site=="Site_34") #1
phyASV.mic.FaviaM.35 <- subset_samples(phyASV.mic.FaviaM,site=="Site_35") #1
phyASV.mic.FaviaM.38 <- subset_samples(phyASV.mic.FaviaM,site=="Site_38") #0
