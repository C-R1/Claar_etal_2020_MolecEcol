# Load necessary packages
library(phyloseq)

# Load necessary data
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")

# Calculate number of samples
phy.f.f.MAaeq <- subset_samples(phy.f.f,field_host_name=="Montipora aequituberculata") 
phy.f.f.MAaeq# n = 67
phy.f.f.Plob <- subset_samples(phy.f.f,field_host_name=="Porites lobata") 
phy.f.f.Plob# n = 56
phy.f.f.Hydno <- subset_samples(phy.f.f,field_host_name=="Hydnophora microconos") 
phy.f.f.Hydno# n = 38
phy.f.f.Platy <- subset_samples(phy.f.f,field_host_name=="Platygyra daedalea") 
phy.f.f.Platy# n = 32
phy.f.f.Fpenta <- subset_samples(phy.f.f,field_host_name=="Favites pentagona") 
phy.f.f.Fpenta# n = 23
phy.f.f.Peyd <- subset_samples(phy.f.f,field_host_name=="Pocillopora grandis") 
phy.f.f.Peyd# n = 15
phy.f.f.FaviaM <- subset_samples(phy.f.f,field_host_name=="Favia matthaii"|field_host_name=="Favia sp") 
phy.f.f.FaviaM# n = 14 (10 if Favia sp excluded)

# Total number of samples
67+56+38+32+23+15+10# = 241
