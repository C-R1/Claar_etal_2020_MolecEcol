# Load necessary packages
library(phyloseq)

# Load necessary data
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData")

# Calculate sample sizes
phyASV.f.c.MAeq <- subset_samples(phyASV.f.c,Coral_Species=="Montipora aequituberculata") 
phyASV.f.c.MAeq# n = 124
phyASV.f.c.Peyd <- subset_samples(phyASV.f.c,Coral_Species=="Pocillopora grandis") 
phyASV.f.c.Peyd# n = 117
phyASV.f.c.Plob <- subset_samples(phyASV.f.c,Coral_Species=="Porites lobata") 
phyASV.f.c.Plob# n = 96
phyASV.f.c.Hydno <- subset_samples(phyASV.f.c,Coral_Species=="Hydnophora microconos") 
phyASV.f.c.Hydno# n = 46
phyASV.f.c.Platy <- subset_samples(phyASV.f.c,Coral_Species=="Platygyra rykyuensis") 
phyASV.f.c.Platy# n = 39
phyASV.f.c.Fpenta <- subset_samples(phyASV.f.c,Coral_Species=="Favites pentagona") 
phyASV.f.c.Fpenta# n = 27
phyASV.f.c.FaviaM <- subset_samples(phyASV.f.c,Coral_Species=="Favia matthaii") 
phyASV.f.c.FaviaM# n = 23

# Total sample size
124+117+96+46+39+27+23 # = 472

# See which sites are included
sample_data(phyASV.f.c.MAeq)$Site
# 3 8 14 15 25 27 30 32 34 35 38

# Calculate sample size by species and site
phyASV.f.c.MAeq.3 <- subset_samples(phyASV.f.c.MAeq,Site==3) #12
phyASV.f.c.MAeq.8 <- subset_samples(phyASV.f.c.MAeq,Site==8) #11
phyASV.f.c.MAeq.14 <- subset_samples(phyASV.f.c.MAeq,Site==14) #13
phyASV.f.c.MAeq.15 <- subset_samples(phyASV.f.c.MAeq,Site==15) #12
phyASV.f.c.MAeq.25 <- subset_samples(phyASV.f.c.MAeq,Site==25) #12
phyASV.f.c.MAeq.27 <- subset_samples(phyASV.f.c.MAeq,Site==27) #13
phyASV.f.c.MAeq.30 <- subset_samples(phyASV.f.c.MAeq,Site==30) #10
phyASV.f.c.MAeq.32 <- subset_samples(phyASV.f.c.MAeq,Site==32) #5
phyASV.f.c.MAeq.34 <- subset_samples(phyASV.f.c.MAeq,Site==34) #12
phyASV.f.c.MAeq.35 <- subset_samples(phyASV.f.c.MAeq,Site==35) #12
phyASV.f.c.MAeq.38 <- subset_samples(phyASV.f.c.MAeq,Site==38) #12

phyASV.f.c.Peyd.3 <- subset_samples(phyASV.f.c.Peyd,Site==3) #12
phyASV.f.c.Peyd.8 <- subset_samples(phyASV.f.c.Peyd,Site==8) #12
phyASV.f.c.Peyd.14 <- subset_samples(phyASV.f.c.Peyd,Site==14) #12
phyASV.f.c.Peyd.15 <- subset_samples(phyASV.f.c.Peyd,Site==15) #9
phyASV.f.c.Peyd.25 <- subset_samples(phyASV.f.c.Peyd,Site==25) #10
phyASV.f.c.Peyd.27 <- subset_samples(phyASV.f.c.Peyd,Site==27) #10
phyASV.f.c.Peyd.30 <- subset_samples(phyASV.f.c.Peyd,Site==30) #12
phyASV.f.c.Peyd.32 <- subset_samples(phyASV.f.c.Peyd,Site==32) #5
phyASV.f.c.Peyd.34 <- subset_samples(phyASV.f.c.Peyd,Site==34) #12
phyASV.f.c.Peyd.35 <- subset_samples(phyASV.f.c.Peyd,Site==35) #11
phyASV.f.c.Peyd.38 <- subset_samples(phyASV.f.c.Peyd,Site==38) #12

phyASV.f.c.Plob.3 <- subset_samples(phyASV.f.c.Plob,Site==3) #9
phyASV.f.c.Plob.8 <- subset_samples(phyASV.f.c.Plob,Site==8) #10
phyASV.f.c.Plob.14 <- subset_samples(phyASV.f.c.Plob,Site==14) #10
phyASV.f.c.Plob.15 <- subset_samples(phyASV.f.c.Plob,Site==15) #3
phyASV.f.c.Plob.25 <- subset_samples(phyASV.f.c.Plob,Site==25) #8
phyASV.f.c.Plob.27 <- subset_samples(phyASV.f.c.Plob,Site==27) #11
phyASV.f.c.Plob.30 <- subset_samples(phyASV.f.c.Plob,Site==30) #10
phyASV.f.c.Plob.32 <- subset_samples(phyASV.f.c.Plob,Site==32) #4
phyASV.f.c.Plob.34 <- subset_samples(phyASV.f.c.Plob,Site==34) #9
phyASV.f.c.Plob.35 <- subset_samples(phyASV.f.c.Plob,Site==35) #11
phyASV.f.c.Plob.38 <- subset_samples(phyASV.f.c.Plob,Site==38) #11

phyASV.f.c.Hydno.3 <- subset_samples(phyASV.f.c.Hydno,Site==3) #5
phyASV.f.c.Hydno.8 <- subset_samples(phyASV.f.c.Hydno,Site==8) #5
phyASV.f.c.Hydno.14 <- subset_samples(phyASV.f.c.Hydno,Site==14) #4
phyASV.f.c.Hydno.15 <- subset_samples(phyASV.f.c.Hydno,Site==15) #4
phyASV.f.c.Hydno.25 <- subset_samples(phyASV.f.c.Hydno,Site==25) #5
phyASV.f.c.Hydno.27 <- subset_samples(phyASV.f.c.Hydno,Site==27) #3
phyASV.f.c.Hydno.30 <- subset_samples(phyASV.f.c.Hydno,Site==30) #4
phyASV.f.c.Hydno.32 <- subset_samples(phyASV.f.c.Hydno,Site==32) #1
phyASV.f.c.Hydno.34 <- subset_samples(phyASV.f.c.Hydno,Site==34) #7
phyASV.f.c.Hydno.35 <- subset_samples(phyASV.f.c.Hydno,Site==35) #4
phyASV.f.c.Hydno.38 <- subset_samples(phyASV.f.c.Hydno,Site==38) #4

phyASV.f.c.Platy.3 <- subset_samples(phyASV.f.c.Platy,Site==3) #6
phyASV.f.c.Platy.8 <- subset_samples(phyASV.f.c.Platy,Site==8) #1
phyASV.f.c.Platy.14 <- subset_samples(phyASV.f.c.Platy,Site==14) #5
phyASV.f.c.Platy.15 <- subset_samples(phyASV.f.c.Platy,Site==15) #0
phyASV.f.c.Platy.25 <- subset_samples(phyASV.f.c.Platy,Site==25) #4
phyASV.f.c.Platy.27 <- subset_samples(phyASV.f.c.Platy,Site==27) #7
phyASV.f.c.Platy.30 <- subset_samples(phyASV.f.c.Platy,Site==30) #3
phyASV.f.c.Platy.32 <- subset_samples(phyASV.f.c.Platy,Site==32) #2
phyASV.f.c.Platy.34 <- subset_samples(phyASV.f.c.Platy,Site==34) #5
phyASV.f.c.Platy.35 <- subset_samples(phyASV.f.c.Platy,Site==35) #1
phyASV.f.c.Platy.38 <- subset_samples(phyASV.f.c.Platy,Site==38) #1

phyASV.f.c.Fpenta.3 <- subset_samples(phyASV.f.c.Fpenta,Site==3) #1
phyASV.f.c.Fpenta.8 <- subset_samples(phyASV.f.c.Fpenta,Site==8) #3
phyASV.f.c.Fpenta.14 <- subset_samples(phyASV.f.c.Fpenta,Site==14) #3
phyASV.f.c.Fpenta.15 <- subset_samples(phyASV.f.c.Fpenta,Site==15) #0
phyASV.f.c.Fpenta.25 <- subset_samples(phyASV.f.c.Fpenta,Site==25) #3
phyASV.f.c.Fpenta.27 <- subset_samples(phyASV.f.c.Fpenta,Site==27) #4
phyASV.f.c.Fpenta.30 <- subset_samples(phyASV.f.c.Fpenta,Site==30) #3
phyASV.f.c.Fpenta.32 <- subset_samples(phyASV.f.c.Fpenta,Site==32) #4
phyASV.f.c.Fpenta.34 <- subset_samples(phyASV.f.c.Fpenta,Site==34) #3
phyASV.f.c.Fpenta.35 <- subset_samples(phyASV.f.c.Fpenta,Site==35) #1
phyASV.f.c.Fpenta.38 <- subset_samples(phyASV.f.c.Fpenta,Site==38) #3

phyASV.f.c.FaviaM.3 <- subset_samples(phyASV.f.c.FaviaM,Site==3) #1
phyASV.f.c.FaviaM.8 <- subset_samples(phyASV.f.c.FaviaM,Site==8) #2
phyASV.f.c.FaviaM.14 <- subset_samples(phyASV.f.c.FaviaM,Site==14) #2
phyASV.f.c.FaviaM.15 <- subset_samples(phyASV.f.c.FaviaM,Site==15) #0
phyASV.f.c.FaviaM.25 <- subset_samples(phyASV.f.c.FaviaM,Site==25) #1
phyASV.f.c.FaviaM.27 <- subset_samples(phyASV.f.c.FaviaM,Site==27) #3
phyASV.f.c.FaviaM.30 <- subset_samples(phyASV.f.c.FaviaM,Site==30) #5
phyASV.f.c.FaviaM.32 <- subset_samples(phyASV.f.c.FaviaM,Site==32) #2
phyASV.f.c.FaviaM.34 <- subset_samples(phyASV.f.c.FaviaM,Site==34) #3
phyASV.f.c.FaviaM.35 <- subset_samples(phyASV.f.c.FaviaM,Site==35) #2
phyASV.f.c.FaviaM.38 <- subset_samples(phyASV.f.c.FaviaM,Site==38) #2
