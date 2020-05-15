library(dichromat)
library(maptools)
library(scales)
library(RColorBrewer)
library(rgdal)
library(ggplot2)

### site data
sites<-read.csv('figures/ki_map_files/ki_sites_SymMicrobes.csv')

###village data
villages<-read.csv("figures/ki_map_files/KI_villagesDCC_2015update.csv", header = TRUE)

villages<-data.frame(c(1.989386333, 2.022090868, 1.983594483, 1.865048333, 
                       1.8508,1.8785,1.906))
villages$lon<-c(-157.4760637, -157.4884092, -157.3683462, -157.5522183,
                -157.257, -157.257, -157.257)
villages$pop<-c(1879, 2311, 955, 441, 1500, 1000, 500)
villages$village<-c("London", "Tabwakea", "Banana", "Poland", 
                    "legend1500", "legend1000", "legend500")
colnames(villages)[1]<-"lat"

### reordering levels for colouring plot
levels(sites$f.pressure)
sites$levels<-as.numeric(factor(sites$f.pressure, levels(factor(sites$f.pressure))[c(4,1,2,3)]))
sites$f.pressure<-factor(sites$f.pressure, levels(factor(sites$f.pressure))[c(4,1,2,3)])

## set palette for fishing pressure
fishing.cols<-c("#01665e", #dark turq
                "#c7eae5", #light turq
                "#5ab4ac", #mid turq
                "#8c510a") #dark brown
sitecols <- c("#8c510a","#c7eae5","#5ab4ac","#01665e")
fishing.cols<-as.data.frame(fishing.cols)
fishing.cols$f.pressure<-levels(sites$f.pressure)
sites$col<-fishing.cols$fishing.cols[match(sites$f.pressure, fishing.cols$f.pressure)]
sites$site.simple <- c("L1","M1","M4","VL1","M5","VH2","VH1","VH3","M3","M2","L2")

levels(sites$f.pressure) = c("Very High", "Medium", "Low", "Very Low")

###
pdf(file="figures/KI_map_KISymMicrobe_sites_villages_inset_bigger.pdf", 
    width = 7.5, height =7)
source("figures/ki_map_files/KI_base_B&W_bigger.R")

# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, 
        add=TRUE,inches=0.3, bg=alpha("red", 0.4))
## legend for village size
text(-157.19, 1.85, "1500 people", cex=0.66)   
text(-157.19, 1.877, "1000 people", cex=0.66)   
text(-157.1945, 1.904, "500 people", cex=0.66)  
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4) 
with(sites, text(lon, lat, label=site.simple, cex=0.3))
legend(-157.265, 2.075,legend=levels(sites$f.pressure), 
       pt.bg=sitecols, pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.275, 2.01, "Human disturbance", srt=90, cex=0.6)
segments(-157.263, 1.968,-157.263, 2.057)
text(-157.315, 1.88, "Bay of\nWrecks", cex = 0.45)
text(-157.524, 1.825, "Vaskess\nBay", cex = 0.45)
par(mar=c(1.3,0.9,0.5,0.5))
source("figures/ki_map_files/KI_base_inset_forbigger.R")

dev.off()


jpeg(filename="figures/KI_map_KISymMicrobe_sites_villages_inset_bigger.jpg",
     width = 7.5, height =7, units="in", res=300)
source("figures/ki_map_files/KI_base_B&W_bigger.R")

# village markers sized by population
symbols(villages$lon, villages$lat, circles=(villages$pop)/10, add=TRUE,inches=0.3, bg=alpha("red", 0.4))
## legend for village size
text(-157.19, 1.85, "1500 people", cex=0.66)   
text(-157.19, 1.877, "1000 people", cex=0.66)   
text(-157.1945, 1.904, "500 people", cex=0.66)  
points(sites$lon, sites$lat, bg=alpha(sites$col,0.8), pch=21, cex=1.4) 
with(sites, text(lon, lat, label=site.simple, cex=0.3))
legend(-157.265, 2.075,legend=levels(sites$f.pressure), 
       pt.bg=sitecols, pch=21, bty="n", pt.cex=1.4, cex=0.6)
text(-157.275, 2.01, "Human disturbance", srt=90, cex=0.6)
segments(-157.263, 1.968,-157.263, 2.057)
text(-157.315, 1.88, "Bay of\nWrecks", cex = 0.45)
text(-157.524, 1.825, "Vaskess\nBay", cex = 0.45)
par(mar=c(1.3,0.9,0.5,0.5))
source("figures/ki_map_files/KI_base_inset_forbigger.R")

dev.off()
