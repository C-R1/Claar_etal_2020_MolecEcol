# Procrustes
# Code modified from http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-Section-01-Main-Lab.html#procrustes-and-ggplot2

# Load necessary libraries
library(plyr)
library(vegan)
library(ade4)
library(ggplot2)

# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Set seed for reproducibility
set.seed(1010)

# Select phyloseq objects
phylo1 <- phyASV.sym.Hydnophora
phylo2 <- phyASV.mic.Hydnophora

# Set filenames 
FILE <- "figures/Figure_5/procrustes_plot_Hydno.jpg"
FILE2 <- "figures/Figure_5/procrustes_plot_Hydno.pdf"

# Ordinate
sym.ord <- ordinate(physeq = phylo1,method = "PCoA",distance = "bray")
mic.ord <- ordinate(physeq = phylo2,method = "PCoA",distance = "bray")

# Extract vectors and incorporate sample data
sym.ord.df <- data.frame(sym.ord$vectors)
mic.ord.df <- data.frame(mic.ord$vectors)

sym.dim1 <- dim(sample_data(phylo1))[2]+1
sym.dim2 <- as.numeric(dim(sample_data(phylo1))[2]+dim(sym.ord.df)[2])
mic.dim1 <- dim(sample_data(phylo2))[2]+1
mic.dim2 <- as.numeric(dim(sample_data(phylo2))[2]+dim(mic.ord.df)[2])

sample_data(phylo1)[,c(sym.dim1:sym.dim2)] <- sym.ord.df
sample_data(phylo2)[,c(mic.dim1:mic.dim2)] <- mic.ord.df

phylo1.ord <- data.frame(sample_data(phylo1))
phylo2.ord <- data.frame(sample_data(phylo2))

# Order both tables by sample ID
phylo1.ord <- phylo1.ord[ order(row.names(phylo1.ord)), ]
phylo2.ord <- phylo2.ord[ order(row.names(phylo2.ord)), ]

rownames(phylo1.ord) == rownames(phylo2.ord) # double check

# Procrustes analysis
pro3 <- protest(phylo1.ord[, c("Axis.1", "Axis.2")],
                phylo2.ord[, c("Axis.1", "Axis.2")],
                scores = "sites", permutations = how(nperm = 9999))
pro3

# Plot
pro.p <- data.frame(rda1=pro3$Yrot[,1],
                    rda2=pro3$Yrot[,2],xrda1=pro3$X[,1],
                    xrda2=pro3$X[,2])
Dist <- data.frame(sample_data(phyASV.sym)[,12])
pro.p <- merge(pro.p,Dist,by="row.names")
rownames(pro.p) <- pro.p$Row.names
pro.p <- pro.p[ -c(1) ]
CoralSp <- data.frame(sample_data(phyASV.sym)[,"Coral_Species"])
pro.p <- merge(pro.p,CoralSp,by="row.names")


pro.plot.Hyd <- ggplot(pro.p) +
  theme(panel.background = element_blank(),
        legend.position = c(0.15,0.1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        axis.title = element_blank()) +
  geom_point(aes(x=rda1, y=rda2, colour=Coral_Species)) +
  geom_point(aes(x=xrda1, y=xrda2, colour=Coral_Species)) +
  geom_segment(aes(x=xrda1,y=xrda2,xend=rda1,yend=rda2,colour=Coral_Species),
               arrow=arrow(length=unit(0.2,"cm"))) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_abline(slope = -1/pro3$rotation[1,2]) +
  geom_abline(slope = pro3$rotation[1,2])+
  scale_color_manual(values="#00FF92",labels=c("H. microconos"))

# Make figure
jpeg(FILE,width=6, height=6,units="in",res=300)
pro.plot.Hyd
dev.off()

pdf(FILE2,width=6, height=6,useDingbats = FALSE)
pro.plot.Hyd
dev.off()

# Save
save(pro.plot.Hyd,file="figures/Figure_5/procrustes_plot_Hydno.RData")