# Procrustes analyses
# Code modified from http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-Section-01-Main-Lab.html#procrustes-and-ggplot2

# Load necessary libraries
library(plyr)
library(vegan)
library(ade4)
library(ggplot2)

# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Set random seed for reproducibility
set.seed(1010)

# Reorder factor levels
sample_data(phyASV.sym)$Coral_Species <- as.factor(sample_data(phyASV.sym)$Coral_Species)
sample_data(phyASV.sym)$Coral_Species <- factor(sample_data(phyASV.sym)$Coral_Species, levels(sample_data(phyASV.sym)$Coral_Species)[c(4,6,7,3,5,2,1)]) 
# Reorder factor levels
sample_data(phyASV.mic)$field_host_name <- as.factor(sample_data(phyASV.mic)$field_host_name)
sample_data(phyASV.mic)$field_host_name <- factor(sample_data(phyASV.mic)$field_host_name, levels(sample_data(phyASV.mic)$field_host_name)[c(4,6,7,3,5,2,1)]) 

# Set phyloseq objects to compare
phylo1 <- phyASV.sym
phylo2 <- phyASV.mic
# Set filenames for figures
FILE <- "figures/Figure_5/procrustes_plot.jpg"
FILE2 <- "figures/Figure_5/procrustes_plot.pdf"

# Ordinate each phyloseq object
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

# Run procrustes analysis
pro3 <- protest(phylo1.ord[, c("Axis.1", "Axis.2")],
        phylo2.ord[, c("Axis.1", "Axis.2")],
        scores = "sites", permutations = how(nperm = 9999))
pro3

pro.p <- data.frame(rda1=pro3$Yrot[,1],
                    rda2=pro3$Yrot[,2],xrda1=pro3$X[,1],
                    xrda2=pro3$X[,2])
Dist <- data.frame(sample_data(phyASV.sym)[,12])
pro.p <- merge(pro.p,Dist,by="row.names")
rownames(pro.p) <- pro.p$Row.names
pro.p <- pro.p[ -c(1) ]
CoralSp <- data.frame(sample_data(phyASV.sym)[,"Coral_Species"])
pro.p <- merge(pro.p,CoralSp,by="row.names")


pro.plot <- ggplot(pro.p) +
  theme(panel.background = element_blank(),
        legend.position = c(0.22,0.25),
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
  geom_abline(slope = pro3$rotation[1,2]) +
  scale_color_manual(values=c("#FF0000", "#FFDB00", 
                              "#49FF00", "#00FF92", 
                              "#0092FF", "#4900FF", 
                              "#FF00DB"), name = "Coral Species")
pro.plot

# Make figures
jpeg(FILE,width=6, height=6,units="in",res=300)
pro.plot
dev.off()

pdf(FILE2,width=6, height=6,useDingbats = FALSE)
pro.plot
dev.off()

# Save
save(pro.plot,file="figures/procrustes_plot_allspp.RData")
