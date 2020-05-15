# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# packageVersion("DESeq2")
library(DESeq2)
library("ggplot2")
library(gridExtra)

# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Set random seed for reproducibility
set.seed(1010)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Assign simple name to object
phylo <- phyASV.mic.raw

#Check
head(sample_data(phyASV.mic.raw))

# Convert phyloseq object to deseq object
phyASV.mic.deseq.C15dom <- phyloseq_to_deseq2(phylo, ~C15_dom)
# Calculate geometric means
geoMeans.C15 <- apply(counts(phyASV.mic.deseq.C15dom), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.C15dom <- estimateSizeFactors(phyASV.mic.deseq.C15dom, geoMeans = geoMeans.C15)
# Run deseq analysis
phyASV.mic.deseq.C15dom <- DESeq(phyASV.mic.deseq.C15dom, fitType = "local")
# Filter results based on alpha and specific contrast
res.C15 <- results(phyASV.mic.deseq.C15dom, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C15_dom","Y","N"))
alpha.C15 <- 0.05
sigtab.C15 <- res.C15[which(res.C15$padj < alpha.C15), ]
sigtab.C15 <- cbind(as(sigtab.C15, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C15), ], "matrix"))
head(sigtab.C15)


# Plot differential abundances
x.C15 <- tapply(sigtab.C15$log2FoldChange, sigtab.C15$Rank2, function(x) max(x))
x.C15 <- sort(x.C15, TRUE)
sigtab.C15$Rank2 <- factor(as.character(sigtab.C15$Rank2), levels=names(x.C15))

x.C15 <- tapply(sigtab.C15$log2FoldChange, sigtab.C15$Rank4, function(x) max(x))
x.C15 <- sort(x.C15, TRUE)
sigtab.C15$Rank4 <- factor(as.character(sigtab.C15$Rank4), levels=names(x.C15))

x.C15 <- tapply(sigtab.C15$log2FoldChange, sigtab.C15$Rank5, function(x) max(x))
x.C15 <- sort(x.C15, TRUE)
sigtab.C15$Rank5 <- factor(as.character(sigtab.C15$Rank5), levels=names(x.C15))

p1 <- ggplot(sigtab.C15, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p1

jpeg("figures/deseq_ASV_C15dom.jpg",height=4,width=4,units="in", res=300)
p1
dev.off()

# Convert phyloseq object to deseq object
phyASV.mic.deseq.C1dom <- phyloseq_to_deseq2(phylo, ~C1_dom)
# Calculate geometric means
geoMeans.C1 <- apply(counts(phyASV.mic.deseq.C1dom), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.C1dom <- estimateSizeFactors(phyASV.mic.deseq.C1dom, geoMeans = geoMeans.C1)
# Run deseq analysis
phyASV.mic.deseq.C1dom <- DESeq(phyASV.mic.deseq.C1dom, fitType = "local")
# Filter results based on alpha and specific contrast
res.C1 <- results(phyASV.mic.deseq.C1dom, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C1_dom","Y","N"))
alpha.C1 <- 0.05
sigtab.C1 <- res.C1[which(res.C1$padj < alpha.C1), ]
sigtab.C1 <- cbind(as(sigtab.C1, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C1), ], "matrix"))
head(sigtab.C1)


# Plot differential abundances
x.C1 <- tapply(sigtab.C1$log2FoldChange, sigtab.C1$Rank2, function(x) max(x))
x.C1 <- sort(x.C1, TRUE)
sigtab.C1$Rank2 <- factor(as.character(sigtab.C1$Rank2), levels=names(x.C1))

x.C1 <- tapply(sigtab.C1$log2FoldChange, sigtab.C1$Rank4, function(x) max(x))
x.C1 <- sort(x.C1, TRUE)
sigtab.C1$Rank4 <- factor(as.character(sigtab.C1$Rank4), levels=names(x.C1))

x.C1 <- tapply(sigtab.C1$log2FoldChange, sigtab.C1$Rank5, function(x) max(x))
x.C1 <- sort(x.C1, TRUE)
sigtab.C1$Rank5 <- factor(as.character(sigtab.C1$Rank5), levels=names(x.C1))

p2 <- ggplot(sigtab.C1, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p2

jpeg("figures/deseq_ASV_C1dom.jpg",height=4,width=4,units="in", res=300)
p2
dev.off()

# Convert phyloseq object to deseq object
phyASV.mic.deseq.C3dom <- phyloseq_to_deseq2(phylo, ~C3_dom)
# Calculate geometric means
geoMeans.C3 <- apply(counts(phyASV.mic.deseq.C3dom), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.C3dom <- estimateSizeFactors(phyASV.mic.deseq.C3dom, geoMeans = geoMeans.C3)
# Run deseq analysis
phyASV.mic.deseq.C3dom <- DESeq(phyASV.mic.deseq.C3dom, fitType = "local")
# Filter results based on alpha and specific contrast
res.C3 <- results(phyASV.mic.deseq.C3dom, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C3_dom","Y","N"))
alpha.C3 <- 0.05
sigtab.C3 <- res.C3[which(res.C3$padj < alpha.C3), ]
sigtab.C3 <- cbind(as(sigtab.C3, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C3), ], "matrix"))
head(sigtab.C3)


# Plot differential abundances
x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank2, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank2 <- factor(as.character(sigtab.C3$Rank2), levels=names(x.C3))

x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank4, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank4 <- factor(as.character(sigtab.C3$Rank4), levels=names(x.C3))

x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank5, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank5 <- factor(as.character(sigtab.C3$Rank5), levels=names(x.C3))

p3 <- ggplot(sigtab.C3, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p3

jpeg("figures/deseq_ASV_C3dom.jpg",height=4,width=4,units="in", res=300)
p3
dev.off()

# Convert phyloseq object to deseq object
phyASV.mic.deseq.C42dom <- phyloseq_to_deseq2(phylo, ~C42_dom)
# Calculate geometric means
geoMeans.C42 <- apply(counts(phyASV.mic.deseq.C42dom), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.C42dom <- estimateSizeFactors(phyASV.mic.deseq.C42dom, geoMeans = geoMeans.C42)
# Run deseq analysis
phyASV.mic.deseq.C42dom <- DESeq(phyASV.mic.deseq.C42dom, fitType = "local")
# Filter results based on alpha and specific contrast
res.C42 <- results(phyASV.mic.deseq.C42dom, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C42_dom","Y","N"))
alpha.C42 <- 0.05
sigtab.C42 <- res.C42[which(res.C42$padj < alpha.C42), ]
sigtab.C42 <- cbind(as(sigtab.C42, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C42), ], "matrix"))
head(sigtab.C42)


# Plot differential abundances
x.C42 <- tapply(sigtab.C42$log2FoldChange, sigtab.C42$Rank2, function(x) max(x))
x.C42 <- sort(x.C42, TRUE)
sigtab.C42$Rank2 <- factor(as.character(sigtab.C42$Rank2), levels=names(x.C42))

x.C42 <- tapply(sigtab.C42$log2FoldChange, sigtab.C42$Rank4, function(x) max(x))
x.C42 <- sort(x.C42, TRUE)
sigtab.C42$Rank4 <- factor(as.character(sigtab.C42$Rank4), levels=names(x.C42))

x.C42 <- tapply(sigtab.C42$log2FoldChange, sigtab.C42$Rank5, function(x) max(x))
x.C42 <- sort(x.C42, TRUE)
sigtab.C42$Rank5 <- factor(as.character(sigtab.C42$Rank5), levels=names(x.C42))

p4 <- ggplot(sigtab.C42, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p4

jpeg("figures/deseq_ASV_C42dom.jpg",height=4,width=4,units="in", res=300)
p4
dev.off()

# Convert phyloseq object to deseq object
phyASV.mic.deseq.C31dom <- phyloseq_to_deseq2(phylo, ~C31_dom)
# Calculate geometric means
geoMeans.C31 <- apply(counts(phyASV.mic.deseq.C31dom), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.C31dom <- estimateSizeFactors(phyASV.mic.deseq.C31dom, geoMeans = geoMeans.C31)
# Run deseq analysis
phyASV.mic.deseq.C31dom <- DESeq(phyASV.mic.deseq.C31dom, fitType = "local")
# Filter results based on alpha and specific contrast
res.C31 <- results(phyASV.mic.deseq.C31dom, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C31_dom","Y","N"))
alpha.C31 <- 0.05
sigtab.C31 <- res.C31[which(res.C31$padj < alpha.C31), ]
sigtab.C31 <- cbind(as(sigtab.C31, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C31), ], "matrix"))
head(sigtab.C31)


# Plot differential abundances
x.C31 <- tapply(sigtab.C31$log2FoldChange, sigtab.C31$Rank2, function(x) max(x))
x.C31 <- sort(x.C31, TRUE)
sigtab.C31$Rank2 <- factor(as.character(sigtab.C31$Rank2), levels=names(x.C31))

x.C31 <- tapply(sigtab.C31$log2FoldChange, sigtab.C31$Rank4, function(x) max(x))
x.C31 <- sort(x.C31, TRUE)
sigtab.C31$Rank4 <- factor(as.character(sigtab.C31$Rank4), levels=names(x.C31))

x.C31 <- tapply(sigtab.C31$log2FoldChange, sigtab.C31$Rank5, function(x) max(x))
x.C31 <- sort(x.C31, TRUE)
sigtab.C31$Rank5 <- factor(as.character(sigtab.C31$Rank5), levels=names(x.C31))

p5 <- ggplot(sigtab.C31, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p5

jpeg("figures/deseq_ASV_C31dom.jpg",height=4,width=4,units="in", res=300)
p5
dev.off()

# Convert phyloseq object to deseq object
phyASV.mic.deseq.D1dom <- phyloseq_to_deseq2(phylo, ~D1_dom)
# Calculate geometric means
geoMeans.D1 <- apply(counts(phyASV.mic.deseq.D1dom), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.D1dom <- estimateSizeFactors(phyASV.mic.deseq.D1dom, geoMeans = geoMeans.D1)
# Run deseq analysis
phyASV.mic.deseq.D1dom <- DESeq(phyASV.mic.deseq.D1dom, fitType = "local")
# Filter results based on alpha and specific contrast
res.D1 <- results(phyASV.mic.deseq.D1dom, cooksCutoff = FALSE, alpha = 0.05,contrast=c("D1_dom","Y","N"))
alpha.D1 <- 0.05
sigtab.D1 <- res.D1[which(res.D1$padj < alpha.D1), ]
sigtab.D1 <- cbind(as(sigtab.D1, "data.frame"), as(tax_table(phylo)[rownames(sigtab.D1), ], "matrix"))
head(sigtab.D1)


# Plot differential abundances
x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank2, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank2 <- factor(as.character(sigtab.D1$Rank2), levels=names(x.D1))

x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank4, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank4 <- factor(as.character(sigtab.D1$Rank4), levels=names(x.D1))

x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank5, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank5 <- factor(as.character(sigtab.D1$Rank5), levels=names(x.D1))

p6 <- ggplot(sigtab.D1, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p6

jpeg("figures/deseq_ASV_D1dom.jpg",height=4,width=4,units="in", res=300)
p6
dev.off()


################
phylo <- phyASV.mic.Favites.raw

# Convert phyloseq object to deseq object
phyASV.mic.deseq.FpentaD1 <- phyloseq_to_deseq2(phylo, ~D1_dom)
# Calculate geometric means
geoMeans.D1 <- apply(counts(phyASV.mic.deseq.FpentaD1), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.FpentaD1 <- estimateSizeFactors(phyASV.mic.deseq.FpentaD1, geoMeans = geoMeans.D1)
# Run deseq analysis
phyASV.mic.deseq.FpentaD1 <- DESeq(phyASV.mic.deseq.FpentaD1, fitType = "local")
# Filter results based on alpha and specific contrast
res.D1 <- results(phyASV.mic.deseq.FpentaD1, cooksCutoff = FALSE, alpha = 0.05,contrast=c("D1_dom","Y","N"))
alpha.D1 <- 0.05
sigtab.D1 <- res.D1[which(res.D1$padj < alpha.D1), ]
sigtab.D1 <- cbind(as(sigtab.D1, "data.frame"), as(tax_table(phylo)[rownames(sigtab.D1), ], "matrix"))
head(sigtab.D1)


# Plot differential abundances
x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank2, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank2 <- factor(as.character(sigtab.D1$Rank2), levels=names(x.D1))

x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank4, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank4 <- factor(as.character(sigtab.D1$Rank4), levels=names(x.D1))

x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank5, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank5 <- factor(as.character(sigtab.D1$Rank5), levels=names(x.D1))

sigtab.D1$Rank2 <- gsub(pattern="p__",replacement="",
                        sigtab.D1$Rank2)

sigtab.D1$Rank4 <- gsub(pattern="o__",replacement="",
                        sigtab.D1$Rank4)


p7 <- ggplot(sigtab.D1, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p7

jpeg("figures/deseq_ASV_FpentaD1.jpg",height=4,width=4,units="in", res=300)
p7 # Fpenta D1
dev.off()


# Convert phyloseq object to deseq object
phyASV.mic.deseq.FpentaC3 <- phyloseq_to_deseq2(phylo, ~C3_dom)
# Calculate geometric means
geoMeans.C3 <- apply(counts(phyASV.mic.deseq.FpentaC3), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.FpentaC3 <- estimateSizeFactors(phyASV.mic.deseq.FpentaC3, geoMeans = geoMeans.C3)
# Run deseq analysis
phyASV.mic.deseq.FpentaC3 <- DESeq(phyASV.mic.deseq.FpentaC3, fitType = "local")
# Filter results based on alpha and specific contrast
res.C3 <- results(phyASV.mic.deseq.FpentaC3, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C3_dom","Y","N"))
alpha.C3 <- 0.05
sigtab.C3 <- res.C3[which(res.C3$padj < alpha.C3), ]
sigtab.C3 <- cbind(as(sigtab.C3, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C3), ], "matrix"))
head(sigtab.C3)


# Plot differential abundances
x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank2, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank2 <- factor(as.character(sigtab.C3$Rank2), levels=names(x.C3))

x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank4, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank4 <- factor(as.character(sigtab.C3$Rank4), levels=names(x.C3))

x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank5, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank5 <- factor(as.character(sigtab.C3$Rank5), levels=names(x.C3))

sigtab.C3$Rank2 <- gsub(pattern="p__",replacement="",
                        sigtab.C3$Rank2)

sigtab.C3$Rank4 <- gsub(pattern="o__",replacement="",
                        sigtab.C3$Rank4)


p8 <- ggplot(sigtab.C3, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p8

jpeg("figures/deseq_ASV_FpentaC3.jpg",height=4,width=4,units="in", res=300)
p8 # Penta C3
dev.off()

################
phylo <- phyASV.mic.Platygyra.raw

# Convert phyloseq object to deseq object
phyASV.mic.deseq.PlatyD1 <- phyloseq_to_deseq2(phylo, ~D1_dom)
# Calculate geometric means
geoMeans.D1 <- apply(counts(phyASV.mic.deseq.PlatyD1), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.PlatyD1 <- estimateSizeFactors(phyASV.mic.deseq.PlatyD1, geoMeans = geoMeans.D1)
# Run deseq analysis
phyASV.mic.deseq.PlatyD1 <- DESeq(phyASV.mic.deseq.PlatyD1, fitType = "local")
# Filter results based on alpha and specific contrast
res.D1 <- results(phyASV.mic.deseq.PlatyD1, cooksCutoff = FALSE, alpha = 0.05,contrast=c("D1_dom","Y","N"))
alpha.D1 <- 0.05
sigtab.D1 <- res.D1[which(res.D1$padj < alpha.D1), ]
sigtab.D1 <- cbind(as(sigtab.D1, "data.frame"), as(tax_table(phylo)[rownames(sigtab.D1), ], "matrix"))
head(sigtab.D1)


# Plot differential abundances
x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank2, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank2 <- factor(as.character(sigtab.D1$Rank2), levels=names(x.D1))

x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank4, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank4 <- factor(as.character(sigtab.D1$Rank4), levels=names(x.D1))

x.D1 <- tapply(sigtab.D1$log2FoldChange, sigtab.D1$Rank5, function(x) max(x))
x.D1 <- sort(x.D1, TRUE)
sigtab.D1$Rank5 <- factor(as.character(sigtab.D1$Rank5), levels=names(x.D1))

sigtab.D1$Rank2 <- gsub(pattern="p__",replacement="",
                        sigtab.D1$Rank2)

sigtab.D1$Rank4 <- gsub(pattern="o__Xanthomonadales",
                        replacement="  Xanthomonadales",
                        sigtab.D1$Rank4)


p9 <- ggplot(sigtab.D1, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p9

jpeg("figures/deseq_ASV_PlatyD1.jpg",height=4,width=4,units="in", res=300)
p9 # Platy D1
dev.off()


# Convert phyloseq object to deseq object
phyASV.mic.deseq.PlatyC3 <- phyloseq_to_deseq2(phylo, ~C3_dom)
# Calculate geometric means
geoMeans.C3 <- apply(counts(phyASV.mic.deseq.PlatyC3), 1, gm_mean)
# Estimate size factors
phyASV.mic.deseq.PlatyC3 <- estimateSizeFactors(phyASV.mic.deseq.PlatyC3, geoMeans = geoMeans.C3)
# Run deseq analysis
phyASV.mic.deseq.PlatyC3 <- DESeq(phyASV.mic.deseq.PlatyC3, fitType = "local")
# Filter results based on alpha and specific contrast
res.C3 <- results(phyASV.mic.deseq.PlatyC3, cooksCutoff = FALSE, alpha = 0.05,contrast=c("C3_dom","Y","N"))
alpha.C3 <- 0.05
sigtab.C3 <- res.C3[which(res.C3$padj < alpha.C3), ]
sigtab.C3 <- cbind(as(sigtab.C3, "data.frame"), as(tax_table(phylo)[rownames(sigtab.C3), ], "matrix"))
head(sigtab.C3)


# Plot differential abundances
x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank2, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank2 <- factor(as.character(sigtab.C3$Rank2), levels=names(x.C3))

x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank4, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank4 <- factor(as.character(sigtab.C3$Rank4), levels=names(x.C3))

x.C3 <- tapply(sigtab.C3$log2FoldChange, sigtab.C3$Rank5, function(x) max(x))
x.C3 <- sort(x.C3, TRUE)
sigtab.C3$Rank5 <- factor(as.character(sigtab.C3$Rank5), levels=names(x.C3))

sigtab.C3$Rank2 <- gsub(pattern="p__",replacement="",
                        sigtab.C3$Rank2)

sigtab.C3$Rank4 <- gsub(pattern="o__Xanthomonadales",
                        replacement="  Xanthomonadales",
                        sigtab.C3$Rank4)

p10 <- ggplot(sigtab.C3, aes(x=Rank4, y=log2FoldChange, color=Rank2)) + 
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.8,0.75)) +
  labs(x = "Order", y = "Log2 Fold Change", color = "")
p10

jpeg("figures/deseq_ASV_PlatyC3.jpg",height=4,width=4,units="in", res=300)
p10 # Platy C3
dev.off()

# Specify which ancestral lineage it is
sigtab.C1$lin <- "C1"
sigtab.C3$lin <-"C3"
sigtab.C15$lin <- "C15"
sigtab.C31$lin <- "C31"
sigtab.C42$lin <- "C42"
sigtab.D1$lin <- "D1"

# Merge all sigtabs
sigtab.allsp.all <- rbind(sigtab.C1,sigtab.C3,sigtab.C15,
                          sigtab.C31,sigtab.C42,sigtab.D1)
# Clean up for plotting
sigtab.allsp.all$Rank2 <- gsub(pattern="p__",replacement="",
                               sigtab.allsp.all$Rank2)
sigtab.allsp.all$Rank4 <- gsub(pattern="o__",replacement="",
                               sigtab.allsp.all$Rank4)

p11 <- ggplot(sigtab.allsp.all, aes(x=Rank4, y=log2FoldChange, 
                             color=Rank2,fill=Rank2)) + 
  geom_point(size=4) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.5),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        legend.key = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        legend.position = c(0.1,0.78),
        axis.ticks.x=element_blank(),
        panel.spacing.x = unit(0.1,"lines")) +
  labs(x = "Order", y = "Log2 Fold Change", color = "") +
  facet_grid(facets = .~lin, drop=TRUE, space= "free",scales="free") +
  geom_hline(yintercept = 0,linetype="dashed")+
  guides(shape=FALSE) + guides(fill=FALSE)+
  scale_color_manual(values=c("#edae49","#2e4057","gray52"))+
  scale_fill_manual(values=c("#edae49","#2e4057","gray52"))
p11


jpeg("figures/Figure_6/Fig_6_new.jpg",height=4.4,width=9.9,units="in", res=300)
p11
dev.off()

pdf("figures/Figure_6/Fig_6_new.pdf",height=4.4,width=9.9,useDingbats = FALSE)
p11
dev.off()


fpentaC3 <- p8 +
  ylim(c(-30,30))+
  scale_color_manual(values = c("#00798c","#edae49","gray52","#d1495b"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.5),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(0.77,0.77),
        panel.border = element_rect(fill=NA))+
  geom_hline(yintercept = 0,linetype="dashed")

fpentaD1 <- p7 +
  ylim(c(-30,30)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill=NA)) +
  scale_color_manual(values = c("#00798c","#edae49","gray52","#d1495b"))+
  geom_hline(yintercept = 0,linetype="dashed")
platyC3 <- p10 +
  ylim(c(-30,30))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill=NA))+
  scale_color_manual(values = c("gray52"))+
  geom_hline(yintercept = 0,linetype="dashed")
platyD1 <- p9 +
  ylim(c(-30,30))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill=NA))+
  scale_color_manual(values = c("gray52"))+
  geom_hline(yintercept = 0,linetype="dashed")

grid.arrange(fpentaC3,fpentaD1,platyC3,platyD1, 
             nrow = 1, widths=c(5,5,1,1))

jpeg("figures/Figure_7/Fig_7.jpg",height=4,width=9,units="in", res=300)
grid.arrange(fpentaC3,fpentaD1,platyC3,platyD1, 
             nrow = 1, widths=c(5,5,1.5,1.5))
dev.off()

pdf("figures/Figure_7/Fig_7.pdf",height=4,width=9,useDingbats = FALSE)
grid.arrange(fpentaC3,fpentaD1,platyC3,platyD1, 
             nrow = 1, widths=c(5,5,1.5,1.5))
dev.off()

# save.image("analyses/16S/deseq_dl.RData")
