# Import Libraries
library(stringr)
library(reshape2)
library(phyloseq)
library(seqinr)
library(phangorn) 
library(caroline)
library(DECIPHER)

# Load necessary data
load("analyses/ITS2/dada2/KI_SymMicrobes_dada.RData")

# Rename phyloseq object
phy.f <- ps

########################## Site Formatting ################################

# Characterize sites by disturbance level
VeryHigh <- c(30,31,32,27)
High <- c(1,6,26,40)
Medium <- c(7,8,12,13,14,22,25,33,34,35)
Low <- c(2,3,4,9,23,24,38)
VeryLow <- c(5,10,11,15,16,17,18,19,20,21,28,29,36,37,39)

sample_data(phy.f)$Site <- factor(sample_data(phy.f)$Site)
sample_data(phy.f)$Dist <- sample_data(phy.f)$Site

for (i in VeryHigh){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$Site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryHigh"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("VeryHigh")))
}

for (i in High){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$Site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("High"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("High")))
}

for (i in Medium){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$Site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Medium"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("Medium")))
}

for (i in Low){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$Site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Low"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("Low")))
}

for (i in VeryLow){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$Site)))
    (sample_data(phy.f)$Dist<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryLow"), as.character(data.frame(sample_data(phy.f))$Dist), as.character("VeryLow")))
}

levels(sample_data(phy.f)$Dist) <- c("VeryHigh","High","Medium","Low","VeryLow")

###################### Physeq formatting and tree #####################################
# Assign new name for clarity
phyASV.f.c <- phy.f
# Remove taxa that have zero reads
phyASV.f.c <- prune_taxa(taxa_sums(phyASV.f.c)>0,phyASV.f.c)

# Subset by Symbiodiniaceae genus
ASV.As <- subset_taxa(phyASV.f.c,Genus=="g__Symbiodinium")
ASV.Bs <- subset_taxa(phyASV.f.c,Genus=="g__Breviolum")
ASV.Cs <- subset_taxa(phyASV.f.c,Genus=="g__Cladocopium")
ASV.Ds <- subset_taxa(phyASV.f.c,Genus=="g__Durusdinium")
# ASV.Es <- subset_taxa(phyASV.f.c,Genus=="g__Effrenium") # no seqs
ASV.Fs <- subset_taxa(phyASV.f.c,Genus=="g__Fugacium")
ASV.Gs <- subset_taxa(phyASV.f.c,Genus=="g__Gerakladium")
# ASV.Hs <- subset_taxa(phyASV.f.c,Genus=="g__cladeH") # no seqs
# ASV.Is <- subset_taxa(phyASV.f.c,Genus=="g__cladeI") # no seqs

# Extract sequences and check for max and min lengths
ASV.seqs <- refseq(phyASV.f.c)
min(width(ASV.seqs))
max(width(ASV.seqs))

# Align sequences within each Symbiodiniaceae genus
ASV.A.seqs <- refseq(ASV.As)
writeXStringSet(ASV.A.seqs, # write to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_A_tree_seqs.fasta")
ASV.A.seqs.aligned <- AlignSeqs(ASV.A.seqs)
writeXStringSet(ASV.A.seqs.aligned, # write the alignment to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_A_tree_seqs_aligned.fasta")
ASV.B.seqs <- refseq(ASV.Bs)
writeXStringSet(ASV.B.seqs, # write to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_B_tree_seqs.fasta")
ASV.B.seqs.aligned <- refseq(ASV.Bs) # There is only one sequence here! Can't align
writeXStringSet(ASV.B.seqs.aligned, # write the alignment to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_B_tree_seqs_aligned.fasta")
ASV.C.seqs <- refseq(ASV.Cs)
writeXStringSet(ASV.C.seqs, # write to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_C_tree_seqs.fasta")
ASV.C.seqs.aligned <- AlignSeqs(ASV.C.seqs)
writeXStringSet(ASV.C.seqs.aligned, # write the alignment to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_C_tree_seqs_aligned.fasta")
ASV.D.seqs <- refseq(ASV.Ds)
writeXStringSet(ASV.D.seqs, # write to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_D_tree_seqs.fasta")
ASV.D.seqs.aligned <- AlignSeqs(ASV.D.seqs)
writeXStringSet(ASV.D.seqs.aligned, # write the alignment to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_D_tree_seqs_aligned.fasta")
ASV.F.seqs <- refseq(ASV.Fs)
writeXStringSet(ASV.F.seqs, # write to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_F_tree_seqs.fasta")
ASV.F.seqs.aligned <- AlignSeqs(ASV.F.seqs)
writeXStringSet(ASV.F.seqs.aligned, # write the alignment to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_F_tree_seqs_aligned.fasta")
ASV.G.seqs <- refseq(ASV.Gs)
writeXStringSet(ASV.G.seqs, # write to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_G_tree_seqs.fasta")
ASV.G.seqs.aligned <- AlignSeqs(ASV.G.seqs)
writeXStringSet(ASV.G.seqs.aligned, # write the alignment to a new FASTA file
                file="data/ITS2/Bioinf/tree/ASV_G_tree_seqs_aligned.fasta")

writeXStringSet(c(ASV.A.seqs.aligned,ASV.B.seqs.aligned,ASV.C.seqs.aligned,ASV.D.seqs.aligned,ASV.F.seqs.aligned,ASV.G.seqs.aligned), 
                file="data/ITS2/Bioinf/tree/ASV_ALL_tree_seqs_aligned.fasta")

#https://rdrr.io/rforge/seqinr/man/dist.alignment.html
#returns sqrt of pairwise genetic distance, then squared the matrices
A.seqs <- read.alignment(file = "data/ITS2/Bioinf/tree/ASV_A_tree_seqs_aligned.fasta", format= "fasta")
A.dis <- (as.matrix(dist.alignment(A.seqs, matrix = "identity" )))^2
write.csv(A.dis, file="data/ITS2/Bioinf/tree/ASV_A_dis_matx.csv")

B.seqs <- read.alignment(file = "data/ITS2/Bioinf/tree/ASV_B_tree_seqs_aligned.fasta", format= "fasta")
B.dis <- (as.matrix(dist.alignment(B.seqs, matrix = "identity" )))^2
write.csv(B.dis, file="data/ITS2/Bioinf/tree/ASV_B_dis_matx.csv")

C.seqs <- read.alignment(file = "data/ITS2/Bioinf/tree/ASV_C_tree_seqs_aligned.fasta", format= "fasta")
C.dis <- (as.matrix(dist.alignment(C.seqs, matrix = "identity" )))^2
write.csv(C.dis, file="data/ITS2/Bioinf/tree/ASV_C_dis_matx.csv")

D.seqs <- read.alignment(file = "data/ITS2/Bioinf/tree/ASV_D_tree_seqs_aligned.fasta", format= "fasta")
D.dis <- (as.matrix(dist.alignment(D.seqs, matrix = "identity" )))^2
write.csv(D.dis, file="data/ITS2/Bioinf/tree/ASV_D_dis_matx.csv")

F.seqs <- read.alignment(file = "data/ITS2/Bioinf/tree/ASV_F_tree_seqs_aligned.fasta", format= "fasta")
F.dis <- (as.matrix(dist.alignment(F.seqs, matrix = "identity" )))^2
write.csv(F.dis, file="data/ITS2/Bioinf/tree/ASV_F_dis_matx.csv")

G.seqs <- read.alignment(file = "data/ITS2/Bioinf/tree/ASV_G_tree_seqs_aligned.fasta", format= "fasta")
G.dis <- (as.matrix(dist.alignment(G.seqs, matrix = "identity" )))^2
write.csv(G.dis, file="data/ITS2/Bioinf/tree/ASV_G_dis_matx.csv")

#give clade distances using average 28s distance from Pochon and Gates 2010
A_B <- matrix(0.219, ncol=ncol(A.dis), nrow=nrow(B.dis), 
              dimnames=list(rownames(B.dis), colnames(A.dis)))
A_C <- matrix(0.1960, ncol=ncol(A.dis), nrow=nrow(C.dis), 
              dimnames=list(rownames(C.dis), colnames(A.dis)))
A_D <- matrix(0.1775, ncol=ncol(A.dis), nrow=nrow(D.dis), 
              dimnames=list(rownames(D.dis), colnames(A.dis)))
A_F <- matrix(0.2085, ncol=ncol(A.dis), nrow=nrow(F.dis), 
              dimnames=list(rownames(F.dis), colnames(A.dis)))
A_G <- matrix(0.216, ncol=ncol(A.dis), nrow=nrow(G.dis), 
              dimnames=list(rownames(G.dis), colnames(A.dis)))
B_C <- matrix(0.114, ncol=ncol(B.dis), nrow=nrow(C.dis), 
              dimnames=list(rownames(C.dis), colnames(B.dis)))
B_D <- matrix(0.1705, ncol=ncol(B.dis), nrow=nrow(D.dis), 
              dimnames=list(rownames(D.dis), colnames(B.dis)))
B_F <- matrix(0.1355, ncol=ncol(B.dis), nrow=nrow(F.dis), 
              dimnames=list(rownames(F.dis), colnames(B.dis)))
B_G <- matrix(0.21, ncol=ncol(B.dis), nrow=nrow(G.dis), 
              dimnames=list(rownames(G.dis), colnames(B.dis)))
C_D <- matrix(0.1520, ncol=ncol(C.dis), nrow=nrow(D.dis), 
              dimnames=list(rownames(D.dis), colnames(C.dis)))
C_F <- matrix(0.0945, ncol=ncol(C.dis), nrow=nrow(F.dis), 
              dimnames=list(rownames(F.dis), colnames(C.dis)))
C_G <- matrix(0.187, ncol=ncol(C.dis), nrow=nrow(G.dis), 
              dimnames=list(rownames(G.dis), colnames(C.dis)))
D_F <- matrix(0.1691, ncol=ncol(D.dis), nrow=nrow(F.dis), 
              dimnames=list(rownames(F.dis), colnames(D.dis)))
D_G <- matrix(0.1795, ncol=ncol(D.dis), nrow=nrow(G.dis), 
              dimnames=list(rownames(G.dis), colnames(D.dis)))
F_G <- matrix(0.2072, ncol=ncol(F.dis), nrow=nrow(G.dis), 
              dimnames=list(rownames(G.dis), colnames(F.dis)))

#build ACDG matrix
col1 <- rbind(A.dis, A_B, A_C, A_D, A_F, A_G)
col2 <- rbind(matrix(NA, nrow=nrow(A.dis), ncol=ncol(B.dis), 
                     dimnames=list(rownames(A.dis), colnames(B.dis))), 
              B.dis, B_C, B_D, B_F, B_G)
col3 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis), ncol=ncol(C.dis),
                     dimnames=list(c(rownames(A.dis), rownames(B.dis)),
                                   colnames(C.dis))), C.dis, C_D, C_F, C_G)
col4 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis), 
                     ncol=ncol(D.dis), dimnames=list(c(rownames(A.dis),
                                                       rownames(B.dis),
                                                       rownames(C.dis)),
                                                     colnames(D.dis))), 
              D.dis, D_F, D_G)
col5 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis)+nrow(D.dis),
                     ncol=ncol(F.dis), dimnames=list(c(rownames(A.dis),
                                                       rownames(B.dis),
                                                       rownames(C.dis),
                                                       rownames(D.dis)),
                                                     colnames(F.dis))), 
              F.dis, F_G)
col6 <- rbind(matrix(NA, nrow=nrow(A.dis)+nrow(B.dis)+nrow(C.dis)+nrow(D.dis)+nrow(F.dis), ncol=ncol(G.dis), dimnames=list(c(rownames(A.dis),                                                               rownames(B.dis), 
                                                 rownames(C.dis), 
                                                 rownames(D.dis),  
                                                 rownames(F.dis)), 
                                                 colnames(G.dis))), G.dis)

ubermatrix <- cbind(col1, col2, col3, col4, col5, col6)
dim(ubermatrix)

#build tree
uber.tree <- phangorn::upgma(ubermatrix)
plot(uber.tree, main="UPGMA")

#write tree to file
write.tree(uber.tree, file="data/ITS2/Bioinf/tree/uber.tre")

class(uber.tree)

# Slot uber tree into the phy_tree slot of the phyloseq object
phy_tree(phyASV.f.c) <- phy_tree(uber.tree)
plot_tree(phyASV.f.c, label.tips = "Genus", color = "Genus")

################# Fix and standardize coral species names #################
sample_data(phyASV.f.c)$Coral_Species <- gsub(pattern="_",replacement=" ", x = sample_data(phyASV.f.c)$Coral_Species)
sample_data(phyASV.f.c)$Coral_Species <- gsub(pattern='matthai$',replacement="matthaii", x = sample_data(phyASV.f.c)$Coral_Species)
sample_data(phyASV.f.c)$Coral_Species <- gsub(pattern='foliosa',replacement="aequituberculata", x = sample_data(phyASV.f.c)$Coral_Species)
sample_data(phyASV.f.c)$Coral_Species <- gsub(pattern='eydouxi',replacement="grandis", x = sample_data(phyASV.f.c)$Coral_Species)
sample_data(phyASV.f.c)$Coral_Species <- gsub(pattern='Platygyra sp',replacement="Platygyra rykyuensis", x = sample_data(phyASV.f.c)$Coral_Species)
sample_data(phyASV.f.c)$Coral_Species <- gsub(pattern='Favia sp',replacement="Favia matthaii", x = sample_data(phyASV.f.c)$Coral_Species)

############ Transform and Calc seqs, subset by coral species ############

# Remove samples with no sequences
phyASV.f.c <- prune_samples(sample_sums(phyASV.f.c)>0,phyASV.f.c) # removed 4 samples

# Transform sample counts to proportional abundance for downstream analyses
phyASV.f.c.p <- transform_sample_counts(phyASV.f.c, function(x) x/sum(x))

# Calculate number of sequences in phyASV.f.c
total_seqs <- sum(taxa_sums(phyASV.f.c))

# Subset coral species
phyASV.f.c.p.Pocillopora <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Pocillopora grandis")
phyASV.f.c.p.Pocillopora <- prune_taxa(taxa_sums(phyASV.f.c.p.Pocillopora)>0,phyASV.f.c.p.Pocillopora)

phyASV.f.c.p.Porites <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Porites lobata")
phyASV.f.c.p.Porites <- prune_taxa(taxa_sums(phyASV.f.c.p.Porites)>0,phyASV.f.c.p.Porites)

phyASV.f.c.p.Montipora <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Montipora aequituberculata")
phyASV.f.c.p.Montipora <- prune_taxa(taxa_sums(phyASV.f.c.p.Montipora)>0,phyASV.f.c.p.Montipora)

phyASV.f.c.p.FaviaM <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Favia matthaii")
phyASV.f.c.p.FaviaM <- prune_taxa(taxa_sums(phyASV.f.c.p.FaviaM)>0,phyASV.f.c.p.FaviaM)

phyASV.f.c.p.Hydnophora <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Hydnophora microconos")
phyASV.f.c.p.Hydnophora <- prune_taxa(taxa_sums(phyASV.f.c.p.Hydnophora)>0,phyASV.f.c.p.Hydnophora)

phyASV.f.c.p.Platygyra <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Platygyra rykyuensis")
phyASV.f.c.p.Platygyra <- prune_taxa(taxa_sums(phyASV.f.c.p.Platygyra)>0,phyASV.f.c.p.Platygyra)

phyASV.f.c.p.Favites <- subset_samples(phyASV.f.c.p,sample_data(phyASV.f.c.p)$Coral_Species=="Favites pentagona")
phyASV.f.c.p.Favites <- prune_taxa(taxa_sums(phyASV.f.c.p.Favites)>0,phyASV.f.c.p.Favites)

#################### Save grouped data as RData file ##########################
save(list=ls(),file="analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData")

save(list=c("phyASV.f.c.p.Pocillopora","phyASV.f.c.p.Porites",
            "phyASV.f.c.p.Montipora","phyASV.f.c.p.FaviaM",
            "phyASV.f.c.p.Hydnophora",
            "phyASV.f.c.p.Platygyra","phyASV.f.c.p.Favites"),
     file="analyses/ITS2/KI_Sym_Microbe_f_coral_grouped_byspecies.RData")
save(list=c("phyASV.f.c.p"),file="data/ITS2/data/KI_Sym_Microbe_phyASVfcp.RData")

