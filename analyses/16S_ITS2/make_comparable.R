# Load necessary libraries
library(phyloseq)

# Load necessary data
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")
load("figures/KI_Sym_Microbes_colors.RData")
load("analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData")

# phy.f.f is microbial data
# phyASV.f.c.p is Symbiodinium data

# # Transform sample counts to relative abundances (McMurdie and Holmes 2014)
# phy.f.f <- transform_sample_counts(phy.f.f,function(x) x / sum(x))

# Make a new column ("ID") with tag number for Microbes data
sample_data(phy.f.f)$ID <- paste("tag_",sample_data(phy.f.f)$local_colony_name,sep="")
# Remove any re-tag information, to ensure that Symbiodinium names are comparable to microbe names 
sample_data(phyASV.f.c.p)$CoralTag <- gsub("_.*","",sample_data(phyASV.f.c.p)$CoralTag)
# Make a new column ("ID") with tag number for Symbiodinium data
sample_data(phyASV.f.c.p)$ID <- paste("tag_",sample_data(phyASV.f.c.p)$CoralTag,sep="")

# Remove samples that we have sequenced twice
phyASV.f.c.p <- subset_samples(phyASV.f.c.p, InputFileName != "KI14FF164.fasta")
phyASV.f.c.p <- subset_samples(phyASV.f.c.p, InputFileName != "KI14FF109.fasta")
phyASV.f.c.p <- subset_samples(phyASV.f.c.p, InputFileName != "KI14FF093.fasta")

# Rename samples using common ID
sample_names(phy.f.f) <- sample_data(phy.f.f)$ID
sample_names(phyASV.f.c.p) <- sample_data(phyASV.f.c.p)$ID

# Keep only Symbiodinium samples that we also have microbe samples for
phyASV.sym <- subset_samples(phyASV.f.c.p, sample_names(phyASV.f.c.p) %in% sample_names(phy.f.f))
# Remove any empty taxa from phyloseq object
phyASV.sym <- subset_taxa(phyASV.sym, taxa_sums(phyASV.sym)>0)

# Keep only microbe samples that we also have Symbiodinium samples for
phyASV.mic <- subset_samples(phy.f.f, sample_names(phy.f.f) %in% sample_names(phyASV.f.c.p))
# Remove any empty taxa from phyloseq object
phyASV.mic <- subset_taxa(phyASV.mic, taxa_sums(phyASV.mic)>0)

#################
tax <- data.frame(tax_table(phyASV.sym))
otu <- data.frame(otu_table(phyASV.sym))
sam <- data.frame(sample_data(phyASV.sym))
mic.sam <- data.frame(sample_data(phyASV.mic))

# Determine whether each sample has any clade A (now Symbiodinium sensu stricto)
A_tax <- rownames(tax)[tax$Genus=="g__Symbiodinium"] # Find otus that are in clade A
A_Y <- otu[(colnames(otu) %in% A_tax)] # Make a list of otus that are clade A
A_Y_sums <- rowSums(A_Y) # Sum instances of clade A being present (a coral may have more than one OTU that is clade A)
A_Y_names <- names(A_Y_sums[A_Y_sums>0]) # Extract the names of samples that have at least some clade A in them
mic.sam$clade_A <- "N" # Default is no clade A
mic.sam$clade_A[rownames(mic.sam) %in% A_Y_names] <- "Y" # Replace as "Y" all samples with clade A in them 

# Determine whether each sample has any clade B (now Breviolum)
B_tax <- rownames(tax)[tax$Genus=="g__Breviolum"]
B_Y <- otu[(colnames(otu) %in% B_tax)]
B_Y_sums <- rowSums(B_Y)
B_Y_names <- names(B_Y_sums[B_Y_sums>0])
mic.sam$clade_B <- "N"
mic.sam$clade_B[rownames(mic.sam) %in% B_Y_names] <- "Y"

# Determine whether each sample has any clade D (now Durusdinium)
D_tax <- rownames(tax)[tax$Genus=="g__Durusdinium"]
D_Y <- otu[(colnames(otu) %in% D_tax)]
D_Y_sums <- rowSums(D_Y)
D_Y_names <- names(D_Y_sums[D_Y_sums>0])
mic.sam$clade_D <- "N"
mic.sam$clade_D[rownames(mic.sam) %in% D_Y_names] <- "Y"

# Determine whether each sample has any clade F (now Fugacium)
F_tax <- rownames(tax)[tax$Genus=="g__Fugacium"]
F_Y <- otu[(colnames(otu) %in% F_tax)]
F_Y_sums <- rowSums(F_Y)
F_Y_names <- names(F_Y_sums[F_Y_sums>0])
mic.sam$clade_F <- "N"
mic.sam$clade_F[rownames(mic.sam) %in% F_Y_names] <- "Y"

# Determine whether each sample has any clade G (now Gerakladium)
G_tax <- rownames(tax)[tax$Genus=="g__Gerakladium"]
G_Y <- otu[(colnames(otu) %in% G_tax)]
G_Y_sums <- rowSums(G_Y)
G_Y_names <- names(G_Y_sums[G_Y_sums>0])
mic.sam$clade_G <- "N"
mic.sam$clade_G[rownames(mic.sam) %in% G_Y_names] <- "Y"

mic.sam
sam_asv_ids <- sam[c("dominant_asv","dominant_asv_id")]

test <- merge(mic.sam, sam_asv_ids,by="row.names")
rownames(test)<-rownames(mic.sam)

test$C15_dom <- test$dominant_asv_id
test$C15_dom[test$C15_dom!="C15"] <- "N"
test$C15_dom[test$C15_dom=="C15"] <- "Y"
test$C15_dom <- unlist(test$C15_dom)

test$C31_dom <- test$dominant_asv_id
test$C31_dom[test$C31_dom!="C31"] <- "N"
test$C31_dom[test$C31_dom=="C31"] <- "Y"
test$C31_dom <- unlist(test$C31_dom)

test$C3_dom <- test$dominant_asv_id
test$C3_dom[test$C3_dom!="C3"] <- "N"
test$C3_dom[test$C3_dom=="C3"] <- "Y"
test$C3_dom <- unlist(test$C3_dom)

test$C1_dom <- test$dominant_asv_id
test$C1_dom[test$C1_dom!="C1"] <- "N"
test$C1_dom[test$C1_dom=="C1"] <- "Y"
test$C1_dom <- unlist(test$C1_dom)

test$C42_dom <- test$dominant_asv_id
test$C42_dom[test$C42_dom!="C42"] <- "N"
test$C42_dom[test$C42_dom=="C42"] <- "Y"
test$C42_dom <- unlist(test$C42_dom)

test$D1_dom <- test$dominant_asv_id
test$D1_dom[test$D1_dom!="D1"] <- "N"
test$D1_dom[test$D1_dom=="D1"] <- "Y"
test$D1_dom <- unlist(test$D1_dom)

# Replace phyloseq sample data with original sample data + clade Y/N appended
sample_data(phyASV.mic) <- test

############## Subset samples by coral species
phyASV.sym.Pocillopora <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Pocillopora grandis")
phyASV.sym.Pocillopora <- subset_taxa(phyASV.sym.Pocillopora,taxa_sums(phyASV.sym.Pocillopora)>0)
phyASV.sym.Porites <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Porites lobata")
phyASV.sym.Porites <- subset_taxa(phyASV.sym.Porites,taxa_sums(phyASV.sym.Porites)>0)
phyASV.sym.Montipora <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Montipora aequituberculata")
phyASV.sym.Montipora <- subset_taxa(phyASV.sym.Montipora,taxa_sums(phyASV.sym.Montipora)>0)
phyASV.sym.FaviaM <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Favia matthaii")
phyASV.sym.FaviaM <- subset_taxa(phyASV.sym.FaviaM,taxa_sums(phyASV.sym.FaviaM)>0)
phyASV.sym.Hydnophora <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Hydnophora microconos")
phyASV.sym.Hydnophora <- subset_taxa(phyASV.sym.Hydnophora,taxa_sums(phyASV.sym.Hydnophora)>0)
phyASV.sym.Platygyra <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Platygyra rykyuensis")
phyASV.sym.Platygyra <- subset_taxa(phyASV.sym.Platygyra,taxa_sums(phyASV.sym.Platygyra)>0)
phyASV.sym.Favites <- subset_samples(phyASV.sym,sample_data(phyASV.sym)$Coral_Species =="Favites pentagona")
phyASV.sym.Favites <- subset_taxa(phyASV.sym.Favites,taxa_sums(phyASV.sym.Favites)>0)


phyASV.mic.Pocillopora <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Pocillopora grandis")
phyASV.mic.Pocillopora <- subset_taxa(phyASV.mic.Pocillopora,taxa_sums(phyASV.mic.Pocillopora)>0)
phyASV.mic.Porites <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Porites lobata")
phyASV.mic.Porites <- subset_taxa(phyASV.mic.Porites,taxa_sums(phyASV.mic.Porites)>0)
phyASV.mic.Montipora <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Montipora aequituberculata")
phyASV.mic.Montipora <- subset_taxa(phyASV.mic.Montipora,taxa_sums(phyASV.mic.Montipora)>0)
phyASV.mic.FaviaM <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Favia matthaii")
phyASV.mic.FaviaM <- subset_taxa(phyASV.mic.FaviaM,taxa_sums(phyASV.mic.FaviaM)>0)
phyASV.mic.Hydnophora <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Hydnophora microconos")
phyASV.mic.Hydnophora <- subset_taxa(phyASV.mic.Hydnophora,taxa_sums(phyASV.mic.Hydnophora)>0)
phyASV.mic.Platygyra <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Platygyra daedalea")
phyASV.mic.Platygyra <- subset_taxa(phyASV.mic.Platygyra,taxa_sums(phyASV.mic.Platygyra)>0)
phyASV.mic.Favites <- subset_samples(phyASV.mic,sample_data(phyASV.mic)$field_host_name =="Favites pentagona")
phyASV.mic.Favites <- subset_taxa(phyASV.mic.Favites,taxa_sums(phyASV.mic.Favites)>0)

phyASV.mic.raw <- phyASV.mic
phyASV.mic.Pocillopora.raw <- phyASV.mic.Pocillopora
phyASV.mic.Porites.raw <- phyASV.mic.Porites
phyASV.mic.Montipora.raw <- phyASV.mic.Montipora
phyASV.mic.FaviaM.raw <- phyASV.mic.FaviaM
phyASV.mic.Hydnophora.raw <- phyASV.mic.Hydnophora
phyASV.mic.Platygyra.raw <- phyASV.mic.Platygyra
phyASV.mic.Favites.raw <- phyASV.mic.Favites

# Transform sample counts to relative abundances (McMurdie and Holmes 2014)
phyASV.mic <- transform_sample_counts(phyASV.mic,function(x) x / sum(x))
phyASV.mic.Pocillopora <- transform_sample_counts(phyASV.mic.Pocillopora,function(x) x / sum(x))
phyASV.mic.Porites <- transform_sample_counts(phyASV.mic.Porites,function(x) x / sum(x))
phyASV.mic.Montipora <- transform_sample_counts(phyASV.mic.Montipora,function(x) x / sum(x))
phyASV.mic.FaviaM <- transform_sample_counts(phyASV.mic.FaviaM,function(x) x / sum(x))
phyASV.mic.Hydnophora <- transform_sample_counts(phyASV.mic.Hydnophora,function(x) x / sum(x))
phyASV.mic.Platygyra <- transform_sample_counts(phyASV.mic.Platygyra,function(x) x / sum(x))
phyASV.mic.Favites <- transform_sample_counts(phyASV.mic.Favites,function(x) x / sum(x))

# Save final phyloseqs for downstream use
save(list=c("phyASV.mic","phyASV.sym",
            "phyASV.sym.Pocillopora","phyASV.sym.Porites",
            "phyASV.sym.Montipora","phyASV.sym.FaviaM",
            "phyASV.sym.Hydnophora","phyASV.sym.Platygyra","phyASV.sym.Favites",
            "phyASV.mic.Pocillopora","phyASV.mic.Porites",
            "phyASV.mic.Montipora","phyASV.mic.FaviaM",
            "phyASV.mic.Hydnophora","phyASV.mic.Platygyra","phyASV.mic.Favites",
            "phyASV.mic.raw",
            "phyASV.mic.Pocillopora.raw","phyASV.mic.Porites.raw",
            "phyASV.mic.Montipora.raw","phyASV.mic.FaviaM.raw",
            "phyASV.mic.Hydnophora.raw","phyASV.mic.Platygyra.raw",
            "phyASV.mic.Favites.raw"),
     file = "analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData")
