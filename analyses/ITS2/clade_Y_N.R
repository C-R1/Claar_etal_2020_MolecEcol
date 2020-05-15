# Load necessary libraries
library(phyloseq)

# Load necessary data
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData")
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")

# Extract each component of the phyloseq object for tree building
tax <- data.frame(tax_table(phyASV.f.c.p))
otu <- data.frame(otu_table(phyASV.f.c.p))
sam <- data.frame(sample_data(phyASV.f.c.p))
  
# Determine whether each sample has any clade A (now Symbiodinium sensu stricto)
A_tax <- rownames(tax)[tax$Genus=="g__Symbiodinium"] # Find otus that are in clade A
A_Y <- otu[(colnames(otu) %in% A_tax)] # Make a list of otus that are clade A
A_Y_sums <- rowSums(A_Y) # Sum instances of clade A being present (a coral may have more than one OTU that is clade A)
A_Y_names <- names(A_Y_sums[A_Y_sums>0]) # Extract the names of samples that have at least some clade A in them
sam$clade_A <- "N" # Default is no clade A
sam$clade_A[rownames(sam) %in% A_Y_names] <- "Y" # Replace as "Y" all samples with clade A in them 

# Determine whether each sample has any clade B (now Breviolum)
B_tax <- rownames(tax)[tax$Genus=="g__Breviolum"]
B_Y <- otu[(colnames(otu) %in% B_tax)]
B_Y_sums <- rowSums(B_Y)
B_Y_names <- names(B_Y_sums[B_Y_sums>0])
sam$clade_B <- "N"
sam$clade_B[rownames(sam) %in% B_Y_names] <- "Y"

# Determine whether each sample has any clade C (now Cladocopium)
C_tax <- rownames(tax)[tax$Genus=="g__Cladocopium"]
C_Y <- otu[(colnames(otu) %in% C_tax)]
C_Y_sums <- rowSums(C_Y)
C_Y_names <- names(C_Y_sums[C_Y_sums>0])
sam$clade_C <- "N"
sam$clade_C[rownames(sam) %in% C_Y_names] <- "Y"

# Determine whether each sample has any clade D (now Durusdinium)
D_tax <- rownames(tax)[tax$Genus=="g__Durusdinium"]
D_Y <- otu[(colnames(otu) %in% D_tax)]
D_Y_sums <- rowSums(D_Y)
D_Y_names <- names(D_Y_sums[D_Y_sums>0])
sam$clade_D <- "N"
sam$clade_D[rownames(sam) %in% D_Y_names] <- "Y"

# Determine whether each sample has any clade F (now Fugacium)
F_tax <- rownames(tax)[tax$Genus=="g__Fugacium"]
F_Y <- otu[(colnames(otu) %in% F_tax)]
F_Y_sums <- rowSums(F_Y)
F_Y_names <- names(F_Y_sums[F_Y_sums>0])
sam$clade_F <- "N"
sam$clade_F[rownames(sam) %in% F_Y_names] <- "Y"

# Determine whether each sample has any clade G (now Gerakladium)
G_tax <- rownames(tax)[tax$Genus=="g__Gerakladium"]
G_Y <- otu[(colnames(otu) %in% G_tax)]
G_Y_sums <- rowSums(G_Y)
G_Y_names <- names(G_Y_sums[G_Y_sums>0])
sam$clade_G <- "N"
sam$clade_G[rownames(sam) %in% G_Y_names] <- "Y"

# Replace phyloseq sample data with original sample data + clade Y/N appended
sample_data(phyASV.f.c.p) <- sam
sample_data(phyASV.f.c.p) # Double check

#############
# Calculate dominant ASVs

# https://github.com/joey711/phyloseq/issues/847
# From https://github.com/cyklee
find.top.asv <- function(x,num){
  require(phyloseq)
  require(magrittr)
  
  otu <- otu_table(x)
  tax <- tax_table(x)
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"ix") # select for index
  
  m <- data.frame(otu@.Data[,unique(unlist(j2))])
  n <- apply(m,1,sort,index.return=T, decreasing=T) %>%
    lapply('[[',"ix") %>%  # Extract index
    lapply(head,n=num) # This to returns the top x tax
  
  p <- list()
  for(i in 1:length(n)){
    p[[i]]<- colnames(m)[n[[i]]]
  }
  m$taxa <- p
  return(m)
}

top_asvs_all <- find.top.asv(phyASV.f.c.p, 1) # Top ASVs per sample

top_asvs <- unique(top_asvs_all$taxa)

for(asv in top_asvs){
  print(tax_table(phyASV.f.c.p)[asv])
}

for(asv in top_asvs){
  write.csv((refseq(phyASV.f.c.p)[asv]),file=paste0("analyses/ITS2/dominant_seqs/refseq_",asv,".csv"))
}

sample_data(phyASV.f.c.p)$dominant_asv <- top_asvs_all$taxa
sample_data(phyASV.f.c.p)$dominant_asv_id <- top_asvs_all$taxa

sample_data(phyASV.f.c.p)$dominant_asv_id[sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV5"] <- "C42"

sample_data(phyASV.f.c.p)$dominant_asv_id[sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV1"] <- "C31"

sample_data(phyASV.f.c.p)$dominant_asv_id[sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV6"] <- "D1"

sample_data(phyASV.f.c.p)$dominant_asv_id[sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV2" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV25" ] <- "C1"

sample_data(phyASV.f.c.p)$dominant_asv_id[sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV3" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV10" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV58" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV162"] <- "C3"

sample_data(phyASV.f.c.p)$dominant_asv_id[sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV4" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV9" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV29" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV34" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV153" | sample_data(phyASV.f.c.p)$dominant_asv_id == "ASV155"] <- "C15"


sample_data(phyASV.f.c.p)[,c("dominant_asv","dominant_asv_id")]

# Save this file for downstream processing
save(phyASV.f.c.p,file="analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData")
