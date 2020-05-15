# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData")
load("figures/KI_Sym_Microbes_colors.RData")

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

sam_data <- data.frame(sample_data(phyASV.f.c.p))

sam_data$dominant_asv_id <- as.factor(unlist(sam_data$dominant_asv_id))

sam_data %>% # Total number of colonies dominated by each Sym lineage
  count(dominant_asv_id)

# Number of colonies dominated by each Sym lineage BY SPECIES
sp_dom <- sam_data %>% 
  count(Coral_Species,dominant_asv_id) %>% 
  as.data.frame()

data.frame(tax_table(phyASV.f.c.p))[data.frame(tax_table(phyASV.f.c.p))$Genus=="g__Symbiodinium",]
unk_as<-c("ASV990","ASV993","ASV1579","ASV1600","ASV1764")
for(asv in unk_as){
  write.csv((refseq(phyASV.f.c.p)[asv]),file=paste0("analyses/ITS2/unk_As/refseq_",asv,".csv"))
}

data.frame(tax_table(phyASV.f.c.p))[data.frame(tax_table(phyASV.f.c.p))$Genus=="g__Breviolum",]

unique(data.frame(tax_table(phyASV.f.c.p))[data.frame(tax_table(phyASV.f.c.p))$Genus=="g__Durusdinium",])

data.frame(tax_table(phyASV.f.c.p))[data.frame(tax_table(phyASV.f.c.p))$Genus=="g__Fugacium",]
unk_fs<-c("ASV1581","ASV1645","ASV1995","ASV3578")
for(asv in unk_fs){
  write.csv((refseq(phyASV.f.c.p)[asv]),file=paste0("analyses/ITS2/unk_Fs/refseq_",asv,".csv"))
}

data.frame(tax_table(phyASV.f.c.p))[data.frame(tax_table(phyASV.f.c.p))$Genus=="g__Gerakladium",]
unk_gs<-c("ASV1812","ASV2181","ASV2209","ASV2212","ASV2649")
for(asv in unk_gs){
  write.csv((refseq(phyASV.f.c.p)[asv]),file=paste0("analyses/ITS2/unk_Gs/refseq_",asv,".csv"))
}
