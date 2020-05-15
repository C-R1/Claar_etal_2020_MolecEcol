# Load necessary libraries
library(phyloseq)
library(biomformat)

# List necessary files
biom_file1 <- "data/16S_qiita_biom/reference-hit-3982-90bp_taxa.biom"
map_file1 <- "data/16S_qiita_biom/11244_prep_3982_qiime_20171026-180525.txt"
ref_file1 <- "data/16S_qiita_biom/reference-hit.seqs-3982-90bp.fa"

# Import biom
phy <- import_biom(BIOMfilename = biom_file1,
                   refseqfilename = ref_file1)

# Import sample data
sampledata0 <- read.delim(map_file1,stringsAsFactors=FALSE,row.names="X.SampleID")
sampledata1 <- sample_data(sampledata0)

# Merge taxa table and sample data and subset
phy.f1.nt <- merge_phyloseq(phy,sampledata1)
phy.f1.nt <- subset_samples(physeq = phy.f1.nt,sample_data(phy.f1.nt)$sample_type=="Whole Coral")
phy.f1.nt <- subset_samples(physeq = phy.f1.nt,sample_data(phy.f1.nt)$sampling_expedition=="KI14")

### Second set of files
biom_file2 <- "data/16S_qiita_biom/reference-hit-3983-90bp_taxa.biom"
map_file2 <- "data/16S_qiita_biom/11244_prep_3983_qiime_20171026-180530.txt"
ref_file2 <- "data/16S_qiita_biom/reference-hit.seqs-3983-90bp.fa"

phy <- import_biom(BIOMfilename = biom_file2,
                   refseqfilename = ref_file2)

sampledata0 <- read.delim(map_file2,stringsAsFactors=FALSE,row.names="X.SampleID")
sampledata1 <- sample_data(sampledata0)

phy.f2.nt <- merge_phyloseq(phy,sampledata1)
phy.f2.nt <- subset_samples(physeq = phy.f2.nt,sample_data(phy.f2.nt)$sample_type=="Whole Coral")
phy.f2.nt <- subset_samples(physeq = phy.f2.nt,sample_data(phy.f2.nt)$sampling_expedition=="KI14")
phy.f2.nt

### Third set of files
biom_file3 <- "data/16S_qiita_biom/reference-hit-4005-90bp_taxa.biom"
map_file3 <- "data/16S_qiita_biom/11244_prep_4005_qiime_20171026-181747.txt"
ref_file3 <- "data/16S_qiita_biom/reference-hit.seqs-4005-90bp.fa"

phy <- import_biom(BIOMfilename = biom_file3,
                   refseqfilename = ref_file3)

sampledata0 <- read.delim(map_file3,stringsAsFactors=FALSE,row.names="X.SampleID")
sampledata1 <- sample_data(sampledata0)

phy.f3.nt <- merge_phyloseq(phy,sampledata1)
phy.f3.nt <- subset_samples(physeq = phy.f3.nt,sample_data(phy.f3.nt)$sample_type=="Whole Coral")
phy.f3.nt <- subset_samples(physeq = phy.f3.nt,sample_data(phy.f3.nt)$sampling_expedition=="KI14")
phy.f3.nt

# Merge all three imported phyloseq objects
phy.f0 <- merge_phyloseq(phy.f1.nt,phy.f2.nt,phy.f3.nt)

# Read in tree for all seqs
all_tree <- read_tree("data/16S_qiita_biom/reference-hit.seqs-all-90bpdrep.tre")
# Add tree to the phyloseq object
phy.f <- merge_phyloseq(phy.f0,all_tree)
# Remove taxa with zero sequences
phy.f <- subset_taxa(phy.f,taxa_sums(phy.f) > 0)

# Characterize sites by disturbance level
sample_data(phy.f)$site <- sample_data(phy.f)$reef_name
VeryHigh <- c("Site_27","Site_30","Site_31","Site_32")
High <- c("Site_1","Site_6","Site_25","Site_26","Site_38","Site_40")
Medium <- c("Site_7","Site_8.5","Site_8","Site_12","Site_13","Site_14","Site_22","Site_33","Site_34","Site_35")
Low <- c("Site_2","Site_3","Site_4","Site_9","Site_23","Site_24")
VeryLow <- c("Site_5","Site_10","Site_11","Site_15","Site_16","Site_17","Site_18","Site_19","Site_20","Site_21","Site_28","Site_29","Site_36","Site_37","Site_39")

# Make site a factor
sample_data(phy.f)$site <- factor(sample_data(phy.f)$site)
# Create a disturbance level column, start as site
sample_data(phy.f)$human_disturbance <- sample_data(phy.f)$site

# Finish creating human_disturbance
for (i in VeryHigh){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$human_disturbance<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryHigh"), as.character(data.frame(sample_data(phy.f))$human_disturbance), as.character("VeryHigh")))
}

for (i in High){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$human_disturbance<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("High"), as.character(data.frame(sample_data(phy.f))$human_disturbance), as.character("High")))
}

for (i in Medium){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$human_disturbance<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Medium"), as.character(data.frame(sample_data(phy.f))$human_disturbance), as.character("Medium")))
}

for (i in Low){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$human_disturbance<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("Low"), as.character(data.frame(sample_data(phy.f))$human_disturbance), as.character("Low")))
}

for (i in VeryLow){
  if(i %in% unique(as.character(data.frame(sample_data(phy.f))$site)))
    (sample_data(phy.f)$human_disturbance<- gsub (paste("\\<",as.character(i),"\\>",sep=""),as.character("VeryLow"), as.character(data.frame(sample_data(phy.f))$human_disturbance), as.character("VeryLow")))
}

# Order levels
levels(sample_data(phy.f)$human_disturbance) <- c("VeryHigh","High","Medium","Low","VeryLow")

# Look at phy.f
phy.f
# Set threshold number
n <- 200
# Subset by minimum number of seqs
phy.f2 <- subset_samples(phy.f,sample_sums(phy.f) > n)
# Look at filtered physeq object
phy.f2


phy.f.f <- subset_samples(phy.f2,sample_names(phy.f2)!="11244.KI14FMD049") # Remove this sample
sample_data(phy.f.f)$human_disturbance <- factor(sample_data(phy.f.f)$human_disturbance,levels=c("Low","Medium","High","VeryHigh")) # Reorder disturbance levels

# Remove site 8.5 (shallow site)
phy.f.f <- subset_samples(phy.f.f,sample_data(phy.f.f)$reef_name!="Site_8.5")

# Save 
save(list=c("phy.f.f"), file="analyses/16S/KI_Sym_Microbes_16S.RData")
