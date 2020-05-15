# Import Libraries
library(stringr)
library(reshape2)
library(phyloseq)
library(seqinr)
library(phangorn) 
library(caroline)

# Load filtered RData object 
load("analyses/16S/KI_Sym_Microbes_16S.RData")

### Fix site 38
data.frame(Dist = sample_data(phy.f.f)$human_disturbance, Site = sample_data(phy.f.f)$reef_name)
sample_data(phy.f.f)$human_disturbance[sample_data(phy.f.f)$reef_name == "Site_38"]<- "Low"

## Fix site 25
sample_data(phy.f.f)$human_disturbance[sample_data(phy.f.f)$reef_name == "Site_25"]<- "Medium"

##################### Fix and standardize coral species names #####################
sample_data(phy.f.f)$field_host_name <- gsub(pattern='matthai$',replacement="matthaii", x = sample_data(phy.f.f)$field_host_name)
sample_data(phy.f.f)$field_host_name <- gsub(pattern='foliosa',replacement="aequituberculata", x = sample_data(phy.f.f)$field_host_name)
sample_data(phy.f.f)$field_host_name <- gsub(pattern='eydouxi',replacement="grandis", x = sample_data(phy.f.f)$field_host_name)
sample_data(phy.f.f)$field_host_name <- gsub(pattern='Platygyra sp',replacement="Platygyra daedalea", x = sample_data(phy.f.f)$field_host_name)

###################### Remove mitochondria and chloroplasts ##############

# Check for unassigned at Kingdom level
get_taxa_unique(phy.f.f,"Rank1") # "Unassigned" should be removed
phy.f.f.1 <- subset_taxa(phy.f.f, tax_table(phy.f.f)[,"Rank1"] != "Unassigned")

# Check for cnidarian mitochondria
"f__cnidarian_mitochondria" %in% get_taxa_unique(phy.f.f.1,"Rank5") # No instances

# Check for chloroplasts
"c__Chloroplast" %in% get_taxa_unique(phy.f.f.1,"Rank3") # Yes, need to remove
temp <- subset_taxa(phy.f.f.1, tax_table(phy.f.f.1)[,"Rank3"] == "c__Chloroplast")
chloroplast_names <- taxa_names(temp)
goodTaxa <- setdiff(taxa_names(phy.f.f.1), chloroplast_names)
phy.f.f.2 <- prune_taxa(goodTaxa, phy.f.f.1)

# Check for mitochondria
temp2 <- subset_taxa(phy.f.f.2, tax_table(phy.f.f.2)[,"Rank5"] == "f__mitochondria") #yes
mito_names <- taxa_names(temp2)
goodTaxa2 <- setdiff(taxa_names(phy.f.f.2), mito_names)
phy.f.f.3 <- prune_taxa(goodTaxa2, phy.f.f.2)

phy.f.f <- phy.f.f.3

###################### 
names(sample_data(phy.f.f))

# Make proportional taxa
phy.f.f.p <- transform_sample_counts(phy.f.f, function(x) x/sum(x))

# Save 
save(list=ls(),file="analyses/16S/KI_Sym_Microbes_16S_f.RData")


##############
# Subset coral species
phy.f.f.p.Pocillopora <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Pocillopora grandis")
phy.f.f.p.Pocillopora <- prune_taxa(taxa_sums(phy.f.f.p.Pocillopora)>0,phy.f.f.p.Pocillopora)

phy.f.f.p.Porites <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Porites lobata")
phy.f.f.p.Porites <- prune_taxa(taxa_sums(phy.f.f.p.Porites)>0,phy.f.f.p.Porites)

phy.f.f.p.Montipora <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Montipora aequituberculata")
phy.f.f.p.Montipora <- prune_taxa(taxa_sums(phy.f.f.p.Montipora)>0,phy.f.f.p.Montipora)

phy.f.f.p.FaviaM <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Favia matthaii")
phy.f.f.p.FaviaM <- prune_taxa(taxa_sums(phy.f.f.p.FaviaM)>0,phy.f.f.p.FaviaM)

phy.f.f.p.Hydnophora <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Hydnophora microconos")
phy.f.f.p.Hydnophora <- prune_taxa(taxa_sums(phy.f.f.p.Hydnophora)>0,phy.f.f.p.Hydnophora)

phy.f.f.p.Faviasp <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Favia sp")
phy.f.f.p.Faviasp <- prune_taxa(taxa_sums(phy.f.f.p.Faviasp)>0,phy.f.f.p.Faviasp)

phy.f.f.p.Platygyra <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Platygyra daedalea")
phy.f.f.p.Platygyra <- prune_taxa(taxa_sums(phy.f.f.p.Platygyra)>0,phy.f.f.p.Platygyra)

phy.f.f.p.Favites <- subset_samples(phy.f.f.p,sample_data(phy.f.f.p)$field_host_name=="Favites pentagona")
phy.f.f.p.Favites <- prune_taxa(taxa_sums(phy.f.f.p.Favites)>0,phy.f.f.p.Favites)

# Save
save(list=c("phy.f.f.p.Pocillopora","phy.f.f.p.Porites",
            "phy.f.f.p.Montipora","phy.f.f.p.FaviaM",
            "phy.f.f.p.Hydnophora","phy.f.f.p.Faviasp",
            "phy.f.f.p.Platygyra","phy.f.f.p.Favites"),
     file="analyses/16S/KI_Sym_Microbes_16S_f_byspecies.RData")
