# Load necessary libraries
library(ggplot2)
library(vegan)
library(lme4)
library(MuMIn)
library(lsmeans)
library(grid)
library(gridExtra)
library(phyloseq)

# Load necessary data
load("analyses/16S/KI_Sym_Microbes_16S_f.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Function to be able to plot the lsmeans results (From J. McDevitt-Irwin)
lsmeans.table <- function(x) {
  slm <- summary(x)
  slm_dat <- as.data.frame(slm[])
  return(slm_dat)
}

# Set random seed for reproducibility
set.seed(1010)

# Rarefy samples to an even depth. Chose 800 sequences for now.
phy.f800 <- rarefy_even_depth(phy.f.f,sample.size = 800,replace = FALSE)
phy.f800 <- subset_samples(phy.f800,sample_sums(phy.f800)>0)

phy.f800

phy.f800 <- subset_samples(phy.f800,sample_data(phy.f800)$field_host_name!="Favia speciosa")
phy.f800 <- subset_samples(phy.f800,sample_data(phy.f800)$field_host_name!="Favia sp")
phy.f800 <- subset_samples(phy.f800,sample_data(phy.f800)$field_host_name!="Acropora sp")
phy.f800 <- subset_samples(phy.f800,sample_data(phy.f800)$field_host_name!="Favites halicora")

phy.f800 <- prune_taxa(taxa_sums(phy.f800) > 0, phy.f800)

# Calculate shannon diversity (Modified from Jamie)
phy.f800.shannon <- data.frame(sample_data(phy.f800), 
                               estimate_richness(phy.f800, measures="Shannon"))
qqnorm(phy.f800.shannon$Shannon) # seems ok

# Pick the best model
# Chose this model structure because we want to know whether coral species and human disturbance influence alpha diversity (so we set them as fixed effects). Reef_name (site) is nested within human disturbance to account for site-by-site variability.
model0 <- lm(Shannon ~ field_host_name*human_disturbance/reef_name, 
             data=phy.f800.shannon)
MuMIn::AICc(model0) ### 637 This is the worst model
model1 <- lm(Shannon ~ field_host_name+human_disturbance/reef_name, 
             data=phy.f800.shannon)
MuMIn::AICc(model1) ### 562 This is the best model
model2 <- lm(Shannon ~ human_disturbance/reef_name, data=phy.f800.shannon)
MuMIn::AICc(model2) ### 589
model3<- lm(Shannon ~ field_host_name, data=phy.f800.shannon)
MuMIn::AICc(model3) ### 573
# so the best model is coral species + human_disturbance/reef_name
a1 <- aov(model1) # to look at the stats
summary(model1)
TukeyHSD(a1)

# Note: lsmeans can't deal with nested structure
model1.1<- lm(Shannon ~ field_host_name + human_disturbance, data=phy.f800.shannon)
model1.1_plot <- lsmeans(model1.1, c("human_disturbance","field_host_name"))
model1.1_plot
plot(model1.1_plot)
pairs(model1.1_plot)
model1.2_plot <- lsmeans(model1.1, c("human_disturbance"))
plot(model1.2_plot)
model1.3_plot <- lsmeans(model1.1, c("field_host_name"))
plot(model1.3_plot)

model1.1_table <- lsmeans.table(model1.1_plot) 
model1.2_table <- lsmeans.table(model1.2_plot) 
model1.3_table <- lsmeans.table(model1.3_plot) 

########################################################
# Make plots
model1.1_table$field_host_name <- gsub("Favia","D.",model1.1_table$field_host_name)
model1.1_table$field_host_name <- gsub("Porites","P.",model1.1_table$field_host_name)
model1.1_table$field_host_name <- gsub("Montipora","M.",model1.1_table$field_host_name)
model1.1_table$field_host_name <- gsub("Pocillopora","P.",model1.1_table$field_host_name)
model1.1_table$field_host_name <- gsub("Hydnophora","H.",model1.1_table$field_host_name)
model1.1_table$field_host_name <- gsub("Platygyra daedalea","P. ryukyuensis",model1.1_table$field_host_name)
model1.1_table$field_host_name <- gsub("Favites","F.",model1.1_table$field_host_name)
model1.1_table$field_host_name <- as.factor(model1.1_table$field_host_name)
# model1.1_table <- model1.1_table[c(4,6,7,3,5,2,1),]
model1.1_table$field_host_name <- factor(model1.1_table$field_host_name, levels=c("M. aequituberculata","P. grandis","P. lobata","H. microconos","P. ryukyuensis","F. pentagona","D. matthaii"))     

p1 <- ggplot(model1.1_table, 
             aes(x=field_host_name, y=lsmean, col=human_disturbance)) + 
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), 
                  size=1, position=position_dodge(width=0.5)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=12),
        panel.background = element_blank(),
        legend.position = 'n',
        axis.text.x = element_text(angle = 45, hjust = 1,face="italic"),
        plot.margin = margin(0,0,0,1,"cm")) + 
  labs(x="") + labs(y="Shannon Fitted Means") + labs(colour="Local Disturbance") +
  scale_colour_manual(values=distcols) + 
  scale_y_continuous(limits = c(0, 4)) +
  NULL
p1

p2 <- ggplot(model1.2_table, 
             aes(x=human_disturbance, y=lsmean, col=human_disturbance)) + 
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), 
                  size=0.5, position=position_dodge(width=0.5)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=12),
        panel.background = element_blank(),
        axis.line = element_line())+
        # axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="Local Disturbance") +
  # labs(x="") +
  labs(y="Shannon Fitted Means") + 
  labs(colour="Local Disturbance") +
  scale_colour_manual(values=distcols) + 
  scale_y_continuous(limits = c(1.4,3.2)) +
  guides(color=FALSE)+
  NULL
p2

model1.3_table$field_host_name <- gsub("Favia","D.",model1.3_table$field_host_name)
model1.3_table$field_host_name <- gsub("Porites","P.",model1.3_table$field_host_name)
model1.3_table$field_host_name <- gsub("Montipora","M.",model1.3_table$field_host_name)
model1.3_table$field_host_name <- gsub("Pocillopora","P.",model1.3_table$field_host_name)
model1.3_table$field_host_name <- gsub("Hydnophora","H.",model1.3_table$field_host_name)
model1.3_table$field_host_name <- gsub("Platygyra daedalea","P. ryukyuensis",model1.3_table$field_host_name)
model1.3_table$field_host_name <- gsub("Favites","F.",model1.3_table$field_host_name)
model1.3_table$field_host_name <- as.factor(model1.3_table$field_host_name)
# model1.3_table <- model1.3_table[c(4,6,7,3,5,2,1),]
model1.3_table$field_host_name <- factor(model1.3_table$field_host_name, levels=c("M. aequituberculata","P. grandis","P. lobata","H. microconos","P. ryukyuensis","F. pentagona","D. matthaii"))           

p3 <- ggplot(model1.3_table, 
             aes(x=field_host_name, y=lsmean, 
                 col=field_host_name)) + 
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), 
                  size=1, position=position_dodge(width=0.5)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=12),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,face="italic"),
        axis.line = element_line(),
        plot.margin = margin(0,0,0,1,"cm")) + 
  labs(x="Coral Species") + labs(y="Shannon Fitted Means") + 
  labs(colour="Coral Species") +
  scale_y_continuous(limits = c(0, 3.5)) +
  guides(color=FALSE)+
  scale_color_manual(values = speccols)+
  NULL
p3

p_alpha_dist <- p2
p_alpha_dist_spp <- p1
p_alpha_spp <- p3

save(p_alpha_dist_spp,p_alpha_dist,p_alpha_spp,file = "figures/mic_alpha_plots.RData")

# Make jpeg of only human disturbance
jpeg(filename = "figures/alpha_diversity_dist.jpeg",width=3,height=4,units="in",res=300)
p2
dev.off()

# Make pdf of only human disturbance
pdf(file = "figures/alpha_diversity_dist.pdf",width=2,height=1.6,useDingbats = FALSE)
p2
dev.off()

# Make jpeg of only coral species 
jpeg(filename = "figures/alpha_diversity_coralsp.jpeg",width=6,height=4,units="in",res=300)
p3
dev.off()

# Make pdf of only coral species 
pdf(file = "figures/alpha_diversity_coralsp.pdf",width=6,height=4)
p3
dev.off()

# Make jpeg of human disturbance + coral species
jpeg(filename = "figures/alpha_diversity.jpeg",width=6,height=6,units="in",res=300)
p1
#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.55, height = 0.32, x = 0.73, y = 0.87)
print(p1)
print(p2, vp = vp)
dev.off()

# Make jpeg of human disturbance + coral species
pdf(file = "figures/alpha_diversity.pdf",width=6,height=6,useDingbats = FALSE)
# p1
#A viewport taking up a fraction of the plot area
p1
dev.off()

