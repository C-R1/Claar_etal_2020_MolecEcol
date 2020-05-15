# Load necessary libraries
library(ggplot2)
library(vegan)
library(lme4)
library(MuMIn)
library(lsmeans)
library(grid)
library(gridExtra)
library(phyloseq)
library(dplyr)

library(devtools)
# install_github("microbiome/microbiome")
library(microbiome)

# Load necessary data
load("analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Function to be able to plot the lsmeans results (From J. McDevitt-Irwin)
lsmeans.table <- function(x) {
  slm <- summary(x)
  slm_dat <- as.data.frame(slm[])
  return(slm_dat)
}

# Set random seed for reproducibility
set.seed(1010)

sample_data(phyASV.f.c)$Dist <- factor(sample_data(phyASV.f.c)$Dist, levels = c("VeryLow","Low","Medium","VeryHigh"))

sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Platygyra rykyuensis"] <- "P. ryukyuensis"
sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Montipora aequituberculata"] <- "M. aequituberculata"
sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Pocillopora grandis"] <- "P. grandis"
sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Porites lobata"] <- "P. lobata"
sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Hydnophora microconos"] <- "H. microconos"
sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Favites pentagona"] <- "F. pentagona"
sample_data(phyASV.f.c)$Coral_Species[sample_data(phyASV.f.c)$Coral_Species=="Favia matthaii"] <- "D. matthaii"

sample_data(phyASV.f.c)$Coral_Species <- factor(sample_data(phyASV.f.c)$Coral_Species, levels = c("M. aequituberculata","P. grandis","P. lobata","H. microconos","P. ryukyuensis","F. pentagona","D. matthaii"))

unique(sample_data(phyASV.f.c)$Coral_Species)

# Rarefy samples to an even depth. Chose 800 sequences.
phyASV.f.c.800 <- rarefy_even_depth(phyASV.f.c,sample.size = 800,
                                    replace = FALSE)
phyASV.f.c.800 <- subset_samples(phyASV.f.c.800,sample_sums(phyASV.f.c.800)>0)

phyASV.f.c.800 <- subset_samples(phyASV.f.c.800,sample_data(phyASV.f.c.800)$Coral_Species!="Favia speciosa")
phyASV.f.c.800 <- subset_samples(phyASV.f.c.800,sample_data(phyASV.f.c.800)$Coral_Species!="Favia sp")
phyASV.f.c.800 <- subset_samples(phyASV.f.c.800,sample_data(phyASV.f.c.800)$Coral_Species!="Acropora sp")
phyASV.f.c.800 <- subset_samples(phyASV.f.c.800,sample_data(phyASV.f.c.800)$Coral_Species!="Favites halicora")

phyASV.f.c.800 <- subset_taxa(phyASV.f.c.800,taxa_sums(phyASV.f.c.800)>0)

# Calculate shannon diversity (Modified from Jamie)
phyASV.f.c.800.shannon <- data.frame(sample_data(phyASV.f.c.800), 
                               estimate_richness(phyASV.f.c.800, 
                                                 measures="Shannon"))
qqnorm(phyASV.f.c.800.shannon$Shannon) # looks terrible
qqline(phyASV.f.c.800.shannon$Shannon) # looks terrible

phyASV.f.c.800.hell <- microbiome::transform(phyASV.f.c.800, transform = "hellinger", target = "OTU")
# Calculate shannon diversity (Modified from Jamie)
phyASV.f.c.800.hell.shannon <- data.frame(sample_data(phyASV.f.c.800.hell), 
                                    estimate_richness(phyASV.f.c.800.hell, measures="Shannon"))
qqnorm(phyASV.f.c.800.hell.shannon$Shannon) # looks terrible
qqline(phyASV.f.c.800.hell.shannon$Shannon) # looks terrible

phyASV.f.c.800.log10 <- microbiome::transform(phyASV.f.c.800, transform = "log10p", target = "OTU")
# Calculate shannon diversity (Modified from Jamie)
phyASV.f.c.800.log10.shannon <- data.frame(sample_data(phyASV.f.c.800.log10), 
                                         estimate_richness(phyASV.f.c.800.log10, measures="Shannon"))
qqnorm(phyASV.f.c.800.log10.shannon$Shannon) # looks the best
qqline(phyASV.f.c.800.log10.shannon$Shannon)

# Pick the best model
# Chose this model structure because we want to know whether coral species and human disturbance influence alpha diversity (so we set them as fixed effects). Reef_name (site) is nested within human disturbance to account for site-by-site variability.
model0 <- lm(Shannon ~ Coral_Species*Dist/Site, 
             data=phyASV.f.c.800.log10.shannon)
MuMIn::AICc(model0) ### 332
model1 <- lm(Shannon ~ Coral_Species+Dist/Site, 
             data=phyASV.f.c.800.log10.shannon)
MuMIn::AICc(model1) ### 298 #best
model2 <- lm(Shannon ~ Dist/Site, data=phyASV.f.c.800.log10.shannon)
MuMIn::AICc(model2) ### 726
model3<- lm(Shannon ~ Coral_Species, data=phyASV.f.c.800.log10.shannon)
MuMIn::AICc(model3) ### 313
model4<- lm(Shannon ~ Coral_Species+Dist, data=phyASV.f.c.800.log10.shannon)
MuMIn::AICc(model4) ### 310 #secondbest

summary(model1) #"Symbiodiniaceae sequence diversity was significantly different among coral species "
a1 <- aov(model1) # to look at the stats
TukeyHSD(a1)

#################################
# Note: lsmeans can't deal with nested structure
model1.1<- lm(Shannon ~ Coral_Species + Dist, data=phyASV.f.c.800.log10.shannon)
model1.1_plot <- lsmeans(model1.1, c("Dist","Coral_Species"))
model1.1_plot
plot(model1.1_plot)
pairs(model1.1_plot)
model1.2_plot <- lsmeans(model1.1, c("Dist"))
plot(model1.2_plot)
model1.3_plot <- lsmeans(model1.1, c("Coral_Species"))
plot(model1.3_plot)

model1.3_plot %>% 
  as.data.frame() %>% 
  arrange(desc(lsmean))

model1.2_plot %>% 
  as.data.frame() %>% 
  arrange(desc(lsmean))

model1.1_table <- lsmeans.table(model1.1_plot) 
model1.2_table <- lsmeans.table(model1.2_plot) 
model1.3_table <- lsmeans.table(model1.3_plot) 

########################################################
# Make plots

p1 <- ggplot(model1.1_table, 
             aes(x=Coral_Species, y=lsmean, col=Dist)) + 
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
  scale_y_continuous(limits = c(2.0, 4)) +
  NULL
p1

p2 <- ggplot(model1.2_table, 
             aes(x=Dist, y=lsmean, col=Dist)) + 
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), 
                  size=0.5, position=position_dodge(width=0.5)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=12),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="") +
  labs(y="Shannon Fitted Means") + 
  labs(colour="Local Disturbance") +
  scale_colour_manual(values=distcols) + 
  scale_y_continuous(limits = c(2,4)) +
  guides(color=FALSE)+
  NULL
p2

p3 <- ggplot(model1.3_table, 
             aes(x=Coral_Species, y=lsmean, 
                 col=Coral_Species)) + 
  geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), 
                  size=1, position=position_dodge(width=0.5)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text=element_text(size=12),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,face="italic"),
        axis.line = element_line()) + 
  labs(x="Coral Species") + labs(y="Shannon Fitted Means") + 
  labs(colour="Coral Species") +
  scale_y_continuous(limits = c(2, 4)) +
  guides(color=FALSE)+
  scale_color_manual(values = speccols)+
  NULL
p3

p_alpha_dist <- p2
p_alpha_dist_spp <- p1
p_alpha_spp <- p3

save(p_alpha_dist_spp,p_alpha_dist,p_alpha_spp,file = "figures/sym_alpha_plots.RData")

# Make jpeg of only human disturbance
jpeg(filename = "figures/alpha_diversity_dist.jpeg",width=3,height=4,units="in",res=300)
p2
dev.off()


# Make jpeg of only human disturbance
jpeg(filename = "figures/sym_alpha_diversity_dist.jpeg",width=3,height=4,units="in",res=300)
p2
dev.off()

# Make pdf of only human disturbance
pdf(file = "figures/sym_alpha_diversity_dist.pdf",width=4,height=3.2,useDingbats = FALSE)
p2
dev.off()

# Make jpeg of only coral species 
jpeg(filename = "figures/sym_alpha_diversity_coralsp.jpeg",width=6,height=4,units="in",res=300)
p3
dev.off()

# Make jpeg of only coral species 
pdf(file = "figures/sym_alpha_diversity_coralsp.pdf",width=6,height=4,useDingbats = FALSE)
p3
dev.off()

# Make jpeg of human disturbance + coral species
jpeg(filename = "figures/sym_alpha_diversity.jpeg",width=6,height=6,units="in",res=300)
p1
#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.55, height = 0.32, x = 0.73, y = 0.87)
print(p1)
print(p2, vp = vp)
dev.off()

# Make jpeg of human disturbance + coral species
pdf(file = "figures/sym_alpha_diversity.pdf",width=6,height=6,useDingbats = FALSE)
p1
dev.off()

