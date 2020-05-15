# Load necessary packages
library(ggplot2)
library(gridExtra)

# Load necessary data
load("figures/sym_alpha_plots.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Plot
p_alpha_dist <- p_alpha_dist + theme(axis.text.x = element_text(size=11),
                                     plot.margin = margin(1,1,30,1)) +
  ylim(c(0,4))+
  scale_x_discrete(labels = c("Very Low","Low","Medium","Very High"))
p_alpha_dist_spp <- p_alpha_dist_spp + 
  theme(axis.title.y = element_text(vjust=2,size=13),
        plot.margin = margin(1,1,1,5),
        axis.line = element_line()) + 
  ylim(c(0,4))
p_alpha_spp <- p_alpha_spp + theme(axis.title.x = element_blank(),
                                   plot.margin = margin(1,1.5,1,1))+
  ylim(c(0,4))

# Make figure
jpeg(filename = "figures/Figure_1/Fig_1.jpg",height = 7, width=6, res=300,units="in")
grid.arrange(p_alpha_dist_spp,p_alpha_dist,p_alpha_spp,
             widths=c(1,1),heights=c(2,1.5),
             nrow=2,ncol=2, layout_matrix = rbind(c(1,1),c(2,3)))
dev.off()

pdf(file = "figures/Figure_1/Fig_1.pdf",height = 7, width=6, useDingbats = FALSE)
grid.arrange(p_alpha_dist_spp,p_alpha_dist,p_alpha_spp,
             widths=c(1,1),heights=c(2,1.5),
             nrow=2,ncol=2, layout_matrix = rbind(c(1,1),c(2,3)))
dev.off()
