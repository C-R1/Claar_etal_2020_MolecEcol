# Load necessary packages
library(ggplot2)
library(gridBase)

# Load necessary data
load("figures/mic_alpha_plots.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Make figure
p_alpha_dist <- p_alpha_dist + theme(axis.text.x = element_text(size=11,
                                                                angle=45,
                                                                hjust = 1),
                                     axis.title.x = element_blank(),
                                     plot.margin = margin(1,1,45,1)) +
  ylim(c(1.25,3))+ 
  scale_x_discrete(labels = c("Low","Medium","Very High"))
p_alpha_dist_spp <- p_alpha_dist_spp + 
  theme(axis.title.y = element_text(vjust=2,size=13),
        plot.margin = margin(1,1,1,15),
        axis.line = element_line())  + 
  ylim(c(0,3.7)) 
p_alpha_spp <- p_alpha_spp + theme(axis.title.x = element_blank(),
                                   plot.margin = margin(1,1.5,1,1)) + 
  ylim(c(1.2,3.05))

jpeg(filename = "figures/Figure_3/Fig_3.jpg",height = 7, width=6, 
     res=300,units="in")
grid.arrange(p_alpha_dist_spp,p_alpha_dist,p_alpha_spp,
             widths=c(1,1),heights=c(2,1.5),
             nrow=2,ncol=2, layout_matrix = rbind(c(1,1),c(2,3)))
dev.off()

pdf(file = "figures/Figure_3/Fig_3.pdf",height = 7, width=6, useDingbats = FALSE)
grid.arrange(p_alpha_dist_spp,p_alpha_dist,p_alpha_spp,
             widths=c(1,1),heights=c(2,1.5),
             nrow=2,ncol=2, layout_matrix = rbind(c(1,1),c(2,3)))
dev.off()
