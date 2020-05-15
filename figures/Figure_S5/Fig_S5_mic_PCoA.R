# Load necessary data
load("analyses/16S/mic_betadisper.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Load necessary packages
library(vegan)
library(phyloseq)

# Make Figure
jpeg(filename="figures/Figure_S5/Fig_S5_mic_PCoA.jpg", 
     width = 8, height = 16, units="in",res = 300)
par(mfrow=c(4,2),mar=c(1,1,3,1))

plot(phy.f.f.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("All corals", line = 1, cex.main=2.2)
ordihull(phy.f.f.bd.braydist, sample_data(phy.f.f)$human_disturbance,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(phy.f.f.Mfol.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("M. aequituberculata")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Mfol.bd.braydist, 
         sample_data(phy.f.f.Mfol)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=c(15,16,17), 
       col=distcols[c(2:3,5)], cex=2, pt.cex=4,
       legend=c("Low","Medium","Very High"))

plot(phy.f.f.Plob.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. lobata")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Plob.bd.braydist, 
         sample_data(phy.f.f.Plob)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Hydno.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("H. microconos")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Hydno.bd.braydist, 
         sample_data(phy.f.f.Hydno)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Platy.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)], pch = c(15,16,17,20,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. daedalea")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Platy.bd.braydist, 
         sample_data(phy.f.f.Platy)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Fpent.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(3,5)],  pch = c(16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. pentagona")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Fpent.bd.braydist, 
         sample_data(phy.f.f.Fpent)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Favmat.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(3,5)],  pch = c(15,16,17,20,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. matthaii")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Favmat.bd.braydist,
         sample_data(phy.f.f.Favmat)$human_disturbance,
         draw = c("polygon"), col = distcols[c(3,5)], alpha=0.2, lwd=0.05)

dev.off()

# Make Figure (pdf)
pdf(file="figures/Figure_S5/Fig_S5_mic_PCoA.pdf", width = 8, height = 16,useDingbats = FALSE)
par(mfrow=c(4,2),mar=c(1,1,3,1))

plot(phy.f.f.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("All corals", line = 1, cex.main=2.2)
ordihull(phy.f.f.bd.braydist, sample_data(phy.f.f)$human_disturbance,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(phy.f.f.Mfol.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("M. aequituberculata")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Mfol.bd.braydist, 
         sample_data(phy.f.f.Mfol)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)
legend("topleft", bty="n", pch=c(15,16,17), 
       col=distcols[c(2:3,5)], cex=2, pt.cex=4,
       legend=c("Low","Medium","Very High"))

plot(phy.f.f.Plob.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. lobata")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Plob.bd.braydist, 
         sample_data(phy.f.f.Plob)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Hydno.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)],  pch = c(15,16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("H. microconos")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Hydno.bd.braydist, 
         sample_data(phy.f.f.Hydno)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Platy.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)], pch = c(15,16,17,20,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. daedalea")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Platy.bd.braydist, 
         sample_data(phy.f.f.Platy)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(2:3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Fpent.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(3,5)],  pch = c(16,17),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. pentagona")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Fpent.bd.braydist, 
         sample_data(phy.f.f.Fpent)$human_disturbance,  
         draw = c("polygon"), col = distcols[c(3,5)], alpha=0.2, lwd=0.05)

plot(phy.f.f.Favmat.bd.braydist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(3,5)],  pch = c(15,16,17,20,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. matthaii")), line = 1, cex.main=2.2)
ordihull(phy.f.f.Favmat.bd.braydist,
         sample_data(phy.f.f.Favmat)$human_disturbance,
         draw = c("polygon"), col = distcols[c(3,5)], alpha=0.2, lwd=0.05)

dev.off()
