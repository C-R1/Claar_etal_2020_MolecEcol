# Load necessary data
load("analyses/ITS2/sym_betadisper.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Load necessary packages
library(vegan)
library(phyloseq)


# Make Figure
jpeg(filename="figures/Figure_S3/sym_PCoA.jpg", 
     width = 8, height = 16, units="in",res = 300)
par(mfrow=c(4,2),mar=c(1,1,3,1))

plot(coral.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("All corals", line = 1, cex.main=2.2)
ordihull(coral.bd.dist, sample_data(phyASV.f.c.p)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Montipora.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("M. aequituberculata")), line = 1, cex.main=2.2)
ordihull(Montipora.bd.dist, sample_data(phyASV.f.c.p.Montipora)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=c(15,16,17,20), 
       col=distcols[c(1:3,5)], cex=2, pt.cex=4,
       legend=c("Very Low","Low","Medium","Very High"))

plot(Pocillopora.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. grandis")), line = 1, cex.main=2.2)
ordihull(Pocillopora.bd.dist, sample_data(phyASV.f.c.p.Pocillopora)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Porites.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. lobata")), line = 1, cex.main=2.2)
ordihull(Porites.bd.dist, sample_data(phyASV.f.c.p.Porites)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Hydnophora.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("H. microconos")), line = 1, cex.main=2.2)
ordihull(Hydnophora.bd.dist, sample_data(phyASV.f.c.p.Hydnophora)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Platygyra.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)], pch = c(16,17,20,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. daedalea")), line = 1, cex.main=2.2)
ordihull(Platygyra.bd.dist, sample_data(phyASV.f.c.p.Platygyra)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Favites.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:2,5)],  pch = c(16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. pentagona")), line = 1, cex.main=2.2)
ordihull(Favites.bd.dist, sample_data(phyASV.f.c.p.Favites)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(FaviaM.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:2,5)],  pch = c(16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. matthaii")), line = 1, cex.main=2.2)
ordihull(FaviaM.bd.dist, sample_data(phyASV.f.c.p.FaviaM)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

dev.off()

# Make Figure (pdf)
pdf(file="figures/Figure_S3/sym_PCoA.pdf",width = 8, height = 16, useDingbats = FALSE)
par(mfrow=c(4,2),mar=c(1,1,3,1))

plot(coral.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title("All corals", line = 1, cex.main=2.2)
ordihull(coral.bd.dist, sample_data(phyASV.f.c.p)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Montipora.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("M. aequituberculata")), line = 1, cex.main=2.2)
ordihull(Montipora.bd.dist, sample_data(phyASV.f.c.p.Montipora)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)
legend("topright", bty="n", pch=c(15,16,17,20), 
       col=distcols[c(1:3,5)], cex=2, pt.cex=4,
       legend=c("Very Low","Low","Medium","Very High"))

plot(Pocillopora.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. grandis")), line = 1, cex.main=2.2)
ordihull(Pocillopora.bd.dist, sample_data(phyASV.f.c.p.Pocillopora)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Porites.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. lobata")), line = 1, cex.main=2.2)
ordihull(Porites.bd.dist, sample_data(phyASV.f.c.p.Porites)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Hydnophora.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:3,5)],  pch = c(15,16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("H. microconos")), line = 1, cex.main=2.2)
ordihull(Hydnophora.bd.dist, sample_data(phyASV.f.c.p.Hydnophora)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Platygyra.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(2:3,5)], pch = c(16,17,20,18),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("P. daedalea")), line = 1, cex.main=2.2)
ordihull(Platygyra.bd.dist, sample_data(phyASV.f.c.p.Platygyra)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(Favites.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:2,5)],  pch = c(16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. pentagona")), line = 1, cex.main=2.2)
ordihull(Favites.bd.dist, sample_data(phyASV.f.c.p.Favites)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)

plot(FaviaM.bd.dist, hull=F, label=F, cex=2,
     main=NULL, col=distcols[c(1:2,5)],  pch = c(16,17,20),
     xaxt="n", yaxt="n", xlab=NULL, ylab=NULL, sub="")
title(expression(italic("F. matthaii")), line = 1, cex.main=2.2)
ordihull(FaviaM.bd.dist, sample_data(phyASV.f.c.p.FaviaM)$Dist,  
         draw = c("polygon"), col = distcols, alpha=0.2, lwd=0.05)
dev.off()