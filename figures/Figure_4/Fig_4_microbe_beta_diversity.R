# Load necessary data
load("figures/mic_beta_plots.RData")
load("figures/KI_Sym_Microbes_colors.RData")

# Make figure
jpeg("figures/Figure_4/Fig_4.jpg",height=6,width=4,units="in",res=300)
par(mfrow=c(2,1),mar=c(1,3,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), 
       heights=c(1,1.55))
boxplot(phy.f.f.bd.braydist,col=distcols[c(2:3,5)],cex.axis=0.9,ylab="",xlab="",
        xaxt='n',yaxt='n',ylim=c(0.4,0.8))
axis(tck=0,side=1,labels = c("Low", "Medium", "Very High"),at = c(1,2,3),padj = -1.5)
axis(tck=0.05,side=2,padj = 1.5)
title(ylab="Distance to Centroid", line=1.5, cex.lab=1.2)

par(mar=c(8,3,1,1))
boxplot(phy.f.f.bd.brayspp,col=speccols,cex.axis=0.9, ylab="",xlab="",
        yaxt='n',xaxt='n',ylim=c(0.15,0.8))
axis(tck=0,side=1,las=2,hadj = 0.9, labels = c("M. aequituberculata","P. grandis","P. lobata","H. microconos","P. ryukyuensis","F. pentagona","D. matthaii"),font=3, at=c(1,2,3,4,5,6,7))
axis(tck=0.05,side=2,padj = 1.5)
title(ylab="Distance to Centroid", line=1.5, cex.lab=1.2)
mtext(side = 1,text = "b",adj = 0.1,padj=-2,col="darkgray")
mtext(side = 1,text = "a,b",adj = 0.22,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.37,padj=-2,col="darkgray")
mtext(side = 1,text = "a",adj = 0.5,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.64,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.78,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.91,padj=-2,col="darkgray")
dev.off()

pdf("figures/Figure_4/Fig_4.pdf",height=6,width=4,useDingbats = FALSE)
par(mfrow=c(2,1),mar=c(1,3,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE), 
       heights=c(1,1.55))
boxplot(phy.f.f.bd.braydist,col=distcols[c(2:3,5)],cex.axis=0.9,ylab="",xlab="",
        xaxt='n',yaxt='n',ylim=c(0.4,0.8))
axis(tck=0,side=1,labels = c("Low", "Medium", "Very High"),at = c(1,2,3),padj = -1.5)
axis(tck=0.05,side=2,padj = 1.5)
title(ylab="Distance to Centroid", line=1.5, cex.lab=1.2)


par(mar=c(8,3,1,1))
boxplot(phy.f.f.bd.brayspp,col=speccols,cex.axis=0.9, ylab="",xlab="",
        yaxt='n',xaxt='n',ylim=c(0.15,0.8))
axis(tck=0,side=1,las=2,hadj = 0.9, labels = c("M. aequituberculata","P. grandis","P. lobata","H. microconos","P. ryukyuensis","F. pentagona","D. matthaii"),font=3, at=c(1,2,3,4,5,6,7))
axis(tck=0.05,side=2,padj = 1.5)
title(ylab="Distance to Centroid", line=1.5, cex.lab=1.2)
mtext(side = 1,text = "b",adj = 0.1,padj=-2,col="darkgray")
mtext(side = 1,text = "a,b",adj = 0.22,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.37,padj=-2,col="darkgray")
mtext(side = 1,text = "a",adj = 0.5,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.64,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.78,padj=-2,col="darkgray")
mtext(side = 1,text = "b",adj = 0.91,padj=-2,col="darkgray")
dev.off()