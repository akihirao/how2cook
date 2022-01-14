#Plot.PCO.ChiopiGRASDi_588indiv_R1R2.R
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(stringr) # handling characters

#setwd("~/work/Chiopi/Chiopi_GRASDi/R_work")

#read a vcf file (588 samples from R1+R2 data)
Chiopi.vcf <- vcfR::read.vcfR("Chiopi_GRASDi_ref1_R1R2.sca_all.snp.maf001.LDpruned.extract.3000.vcf.gz")
Stock.data <- read.table("ID_Stock_lab.txt", sep = "\t", header = TRUE)
outlier.indiv <- c("C0154","C0185","C0192","C0234","C0318","C0450")
#outlier.indiv.2nd <- c("C0086","C0087","C0154","C0185","C0192","C0234","C0318")


Chiopi.genind <- vcfR2genind(Chiopi.vcf)

pop(Chiopi.genind) <- Stock.data$Stock
no.pop <- length(unique(pop(Chiopi.genind)))

Chiopi.id <- indNames(Chiopi.genind)
Chiopi.id.non.outlier <- Chiopi.id[-which(Chiopi.id %in% outlier.indiv)]
Chiopi.non.outlier.genind <- Chiopi.genind[Chiopi.id.non.outlier,]
#Chiopi.id.non.outlier.2nd <- Chiopi.id[-which(Chiopi.id %in% outlier.indiv.2nd)]
#Chiopi.non.outlier.2nd.genind <- Chiopi.genind[Chiopi.id.non.outlier.2nd,]


Chiopi.gl <- vcfR2genlight(Chiopi.vcf)
Chiopi.pca <- glPca(Chiopi.gl, nf = 3, n.cores= 2)


X <- tab(Chiopi.genind, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)

X.non.outlier <- tab(Chiopi.non.outlier.genind, freq = TRUE, NA.method = "mean")
pca1.non.outlier <- dudi.pca(X.non.outlier, scale = FALSE, scannf = FALSE, nf = 3)

#X.non.outlier.2nd <- tab(Chiopi.non.outlier.2nd.genind, freq = TRUE, NA.method = "mean")
#pca1.non.outlier.2nd <- dudi.pca(X.non.outlier.2nd, scale = FALSE, scannf = FALSE, nf = 3)


no.eig <- length(Chiopi.pca$eig)
PCA.1.contribute <- round(100*Chiopi.pca$eig[1]/sum(Chiopi.pca$eig),digits=3)
PCA.2.contribute <- round(100*Chiopi.pca$eig[2]/sum(Chiopi.pca$eig),digits=3)
PCA.1.contribute.lab <- str_c("PC1"," (",PCA.1.contribute,"%)")
PCA.2.contribute.lab <- str_c("PC2"," (",PCA.2.contribute,"%)")



#s.class(pca1$li, pop(Traja.genind))
#par(mfrow=c(1,3))
par(mfrow=c(1,1))
col.pop <- funky(no.pop)
s.class(pca1$li, pop(Chiopi.genind), xax=1, yax=2, col=transp(col.pop,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
#title("PCA plot for Chionoecetes opilio")

#s.class(pca1.non.outlier$li, pop(Chiopi.non.outlier.genind), xax=1, yax=2, col=transp(col.pop,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
#title("PCA plot for Chionoecetes opilio")

#s.class(pca1.non.outlier.2nd$li, pop(Chiopi.non.outlier.2nd.genind), xax=1, yax=2, col=transp(col.pop,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)
#title("PCA plot for Chionoecetes opilio")


Chiopi.gl.tre <- nj(dist(as.matrix(Chiopi.gl)))
plot(Chiopi.gl.tre)

pop(Chiopi.gl) <- Stock.data$Stock


#barplot(100*Traja.pca$eig/sum(Traja.pca$eig), col = heat.colors(50), main = "PCA Eigenevalues")



Chiopi.pca.scores <- as.data.frame(Chiopi.pca$scores)
Chiopi.pca.scores$pop <- pop(Chiopi.gl)

library(ggplot2)
set.seed(9)
p <- ggplot(Chiopi.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + ggtitle("PCA plot for Chionoecetes opilio") + xlab(PCA.1.contribute.lab) + ylab(PCA.2.contribute.lab)
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
#p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()




