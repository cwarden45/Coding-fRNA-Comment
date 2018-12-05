dNadd = 0.001
#not previously added for dS or dSadj, but go ahead and add it this time
#dN has more missing values than dS
#adding before scaling (so, scale log values *before* PCA, not just for linear regression)

### uncomment for 4-species divergence, 169 gene set (with ribosomal genes removed) ###
#input.file = "../fRNA_wo_ribosome_169_genes_expanded.txt"
#depVars = c("dN","dS","dSadj","dNdS","dNdSadj")
#yeast = read.table(input.file, head=T, sep="\t")
#d = data.frame(dN=yeast$dN+dNadd, dS=yeast$dS+dNadd, dSadj=yeast$dSadj+dNadd, dNdSadj=yeast$dNdSadj, dNdS=yeast$dNdS+dNadd,
#					cai=yeast$cai, len=yeast$len, expr=yeast$expr, disp=yeast$disp+dNadd,cov=yeast$coverage)
#rownames(d)=yeast$orf
#d = na.omit(d)
#d = data.frame(dN=scale(log(d$dN)), dS=scale(log(d$dS)), dSadj=scale(log(d$dSadj)), dNdSadj=scale(log(d$dNdSadj)), dNdS=scale(log(d$dNdS)),
#					cai=scale(log(d$cai)), len=scale(log(d$len)), expr=scale(log(d$expr)), disp=scale(log(d$disp)), cov=scale(log(d$cov)))
#output.prefix = "169gene"
#legend.inset = -0.4
#legend.cex = 1

### uncomment for 4-species divergence, 195 gene set (all genes with coding fRNAs)###
#input.file = "../fRNA_WITH_ribosome_195_genes_expanded.txt"
#depVars = c("dN","dS","dSadj","dNdS","dNdSadj")
#yeast = read.table(input.file, head=T, sep="\t")
#d = data.frame( dN=yeast$dN+dNadd, dS=yeast$dS+dNadd, dSadj=yeast$dSadj+dNadd, dNdSadj=yeast$dNdSadj+dNadd, dNdS=yeast$dNdS+dNadd,
#				cai=yeast$cai, len=yeast$len, expr=yeast$expr, disp=yeast$disp+dNadd,cov=yeast$coverage)
#rownames(d)=yeast$orf
#d = na.omit(d)
#d = data.frame(dN=scale(log(d$dN)), dS=scale(log(d$dS)), dSadj=scale(log(d$dSadj)), dNdSadj=scale(log(d$dNdSadj)), dNdS=scale(log(d$dNdS)),
#					cai=scale(log(d$cai)), len=scale(log(d$len)), expr=scale(log(d$expr)), disp=scale(log(d$disp)), cov=scale(log(d$cov)))
#output.prefix = "195gene_WITH_ribosome"
#legend.inset = -0.4
#legend.cex = 1

### uncomment for re-calculated 2-species divergence, 169 gene set (with ribosomal genes removed) ###
#input.file = "../fRNA_wo_ribosome_169_genes_expanded.txt"
#depVars = c("small_dN","small_dS","small_dSadj","small_dNdS","small_dNdSadj")
#yeast = read.table(input.file, head=T, sep="\t")
#d = data.frame(small_dSadj=yeast$small_dSadj+dNadd, small_dNdSadj=yeast$small_dNdSadj+dNadd, small_dNdS=yeast$small_dNdS+dNadd, small_dN=yeast$sacCer_sacPar_dN+dNadd, small_dS=yeast$sacCer_sacPar_dS+dNadd,
#				cai=yeast$cai, len=yeast$len, expr=yeast$expr, disp=yeast$disp+dNadd,cov=yeast$coverage)
#rownames(d)=yeast$orf
#d = na.omit(d)
#d = data.frame(small_dN=scale(log(d$small_dN)), small_dS=scale(log(d$small_dS)), small_dSadj=scale(log(d$small_dSadj)), small_dNdSadj=scale(log(d$small_dNdSadj)), small_dNdS=scale(log(d$small_dNdS)),
#					cai=scale(log(d$cai)), len=scale(log(d$len)), expr=scale(log(d$expr)), disp=scale(log(d$disp)), cov=scale(log(d$cov)))
#output.prefix = "169gene_small_divergence"
#legend.inset = -0.45
#legend.cex = 0.8

### uncomment for re-calculated 2-species divergence, 195 gene set (all genes with coding fRNAs) ###
input.file = "../fRNA_WITH_ribosome_195_genes_expanded.txt"
depVars = c("small_dN","small_dS","small_dSadj","small_dNdS","small_dNdSadj")
yeast = read.table(input.file, head=T, sep="\t")
d = data.frame( small_dSadj=yeast$small_dSadj+dNadd, small_dNdSadj=yeast$small_dNdSadj+dNadd, small_dNdS=yeast$small_dNdS+dNadd, small_dN=yeast$sacCer_sacPar_dN+dNadd, small_dS=yeast$sacCer_sacPar_dS+dNadd,
				cai=yeast$cai, len=yeast$len, expr=yeast$expr, disp=yeast$disp+dNadd,cov=yeast$coverage)
rownames(d)=yeast$orf
d = na.omit(d)
d = data.frame(small_dN=scale(log(d$small_dN)), small_dS=scale(log(d$small_dS)), small_dSadj=scale(log(d$small_dSadj)), small_dNdSadj=scale(log(d$small_dNdSadj)), small_dNdS=scale(log(d$small_dNdS)),
					cai=scale(log(d$cai)), len=scale(log(d$len)), expr=scale(log(d$expr)), disp=scale(log(d$disp)), cov=scale(log(d$cov)))
output.prefix = "195gene_WITH_ribosome_small_divergence"
legend.inset = -0.45
legend.cex = 0.8

#################################################################
### create different output files by changing variables above ###
#################################################################

print(dim(yeast))
print(dim(d))

#used in more recent COH QC scripts
#usually I would reduce the genes as features (to number of samples)...however, the PCA Regression reduced number of variables (to number of genes)
pca.values = prcomp(as.matrix(t(d)))
pc.values = data.frame(pca.values$rotation)
variance.explained = (pca.values$sdev)^2 / sum(pca.values$sdev^2)
pca.table = data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

pc.mat = t(pca.table[,3:ncol(pca.table)])

compID = c()
cor.coefficient = c()
pcr.glm.pvalue = c()

dir.create("R-base_PCA_plots_Scaled")
dir.create("R-base_Clustering_plots_Scaled")

for (i in 1:length(depVars)){
	rate = depVars[i]
	
	for (j in 1:ncol(pc.mat)){
		PC = pc.mat[,j]

		tempID = paste(rate,"_versus_PC",j,sep="")
		#print(tempID)
		#print(tempID)
		compID = c(compID, tempID)
		
		Var2 = d[,rate]+dNadd
		Var1 = pc.values[,j] #don't log transform...often has negative numbers
		
		#I saw some similar code with all varibles (but that wouldn't be the main results, and would have been over-fit if it was)
		cor.coeff = cor(Var1,Var2,use="pairwise.complete.obs")
		cor.coefficient = c(cor.coefficient, cor.coeff)
		fit = glm(Var2 ~ Var1)
		result = summary(fit)$coefficients
		pcr.glm.pvalue = c(pcr.glm.pvalue, result[2,4])
		
		plot.file = paste("R-base_PCA_plots_Scaled/",output.prefix,"_",tempID,".png",sep="")
		png(plot.file)
		plot(Var1,Var2,
			xlab=paste("scale(ln(",rate,"))",sep=""), ylab = paste("PC",j,sep=""),
			pch=16, main = paste("r = ",signif(cor.coeff, digits=3)," , glm p = ",round(result[2,4], digits=3),sep=""))
		abline(fit, col="red")
		dev.off()
	}#end for (j in 1:ncol(g$scores))
	
	### smaller divergence
}#end for (i in 1:length(depVars))

pcr.glm.fdr = p.adjust(pcr.glm.pvalue, "fdr")

stat.table = data.frame(compID, cor.coefficient, pcr.glm.pvalue, pcr.glm.fdr)
stat.file = paste("R-base_PCA_plots_Scaled/",output.prefix,"_PCA_stats.txt",sep="")
write.table(stat.table, stat.file, quote=F, sep="\t", row.names=F)

###try "normal" PCA plot (to see relationship between variables)
#again, use code similar to what I used as template for COH QC analysis

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

pca.values = prcomp(as.matrix(d))
pc.values = data.frame(pca.values$rotation)
variance.explained = (pca.values$sdev)^2 / sum(pca.values$sdev^2)
pca.table = data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

variance.file = paste("R-base_PCA_plots_Scaled/",output.prefix,"_R-Base_prcomp_percent_variance_explained.png",sep="")
png(variance.file)
bar.labels = paste("PC",1:length(variance.explained),sep="")
bar.obj = barplot(100 * variance.explained, names.arg = bar.labels, ylim=c(0,100),
			col="blue", main="Percent Variance Explained", las=2)
#from https://stackoverflow.com/questions/12481430/how-to-display-the-frequency-at-the-top-of-each-factor-in-a-barplot-in-r
text(x=bar.obj, y=100 * variance.explained,
	label=paste(signif(100 * variance.explained,1),"%",sep=""),
	xpd=T, pos = 3, cex=0.7)
dev.off()

pca.text.file = paste("R-base_PCA_plots_Scaled/",output.prefix,"_R-Base_prcomp_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

pc.mat = data.frame(t(pca.table[,3:ncol(pca.table)]))

variables = rownames(pc.mat)
color.palette=fixed.color.palatte[1:nrow(pc.mat)]

png(paste("R-base_PCA_plots_Scaled/",output.prefix,"_PC1_vs_PC2.png",sep=""))
par(mar=c(5,5,3,10))
plot(pc.mat$PC1, pc.mat$PC2, col = color.palette, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
			ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19,
	    		main = "PCA Gene Reduction")
text(pc.mat$PC1, pc.mat$PC2, labels = variables, pos = 4, xpd=T)
legend("right",legend=variables,col=color.palette, 
		xpd=T, inset = legend.inset, pch=19, cex=legend.cex)
dev.off()

colLab = function(n, labelColors, clusMember) { 
   if(is.leaf(n)) { 
       a <- attributes(n) 
	   #print(a)
       # clusMember - vector of sample names (ordered to match label color.palette)
       # labelColors - a vector of color.palette for the above grouping 
       labCol <- labelColors[clusMember == a$label]
	   #print(labCol)
       attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol) 
   } 
   n 
}

	cluster.file = paste("R-base_Clustering_plots_Scaled/Euclidian_Distance_Complete_",output.prefix,"_Hierarchical_Clustering_na_omit.png",sep="")
		#cor.mat = cor(as.matrix(d), use="pairwise.complete.obs")
		#dis.mat = 1 - cor.mat
		#dist1 = as.dist(dis.mat)
		
		dist1 = dist(as.matrix(t(d)))
		
hc = hclust(dist1)
dend1 = as.dendrogram(hc)
png(file = cluster.file)
par(mar=c(5,5,3,10))
dend1 = dendrapply(dend1, colLab, labelColors=color.palette, clusMember=variables) 
a = attributes(dend1) 
attr(dend1, "nodePar") = c(a$nodePar, lab.col = color.palette) 
plot(dend1, horiz=T, main = "Euclidian Distance Hierarchical Clustering")
dev.off()