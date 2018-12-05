### uncomment for 169 gene set (with ribosomal genes removed) ###
input.file = "../fRNA_wo_ribosome_169_genes_expanded.txt"
output.prefix = "169gene"

### uncomment for 195 gene set (all genes with coding fRNAs) ###
#input.file = "../fRNA_WITH_ribosome_195_genes_expanded.txt"
#output.prefix = "195gene_WITH_ribosome"

#################################################################
### create different output files by changing variables above ###
#################################################################


input.table = read.table(input.file, head=T, sep="\t")
inv.length = 1/input.table$gene_length
#if log-transformed, need to add something above zero
input.table$disp = input.table$disp + 0.1
essential.num = as.numeric(input.table$essential)#essential boolean alterantive to disp
essential.num = essential.num + 0.1
sum.length= function(char.value){
	lengths = unlist(strsplit(as.character(char.value),split=","))
	lengths = as.numeric(lengths)
	return(sum(lengths))
}#end def sum.length
fold.length.sum = sapply(input.table$fold_length, sum.length)
input.table = data.frame(input.table, inv.length, essential.num, fold.length.sum)

compID = c()
cor.coefficient = c()
glm.pvalue = c()
altVars = c("gene_length","inv.length","expr","abund","cai","disp","essential.num","fold.length.sum")

dir.create("plots_non-divergence")

plotCol=rep("black",nrow(input.table))
plotCol[is.na(input.table$dN)]="gray"

#1-variable comparison
for (i in 1:length(altVars)){
	testVar = altVars[i]
	tempID = paste("fRNAcov_versus_",testVar,sep="")
	compID = c(compID, tempID)
	
	logVar2 = log(input.table[,testVar])
	logVar1 = log(input.table[,"coverage"])
	
	#I saw some similar code with all varibles (but that wouldn't be the main results, and would have been over-fit if it was)
	cor.coeff = cor(logVar1,logVar2,use="pairwise.complete.obs")
	cor.coefficient = c(cor.coefficient, cor.coeff)
	fit = glm(logVar2 ~ logVar1)
	result = summary(fit)$coefficients
	glm.pvalue = c(glm.pvalue, result[2,4])
	
	plot.file = paste("plots_non-divergence/",output.prefix,"_",tempID,".png",sep="")
	png(plot.file)
	par(mar=c(8,5,3,2))
	plot(logVar1,logVar2,
		xlab=paste("Scaled-ln ","fRNA coverage",sep=""), ylab = paste("ln(",testVar,")",sep=""), col=plotCol,
		pch=16, main = paste("r = ",signif(cor.coeff, digits=3)," , glm p = ",round(result[2,4], digits=3),sep=""))
	abline(fit, col="red")
	legend("bottom",legend=c("Wall dN","Missing Wall dN"),col=c("black","gray"),
				pch=16, ncol=2, xpd=T, inset=-0.3)
	dev.off()	

}#end for (i in 1:length(altVars))

glm.fdr = p.adjust(glm.pvalue, "fdr")

stat.table = data.frame(compID, cor.coefficient, glm.pvalue, glm.fdr)
stat.file = paste(output.prefix,"_stats.txt",sep="")
write.table(stat.table, stat.file, quote=F, sep="\t", row.names=F)