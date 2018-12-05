source("Drummond_functions.R")

### uncomment for 4-species divergence, 169 gene set (with ribosomal genes removed) ###
input.file = "../fRNA_wo_ribosome_169_genes_expanded.txt"
depVars = c("dN","dS","dSadj","dNdS","dNdSadj")
output.prefix = "169gene"

### uncomment for 4-species divergence, 195 gene set (all genes with coding fRNAs)###
#input.file = "../fRNA_WITH_ribosome_195_genes_expanded.txt"
#depVars = c("dN","dS","dSadj","dNdS","dNdSadj")
#output.prefix = "195gene_WITH_ribosome"

### uncomment for re-calculated 2-species divergence, 169 gene set (with ribosomal genes removed) ###
#input.file = "../fRNA_wo_ribosome_169_genes_expanded.txt"
#depVars = c("sacCer_sacPar_dN","sacCer_sacPar_dS","small_dSadj","small_dNdS","small_dNdSadj")
#output.prefix = "169gene_small_divergence"

### uncomment for re-calculated 2-species divergence, 195 gene set (all genes with coding fRNAs) ###
#input.file = "../fRNA_WITH_ribosome_195_genes_expanded.txt"
#depVars = c("sacCer_sacPar_dN","sacCer_sacPar_dS","small_dSadj","small_dNdS","small_dNdSadj")
#output.prefix = "195gene_WITH_ribosome_small_divergence"

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

dNadd = 0.001
#not previously added for dS or dSadj, but go ahead and add it this time
#dN has more missing values than dS

dir.create("plots_divergence")

#1-variable comparison
for (i in 1:length(depVars)){
	rate = depVars[i]
	tempID = paste("fRNAcov_versus_",rate,sep="")
	compID = c(compID, tempID)
	
	logVar2 = log(input.table[,rate]+dNadd)
	logVar1 = log(input.table[,"coverage"])
	
	#I saw some similar code with all varibles (but that wouldn't be the main results, and would have been over-fit if it was)
	cor.coeff = cor(logVar1,logVar2,use="pairwise.complete.obs")
	cor.coefficient = c(cor.coefficient, cor.coeff)
	fit = glm(logVar2 ~ logVar1)
	result = summary(fit)$coefficients
	glm.pvalue = c(glm.pvalue, result[2,4])
	
	plot.file = paste("plots_divergence/",output.prefix,"_",tempID,".png",sep="")
	png(plot.file)
	plot(logVar1,logVar2,
		xlab=paste("Scaled-ln ","fRNA coverage",sep=""), ylab = paste("ln(",rate,"+",dNadd,")",sep=""),
		pch=16, main = paste("r = ",signif(cor.coeff, digits=3)," , glm p = ",round(result[2,4], digits=3),sep=""))
	abline(fit, col="red")
	dev.off()
	
	for(j in 1:length(altVars)){
		compVar = altVars[j]
		
		tempID = paste(compVar,"_versus_",rate,sep="")
		compID = c(compID, tempID)
		
		logVar2 = log(input.table[,rate]+dNadd)
		#skip scaling, even though I saw some code doing that for independent variables
		logVar1 = log(input.table[,compVar])
		
		#I saw some similar code with all varibles (but that wouldn't be the main results, and would have been over-fit if it was)
		cor.coeff = cor(logVar1,logVar2,use="pairwise.complete.obs")
		cor.coefficient = c(cor.coefficient, cor.coeff)
		fit = glm(logVar2 ~ logVar1)
		result = summary(fit)$coefficients
		glm.pvalue = c(glm.pvalue, result[2,4])
		
		plot.file = paste("plots_divergence/",output.prefix,"_",tempID,".png",sep="")
		png(plot.file)
		plot(logVar1,logVar2,
			xlab=paste("Scaled-ln ",compVar,sep=""), ylab = paste("ln(",rate,"+",dNadd,")",sep=""),
			pch=16, main = paste("r = ",signif(cor.coeff, digits=3)," , glm p = ",round(result[2,4], digits=3),sep=""))
		abline(fit, col="red")
		dev.off()
		}#end 	for(i in 1:length(altVars))
}#end for (i in 1:length(depVars))

#adjust for gene length (2 variable comparison)
logInvLength = log(input.table[,"inv.length"])

for (i in 1:length(depVars)){
	rate = depVars[i]
	tempID = paste("fRNAcov_versus_",rate,"_adj_length",sep="")
	compID = c(compID, tempID)
	
	logVar2 = log(input.table[,rate]+dNadd)
	logVar1 = log(input.table[,"coverage"])
	
	#I saw some similar code with all varibles (but that wouldn't be the main results, and would have been over-fit if it was)
	cor.result = partial.cor.test( logVar2, logVar1, logInvLength)
	cor.coeff = cor.result$estimate
	cor.coefficient = c(cor.coefficient, cor.coeff)
	fit = glm(logVar2 ~ logVar1 + logInvLength)
	result = summary(fit)$coefficients
	glm.pvalue = c(glm.pvalue, result[2,4])
	
	for(j in 1:length(altVars)){
		compVar = altVars[j]
		
		tempID = paste(compVar,"_versus_",rate,"_adj_length",sep="")
		compID = c(compID, tempID)
		
		logVar2 = log(input.table[,rate]+dNadd)
		#skip scaling, even though I saw some code doing that for independent variables
		logVar1 = log(input.table[,compVar])
		
		#I saw some similar code with all varibles (but that wouldn't be the main results, and would have been over-fit if it was)
		cor.result = partial.cor.test( logVar2, logVar1, logInvLength)
		cor.coeff = cor.result$estimate
		cor.coefficient = c(cor.coefficient, cor.coeff)
		fit = glm(logVar2 ~ logVar1 + logInvLength)
		fit = glm(logVar2 ~ logVar1)
		result = summary(fit)$coefficients
		glm.pvalue = c(glm.pvalue, result[2,4])
		}#end 	for(i in 1:length(altVars))
}#end for (i in 1:length(depVars))

glm.fdr = p.adjust(glm.pvalue, "fdr")

stat.table = data.frame(compID, cor.coefficient, glm.pvalue, glm.fdr)
stat.file = paste(output.prefix,"_stats.txt",sep="")
write.table(stat.table, stat.file, quote=F, sep="\t", row.names=F)