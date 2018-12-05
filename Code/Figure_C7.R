
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
print(dim(input.table))
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
print(dim(input.table))

otherVars = c("cai","len","expr","disp","coverage","abund","fold.length.sum")

variables = c(depVars, otherVars)
var.mat = input.table[,match(variables,names(input.table))]
rownames(var.mat)=input.table$orf

#use code similar to what I used as template for COH QC analysis
dir.create("R-base_Clustering_plots")

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())
color.palette=fixed.color.palatte[1:ncol(var.mat)]

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

	cluster.file = paste("R-base_Clustering_plots/Pearson_Dissimilarity_Single_",output.prefix,"_Hierarchical_Clustering_WITH_NA.png",sep="")
		cor.mat = cor(as.matrix(var.mat), use="pairwise.complete.obs")
		dis.mat = 1 - cor.mat
		dist1 = as.dist(dis.mat)
		
hc = hclust(dist1, method="single")
dend1 = as.dendrogram(hc)
png(file = cluster.file)
par(mar=c(5,5,3,10))
dend1 = dendrapply(dend1, colLab, labelColors=color.palette, clusMember=variables) 
a = attributes(dend1) 
attr(dend1, "nodePar") = c(a$nodePar, lab.col = color.palette) 
plot(dend1, horiz=T, main = "Pearson Dissimilarity Hierarchical Clustering")
dev.off()