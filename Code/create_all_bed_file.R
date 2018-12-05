#add blue color to non-ribosomal coding fRNAs, red for ribosomal coding fRNAs, green for other coding fRNAs (such as genes with introns), black for other fRNAs

coding.fRNA.file = "../fRNA_WITH_ribosome_195_genes_expanded.txt"
coding.fRNA.file2 = "../fRNA_wo_ribosome_169_genes_expanded.txt"
all.fRNA.file = "fold_positions_confirmed_with_annotations.txt"
bed.output = "../all_fRNA_track_sacCer1.bed"

#filtered coding fRNAs
coding.fRNA.table = read.table(coding.fRNA.file, head=T, sep="\t")
print(dim(coding.fRNA.table))
folds = c()
for (i in 1:nrow(coding.fRNA.table)){
	temp.folds = unlist(strsplit(as.character(coding.fRNA.table$folds[i]),split=","))
	folds=c(folds,temp.folds)
}#end for (i in 1:length(coding.fRNA.table))
print(length(folds))
folds = unique(folds)
print(length(folds))

#unfiltered coding fRNAs
coding.fRNA.table2 = read.table(coding.fRNA.file2, head=T, sep="\t")
print(dim(coding.fRNA.table2))
folds2 = c()
for (i in 1:nrow(coding.fRNA.table2)){
	temp.folds = unlist(strsplit(as.character(coding.fRNA.table2$folds[i]),split=","))
	folds2=c(folds2,temp.folds)
}#end for (i in 1:length(coding.fRNA.table))
print(length(folds2))
folds2 = unique(folds2)
print(length(folds2))

#all fRNAs?
all.fRNA.table = read.table(all.fRNA.file, head=T, sep="\t")
print(dim(all.fRNA.table))

matched.fRNA = all.fRNA.table$name[match(folds,all.fRNA.table$name)]
print(length(matched.fRNA[is.na(matched.fRNA)]))

matched.fRNA2 = all.fRNA.table$name[match(folds2,all.fRNA.table$name)]
print(length(matched.fRNA2[is.na(matched.fRNA2)]))

#other coding fRNAs (such as genes with introns)
other.coding.fRNA = all.fRNA.table$name[all.fRNA.table$type == "Coding"]
print(length(other.coding.fRNA))
other.coding.fRNA=other.coding.fRNA[-match(folds,other.coding.fRNA)]
print(length(other.coding.fRNA))
other.coding.fRNA=as.character(other.coding.fRNA)

#no unmapped folds ... good!

output.lines = c("browser position chr6:75804-75850")
output.lines=c(output.lines,"track name=EvoFold fRNAs description=\"total EvoFold fRNAs\" itemRgb=\"On\"")
#output.lines=c(output.lines,"track name=EvoFold description=\"non-ribosomal coding fRNA\" color=0,0,255,")
folds2 = sort(folds2)
for (i in 1:length(folds2)){
	fold.arr = unlist(strsplit(folds2[i],split="_"))
	fold.length = all.fRNA.table$length[match(folds2[i],all.fRNA.table$name)]
	fold.arr = c(fold.arr,folds2[i],fold.length,".",fold.arr[2:3],"0,0,255,")
	output.lines = c(output.lines,paste(fold.arr,collapse="\t"))
}#end for (i in 1:length(folds2)

#output.lines=c(output.lines,"track name=EvoFold description=\"ribosomal coding fRNA\" color=255,0,0,")
ribo.folds = folds[-match(folds2,folds)]
ribo.folds = sort(ribo.folds)
for (i in 1:length(ribo.folds)){
	fold.arr = unlist(strsplit(ribo.folds[i],split="_"))
	fold.length = all.fRNA.table$length[match(ribo.folds[i],all.fRNA.table$name)]
	fold.arr = c(fold.arr,ribo.folds[i],fold.length,".",fold.arr[2:3],"255,0,0,")
	output.lines = c(output.lines,paste(fold.arr,collapse="\t"))
}#end for (i in 1:length(ribo.folds)

#output.lines=c(output.lines,"track name=EvoFold description=\"other coding fRNA\" color=0,255,0,")
other.coding.fRNA = sort(other.coding.fRNA)
for (i in 1:length(other.coding.fRNA)){
	fold.arr = unlist(strsplit(other.coding.fRNA[i],split="_"))
	fold.length = all.fRNA.table$length[match(other.coding.fRNA[i],all.fRNA.table$name)]
	fold.arr = c(fold.arr,other.coding.fRNA[i],fold.length,".",fold.arr[2:3],"0,255,0,")
	output.lines = c(output.lines,paste(fold.arr,collapse="\t"))
}#end for (i in 1:length(other.coding.fRNA)

#output.lines=c(output.lines,"track name=EvoFold description=\"non-coding fRNA\" color=0,0,0,")
noncoding.folds = as.character(all.fRNA.table$name[-match(c(folds,other.coding.fRNA),all.fRNA.table$name)])
noncoding.folds = sort(noncoding.folds)
for (i in 1:length(noncoding.folds)){
	fold.arr = unlist(strsplit(noncoding.folds[i],split="_"))
	fold.length = all.fRNA.table$length[match(noncoding.folds[i],all.fRNA.table$name)]
	fold.arr = c(fold.arr,noncoding.folds[i],fold.length,".",fold.arr[2:3],"0,0,0,")
	output.lines = c(output.lines,paste(fold.arr,collapse="\t"))
}#end for (i in 1:length(noncoding.folds)

fileConn=file(bed.output)
writeLines(output.lines, fileConn)
close(fileConn)
