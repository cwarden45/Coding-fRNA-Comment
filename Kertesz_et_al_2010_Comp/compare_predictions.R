EvoFold.file = "EvoFold_liftOver_sacCer1-to-sacCer2.bed"
PARS.file = "../Kertesz_et_al_2010/GSE22393_processed_merged_PARS_sacCer2_1.bed"
genes.file = "../Kertesz_et_al_2010/sce_transcriptome_global.tab"

EvoFold.table = read.table(EvoFold.file, head=F, sep="\t")
PARS.table = read.table(PARS.file, head=F, sep="\t")
genes.table = read.table(genes.file, head=F, sep="\t")

exon.table = genes.table[genes.table$V5 == "Exon",]
Roman.chr = rep(NA,nrow(exon.table))
Roman.chr[exon.table$V1 == 1]="chrI"
Roman.chr[exon.table$V1 == 2]="chrII"
Roman.chr[exon.table$V1 == 3]="chrIII"
Roman.chr[exon.table$V1 == 4]="chrIV"
Roman.chr[exon.table$V1 == 5]="chrV"
Roman.chr[exon.table$V1 == 6]="chrVI"
Roman.chr[exon.table$V1 == 7]="chrVII"
Roman.chr[exon.table$V1 == 8]="chrVIII"
Roman.chr[exon.table$V1 == 9]="chrIX"
Roman.chr[exon.table$V1 == 10]="chrX"
Roman.chr[exon.table$V1 == 11]="chrXI"
Roman.chr[exon.table$V1 == 12]="chrXII"
Roman.chr[exon.table$V1 == 13]="chrXIII"
Roman.chr[exon.table$V1 == 14]="chrXIV"
Roman.chr[exon.table$V1 == 15]="chrXV"
Roman.chr[exon.table$V1 == 16]="chrXVI"

unstranded.start = rep(NA, nrow(exon.table))
unstranded.start[exon.table$V4 >= exon.table$V3] = exon.table$V3[exon.table$V4 >= exon.table$V3]
unstranded.start[exon.table$V3 > exon.table$V4] = exon.table$V4[exon.table$V3 > exon.table$V4]

unstranded.stop = rep(NA, nrow(exon.table))
unstranded.stop[exon.table$V4 >= exon.table$V3] = exon.table$V4[exon.table$V4 >= exon.table$V3]
unstranded.stop[exon.table$V3 > exon.table$V4] = exon.table$V3[exon.table$V3 > exon.table$V4]

exon.table = data.frame(Roman.chr, unstranded.start, unstranded.stop, exon.table)

library("GenomicRanges")

EvoFold_gr = GRanges(Rle(EvoFold.table$V1),
    		IRanges(start=EvoFold.table$V2, end=EvoFold.table$V3))
    		
PARS_gr = GRanges(Rle(PARS.table$V1),
    		IRanges(start=PARS.table$V2, end=PARS.table$V3))
    		
exon_gr = GRanges(Rle(exon.table$Roman.chr),
    		IRanges(start=exon.table$unstranded.start, end=exon.table$unstranded.stop))
    		
EvoFold_PARS_gr = intersect(EvoFold_gr, PARS_gr)
EvoFold_PARS.table = data.frame(EvoFold_PARS_gr)
EvoFold.All.ID = paste(EvoFold_PARS.table$seqnames, EvoFold_PARS.table$start, sep=":")

#not really useful, if you don't have a PARS score
#EvoFold_exon_gr = intersect(EvoFold_gr, exon_gr)
#EvoFold_exon.table = data.frame(EvoFold_exon_gr)
#EvoFold.exon.ID = paste(EvoFold_exon.table$seqnames, EvoFold_exon.table$start, sep=":")

PARS_exon_gr = intersect(PARS_gr, exon_gr)
PARS_exon.table = data.frame(PARS_exon_gr)
PARS.exon.ID = paste(PARS_exon.table$seqnames, PARS_exon.table$start, sep=":")

EvoFold_exons_PARS_gr = intersect(EvoFold_gr, PARS_exon_gr)
EvoFold_exons_PARS.table = data.frame(EvoFold_exons_PARS_gr)
EvoFold_PARS.exon.ID = paste(EvoFold_exons_PARS.table$seqnames, EvoFold_exons_PARS.table$start, sep=":")

PARS.ID = paste(PARS.table$V1, PARS.table$V2, sep=":")

png("PARS_density.png")
par(mfcol=c(1,2))

#################
### all sites ###
#################

background = PARS.table$V5
print(length(background))
print(quantile(background))

EvoFold = PARS.table$V5[match(EvoFold.All.ID,PARS.ID,nomatch=0)]
print(length(EvoFold))
print(quantile(EvoFold))

den = density(background, from=-8, to=8)
plot(den$x, den$y, type="l", xlab = "PARS Score", ylab = "Density",
		xlim=c(-8,8), col="gray", main = "All Genomic Sites")
mtext(paste("n = ",length(background)," versus n= ",length(EvoFold),sep=""), side=3, line=0.25)
den = density(EvoFold, na.rm=T, from=-8, to=8)
lines(den$x, den$y, type = "l", col="blue")
legend("topleft",legend=c("All","EvoFold"),col=c("gray","blue"),
		lwd=2, ncol=1, cex=0.8)

##################
### Exon sites ###
##################

background = PARS.table$V5[match(PARS.exon.ID,PARS.ID,nomatch=0)]
print(length(background))
print(quantile(background))

EvoFold = PARS.table$V5[match(EvoFold_PARS.exon.ID,PARS.ID,nomatch=0)]
print(length(EvoFold))
print(quantile(EvoFold))

den = density(background, from=-8, to=8)
plot(den$x, den$y, type="l", xlab = "PARS Score", ylab = "Density",
		xlim=c(-8,8), col="gray", main = "Exonic Sites")
mtext(paste("n = ",length(background)," versus n= ",length(EvoFold),sep=""), side=3, line=0.25)
den = density(EvoFold, na.rm=T, from=-8, to=8)
lines(den$x, den$y, type = "l", col="blue")
legend("topleft",legend=c("All","EvoFold"),col=c("gray","blue"),
		lwd=2, ncol=1, cex=0.8)
dev.off()
