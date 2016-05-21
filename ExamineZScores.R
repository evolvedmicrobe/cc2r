library(cc2r)
library(pbbamr)
library(ggplot2)
library(reshape2)
setwd("/Users/nigel/pacbio/echidnaTrain")
bam_name = "3150194-0001.aligned.bam"
fasta_name = "/Users/nigel/pacbio/AllReferences.fna"
ind = loadPBI(bam_name, loadSNR = TRUE)
ind$alnLength = ind$tend - ind$tstart
head(ind)
ind$acc = 1 - (ind$mismatches + ind$inserts + ind$dels) / ind$alnLength
indf = ind[ind$alnLength > 500,]
ind2 = ind[ind$alnLength > 500 & ind$alnLength < 1200,]
inds = ind2[sample(nrow(ind2), 5000),]

alns = loadDataAtOffsets(inds$offset, bam_name, fasta_name)

mkString <- function(x) {
  paste(as.character(x[x!="-"]), sep="", collapse ="")
}
snrs = inds[,grep("snr", colnames(inds))]
getScores <- function(i) {
  aln = alns[[i]]
  tpl = mkString(aln$ref)
  read = mkString(aln$read)
  pw = aln$pw[aln$pw!=0]
  pw[pw > 3] = 3
  snr = as.numeric(snrs[i,])
  mdl = "S/P1-C1"
  getScore(read, tpl, mdl, pw, snr)
}
data = lapply(1:nrow(inds), getScores)
data = do.call(rbind, data)
data$Z = (data$ll - data$mean) / sqrt(data$var)
ggplot(data, aes(x=Z)) + geom_density(fill="blue") + theme_bw() + labs(x="Z Score", title="Sequel PW Model Z Scores")

d2 = cbind(data, inds)
ggplot(d2, aes(x=acc, y=Z)) + geom_point() + geom_smooth() + labs(x="Aligned Accuracy", y="Z Score", title="Z-Score vs. Accuracy")

tpl = "AAACAGATCACCCGCTGAGCGGGTTATCTGTT"
mdl = "P6-C4"
res = getTransitionParameters(tpl, mdl, snrs=rep(6,4))
res$bp2 = factor(strsplit(tpl, "")[[1]])
head(res)

tpl = "CCCCC"
mdl = "P6-C4"
snrs = rep(6,4)
getVal <- function(snr) {
  lsnr = snrs
  lsnr[2] = snr
  tp = getTransitionParameters(tpl, mdl, snrs=lsnr)
  tp$snr = snr
  tp[3,]
}
res = lapply(seq(4,21, .5), getVal)
vals = do.call(rbind, res)
tp = melt(vals, id.vars = c("bp", "snr"))
ggplot(tp, aes(x=snr, y=value)) + geom_line() + facet_wrap(~variable, scales="free") + labs(y="Probability", x="SNR")
