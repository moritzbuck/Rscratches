
library(plyr)
library(reshape2)
library(vegan)

annot_data=read.table("/pica/v3/b2011138/TroutBogHypo/merged_assemblies/year_wise/bins/annotations")

data = read.table("all_contig_coverages.csv",sep="\t",h=T,row.names=1)
data$clusts = read.table("concoct_2000_clustering_gt2000.csv",h=F, sep=",", row.names=1)[row.names(data),]
val.cols=grep("IH",colnames(data))

raw.reads = t(apply(data, 1,function(x) x[val.cols]*x[1]))
data.raw.reads = data
data.raw.reads[,val.cols] = raw.reads

filtered.data.raw.reads = data.raw.reads[!is.na(data.raw.reads$clusts),]
normed.data.raw.reads = filtered.data.raw.reads
normed.data.raw.reads[,val.cols] = apply(normed.data.raw.reads[,val.cols], 2, function(x) x/sum(x))
                         
tsub = ddply(normed.data.raw.reads, .(clusts), summarize, GC=mean(GC), ass_size = sum(length)/1000000, ncontigs = length(length))
tdat = ddply(normed.data.raw.reads, .(clusts), function(x) colSums(x[,val.cols]))

bin.stats = tsub
row.names(bin.stats) = bin.stats$clusts

normed.data.tseries = cbind(tsub[,3:2], tdat[,2:ncol(tdat)], clusts=tdat[,1])
row.names(normed.data.tseries) = normed.data.tseries$clusts

meta.data = read.table("metadata_trout.csv", sep=",", row.names=1, h=T)
meta.data$data = as.Date(meta.data$Date.MM.DD.YY)
meta.data$year=format(meta.data$data, "%Y")
meta.data$day=as.Date(format(meta.data$data, "1970-%m-%d"))

correction.factors = colSums(raw.reads)/mean(colSums(raw.reads))

normed.data.covs = data
normed.data.covs[,val.cols] = t(apply(normed.data.covs[,val.cols], 1, function(x) x/correction.factors))

iqrs = apply(normed.data.covs [, val.cols],1, IQR)
tseries.iqrs = apply(normed.data.tseries[, val.cols],1, IQR)

acIs = c("127","176","121","86","103","205","91","252","122")
od1s = c("173","47","158","178","101","98","153")

len = sort(data$length)
xs = seq(0,10000, 200)
ys = sapply(xs, function(x) sum(len[len > x])/sum(len) )
plot(xs, ys, type="l")
