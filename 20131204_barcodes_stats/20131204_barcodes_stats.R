library(Biotrings)

#reads a file with 16 nucleotide adapters (one line, one adapter)
adapters=read.table("adapters.txt")$V1
adapters=as.vector(adapters)

#clean-up
adapters=adapters[nchar(adapters)==16]

#count all adapters
counts=table(as.factor(adapters))

#more clean up
counts=counts[sapply(names(counts),function(x) !grepl("[^NTGCA]",x))]

#counts both halves of adapters
counts.right=table(substr(adapters,9,16))
counts.left=table(substr(adapters,1,8))

#our adapters of interest		  
rights=c("TAGATCGC","CTCTCTAT","TATCCTCT","AGAGTAGA","GTAAGGAG","ACTGCATA","AAGGAGTA","CTAAGCCT")
lefts=c("TAAGGCGA","CGTACTAG","AGGCAGAA","TCCTGAGC","GGACTCCT","TAGGCATG","CTCTCTAC","CAGAGAGG")
lefts=as.vector(sapply(lefts,function(x) toString(reverseComplement(DNAString(x)))))

lefts.steph=c("GTGGCCTT","GTTTCGGA","CGTACGTA","GAGTGGAT","ACTGATAT","ATTCCTTT")
rights.steph="TCTTTCCC"

#adapters where both sides match
goods.counts=counts[paste(lefts,rights,sep="")]

#adapters where one side matches
goods.counts.left=counts.left[lefts]
goods.counts.right=counts.right[rights]

#adapters where at least one half matches
one.half.good=sapply(1:8,function(x) {
sum(counts[grepl(rights[x],names(counts)) | grepl(lefts[x],names(counts))])
})

#adapters with an n in the third position of the left and the rest matches perfectly 
n.ed=counts[sapply(paste(lefts,rights,sep=""),function(x) {substr(x,3,3)<-"N";x})]

#the right or left side has been swapped with an other of our adapters
swaps.right=sapply(1:8,function(x) sum(counts[paste(lefts[x],rights[setdiff(1:8,x)],sep="")],na.rm=T))
swaps.left=sapply(1:8,function(x) sum(counts[paste(lefts[setdiff(1:8,x)],rights[x],sep="")],na.rm=T))

#the right or left side has been swapped with one of Steph's adapters
swaps.right.steph=sapply(1:8,function(x) sum(counts[paste(lefts[x],rights.steph,sep="")],na.rm=T))
swaps.left.steph=sapply(1:8,function(x) sum(counts[paste(lefts.steph[1:6],rights[x],sep="")],na.rm=T))

#counts all adapters sequences with one and only one single nucleotide mismatch (include the N on third base)
smm.counts=sapply(paste(lefts,rights,sep=""),function(x)
  {
    sum(sapply(1:16,function(y)
               {
                 prefix=if(y!=0) substr(x,1,y-1) else ""
                 suffix=if(y!=16)substr(x,y+1,16) else ""
                 rege=paste(prefix,"[^",substr(x,y,y),"]",suffix,sep="")
                 sum(counts[sapply(names(counts),function(x) grepl(rege,x))])
               }))
  })

#summarize table
results=data.frame(perfect.match=goods.counts,single.nucleotide.mismatch=smm.counts-n.ed,with.n.in.third=n.ed,barcode.swaps=(swaps.left+swaps.right)/2,steph.swaps=(swaps.left.steph+swaps.right.steph)/2,others=one.half.good-goods.counts-smm.counts-(swaps.left+swaps.right)/2-(swaps.left.steph+swaps.right.steph)/2,one.half.good)
row.names(results)=c("Alex_1","Alex_2","Alex_3","Alex_4","Alex_5","Alex_6","Alex_7","Alex_8")
results$sequence=paste(lefts,rights,sep="")

#normalize to the number of adapters where at least one half matches perfectly
results[,1:6]=results[,1:6]/results[,7]

#plot and outp table
pdf("20131204_barcodes_stats.pdf")
barplot(t(results[,1:6]),horiz=F,legend.text=T,las=2,args.legend=list(x=18,y=0.5,bty="n"),xlim=c(0,18),col=3:11)
dev.off()
write.table(results,file="20131204_barcodes_stats.csv",sep=",")

