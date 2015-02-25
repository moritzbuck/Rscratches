library(ape)
library(pvclust)
library(snow)
library(Rmpi)
#library(MASS) # just needed for test dataset 'Boston'

# parallel version of pvclust
#cl <- makeCluster(10, type="MPI")  # 10 number of nodes

#example
#data(Boston)
#boston.pv <- parPvclust(cl, Boston, nboot=1000, weight=TRUE)

frame<-read.csv(file="humic.framme.csv", row.names=1)

tetras=as.matrix(frame[,6:77])
row.names(tetras)=row.names(frame)

bad.comp=apply(tetras,1,function(x) sum(x>0.05))!=0
bad.len=frame$length<1000
tetras.trimmed=tetras[!bad.comp & !bad.len,]
                                        #frame.pv <- parPvclust(cl, t(tetras), nboot=1000, weight=TRUE)

#result_phylo <- as.phylo(frame.pv$hclust)
#write.tree(result_phylo)
deeds=apply(tetras.trimmed,1,function(x) apply(tetras.trimmed,1,function(y) cor(x,y)))
tetras.logged=sapply(tetras.trimmed,function(x) -log10((tetras.trimmed[x,]*(frame[x,"length"]-3)+1)/(frame[x,"length"]-3)))
deeds.loged=apply(tetras.logged,1,function(x) apply(tetras.logged,1,function(y) cor(x,y)))
