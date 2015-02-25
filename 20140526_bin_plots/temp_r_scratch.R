library(rgl)
library(mclust)
library(lattice)

binning.path = "clustering_gt1000.csv"

binner = read.table(binning.path,head=F,as.is=T,sep=",")
frame = read.table("frame.csv",sep=",",h=T,row.names=1)
sccs=read.table("scc_counts.csv",row.names=1, col.names=c("NA","abscent","single","multiple"),h=F)
sccs[apply(sccs,1,sum) == 0,1]=40
t=levels(as.factor(binner$V2))
names(t)=t
clusters=lapply(t,function(x) binner$V1[binner$V2==x])

plot.sccs=function()
    {
        bin.sizes = sapply(clusters,function(x) sum(frame[x,"length"]))[row.names(sccs)]
        total.coverage= sapply(clusters,function(x) sum(apply(frame[x,3:5],1,sum)*frame[x,"length"]))[row.names(sccs)]
        
        plot(40-sccs$abscent~total.coverage/bin.sizes)
    }


plot.clusters=function(x,name=NA)
    {
        if(!is.na(name)) pdf(name)
        for(b in x)
        {
                                        #bin.name = 42
            bin = binner$V1[binner$V2 == b]
            frame = read.table("frame.csv",sep=",",h=T,row.names=1)
            tmp=data.frame(log2(frame[,c(3:5)]+1),lens=frame$length)
            tmp=tmp[apply(tmp[,1:3],1,sum)!=-Inf,]
            cat(b)
            cat(" --- ")
                                        #        clust = Mclust(tmp[row.names(tmp) %in% bin,1:3],G=1:10)
                                        #        heavi=tmp$lens[row.names(tmp) %in% bin]/max(tmp$lens[row.names(tmp) %in% bin])
                                        #        fitnew = do.call("me.weighted",c(list(weights=heavi),clust))
                                        #            for(n in intersect(names(fitnew),names(clust))) clust[[n]] = fitnew[[n]]
            tmp$bin = 0
            tmp[bin,"bin"] = 1
                                        #        tmp[names(predict(clust)$class),"bin"] = predict(clust)$class
            
                                        #        plot3d(tmp[tmp$bin==0,1],tmp[tmp$bin==0,2],tmp[tmp$bin==0,3])
                                        #        plot3d(tmp[tmp$bin!=0,1],tmp[tmp$bin!=0,2],tmp[tmp$bin!=0,3],col=tmp$bin[tmp$bin!=0]+1,add=T,size=5)
            
            bin.sizes=sapply(1:max(tmp$bin),function(x) sum(tmp$lens[tmp$bin ==x]))
            
            ces=rep(0.4,nrow(tmp))
            ces[tmp$bin!=0]=1
            
            pairs(tmp[1:3], col=tmp$bin+1,cex=ces,main=b,pch=20)

            
        }
        if(!is.na(name)) dev.off()
        cat("\n")
    }

recluster=function(x,name=NA,weighted = TRUE, logp1=FALSE,log=TRUE, reclust=TRUE,pca=FALSE,d3=FALSE)
    {
        b=x
        if(!is.na(name)) pdf(name)
        bin = binner$V1[binner$V2 %in% b]
        tmp=data.frame(frame[,c(3:5)],lens=frame$length)
        if(log) tmp=data.frame(log10(frame[,c(3:5)]),lens=frame$length)
        if(logp1) tmp=data.frame(log10(frame[,c(3:5)]+1),lens=frame$length)
        tmp=tmp[apply(tmp[,1:3],1,sum)!=-Inf,]
        tmp$bin = 0
        if(reclust)
            {
                clust = Mclust(tmp[row.names(tmp) %in% bin,1:3],G=1:10)
                if(weighted)
                    {
                        heavi=tmp$lens[row.names(tmp) %in% bin]/max(tmp$lens[row.names(tmp) %in% bin])
                        fitnew = do.call("me.weighted",c(list(weights=heavi),clust))
                        for(n in intersect(names(fitnew),names(clust))) clust[[n]] = fitnew[[n]]
                    }
                tmp[names(predict(clust)$class),"bin"] = predict(clust)$class
            }
        else
            if(!is.null(x))
                tmp[bin,"bin"] = 1


        if(d3)
            {
                plot3d(tmp[tmp$bin==0,1],tmp[tmp$bin==0,2],tmp[tmp$bin==0,3])
                if(!is.null(x))            
                    plot3d(tmp[tmp$bin!=0,1],tmp[tmp$bin!=0,2],tmp[tmp$bin!=0,3],col=tmp$bin[tmp$bin!=0]+1,add=T,size=5)
            }
        
        if(!is.null(x))            
            bin.sizes=sapply(1:max(tmp$bin),function(x) sum(tmp$lens[tmp$bin ==x]))
        else bin.sizes = NULL
        
        ces=rep(0.4,nrow(tmp))
        ces[tmp$bin!=0]=1
        if(pca)
            {
                pcap=prcomp(tmp[1:3])
                tmp[1:3]=predict(pcap)
            }
        
        if(!is.null(x))            
            title = paste("log10 of coverages (bin ",b,")",sep="")
        else title=""

        
        pairs(tmp[1:3], col=tmp$bin+1,cex=ces,main=title ,pch=20)
                 
        if(!is.na(name)) dev.off()
        tmp
    }

coverage.plot=function(bin,samples=3:5,maxi=NA )
{
    par(mfrow = c(length(samples),1))
    if(!is.na(maxi)) ylim=c(0, maxi)
    else ylim = c(0,max(frame[binner[binner$V2==bin,"V1"],samples]))

    for(i in samples)
    {
        plot(frame[binner[binner$V2==bin,"V1"],1:5]$length,frame[binner[binner$V2==bin,"V1"],i],ylim=ylim)
    }
    
}

ff8.plot=function(sample=5,minl=1000, bin.list=NULL, min.bin.size=1000000, comp=NA,crop=TRUE,transform=rank,rep.filter=0.15,ctransf=log10)
    {
        pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
"#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
"#FFFF99", "#B15928")
        tmp=frame[,c(1,2,sample)]
        tmp=tmp[apply(frame[,grep("freq",colnames(frame))],1,max)<rep.filter, ]
        binnerido=binner$V2
        names(binnerido)=binner$V1
        tmp$bin=binnerido[row.names(tmp)]
        

        
        g.bins=names(which(sapply(levels(as.factor(binnerido)), function(x) sum(tmp[tmp$bin == x,]$length)) > min.bin.size ))
        if(!is.null(bin.list)) g.bins=intersect(g.bins,bin.list)
        tmp=tmp[tmp$bin %in% g.bins,]
        if(crop) tmp=tmp[tmp[,3]>1,]
        tmp=tmp[tmp$length > minl,]
        binnerido = as.factor(binnerido[row.names(tmp)])
        palet=colorRampPalette(pal12)(length(levels(binnerido)))
        cols=palet[as.numeric(binnerido[row.names(tmp)])]

        pca=prcomp(frame[row.names(tmp),grep("freq",colnames(frame))])
        sc=predict(pca)
#        print(summary(pca))
        if(is.na(comp)) x = tmp[,"gc_content"]
        else x = sc[row.names(tmp),comp]
        par(mfrow = c(1,2))                
        plot(transform(x),ctransf(tmp[,3]),pch=20,col=cols,cex=log10(tmp[,"length"])-2,xlab="GC-content",ylab="log10 of coverage")
        if(is.na(comp)) legend(x=min(transform(x)) ,y=max(ctransf(tmp[,3])),levels(binnerido),pch=19,col=palet)
        else legend(x=min(transform(x)),y=max(ctransf(tmp[,3])),levels(binnerido),pch=19,col=palet)
        plot(log10(tmp[,"length"]),ctransf(tmp[,3]),pch=20,col=cols,cex=log10(tmp[,"length"])-2,xlab="log10 of length",ylab="log10 of length")
        if(is.na(comp)) legend(x=min(transform(x)) ,y=max(ctransf(tmp[,3])),levels(binnerido),pch=19,col=palet)
        else legend(x=min(transform(x)),y=max(ctransf(tmp[,3])),levels(binnerido),pch=19,col=palet)

    }
 


cluster.contigs=function(lengths,covs)
    {
        clust.x= predict(Mclust(lengths),G=2)$class
        clust.y= predict(Mclust(covs),G=2)$class
        cls=interaction(clust.x,clust.y)
        lens=sapply(levels(cls),function(x) sum(dat$V1[cls==x]))
        cls==names(which.max(lens))
    }

