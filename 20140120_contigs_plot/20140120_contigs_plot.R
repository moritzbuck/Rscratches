#frame=read.csv("/home/moritz/GEFES/views/projects/acI/binning/frame.csv")
require(gplots)

contigs.dists=function(pool,frame,cutoff.length=5000,cutoff.cov=1, max.tet.freq=0.15,scheme="gc",xlog=TRUE,ylog=TRUE)
    {
        if(length(pool)==1)
            coverage=frame[,pool]
        else
            coverage=apply(frame[,pool],1,sum)
        
        freqs=grep("freq", colnames(frame))
        
        goods=frame$length > cutoff.length & apply(frame[,freqs],1,function(x) sum(x > max.tet.freq) ==0 ) & coverage > cutoff.cov
        if(scheme=="gc")
            {
                ramp=rev(rich.colors(n=100,"blues"))

                factor=frame$gc_content[goods]
            }
        if(scheme=="pca")
            {
                ramp=rich.colors(n=100)

                pca.tetras=prcomp(as.matrix(frame[goods,freqs]))
                pred.tetras=predict(pca.tetras)
                factor=pred.tetras[,1]
            }

        quants=quantile(factor,seq(0,1,0.01))
        cols=ramp[sapply(factor,function(x) max(which(x >= quants)))]

        x=if(xlog) log10(frame$length[goods]) else frame$length[goods]
        y=if(ylog) log10(coverage[goods]) else coverage[goods]
        
        plot(x,y,pch=20,cex=0.8,xlab="length of contig",ylab="sum of coverage",col=cols)
   

    }



cov.pca=function(pool,frame,cutoff.length=5000,max.tet.freq=0.15,scheme="cov")
    {
        if(length(pool)==1)
            coverage=frame[,pool]
        else
            coverage=apply(frame[,pool],1,sum)

        freqs=grep("freq", colnames(frame))

        goods=frame$length > cutoff.length & apply(frame[,freqs],1,function(x) sum(x > max.tet.freq) ==0 )
        coverage=coverage[goods]
        gccontent=frame[goods,"gc_content"]
 
        freqs=grep("freq", colnames(frame))
        pca.tetras=prcomp(as.matrix(frame[goods, freqs]))
        pred.tetras=predict(pca.tetras)
        vars=paste(summary(pca.tetras)$importance[2,c(1,2)]*100,col="%",sep="")

        main=if(length(pool)==1) colnames(frame)[pool] else paste(colnames(frame[,pool]),collapse="/")
        main=paste(scheme,main,collapse=": ", sep=": ")
        if(scheme=="cov")
            {
                ramp=rev(rich.colors(n=100,palette="blues"))
                factor=coverage

            }
        if(scheme=="gc")
            {
                ramp=rev(rich.colors(n=100))
                factor=gccontent
            }

        cov.quants=quantile(factor,seq(0,1,0.01))
        cols=ramp[sapply(factor, function(x) max(which(x >= cov.quants)))]        

        plot(pred.tetras[,1],pred.tetras[,2], pch=20, cex=0.8,col=cols, xlab=paste("First PC (",vars[1],")",sep=""),ylab=paste("Second PC (",vars[2],")",sep=""),main=main )
    }


all.plots=function(CoL,MtF)
{
    cov.pca(pool=c(4:6),frame=frame,cutoff.length=CoL, max.tet.freq=MtF)
    cov.pca(pool=4,frame=frame,cutoff.length=CoL, max.tet.freq=MtF)
    cov.pca(pool=5,frame=frame,cutoff.length=CoL, max.tet.freq=MtF)
    cov.pca(pool=6,frame=frame,cutoff.length=CoL, max.tet.freq=MtF)
    cov.pca(pool=c(4:6),frame=frame,cutoff.length=CoL, max.tet.freq=MtF, scheme="gc")
    cov.pca(pool=4,frame=frame,cutoff.length=CoL, max.tet.freq=MtF, scheme="gc")
    cov.pca(pool=5,frame=frame,cutoff.length=CoL, max.tet.freq=MtF, scheme="gc")
    cov.pca(pool=6,frame=frame,cutoff.length=CoL, max.tet.freq=MtF, scheme="gc")

}
