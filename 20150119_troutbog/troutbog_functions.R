library(ggplot2)
library(reshape2)

plot.tseries = function(bins)
    {
        p.data = melt(normed.data.tseries, c("ass_size","GC", "clusts"))
        p.data = p.data[p.data$clusts %in% bins,]
        p.data$variable = as.character(p.data$variable)
        p.data=cbind(p.data, meta.data[p.data$variable,])
        ggplot(p.data, aes(x=day, y=value, col=as.factor(clusts)))+
            geom_line()+
                facet_wrap("year")+
                    scale_x_date()+
                        scale_y_log10()
    }

plot.tseries.2 = function(bins)
    {
        p.data = melt(cbind(contig=row.names(normed.data.covs),normed.data.covs)[normed.data.covs$clusts %in% bins, ], c("contig","GC","length", "clusts"))
        p.data$variable = as.character(p.data$variable)
        p.data=cbind(p.data, meta.data[p.data$variable,])
        p.data = p.data [ p.data$value != 0,]
        ggplot(p.data, aes(x=day, y=value, col=contig, ymin=mean(value)-1.96*sd(value), ymax=mean(value)+1.96*sd(value)))+
            geom_line()+
                facet_grid(clusts~year)+
                    scale_x_date()+
                        scale_y_log10()+ theme(legend.position="none") #+geom_ribbon()
    }



plot.mds = function(dd= normed.data.covs , c.list = 1:nrow(dd), method="euclidean")
    {
        mds = metaMDS(t(dd[c.list,val.cols]), distance=method)$points 
        mds = cbind(mds,meta.data[row.names(mds),])
        mds = mds[order(mds$data),]
        ggplot(mds, aes(x=MDS1,y=MDS2,col=year, label=format(day,"%m-%d")))+geom_path()+geom_text()
    }


plot.cov.pcas = function(dd = normed.data.covs , b.list = levels(as.factor(dd$clusts)) , c.list = !is.na(dd$clusts), fact = 1, a.e.s = aes(x=PC1, y=PC2, col=as.factor(clusts), size=length), bins = NULL)
    {
        dd=dd[c.list,]
        dd = dd[dd$clusts %in% b.list,]
        dd$tot = apply(dd[,val.cols], 1,sum)

        pca = predict(prcomp(asinh(dd[,val.cols]/fact)))
        dd[,val.cols] = pca
        colnames(dd)[val.cols]=colnames(pca)
        if(!is.null(bins))
            {
                dd$clusts[!dd$clusts %in% bins] = -1
            }
        ggplot(dd, a.e.s)+geom_point()+scale_size_continuous(trans="log10")

    }



