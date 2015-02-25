library(Biobase)

bit_filt = function(dat,bit){
    apply(dat,2,function(x) sapply(strsplit(x,","),function(y) sum(as.numeric(y[-1]) > bit)))
}

plot(sapply(1:500,function(x) sum(apply(bit_filt(data_good,x)==1,1,sum)==ncol(data_good))), type = "l")

var_filt = function(dat, cutof = 27){
    bitscores = apply(dat,1,function(x) unlist(sapply(strsplit(x,","),function(y) as.numeric(y[-1]))))
    bitscores = sapply(bitscores,function(l) {names(l) = as.vector(sapply(names(l),substr,1,7)); sort(l)})
    km = lapply(bitscores, function(x) kmeans(log10(x),2))
    hits = sapply(names(km), function(k) names(which(km[[k]]$cluster == which.max(km[[k]]$centers) | bitscores[[k]] > cutof)))
    sapply(hits,function(h) sum(isUnique(h))) == ncol(dat)
}


plotage =function(se)
    {
        for(i in se)
            {
                cols =(km[[i]]$cluster == which.max(km[[i]]$centers) | bitscores[[i]] > 27) +1
                cols[km[[i]]$cluster != which.max(km[[i]]$centers)] = 3
                plot(sort(log10(bitscores[[i]])), col =  cols, pch = 19, ylab="", main=names(km)[i], xlab=paste(length(hits[[i]]), length(unique(hits[[i]])), sep="/"))
                abline(h=log10(27))
            }
    }
