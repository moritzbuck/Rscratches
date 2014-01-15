table=read.table("minputlength.csv",sep=c("\t"),h=F, skip=1)
samples=apply(table,1,function(x) sub(" ","",grep("[A-Z][0-9]*",x,val=T)))
samples=sapply(samples,function(x) as.vector(x[x!=""]))
things=apply(table,1,function(x) sub(" ", "", grep("[A-Z][0-9]*",x,val=T,invert=T)))
things=sapply(things,function(x) as.vector(x[x!=""]))

label1=sapply(things,function(x) paste(sub(" ","",x),collapse="_"))
samples.levels=levels(as.factor(unlist(samples)))
names(samples)=label1
df.1=as.data.frame(matrix(data=FALSE,nrow=length(label1),ncol=length(samples.levels)),row.names=label1,col.names=samples.levels)
colnames(df.1)=samples.levels

for(l in label1)
    {
        df.1[l,samples[[l]]]=TRUE
    }


label2=levels(as.factor(unlist(things)))
df.2=as.data.frame(matrix(data=FALSE,nrow=length(label2),ncol=length(samples.levels)),row.names=label2,col.names=samples.levels)
colnames(df.2)=samples.levels

for(l in label2)
    {
        samps=levels(as.factor(unlist(samples[sapply(things,function(x) l %in% x)])))
        df.2[l,samps]=TRUE
    }

save(df.2,df.1,file="matrixified_data.Rdata")
write.table(df.2,file="matrixified_data.csv",sep=",")
