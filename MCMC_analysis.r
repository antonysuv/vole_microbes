library("mcmcr")
library("ggplot2")
library('reshape')
library("phateR")

load("micro_matrix_for_anton_Oct2.RData")
y=rbind(y.n,y.cs,y.s)
K=apply(y,1,sum)


mcmcsub=function(input,n_dim,parameter)
{
    input_sub=subset(input,parameters=parameter)
    input_sub=as.data.frame(collapse_chains(input_sub))
    coln=colnames(input_sub)
    mcmclist=list()
    for (i in 1:n_dim)
    {
        nam=paste(parameter,"dim",i,sep="")
        mcmclist[[i]]=input_sub[,grep(coln,pattern=paste(i,",",sep=""))]
    
    }    
    return(mcmclist)
}    


ppred=function(input,sim)
{
    
    posterior_pred=rbind()
    id=1
    for (i in input)
    {
        vec=rbind()
        for (ns in 1:sim)
        {
            n=rnbinom(1,size=9,prob=0.000263)
            ps=sample(1:nrow(i),1)
            vec_sample=cbind(id,t(rmultinom(1,n,i[ps,])))
            vec=rbind(vec,vec_sample)
        }
        posterior_pred=rbind(posterior_pred,vec)
        id=id+1
    }
    return(posterior_pred)
}    



cols=c("royalblue1","darkorange","gray55","darkred")

### MCMC mtDNA
load("MCMC_mtDNA.RData")
ch=fit$mcmc
ppar=mcmcsub(ch,4,"p")
p_ppred=ppred(ppar,2000)

par(mfrow=c(1,2))
#t-SNE
ts_out=Rtsne(p_ppred[,2:ncol(p_ppred)])
colv=ifelse(p_ppred[,1]==1,"royalblue1",ifelse(p_ppred[,1]==2,"darkorange",ifelse(p_ppred[,1]==3,"gray55","darkred")))
plot(ts_out$Y,col=alpha(colv,0.2),pch=16,cex=1,main="t-SNE",xlab="tSNE1",ylab="tSNE2")
legend("topleft",c("N","CN","CS","S"),fill=cols)
#PHATE
ph_out=phate(p_ppred[,2:ncol(p_ppred)])
plot(ph_out$embedding,col=alpha(colv,0.2),pch=16,cex=1,main="PHATE")
legend("topleft",c("N","CN","CS","S"),fill=cols)


cols=c("royalblue1","darkorange","gray55","darkred","darkgreen")
### MCMC lineage
load("MCMC_lineage.RData")
ch=fit$mcmc
ppar=mcmcsub(ch,5,"p")
p_ppred=ppred(ppar,2000)

par(mfrow=c(1,2))
#t-SNE
ts_out=Rtsne(p_ppred[,2:ncol(p_ppred)])
colv=ifelse(p_ppred[,1]==1,"royalblue1",ifelse(p_ppred[,1]==2,"darkorange",ifelse(p_ppred[,1]==3,"gray55",ifelse(p_ppred[,1]==4,"darkred","darkgreen"))))
plot(ts_out$Y,col=alpha(colv,0.2),pch=16,cex=1,main="t-SNE",xlab="tSNE1",ylab="tSNE2")
legend("topleft",c("N","C1","C2","C3","S"),fill=cols)
#PHATE
ph_out=phate(p_ppred[,2:ncol(p_ppred)])
plot(ph_out$embedding,col=alpha(colv,0.2),pch=16,cex=1,main="PHATE")
