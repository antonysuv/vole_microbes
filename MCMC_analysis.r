library("mcmcr")
library("ggplot2")
library('reshape')
library("phateR")
library("abind")
library("dplyr")
library("Rtsne")

load("micro_matrix_for_anton_Oct2.RData")
y=rbind(y.n,y.cs,y.s)
K=apply(y,1,sum)

#Parameter selection   
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

#Subsample from posterior
possub=function(input,n_samples)
{
    subs=lapply(input,sample_n,n_samples)
    subs_concat=abind(subs,along=1)
    return(subs_concat)
}    

#Generate OTU counts from posterior predictive 
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
#Simulating from posterior predictive 
ppar=mcmcsub(ch,4,"p")
p_ppred=ppred(ppar,2000)
p_sub=possub(ppar,2000)

#Climate parameter
bpar=mcmcsub(ch,4,"b")
b_sub=possub(bpar,2000)

#Genetic distance parameter
dpar=mcmcsub(ch,4,"d")
d_sub=possub(dpar,2000)

colv=c(rep("royalblue1",2000),rep("darkorange",2000),rep("gray55",2000),rep("darkred",2000))


#t-SNE
ts_out=Rtsne(p_ppred[,2:ncol(p_ppred)])
ts_outp=Rtsne(p_sub)
ts_outb=Rtsne(b_sub)
ts_outd=Rtsne(d_sub)

#PHATE
ph_out=phate(p_ppred[,2:ncol(p_ppred)])
ph_outp=phate(p_sub)
ph_outb=phate(b_sub)
ph_outd=phate(d_sub)

#Plotting
quartz(width=6, height=11)
par(mfrow=c(4,2),mar = c(4,5,3,2))
plot(ts_out$Y,col=alpha(colv,0.2),pch=16,cex=1,main="t-SNE",xlab="tSNE1",ylab="tSNE2")
title("Posterior predictive OTU", line = 0.2)
legend("topleft",c("N","CN","CS","S"),fill=cols,bty = "n")
plot(ph_out$embedding,col=alpha(colv,0.2),pch=16,cex=1,main="PHATE")

plot(ts_outp$Y,col=alpha(colv,0.2),pch=16,cex=1,,xlab="",ylab="")
title("Probability OTU", line = 0.3)
plot(ph_outp$embedding,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")


plot(ts_outb$Y,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")
title("Climate", line = 0.3)
plot(ph_outb$embedding,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")


plot(ts_outd$Y,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")
title("Genetic distance", line = 0.3)
plot(ph_outd$embedding,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")

quartz.save("mtDNA.pdf", type = "pdf")


cols=c("royalblue1","darkorange","gray55","darkred","darkgreen")
### MCMC lineage
load("MCMC_lineage.RData")
ch=fit$mcmc
#Simulating from posterior predictive 
ppar=mcmcsub(ch,5,"p")
p_ppred=ppred(ppar,2000)
p_sub=possub(ppar,2000)

#Climate parameter
bpar=mcmcsub(ch,5,"b")
b_sub=possub(bpar,2000)

#Genetic distance parameter
dpar=mcmcsub(ch,5,"d")
d_sub=possub(dpar,2000)

colv=c(rep("royalblue1",2000),rep("darkorange",2000),rep("gray55",2000),rep("darkred",2000),rep("darkgreen",2000))

#t-SNE
ts_out=Rtsne(p_ppred[,2:ncol(p_ppred)])
ts_outp=Rtsne(p_sub)
ts_outb=Rtsne(b_sub)
ts_outd=Rtsne(d_sub)

#PHATE
ph_out=phate(p_ppred[,2:ncol(p_ppred)])
ph_outp=phate(p_sub)
ph_outb=phate(b_sub)
ph_outd=phate(d_sub)

#Plotting
quartz(width=6, height=11)
par(mfrow=c(4,2),mar = c(4,5,3,2))
plot(ts_out$Y,col=alpha(colv,0.2),pch=16,cex=1,main="t-SNE",xlab="tSNE1",ylab="tSNE2")
title("Posterior predictive OTU", line = 0.2)
legend("topleft",c("N","C1","C2","C3","S"),fill=cols,bty = "n")
plot(ph_out$embedding,col=alpha(colv,0.2),pch=16,cex=1,main="PHATE")

plot(ts_outp$Y,col=alpha(colv,0.2),pch=16,cex=1,,xlab="",ylab="")
title("Probability OTU", line = 0.3)
plot(ph_outp$embedding,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")


plot(ts_outb$Y,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")
title("Climate", line = 0.3)
plot(ph_outb$embedding,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")


plot(ts_outd$Y,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")
title("Genetic distance", line = 0.3)
plot(ph_outd$embedding,col=alpha(colv,0.2),pch=16,cex=1,xlab="",ylab="")
quartz.save("Lineage.pdf", type = "pdf")




#Simulation from real data


sim_real=function(input,size)
{
    real_sim=rbind()
    for ( i in (1:size))
    {
        s=sample_n(data.frame(input), 1)
        K=sum(s)
        m=t(rmultinom(1,K,s/K))
        real_sim=rbind(real_sim,m)
    }
    return(real_sim)
}

zz=rbind(sim_real(y.n,2000),sim_real(y.cs,2000),sim_real(y.s,2000))
colv=c(rep("royalblue1",2000),rep("darkorange",2000),rep("darkred",2000))
ts_out=Rtsne(zz)
plot(ts_out$Y,col=alpha(colv,0.2),pch=16,cex=1,main="t-SNE",xlab="tSNE1",ylab="tSNE2")

zz=rbind(cbind(y.n,meta.n[,c("ClimPC1","mtDNA","Nuclear","Lineage")]),cbind(y.cs,meta.cs[,c("ClimPC1","mtDNA","Nuclear","Lineage")]),cbind(y.s,meta.s[,c("ClimPC1","mtDNA","Nuclear","Lineage")]))


