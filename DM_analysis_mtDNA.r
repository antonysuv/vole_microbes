#!/usr/bin/env Rscript
library('rjags')
library('runjags') 


load("micro_matrix_for_anton_Oct2.RData")
y=rbind(y.n,y.cs,y.s)
K=apply(y,1,sum)
Nsample=nrow(y)
Notu=ncol(y)
regions=c(ifelse(meta.n$mtDNA=="N",1,0),ifelse(meta.cs$mtDNA=="N",2,3),ifelse(meta.s$mtDNA=="S",4,0))
Nregion=length(unique(regions))

#Climate data prep
clim=c(meta.n$ClimPC1,meta.cs$ClimPC1,meta.s$ClimPC1)
#Normalization of climate
clim_norm=(clim-mean(clim))/sd(clim)
#Mu for empirical Bayes prior constraction  
clim_mu=as.vector(tapply(clim_norm,regions,mean))
#Sigma for empirical Bayes prior constraction
clim_sigma=as.vector(tapply(clim_norm,regions,sd))
#Pointer for regions
v_fac=rep(c(1,2,3),times=c(length(n.dist.vec),length(cs.dist.vec),length((s.dist.vec))))

#Genetic distance data prep
d=c(n.dist.vec,cs.dist.vec,s.dist.vec)
d_norm=(d-mean(d))/sd(d)
d_mu=as.vector(tapply(d_norm,v_fac,mean))
d_sigma=as.vector(tapply(d_norm,v_fac,sd))

# Nsample = number of samples (i.e. OTU rows)
# Notu = number of OTUs (i.e. OTU columns)
# regions = numerical vector corresponding to region assignments 
# a,b,d = regression parameters

rjags.model="
model
{
    for (i in 1:Nsample)
    {

        y[i,1:Notu] ~ dmulti(p[regions[i],1:Notu],K[i])
        
    }
    for (r in 1:Nregion)
    {
        p[r,1:Notu] ~ ddirch(alpha[r,1:Notu])
        for (otu_index in 1:Notu)
        {
            alpha[r,otu_index] = exp(a[r,otu_index]+b[r,otu_index]*climate[r,otu_index]+d[r,otu_index]*distance[r,otu_index]+err)+0.01
            a[r,otu_index] ~ dnorm(0,5)
            b[r,otu_index] ~ dnorm(0,5)
            d[r,otu_index] ~ dnorm(0,5) 
            climate[r,otu_index]~dnorm(mu_climate[r],sigma_climate[r])
            distance[r,otu_index]~dnorm(mu_distance[r],sigma_distance[r])
            
        }

    } 
   
    mu_climate=c(-1.2191640,0.8361488,0.7075440,0.4595694)
    sigma_climate=c(0.4465319,0.3590162,0.4321223,0.6069088)
    mu_distance=c(0.02421556,0.54637980,0.54637980,-0.77520118)
    sigma_distance=c(0.2995745,1.2992233,1.2992233,0.3495740)
    err ~ dnorm(0,3)
}
" 



mydata=list(y=y,K=K,Nsample=Nsample,Notu=Notu,regions=regions,Nregion=Nregion)


fit=run.jags(
                model = rjags.model,
                data = mydata,
                n.chains = 1,
                adapt =   50,
                burnin = 50,
                sample =  500,
                monitor = c('p','alpha','a','b','d') 
            )
save(fit, file = "MCMC_mtDNA.RData")


sample(K,size=1)

all=rbind(as.data.frame(fit$mcmc[[1]][1:100,]),as.data.frame(fit$mcmc[[2]][1:100,]),as.data.frame(fit$mcmc[[3]][1:100,]))
reg1N=all[,seq(1,1237,4)]
reg2NC=all[,seq(2,1238,4)]
reg3SC=all[,seq(3,1239,4)]
reg4S=all[,seq(4,1240,4)] 


reg1N=all[,seq(3721,4957,4)]
reg2NC=all[,seq(3722,4958,4)]
reg3SC=all[,seq(3723,4959,4)]
reg4S=all[,seq(3724,4960,4)] 

reg1N=all[,seq(4961,6197,4)]
reg2NC=all[,seq(4962,6198,4)]
reg3SC=all[,seq(4963,6199,4)]
reg4S=all[,seq(4964,6200,4)] 


reg1N=all[,seq(2877,3717,4)]
reg2NC=all[,seq(2878,3718,4)]
reg3SC=all[,seq(2879,3719,4)]
reg4S=all[,seq(2880,3720,4)] 



all_postpred=c()
for (i in 1:1000)
{
    zran=sample(K,1)
    r1=t(rmultinom(1,zran,reg1N[1,]))
    zran=sample(K,1)
    r2=t(rmultinom(1,zran,reg1N[1,]))
    zran=sample(K,1)
    r3=t(rmultinom(1,zran,reg1N[1,]))
    zran=sample(K,1)
    r4=t(rmultinom(1,zran,reg1N[1,]))
    all_postpred=rbind(all_postpred,r1,r2,r3,r4)
}    
co=rep(c(1,2,3,4),times=1000)    
          