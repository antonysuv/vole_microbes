#!/usr/bin/env Rscript
library('rjags')
library('runjags')


load("micro_matrix_for_anton_Oct2.RData")
y=rbind(y.n,y.cs,y.s)
K=apply(y,1,sum)
Nsample=nrow(y)
Notu=ncol(y)
n_mtnuc=as.vector(apply(meta.n[ ,c("mtDNA","Nuclear")],1,paste,collapse=""))
ns_mtnuc=as.vector(apply(meta.cs[ ,c("mtDNA","Nuclear")],1,paste,collapse=""))
s_mtnuc=as.vector(apply(meta.s[ ,c("mtDNA","Nuclear")],1,paste,collapse=""))

regions=c(ifelse(n_mtnuc=="NN",1),ifelse(ns_mtnuc=="NN",2,ifelse(ns_mtnuc=="NS",3,ifelse(ns_mtnuc=="Shybrid",4,ifelse(ns_mtnuc=="SN",5,ifelse(ns_mtnuc=="SS",6,0))))),ifelse(s_mtnuc=="SS",7))
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
   
    mu_climate=c(-1.2191640,0.6317949,1.1426795,0.8426714,0.3857724,1.1426795,0.4595694)
    sigma_climate=c(0.4465319,0.3266884,0.000001,0.000001,0.2915717,0.000001,0.6069088)
    mu_distance=c(0.02421556,0.54637980,0.54637980,0.54637980,0.54637980,0.54637980,-0.77520118)
    sigma_distance=c(0.2995745,1.2992233,1.2992233,1.2992233,1.2992233,1.2992233,0.3495740)
    err ~ dnorm(0,3)
}
" 



mydata=list(y=y,K=K,Nsample=Nsample,Notu=Notu,regions=regions,Nregion=Nregion)


fit=run.jags(
                model = rjags.model,
                data = mydata,
                n.chains = 3,
                adapt =   500,
                burnin = 1000,
                sample =  10000,
                monitor = c('p','alpha','a','b','d') 
            )
save(fit, file = "MCMC_mtDNAntDNA.RData")


#Test simulation
mydata=list(
                y=t(cbind(rmultinom(10,size=rpois(10,100),prob=c(0.15,0.1,0.7,0.05)),
                    rmultinom(10,size=rpois(10,100),prob=c(0.1,0.1,0.1,0.7)),
                    rmultinom(10,size=rpois(10,100),prob=c(0.6,0.3,0.05,0.05)))),
                K=apply(y,1,sum),
                Nsample=nrow(y),
                Notu=ncol(y),
                regions=rep(c(1,2,3),each=10),
                Nregion=length(unique(regions))
                #mtdna=sample(c(1,2),size=30,replace=T),
                #Nmtdna=length(unique(mtdna))
    )