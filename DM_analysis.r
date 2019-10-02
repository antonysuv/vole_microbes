library('rjags')
library('runjags')



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
        a[r] ~ dnorm(0,3)
        b[r] ~ dnorm(0,3)
        d[r] ~ dnorm(0,3) 
        climate[r]~dnorm(mu_climate[r],sigma_climate[r])
        distance[r]~dnorm(mu_distance[r],sigma_distance[r])
        for (otu_index in 1:Notu)
        { 
           alpha[r,otu_index] = exp(a[r]+b[r]*climate[r]+d[r]*distance[r]+err) 
        }    

    } 
   
    mu_climate=c(1,2,3)
    sigma_climate=c(1,2,10)
    mu_distance=c(1,2,3)
    sigma_distance=c(1,2,10)
    err ~ dnorm(0,3)
}
" 




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

fit=run.jags(
                model = rjags.model,
                data = mydata,
                n.chains = 3,
                adapt =   500,
                burnin = 300,
                sample =  20000,
                monitor = c('p','alpha','a','b') 
            )
