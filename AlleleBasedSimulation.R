#install.packages('mvtnorm')# for multivariate normal distribution
#library(mvtnorm)# for multivariate normal distribution
rm(list=ls())
n<-10 #loci
a<-5 #alleles per loci
Htheta<-0 #mean of fitness function
omega2<-1 #variance of fitness function
K<-100 #carrying capacity
Tmax<-1000 #time steps

#BEGIN Initialize a population
MO<-8 #maximum offspring per breeding pair.  Binomial at line 53
B<-8 # for mutation rate below
Ne<- 2*B/(2*B-1)*K
mu<-0.001 # mutation rate As per Holt et al.
alpha2<-0.05
Vs<-2 #See Burger and Lynch 1995, p. 152
sig_SHC<-4*n*mu*alpha2*Ne/(1+ (alpha2*Ne/Vs))
sigH<-rep(sig_SHC,a)
muH<-matrix(data=0,nrow=n,ncol=a)
pool<-rmvnorm(n=n,mean=muH[1,],sigma=diag(sigH))
sigp<-1 
Hg<-array(NA,c(n,2,K))#will hold alleles for all individuals in the population.
Hp<-vector()#will hold phenotypes for all individuals in the population.
for (i in 1:K){
  j<-cbind(sample.int(n=5,size=10,replace=TRUE),(sample.int(n=5,size=10,replace=TRUE)))
  for (k in 1:n){
    Hg[k,,i]<-pool[k,j[k,]]#alleles for 10 loci, N individuals
  }
  Hp[i]<-rnorm(1,sum(Hg[,,i]),sigp)#Phenotype for all individuals
}
#END Initialize a population

#BEGIN simulation
meanP<-vector()#To hold mean phenotype through time. Output
N<-vector()# to hold abundance through time.  Output
Het<-vector()# to hold population-level genetic heterozygocity.  Output
N[1]<-K# 
for (t in 1:Tmax){
  moms<-sample.int(n=N[t],size=floor(N[t]/2))
  dads<-seq(1:N[t])[-moms]
  dads<-dads[1:length(moms)]#if odd number population size
  parents<-cbind(moms, dads)
  offspring<-0
  Hg2<-array(NA,c(n,2,K*MO))# '2' because diploid
  Hp2<-vector()
  Hsurv<-vector()
  if (t>500){
    Htheta<-1.5
  }
  
  for (i in 1:dim(parents)[1]){
    o<-rbinom(n=1,size=MO,prob=0.5)#number of offspring
    if (o>0) {
      for (j in 1:o){
        offspring<-offspring+1
        for (k in 1:2){
          seg<-sample.int(2,n,replace=TRUE)#inherit from mom or dad
          Hg2[,k,offspring]<- rowSums(Hg[,,parents[i,k]]*cbind(1-(seg-1),seg-1))
          #mutation
          if (as.numeric(runif(1)<(n*mu))==1) {
            u<-sample.int(5,1)
            Hg2[u,k,offspring]<-Hg2[u,k,offspring]+rnorm(1,0,sqrt(alpha2))
          }
        }
      }
    }
  }
  Hg2<-Hg2[,,!is.na(Hg2[1,1,])]#Remove unborn offspring
  Hp2<-rnorm(n=dim(Hg2)[3],mean=colSums((colSums(Hg2))),sd=sigp)
  Hsurv<-rbinom(n=dim(Hg2)[3],size=1,exp(-(Hp2-Htheta)^2/(2*omega2)))
  Hg<-Hg2[,,Hsurv==1]
  Hp<-Hp2[Hsurv==1]
  
  #Apply density-dependence
  if (dim(Hg)[3]>K){
    surv<-sample.int(n=dim(Hg)[3],size=K)
    Hg<-Hg[,,surv]
    Hp<-Hp[surv]
  }
  meanP[t+1]<-mean(Hp)
  N[t+1]<-dim(Hg)[3]
  
  #BEGIN calculate heterozygocity
  H<-matrix(NA,dim(Hg)[3],n)
  for (i in 1:dim(Hg)[3]){
    for (j in 1:n){
      H[i,j]<-Hg[j,1,i]==Hg[j,2,i]
    }
  }
  Het[t+1]<-1-mean(H)
  #End calculate heterozygocity
  
}#END simulation (time)

par(mfrow = c(3, 1),     #
    oma = c(2, 1, 0, 1), # two rows of text at the outer left and bottom margin
    mar = c(3, 2, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA) 
plot(seq(2:t),meanP[2:t],'l',xlab='',ylab='Mean Phenotype',xaxt='n',col=2)
plot(seq(2:t),Het[2:t],'l',xlab='',ylab='Heterozygocity',xaxt='n',col=3,ylim=c(0,1))
plot(seq(2:t),N[2:t],'l',xlab='Time',ylab='Abundance',col=4,ylim=c(0,K*1.1))

