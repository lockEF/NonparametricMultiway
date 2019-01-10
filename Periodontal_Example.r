#install.packages('msm')
library(msm)
library(MASS)

###Load data and functions
Data <- read.csv('PeriodontalData.csv')
source('psb.cp.binom.r')
#see DataDictionary.docx for description of variables

M=1160
nG = Data$N 
kG <-Data$Count
UniqueID = unique(Data$id)
##put data into array form
Mi <- 290
Mt <- 4
nGa <- array(dim = c(Mi,Mt))
kGa <- array(dim = c(Mi,Mt))
for(i in 1:Mi){for(j in 1:Mt){
  nGa[i,j] <- nG[Data$id==UniqueID[i]&Data$tooth==j]
  kGa[i,j] <- kG[Data$id==UniqueID[i]&Data$tooth==j]
}}

#organize subject-level covariates
gend.vec<-c();age.vec<-c(); a1c.vec <- c(); smoker.vec <- c(); bmi.vec <- c()

for(i in 1:Mi){
  gend.vec[i] <- Data$gender[Data$id==UniqueID[i]][1]
  age.vec[i] <- Data$age[Data$id==UniqueID[i]][1]
  a1c.vec[i] <- Data$hba1c[Data$id==UniqueID[i]][1]  
  smoker.vec[i] <- Data$smoker[Data$id==UniqueID[i]][1] 
  bmi.vec[i] <- Data$bmi[Data$id==UniqueID[i]][1]  
}
Cov <- cbind(gend.vec,age.vec,a1c.vec,smoker.vec,bmi.vec)
Cov <- scale(Cov,scale=TRUE)
X <- Cov

#run rank-1 probit stick-breaking model 
#(this can take substantial time; the results from a previous run are saved as 'Res.rda'
#and may be loaded below)
Res <- psb.cp.binom(nGa,kGa,X,1,1,draws=20000,burn=10000,Thin=20)
#save(Res, file='Res.rda')
load(file='Res.rda')

#covariate weight posterior means:
rowMeans(Res[[1]])
#region weight posterior means:
rowMeans(Res[[1]])
#95% credible intervals:
apply(Res[[1]], MARGIN=c(1),FUN='quantile',probs=c(0.025,0.975))
apply(Res[[2]], MARGIN=c(1),FUN='quantile',probs=c(0.025,0.975))



#plot results
#install.packages('ggplot2')
#install.packages('gridExtra')
library(ggplot2)
library(gridExtra)
predNames <- c("Sex (F)", "Age", "A1C", "Smoker", "BMI")
predCoeffs <- rowMeans(Res[[1]])
estimates.upper.pred = apply(Res[[1]], MARGIN=c(1),FUN='quantile',probs=c(0.025,0.975))[2,]
estimates.lower.pred = apply(Res[[1]], MARGIN=c(1),FUN='quantile',probs=c(0.025,0.975))[1,]
Summary = data.frame(predNames,predCoeffs)

Summary$predNames <- factor(Summary$predNames, levels = Summary$predNames)

typeNames <- c("Molar", "Premolar", "Canine", "Incisor")
typeCoeffs <- rowMeans(Res[[2]])
estimates.upper.type = apply(Res[[2]], MARGIN=c(1),FUN='quantile',probs=c(0.025,0.975))[2,]
estimates.lower.type = apply(Res[[2]], MARGIN=c(1),FUN='quantile',probs=c(0.025,0.975))[1,]
SumType = data.frame(typeNames,typeCoeffs)
SumType$typeNames <- factor(SumType$typeNames, levels = SumType$typeNames)

png(file='Rank1weights.png',width=700,height=350)
plot1 <- ggplot(Summary, aes(x = predNames, y = predCoeffs)) +  
  labs(x="Covariate",y = "Coefficients")+
  geom_bar(position = position_dodge(), stat="identity", fill="gray") + 
  geom_errorbar(aes(ymin=estimates.lower.pred, ymax=estimates.upper.pred), width=0.50) +
  ggtitle("Covariate weights") + # plot title
  theme_bw() + # remove grey background (because Tufte said so)
  theme(panel.grid.major = element_blank()) +# r
  geom_hline(yintercept = 0)
plot2 <- ggplot(SumType, aes(x = typeNames, y = typeCoeffs)) +  
  labs(x="Tooth Type",y = "Coefficients")+
  geom_bar(position = position_dodge(), stat="identity", fill="gray") + 
  geom_errorbar(aes(ymin=estimates.lower.type, ymax=estimates.upper.type), width=0.50) +
  ggtitle("Tooth type weights") + # plot title
  theme_bw() + # remove grey background (because Tufte said so)
  theme(panel.grid.major = element_blank()) +# r
  geom_hline(yintercept = 0)
grid.arrange(plot1,plot2,ncol=2)
dev.off()

##Now show continuum of distributions
res <- Res
it=500
plot(Res$ppVec[1:10,500],Res$Vprobs[1:10,1,500])
res$Int[,,500]
#save results for latent = -3,0, and 3
latentvals <- c(-3,-1,1,3)
Mi.test <- 2
props.sim <- list()
nGa.test <- matrix(rep(c(48,48,24,48),Mi.test),nrow=Mi.test,ncol=Mt,byrow=TRUE)
for(jj in 1:length(latentvals)){ 
  props.sim[[jj]] <- matrix(nrow=1000,ncol=4) 
  for(it in 1:1000){
    latent <- latentvals[jj]
    #med=0
    V_new <-array(dim=c(Mi.test,4,K))
    Z_new <- array(dim=c(Mi.test,4,K-1))
    ZZ_new <- list(c(1),t(res$VregionsZ[,,it]),t(res$VprobsZ[,,it]))
    ZZ_new[[1]] = matrix(rnorm(Mi.test*Rz,0,sd=sqrt(res$sigmaVec[it])),ncol=Rz)
    for(r in 1:Rz){
      Zr = lapply(ZZ_new, function(x) x[,r])
      ZZarray_new = array(apply(expand.grid(Zr), 1, prod), dim=c(Mi.test,Mt,K-1))
    }
    covvals <- latent*res$Vregions[,1,it] %*% t(res$Vprobs[,1,it])
    for(k in 1:(K-1)){ for(j in 1:4){
      #Z_new[,j,k] = t(Xtest%*%B[,j,k])+ZZarray_new[,j,k]}}
      Z_new[,j,k] = covvals[j,k] + res$IntVec[j,k,it]+ZZarray_new[,j,k]}}
    V_new[,,1:(K-1)] <- pnorm(Z_new)	
    V_new[,,K]= 1
    Pi_new=V_new
    Pi_new[,,2] =  V_new[,,2]*(1-V_new[,,1])
    for(k in 3:K){ 
      Pi_new[,,k] = V_new[,,k]*apply(1-V_new[,,1:(k-1)],FUN=prod,MARGIN=c(1,2))
    }
    Probs <- matrix(nrow=Mi.test,ncol=Mt)
    for(i in 1:Mi.test){for(j in 1:4){ 
      Probs[i,j] <- res$ppVec[rmultinom(1,1,prob=Pi_new[i,j,])==1,it]
    }}
    
    kGa.sim <- matrix(nrow=Mi.test,ncol=Mt)
    for(j in 1:4){kGa.sim[,j] <- rbinom(Mi.test,nGa.test[,j],prob=Probs[,j])}
    #props.sim_m3 <- kGa.sim/nGa.test
    #props.sim_0 <- kGa.sim/nGa.test
    props.sim[[jj]][it,] <- kGa.sim[1,]/nGa.test[1,]
  }}

hist(props.sim[[1]])
hist(props.sim[[2]])
hist(props.sim[[3]])
freqs <- matrix(nrow=length(latentvals),ncol=25)
for(i in 0:24){ for(j in 1:length(latentvals)){
  freqs[j,i+1] = mean(props.sim[[j]]> i/24 - 0.001 & props.sim[[j]]< i/24 +1/48+ 0.001)
}}
png(file='DistributionContinuum.png',width=600,height=350)
plot(c(0,1/48,c(1:24)/24),c(freqs[1,1], freqs[1,]), type='S', lwd=3, col = 'gray90', ylim=c(0,0.2), xlab = 'Proportion diseased',ylab = 'Frequency')
points(c(0,1/48,c(1:24)/24),c(freqs[2,1], freqs[2,]), type='S', lwd=3, col = 'gray65')
points(c(0,1/48,c(1:24)/24),c(freqs[4,1], freqs[4,]), type='S', lwd=3, col = 'black')
points(c(0,1/48,c(1:24)/24),c(freqs[3,1], freqs[3,]), type='S', lwd=3, col = 'gray40')
legend(x=0.1,y=0.18,legend=c(expression(paste('XB'[1], '= -3')),expression(paste('XB'[1], '= -1')),expression(paste('XB'[1], '= 1')),expression(paste('XB'[1], '= 3'))), lwd =3, col=c('gray90','gray65','gray40','black'))
dev.off()



