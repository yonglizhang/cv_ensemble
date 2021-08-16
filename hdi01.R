
options(scipen=999)
library(MASS)



.libPaths("/ibrix/home8/yongli/RCODE")
library(ncvreg)

.libPaths("/ibrix/home8/yongli/RCODE")
library(glmnet)

.libPaths("/ibrix/home8/yongli/RCODE")
library(hdi)




np<-1; 
nla<-100;
B2_seq<-c(2, 5, 10);


penseq1<-c("lasso","MCP","SCAD")
penseq2<-c("MCP")
penlist<-list(penseq1,penseq2)
penlist[[1]]
penlist[[2]]





data(riboflavin)





XX <- as.matrix(riboflavin[,-1])
round(XX[1:20,1:8],3)
genelistname<-colnames(XX)

XX<-as.matrix(riboflavin[,-1],dimnames = NULL)
XX<-matrix(XX,nrow(XX), ncol(XX),dimnames = NULL)


YY <- riboflavin[,1]

N<-length(YY);


p<-ncol(XX); #genes<-colnames(data.frame(cbind(XX,YY),row.names = NULL))[1:p];

mre_matrix<-NULL

val_ind<-sample(1:N,np,replace=FALSE)
val_ind

#val_ind<-1
fit_ind<-setdiff((1:N),val_ind);
X<-XX[fit_ind,]; dim(X); 
y<-YY[fit_ind]; length(y);
yp<-YY[val_ind]; length(yp);
dat<-data.frame(cbind(X,y),row.names=NULL);dim(dat);
if (np==1) {Xp<-t(as.matrix(XX[val_ind,]))} else {Xp<-(as.matrix(XX[val_ind,]))};  dim(Xp)
dim(X)
dim(Xp)


n<-nrow(X); 
n

la_seq<-ncvreg(X, y, family=c("gaussian"), penalty="lasso",  nlambda=nla )$lambda

B1_seq<-c(n);
###########################Ker Chau Li##########################################
###########################Ker Chau Li##########################################                  
###########################Ker Chau Li##########################################
###########################Ker Chau Li##########################################

###############################################################################z#
# This file provides R code for implementing MCV2 in the following paper
# (Search space of $M$ and $p_k$ is more flexible than that in the paper)
#
# A model averaging approach for high-dimensional regression
#
# Author : Ando, T. and Li, K.C.
# 
# Journal: Journal of the American Statistical Association, 2014
#
################################################################################
#Data generation

#n <- 100
#p <- 1000
#s <- 20
#X <- matrix(runif(n*p,-1,1),nrow=n,ncol=p)
#b <- rep(0,len=p)
#b[1:s] <- runif(s,-1,1)
#z <- X%*%b
#y <- z+rnorm(n,0,0.1)


#Correlation test

COR <- rep(0,len=p)
for(i in 1:p){COR[i] <- cor.test(X[,i],y,"two.sided")$p.value}

Index <- cbind(1:p,COR)
Index <- Index[order(Index[,2],decreasing=F),1:2]


#maximum number of regressors to be used

Maxp <- sum(Index[,2]<=0.05) 
Maxp
Maxp <- n
Maxp
#Candidates of the number of regressors
#To be specified by users

Q <- seq(0.05,0.40,len=8)*n


#Model search

minCV.score <- 10^10

for(Qind in 1:length(Q)){

q <- Q[Qind]

if(trunc(Maxp/q)-abs(Maxp/q)==0){MaxM <- trunc(Maxp/q)}
if(trunc(Maxp/q)-abs(Maxp/q)!=0){MaxM <- trunc(Maxp/q)+1}

for(M in 2:MaxM){

MU <- matrix(0,n,M)

for(k in 1:M){

if(k<MaxM){
USE <- Index[(q*(k-1)+1):(q*k)]
Zk <- X[,USE]
Hk <- Zk%*%solve(t(Zk)%*%Zk)%*%t(Zk)
Dk <- diag(1/as.vector(diag(diag(1,n)-Hk)))
Hk <- Dk%*%(Hk-diag(1,n))+diag(1,n)
MU[,k] <- Hk%*%y
}
if(k==MaxM){
USE <- Index[(q*(k-1)+1):Maxp]
Zk <- X[,USE]
Hk <- Zk%*%solve(t(Zk)%*%Zk)%*%t(Zk)
Dk <- diag(1/as.vector(diag(diag(1,n)-Hk)))
Hk <- Dk%*%(Hk-diag(1,n))+diag(1,n)
MU[,k] <- Hk%*%y
}

CV <- function(w){
sum( (MU%*%w-y)^2 )
}
w <- rep(0.1,len=M)
fit <- optim(w,fn=CV,method="L-BFGS-B",lower=rep(0,len=M),
upper=rep(1,len=M))
w <- fit$par
CV.score <- fit$value

if(CV.score<=minCV.score){
minCV.score <- CV.score
Opt.M <- M
Opt.w <- w
Opt.q <- q
}

}}}


#Construction of the final estimator

M <- Opt.M
q <- Opt.q
w <- Opt.w

if(trunc(Maxp/q)-abs(Maxp/q)==0){MaxM <- trunc(Maxp/q)}
if(trunc(Maxp/q)-abs(Maxp/q)!=0){MaxM <- trunc(Maxp/q)+1}

MU <- matrix(0,np,M)
for(k in 1:M){
if(k<MaxM){
USE <- Index[(q*(k-1)+1):(q*k)]
Zk <- X[,USE]
Zkp <- Xp[,USE]
Hk <- Zkp%*%solve(t(Zk)%*%Zk)%*%t(Zk)
MU[,k] <- Hk%*%y
}
if(k==MaxM){
USE <- Index[(q*(k-1)+1):Maxp]
Zk <- X[,USE]
Zkp <- Xp[,USE]
Hk <- Zkp%*%solve(t(Zk)%*%Zk)%*%t(Zk)
MU[,k] <- Hk%*%y
}
}

#The final estimator

pred <- MU%*%w


#The final estimator
pred <- MU%*%w

pe_mv<-  mean ((pred-yp)^2)
pe_mv
x1<-c("kerchauli",pe_mv)
x1




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
for (youok in 1:length(penlist))

{

print(youok)
penseq<-penlist[[youok]]
print(penseq)



NLA<- nla * length( penseq);
print(NLA)
candpen<-paste0((penseq),collapse="")
print(candpen)


for (B2 in B2_seq)
{



ama<-split(sample(1:n, n,replace=F),1:B2)

tm<-matrix(0,B2,NLA); 
                
for (kb2 in 1:B2)
                {
                

                    
                #cvindtr<-sample(1:n,m,replace=F); cvindvl<-setdiff(1:n,cvindtr) 
                cvindvl<-ama[[kb2]]; cvindtr<-setdiff(1:n,cvindvl)
                
                Xtr<-X[cvindtr,]; ytr<-y[cvindtr]; dim(Xtr)
                Xvl<-X[cvindvl,]; yvl<-y[cvindvl];
                for (penk in 1:length(penseq))

                                {
                                pen<-penseq[penk]; 
                                pen.beta.hat3<-ncvreg(Xtr, ytr, family=c("gaussian"), penalty=pen, lambda=la_seq )$beta; 
                                tm1<-apply(((cbind(1,Xvl)%*% pen.beta.hat3-yvl)^2),2,mean)
                                length(tm1)
                                tm[kb2,(nla*penk-nla+1):(nla*penk) ]<-tm1
                                pen.beta.hat3<-NULL
                                }
                
                }
tm_mean<-apply(tm,2,mean)
tm_mean
opt_la<-min((1:NLA)[tm_mean==min(tm_mean)])
opt_la

penk_sel<-ceiling(opt_la/nla)
penk_sel
pen_sel<-penseq[penk_sel]
pen_sel
la_sel<-opt_la-(penk_sel-1)*nla
la_sel

                

la_o<-la_seq[la_sel]
la_o

sel_beta_hat<-ncvreg(X, y, family=c("gaussian"), penalty=pen_sel, lambda=la_o)$beta
pred_sel<-cbind(1, Xp)%*% sel_beta_hat

sel_pe<-mean((pred_sel-yp)^2)

         
x2<-c("SEL", sel_pe)
x2

########################Bootstrapping#####################
########################Bootstrapping#####################
########################Bootstrapping#####################
########################Bootstrapping#####################

for ( B1 in B1_seq)

{


pen.beta.hat0<-matrix(0,p+1,NLA)
for (penk in 1:length(penseq))

                 {
                  pen<-penseq[penk]; 
                  pen.beta.hat0[, (nla*penk-nla+1):(nla*penk )]<-ncvreg(X, y, family=c("gaussian"), penalty=pen, lambda=la_seq )$beta;

                 }
                 
                  
dim(pen.beta.hat0)
can_pen0<-as.matrix(t(pen.beta.hat0[-1,]!=0))
dim(can_pen0)


####################################################

bagg_betahat_matrix<-matrix(0, p+1, B1)

opt_la_seq<-rep(0,B1)


for (kb1 in 1:B1)
                {          
                #kb1<-1
                bind<-sample(1:n,n,replace=T);      
                Xb<-X[bind,]; yb<-y[bind];

               
                tm<-matrix(0,B2,NLA); 

                for (kb2 in 1:B2)
                                {
                                #kb2<-1
                                cvindvl<-ama[[kb2]]; cvindtr<-setdiff(1:n,cvindvl)
                                
                                Xbtr<-Xb[cvindtr,]; ybtr<-yb[cvindtr];
                                Xbvl<-Xb[cvindvl,]; ybvl<-yb[cvindvl];
                                
                                for (penk in 1:length(penseq))

                                                {
                                                pen<-penseq[penk];  
                                                pen.beta.hat33<-ncvreg(Xbtr, ybtr, family=c("gaussian"), penalty=pen, lambda=la_seq )$beta; 
                                                tm1<-apply(((cbind(1,Xbvl)%*% pen.beta.hat33-ybvl)^2),2,mean)
                                                length(tm1)
                                                
                                                tm[kb2,(nla*penk-nla+1):(nla*penk ) ]<-tm1
                                                pen.beta.hat33<-NULL
                                                }
                                
                                }
                                
                                
                                
                                
                tm_mean<-apply(tm,2,mean)
                tm_mean
                opt_la<-min((1:NLA)[tm_mean==min(tm_mean)])
                opt_la
                opt_la_seq[kb1]<-opt_la

                penk_sel<-ceiling(opt_la/nla)
                penk_sel
                pen_sel<-penseq[penk_sel]
                pen_sel
                la_sel<-opt_la-(penk_sel-1)*nla
                la_sel
 
                la_o<-la_seq[la_sel]
                la_o

                bagg_betahat_matrix[,kb1]<-ncvreg(Xb, yb, family=c("gaussian"),  penalty=pen_sel, lambda=la_o)$beta

                
                Xb<-NULL; yb<-NULL;
                }

print(opt_la_seq)

################################################################################
################################################################################
##########################sma##################################################
################################################################################
################################################################################
sma_betahat_matrix<-pen.beta.hat0[,opt_la_seq]
dim(sma_betahat_matrix)


sma_m<-(sma_betahat_matrix!=0)
dim(sma_m)

sma_uniqueproc<-length(unique( opt_la_seq))
sma_uniqueproc
sma_freq1<-max( data.frame(table(opt_la_seq))[,"Freq"])/B1
sma_freq1 
 
sma_m_s <- data.frame(table(apply(sma_m, 2, paste, collapse = ", ")))
sma_uniquemodel<-length(sma_m_s[,"Freq"])
sma_uniquemodel
sma_freq2<-max(sma_m_s[,"Freq"])/B1
sma_freq2
###########################################################

beta_hat_mean<-apply(sma_betahat_matrix,1,mean)
length(beta_hat_mean)


sma_pred_mean<-cbind(1,Xp)%*% beta_hat_mean


sma_pe_mean<-mean((sma_pred_mean-yp)^2)
(sma_pe_mean)




x31<-c("sma_sel",sma_uniqueproc, sma_freq1 , sma_uniquemodel, sma_freq2,  "sma_pred",sma_pe_mean)
x31

################################################################################
################################################################################
##########################bagg##################################################
################################################################################
################################################################################


dim(bagg_betahat_matrix)
bagg_m<-(bagg_betahat_matrix!=0)
dim(bagg_m)

bagg_uniqueproc<-length(unique( opt_la_seq))
bagg_uniqueproc
bagg_freq1<-max( data.frame(table(opt_la_seq))[,"Freq"])/B1
bagg_freq1 
 
bagg_m_s <- data.frame(table(apply(bagg_m, 2, paste, collapse = ", ")))
bagg_uniquemodel<-length(bagg_m_s[,"Freq"])
bagg_uniquemodel
bagg_freq2<-max(bagg_m_s[,"Freq"])/B1
bagg_freq2
###########################################################

beta_hat_mean<-apply(bagg_betahat_matrix,1,mean)
length(beta_hat_mean)

bagg_pred_mean<-cbind(1,Xp)%*% beta_hat_mean


bagg_pe_mean<-mean((bagg_pred_mean-yp)^2)
(bagg_pe_mean)





x32<-c("bagg_sel",bagg_uniqueproc, bagg_freq1 , bagg_uniquemodel, bagg_freq2,  "bagg_pred",bagg_pe_mean)
x32


################################################################################
################################################################################
################################################################################
################################################################################

mre<-t(c(B1,B2,n,np, candpen,x1,x2,x31,x32))
print(mre)

vvv<-round(abs(rnorm(1)*100000000))
filename<-paste(vvv,".",np,"txt",sep="")
write.table(mre, file=filename, row.names=FALSE,col.names=FALSE)
print(filename)


}
}
}