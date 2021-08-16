
unlink(".RData")
options(scipen=999)

library(MASS)

.libPaths("/ibrix/home8/yongli/RCODE")
library(ncvreg)


s0<-3;
rou<-(0.9);  

tc<-rnorm(s0)*0.5; 
#tc<-rep(1,s0)

n <-100;  p <- 2000;   np <- 10000;


b <- rep(0,len=p)
b[sample(1:p,s0,replace=FALSE)] <- tc
#b[1:s0] <- tc



B1_seq<-c(n);   B2_seq<-c(10);  nla<-n;

penseq1<-c("lasso","MCP","SCAD")
penseq2<-c("MCP")
penlist<-list(penseq1,penseq2)
penlist[[1]]
penlist[[2]]


Cov.Matrix<-matrix(0,p,p);  columnmatrix<-col(Cov.Matrix);  rowmatrix<-row(Cov.Matrix)
Cov.Matrix<-(rou)^abs(columnmatrix-rowmatrix); 
X<-mvrnorm(n, rep(0,p), Cov.Matrix); 
Xp<-mvrnorm(np, rep(0,p), Cov.Matrix);





      
z <- X%*%b;    zp <- Xp%*%b;

devi<- 0.2
devi

y <- z+rnorm(n,0,devi);    yp <- zp+rnorm(np,0,devi);


snratio<-sqrt(sum(z^2))/(devi*sqrt(n));             
true_mse<-mean((lm(y~X[,b!=0])$fitted-z)^2);
la_seq<-ncvreg(X, y, family=c("gaussian"),  penalty="lasso", nlambda=nla )$lambda
#max(la_seq)
#min(la_seq)
#la_max<-max(abs(t(scale(X))%*%y)/(n)) *sqrt(n/(n-1));
#la_max
#la_min<-max(abs(t(scale(X))%*%y)/(n)) *sqrt(n/(n-1))*ifelse(n>p,.001,.05);
#la_min
x0<-c(n,p,s0,rou, devi,snratio, true_mse);  
x0





###########################Rigollet Tsabakov####################################
###########################Rigollet Tsabakov####################################                  
###########################Rigollet Tsabakov####################################
###########################Rigollet Tsabakov####################################


T0<-50000;
T1<-7000;

R<-min(n,p)
H<-2*sum(choose(p,0:R)* ((0:R)/(2*p*exp(1)))^(0:R))
H
beta_matrix<-matrix(0,p+1,T0+T1)
be<-4*devi^2

p1<-rep(0,p)

for (tt in 1:(T0+T1))

                    {
                    q1<-p1
                    ind<-sample(1:p,1)
                    q1[ind]<-1-p1[ind]
                    
                    if (sum(q1)<R) {pi_q<-1/H*(sum(q1)/(2*p*exp(1)))^(sum(q1)) } else{ 
                    if (sum(q1)==p) {pi_q<-1/2 } else {pi_q<-0} }
                    pi_q

                    if (sum(p1)<R) {pi_p<-1/H*(sum(p1)/(2*p*exp(1)))^(sum(p1)) } else{ 
                    if (sum(p1)==p) {pi_p<-1/2 } else {pi_p<-0} }
                    pi_p
                    
                                        
                    pq1<-rbind(q1,p1)
                    
                    
                    q11<-c(1,q1)
                    M_q<-min(sum(q11!=0),n)                    
                    p11<-c(1,p1)
                    M_p<-min(sum(p11!=0),n)
                   
                    XX<-cbind(1,X)
                    model_q<-lm(y~XX[,q11!=0]-1)
                    rss_q<-deviance(model_q)                  
                    model_p<-lm(y~XX[,p11!=0]-1)
                    rss_p<-deviance(model_p)
                    
                    
                    risk_q<-rss_q/n +2*devi^2*M_q/n -devi^2
                    risk_q
                    
                    risk_p<-rss_p/n +2*devi^2*M_p/n -devi^2
                    risk_p

                    v_q<- log(pi_q) + (-n *risk_q/be) 
                    v_q
                    v_p<- log(pi_p) + (-n *risk_p/be)
                    v_p
                    r_pq<-min(exp(v_q-v_p),1)
                    r_pq
                    #r_pq<-0.5
                    #print(r_pq)
                    
                    
                    pq_sel<-sample(c(1,2),size=1,prob=c(r_pq,1-r_pq))
                    pq_sel
                    
                    p1<-pq1[pq_sel,]
                    Xy<-data.frame(X,y)
                    zz<-(names(Xy)[1:p])[p1!=0]
                    ss<-as.formula(paste("y~",paste(c("1",zz),collapse="+"),sep="") )
                    model_p<-lm(ss,Xy)
                    summary(model_p)
                    beta_matrix[c(1,p1)!=0,tt]<-(model_p)$coef
                    
                    
                    XX<-NULL; model_p<-NULL;
                    }
                    
                    
                    
#beta_matrix[c(1,b)!=0,]

rt_beta_hat<-apply(beta_matrix[,(T0+1):(T0+T1)],1,mean,na.rm=TRUE)


pred_rt<-cbind(1, Xp)%*% rt_beta_hat
est_rt<-cbind(1,X)%*% rt_beta_hat


mse_rt<-mean((est_rt-z)^2)
mse_rt
pmse_rt<-mean((pred_rt-zp)^2)
pmse_rt
pe_rt<-mean((pred_rt-yp)^2)
pe_rt

x01<-c("rt",mse_rt,pmse_rt,pe_rt)
x01





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




MUHAT <- matrix(0,n,M)
MU <- matrix(0,np,M)
for(k in 1:M){
if(k<MaxM){
USE <- Index[(q*(k-1)+1):(q*k)]
Zk <- cbind(1,X[,USE])
Zkp <- cbind(1,Xp[,USE])

##########################
Hk <- Zkp%*%solve(t(Zk)%*%Zk)%*%t(Zk)
Ek <- Zk%*%solve(t(Zk)%*%Zk)%*%t(Zk)

#########################
MU[,k] <- Hk%*%y
MUHAT[,k] <- Ek%*%y


}
if(k==MaxM){
USE <- Index[(q*(k-1)+1):Maxp]
Zk <- cbind(1,X[,USE])
Zkp <- cbind(1,Xp[,USE])

##########################
Hk <- Zkp%*%solve(t(Zk)%*%Zk)%*%t(Zk)
Ek <- Zk%*%solve(t(Zk)%*%Zk)%*%t(Zk)

#########################
MU[,k] <- Hk%*%y
MUHAT[,k] <- Ek%*%y
}
}



#The final estimator

pred_al <- MU%*%w
est_al <- MUHAT%*%w


mse_al<-mean((est_al-z)^2)
mse_al
pmse_al<-  mean ((pred_al-zp)^2)
pmse_al
pe_al<-  mean ((pred_al-yp)^2)
pe_al
x1<-c("kerchauli",mse_al, pmse_al, pe_al)
x1









########################Selection#########################
########################Selection#########################
########################Selection#########################
########################Selection#########################

for (youok in 1:length(penlist))

{

penseq<-penlist[[youok]]

NLA<- nla * length(penseq);

candpen<-paste0((penseq),collapse="")


for (B2 in B2_seq)
{
#B2<-10; 
ama<-split(sample(1:n, n,replace=F),1:B2)

tm<-matrix(0,B2,NLA); 
                
for (kb2 in 1:B2)
                {

                cvindvl<-ama[[kb2]]; cvindtr<-setdiff(1:n,cvindvl)
                
                Xtr<-X[cvindtr,]; ytr<-y[cvindtr]; 
                Xvl<-X[cvindvl,]; yvl<-y[cvindvl];
                
                for (penk in 1:length(penseq))

                                              {
                                              pen<-penseq[penk];
                                              pen.beta.hat3<-ncvreg(Xtr, ytr, family=c("gaussian"), penalty=pen, lambda=la_seq )$beta; 
                                              tm1<-apply(((cbind(1,Xvl)%*% pen.beta.hat3-yvl)^2),2,mean)
                                              length(tm1)
                                              tm[kb2,(nla*penk-nla+1):(nla*penk)]<-tm1
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


la_o<-la_seq[la_sel]

sel_beta_hat<-ncvreg(X, y, family=c("gaussian"), penalty=pen_sel, lambda=la_o )$beta
pred_sel<-cbind(1, Xp)%*% sel_beta_hat
est_sel<-cbind(1, X)%*% sel_beta_hat

sel_mse<-mean((est_sel-z)^2)
sel_mse
sel_pmse<-mean((pred_sel-zp)^2)
sel_pe<-mean((pred_sel-yp)^2)
sel_pmse
sel_pe


sel_fn_mean<-length(setdiff((1:p)[b!=0],(1:p)[sel_beta_hat[-1]!=0]))
sel_fn_mean
sel_fp_mean<-length(setdiff((1:p)[sel_beta_hat[-1]!=0], (1:p)[b!=0]))
sel_fp_mean


         
x2<-c("SEL", sel_fn_mean,sel_fp_mean,sel_mse, sel_pmse,sel_pe)
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

                bagg_betahat_matrix[,kb1]<-ncvreg(Xb, yb, family=c("gaussian"),  penalty=pen_sel, lambda=la_o )$beta

                
                Xb<-NULL; yb<-NULL;
                }



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
sma_est_mean<-cbind(1,X)%*% beta_hat_mean

sma_mse_mean<-mean((sma_est_mean-z)^2)
sma_mse_mean
sma_pmse_mean<-mean((sma_pred_mean-zp)^2)
sma_pmse_mean
sma_pe_mean<-mean((sma_pred_mean-yp)^2)
sma_pe_mean


############################quantile

sma_fnfp<-NULL

for (thr in c(0.5,0.1,0.05,0.01,0))
{
beta_hat_quan1<-apply(abs(sma_betahat_matrix),1,quantile, probs=c(thr), na.rm=TRUE)
beta_hat_quan1[beta_hat_quan1!=0]
sma_fnfp<-c(sma_fnfp, length(setdiff((1:p)[b!=0],(1:p)[beta_hat_quan1[-1]!=0])))

sma_fnfp<-c(sma_fnfp, length(setdiff((1:p)[beta_hat_quan1[-1]!=0], (1:p)[b!=0])))

}
sma_fnfp


x31<-c("sma_sel",sma_uniqueproc, sma_freq1 , sma_uniquemodel, sma_freq2, sma_fnfp, "sma_pred",sma_mse_mean, sma_pmse_mean, sma_pe_mean)


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
bagg_est_mean<-cbind(1,X)%*% beta_hat_mean


bagg_mse_mean<-mean((bagg_est_mean-z)^2)
bagg_mse_mean
bagg_pmse_mean<-mean((bagg_pred_mean-zp)^2)
bagg_pmse_mean
bagg_pe_mean<-mean((bagg_pred_mean-yp)^2)
bagg_pe_mean


############################quantile

bagg_fnfp<-NULL

for (thr in c(0.5,0.1,0.05,0.01,0))
{
beta_hat_quan1<-apply(abs(bagg_betahat_matrix),1,quantile, probs=c(thr), na.rm=TRUE)
beta_hat_quan1[beta_hat_quan1!=0]
bagg_fnfp<-c(bagg_fnfp, length(setdiff((1:p)[b!=0],(1:p)[beta_hat_quan1[-1]!=0])))

bagg_fnfp<-c(bagg_fnfp, length(setdiff((1:p)[beta_hat_quan1[-1]!=0], (1:p)[b!=0])))

}
bagg_fnfp



x32<-c("bagg_sel",bagg_uniqueproc, bagg_freq1 , bagg_uniquemodel, bagg_freq2, bagg_fnfp, "bagg_pred",bagg_mse_mean, bagg_pmse_mean, bagg_pe_mean)


################################################################################
################################################################################
################################################################################
################################################################################

mre<-t(c(B1,B2,candpen,x0,x01, x1,x2,x31,x32))
print(mre)

vvv<-round(abs(rnorm(1)*100000000))
filename<-paste(vvv,".",s0,"txt",sep="")
write.table(mre, file=filename, row.names=FALSE,col.names=FALSE)
filename


}
}
}