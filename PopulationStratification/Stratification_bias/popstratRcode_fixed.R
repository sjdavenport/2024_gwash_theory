m=200
n=100

h2=0.3

sigm_s=0.3



rho=0.3


Sigm = rho^abs(outer(1:(m), 1:(m), "-"))

out.Sigm=eigen(Sigm)



out.Sigm=eigen(Sigm)
#v1=out.Sigm$vectors[,1]
#v1

Sigm_half=out.Sigm$vectors %*% diag(out.Sigm$values^{1/2}) %*% t(out.Sigm$vectors)




n_sim=10^3
res_reg_free=rep(NA,n_sim)
res_reg_free_inter=rep(NA,n_sim)

res_reg_fixed=rep(NA,n_sim)

res_reg_GWASH=rep(NA,n_sim)

do.scale=FALSE

F_ST_=c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5)

#F_ST_=0
n_F_ST=length(F_ST_)

res_m_sim=rep(NA,n_F_ST)
res_s_sim=rep(NA,n_F_ST)

res_m_fixed_sim=rep(NA,n_F_ST)
res_s_fixed_sim=rep(NA,n_F_ST)

res_m_GWASH_sim=rep(NA,n_F_ST)
res_s_GWASH_sim=rep(NA,n_F_ST)


res_m_intercept=rep(NA,n_F_ST)
res_s_intercept=rep(NA,n_F_ST)


res_theory=rep(NA,n_F_ST)
res_Cf=rep(NA,n_F_ST)


for (v in 1:n_F_ST)
{
  print(v)  
  F_ST=F_ST_[v]
  #temp=rnorm(m)
  #temp=temp/sqrt(sum(temp^2))
  f_vec=rnorm(m)*sqrt(F_ST)
  
  
  for (r in 1:n_sim)
  {
    
    
    X=array(dim=c(n,m))
    
    X=array(dim=c(n,m))
    Z=array(rnorm(n*m),dim=c(n,m))
    for (i in 1:n)
      if (i<(n/2 +1)) X[i,]=(Sigm_half %*% Z[i,]-f_vec)/(sqrt(1+f_vec^2)) else X[i,]=(Sigm_half %*% Z[i,]+f_vec)/(sqrt(1+f_vec^2))
    
    bet=c(rnorm(m,mean=0,sd=sqrt((h2)/m)))
    
    xi=c(rep(-sigm_s,n/2),rep(sigm_s,n/2))
    
    Y=X %*% bet+xi+rnorm(n,mean=0,sd=sqrt(1-h2-sigm_s^2))
    
    V_y=var(Y)
    
    if (do.scale)
    {
      
      X=scale(X,center=FALSE)
      #Y=scale(Y)
      Y=Y/as.numeric(sqrt(V_y))
      
      u2=(t(X) %*% Y)^2/(n-1)
      R=t(X) %*% X/(n-1)
      #u2=(t(X) %*% Y)^2/(n)
      #R=t(X) %*% X/(n)
    } else {
      u2=(t(X) %*% Y)^2/n
      R=t(X) %*% X/n
      
    }
    
    d2=diag(R)
    hatell=apply(R^2,2,sum)
    
    
    res_reg_free[r]=mean(u2*(n/m*hatell-mean(n/m*hatell))) /mean((n/m*hatell-mean(n/m*hatell))^2)
    
    res_reg_fixed[r]=mean((u2-1)*(n/m*hatell-1)) /mean((n/m*hatell-1)^2)
    
    res_reg_GWASH[r]=mean(u2-1)/mean(n/m*hatell-1)
    
    res_reg_free_inter[r]=mean(u2)-res_reg_free[r]*(n/m*mean(hatell)-1)
    #res_reg_free_inter[r]=(n/m*mean(hatell/m)-1/m)
    
    
    #res_den_free[r]=mean((hatell-mean(hatell))^2)
    
  }
  
  res_m_sim[v]=mean(res_reg_free)-h2
  res_s_sim[v]=sd(res_reg_free)/sqrt(n_sim)
  
  res_m_fixed_sim[v]=mean(res_reg_fixed)-h2
  res_s_fixed_sim[v]=sd(res_reg_fixed)/sqrt(n_sim)
  
  res_m_GWASH_sim[v]=mean(res_reg_GWASH)-h2
  res_s_GWASH_sim[v]=sd(res_reg_GWASH)/sqrt(n_sim)
  
  
  res_m_intercept[v]=mean(res_reg_free_inter)
  res_s_intercept[v]=sd(res_reg_free_inter)/sqrt(n_sim)
  
  
  C_f=mean(f_vec^2/(1+f_vec^2))
  res_theory[v]=sigm_s^2/C_f
  res_Cf[v]=C_f
}




res_m_sim
res_s_sim

res_theory


if (do.scale==TRUE) plot(F_ST_,res_m_sim,xlab='F_ST',ylab='bias',lwd=2,ylim=c(min(res_m_sim-2*res_s_sim,res_m_fixed_sim),max(res_m_sim+2*res_s_sim,res_m_fixed_sim)),main=
                           paste('m=',m,',n=',n,',Stand.',sep='')) else plot(F_ST_,res_m_sim,xlab='F_ST',ylab='bias',lwd=2,ylim=c(min(res_m_sim-2*res_s_sim,res_m_fixed_sim,res_m_GWASH_sim),max(res_m_sim+2*res_s_sim,res_m_fixed_sim)),main=paste('m=',m,',n=',n,',Not Stand.',sep=''))

points(F_ST_+0.01,res_m_fixed_sim,lwd=2,col='red')

points(F_ST_-0.01,res_m_GWASH_sim,lwd=2,col='blue')


for (j in 1:n_sim)
{  
  lines(c(F_ST_[j],F_ST_[j]),c(res_m_sim[j]-2*res_s_sim[j],res_m_sim[j]+2*res_s_sim[j]),lty='dashed')
  lines(c(F_ST_[j]+0.01,F_ST_[j]+0.01),c(res_m_fixed_sim[j]-2*res_s_fixed_sim[j],res_m_fixed_sim[j]+2*res_s_fixed_sim[j]),col='red',lty='dashed')
  lines(c(F_ST_[j]-0.01,F_ST_[j]-0.01),c(res_m_GWASH_sim[j]-2*res_s_GWASH_sim[j],res_m_GWASH_sim[j]+2*res_s_GWASH_sim[j]),col='blue',lty='dashed')
}  
lines(F_ST_,res_theory)
