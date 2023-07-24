library(nleqslv)
library(nloptr)
library(MASS)
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}

set.seed(100)
mean_z_1=0
mean_z_2=0
mean_w_p=0
#corr1=0.5
# W and Z1
#corr2=0.5
# Wp and Z2
#corr3=0.5
corr=0.5
# Z2 and Z1
#sigma_1=1
#sigma_2=1
sigma_3=1
sigma_1=1
sigma_2=1
#var_z_w_p=matrix(c(sigma_1^2,corr1*sigma_1*sigma_2,corr2*sigma_1*sigma_3,corr1*sigma_1*sigma_2,sigma_2^2,
#corr3*sigma_2*sigma_3,
#corr2*sigma_1*sigma_3,corr3*sigma_2*sigma_3,sigma_3^2),nrow=3,ncol=3)
var_z_w_p=matrix(c(sigma_1^2,corr*sigma_1*sigma_2,corr*sigma_1*sigma_2,sigma_2^2),nrow=2,ncol=2)
theta=c(-2,0.5,0.5)
N=5e4
dw=1
dwz1=1
dwz2=1
gamma_int = c(-0.8,1,0.3,0.7)#0.7




discretize<-function(data){
  quant15z2=quantile(data$Z2,probs = 0.15)
  quant85z2=quantile(data$Z2,probs = 0.85)
  data$Z2d=rep(0,times=(nrow(data)))
  data$Z2d[which(data$Z2<(quant15z2))]=0
  data$Z2d[which(data$Z2>=(quant15z2) & data$Z2<=quant85z2)]=1
  data$Z2d[which(data$Z2>quant85z2)]=2
  
  
  quant15wp=quantile(data$W_p,probs = 0.15)
  quant85wp=quantile(data$W_p,probs = 0.85)
  data$Wd=rep(0,times=(nrow(data)))
  data$Wd[which(data$W_p<(quant15wp))]=0
  data$Wd[which(data$W_p>=(quant15wp) & data$W_p<=quant85wp)]=1
  data$Wd[which(data$W_p>quant85wp)]=2
  return(data)
}


joints<-function(data){
  mat0=matrix(0,nrow=3,ncol=3)
  mat1=matrix(0,nrow=3,ncol=3)
  data=discretize(data)
  for(i in 1:3){
    for(j in 1:3){
      mat0[i,j]=sum((data$Z2d==(i-1) & data$Wd==(j-1) & data$D==0))
    }
  }
  
  for(i in 1:3){
    for(j in 1:3){
      mat1[i,j]=sum((data$Z2d==(i-1) & data$Wd==(j-1) & data$D==1))
    }
  }
  mat0=mat0/nrow(data)
  mat1=mat1/nrow(data)
  return(list(mat0,mat1))
}


totals<-function(data){
  total_D=sum(data$D)
  total_W=sum(data$W_p)
  total_Z2=sum(data$Z2)
  return(c(nrow(data),total_D,total_W,total_Z2))
}

simu_popu<-function(N,mean_w_p,mean_z_1,mean_z_2,var_z_w_p,theta,dw,dwz1,dwz2,sigma_3){
  cov<- mvrnorm(n = N, mu = c(mean_z_1,mean_z_2), Sigma = var_z_w_p)
  data <- data.frame(Z1 = cov[, 1], Z2 = cov[, 2])
  # Generate random uniforms
  #set.seed(5678)
  U1 <- runif(N)
  #set.seed(4321)
  # Generate Disease Status
  DISEASE <- expit(theta[1] + theta[2] * data$Z1 + theta[3]*data$Z2)
  data$D   <- ifelse(DISEASE > U1, 1, 0)
  # Relate W_p and D
  data$W_p <- rnorm(n=N,mean=mean_w_p,sd=sigma_3)+dw* data$D +dwz1*data$Z1 + dwz2*data$Z2 
  data$id=c(1:N)
  # Generate Sampling Status
  return(data)
}






simu_int<-function(data,gamma_int){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(gamma_int[1] +  gamma_int[2]* data$D + gamma_int[3] * data$W_p +  gamma_int[4] * data$Z2)
  S_i  <- ifelse(SELECT > U2i, T, F)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}

calibration<-function(intdata,total){
  prop<-function(gamma){
    z=c(0,0,0,0)
    for(i in 1:nrow(intdata)){
      vec=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
      z = z + 1/(as.vector(expit(gamma %*% vec))) * vec 
    }
    y= z - total
    y
  }
  
  start=c(0,0,0,0)
  z = nleqslv(x=start, fn=prop,method="Newton",global = "dbldog",control=list(trace=1,allowSingular=TRUE))
  estweights = 1/expit(cbind(1,intdata$D,intdata$W_p,intdata$Z2) %*% z$x)
  return(estweights)
}





weighted<-function(intdata,estweights){
  modelinit=glm(D ~ Z1 +Z2, family = stats::quasibinomial(),data=intdata)
  start <- coef(modelinit)
  design <- survey::svydesign(data = intdata,ids = 1:length(intdata$D), strata = NULL,
                              weights = estweights)
  mod=survey::svyglm(data=data_i,stats::formula(D ~ Z1 + Z2), design = design,
                     family = quasibinomial(),
                     start = as.numeric(start))
  final<-coef(mod)
  mod<-summary(mod)
  sd<-mod$coefficients[,2]
  return(list(final,sd))
}


unweighted<-function(intdata){
  modelinit=glm(D ~ Z1 + Z2, family = stats::quasibinomial(),data=intdata)
  start <- coef(modelinit)
  return(start)
}

rep=1000

mat_init=matrix(0,nrow=rep,ncol=3)
mat_ps=matrix(0,nrow=rep,ncol=3)
mat_cal=matrix(0,nrow=rep,ncol=3)



weights_ps<-function(data,intdata){
  data=discretize(data)
  intdata=discretize(intdata)
  intdata$ps_weights=rep(0,times=nrow(intdata))
  prob_int=joints(intdata)
  prob_tot=joints(data)
  prob_int[[1]][prob_int[[1]]==0]=min(prob_int[[1]][prob_int[[1]]!=0])
  prob_int[[2]][prob_int[[2]]==0]=min(prob_int[[2]][prob_int[[2]]!=0])
  prob_tot[[1]][prob_tot[[1]]==0]=min(prob_tot[[1]][prob_tot[[1]]!=0])
  prob_tot[[2]][prob_tot[[2]]==0]=min(prob_tot[[2]][prob_tot[[2]]!=0])
  for(i in 1:3){
    for(j in 1:3){
      intdata[(intdata$Z2d == (i-1) & intdata$Wd==(j-1) & intdata$D==0),]$ps_weights=prob_tot[[1]][i,j]/prob_int[[1]][i,j]
    }
  }
  for(i in 1:3){
    for(j in 2:3){
      intdata[(intdata$Z2d == (i-1) & intdata$Wd==(j-1) & intdata$D==1),]$ps_weights=prob_tot[[2]][i,j]/prob_int[[2]][i,j]
    }
  }
  
  return(intdata$ps_weights)
}


variance_cali<-function(intdata,estweights,theta_hat){
  N_est=N
  ## g_theta
  z=c(1,intdata$Z1[1],intdata$Z2[1])
  g_theta = matrix(0,nrow=(length(z)),ncol=(length(z)))
  intprob = 1/estweights
  for(i in 1:nrow(intdata)){
    z=c(1,intdata$Z1[i],intdata$Z2[i])
    g_theta = g_theta + as.numeric(1/intprob[i] * exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z))^2) * 
      z%*%t(z)
  }
  g_theta=g_theta/N_est
  
  ## g_alpha
  z=c(1,intdata$Z1[1],intdata$Z2[1])
  x=c(1, intdata$D[1], intdata$W_p[1], intdata$Z2[1])
  g_alpha = matrix(0,nrow=(length(z)),ncol=(length(x)))
  
  for(i in 1:nrow(intdata)){
    z=c(1,intdata$Z1[i],intdata$Z2[i])
    x=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
    g_alpha = g_alpha + as.numeric(1/intprob[i] *(1-intprob[i]) *(intdata$D[i] -exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z))))*
      z%*%t(x)
  }
  g_alpha = -g_alpha/N_est
  
  ## H
  x=c(1, intdata$D[1], intdata$W_p[1], intdata$Z2[1])
  H = matrix(0,nrow=(length(x)),ncol=(length(x)))
  
  
  for(i in 1:nrow(intdata)){
    x=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
    H = H + as.numeric(1/intprob[i]*(1-intprob[i]))* x%*%t(x)
  }
  H = -H/N_est
  
  ## E1
  z=c(1,intdata$Z1[1],intdata$Z2[1])
  E1 = matrix(0,nrow=(length(z)),ncol=(length(z)))
  
  for(i in 1:nrow(intdata)){
    z=c(1,intdata$Z1[i],intdata$Z2[i])
    g = as.numeric(1/intprob[i] * (intdata$D[i] - exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * z
    E1 = E1 + g %*% t(g)
  }
  E1 = E1/N_est
  
  ## E2
  x=c(1, intdata$D[1], intdata$W_p[1], intdata$Z2[1])
  z=c(1,intdata$Z1[1],intdata$Z2[1])
  E2 = matrix(0,nrow=(length(x)),ncol=(length(z)))
  for(i in 1:nrow(intdata)){
    z=c(1,intdata$Z1[i],intdata$Z2[i])
    x=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
    E2= E2 + as.numeric(1/intprob[i]^2 *(1-intprob[i])*
                          (intdata$D[i] - exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * x %*% t(z)
  }
  E2= (g_alpha) %*%  (solve(H)) %*% (E2)
  E2= E2/N_est
  
  ## E3
  E3 = t(E2)
  
  # E4
  x=c(1, intdata$D[1], intdata$W_p[1], intdata$Z2[1])
  z=c(1,intdata$Z1[1],intdata$Z2[1])
  E4 = matrix(0,nrow=(length(x)),ncol=(length(x)))
  for(i in 1:nrow(intdata)){
    x=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
    E4 = E4 + 1/intprob[i]^2 *(1-intprob[i])*x %*% t(x)
  }
  
  
  
  E4 =  (g_alpha) %*%  (solve(H)) %*% (E4) %*% t((solve(H))) %*% t(g_alpha) 
  E4= E4/N_est
  
  E= E1-E2-E3+E4
  
  V = solve(g_theta) %*% E %*% t(solve(g_theta))
  V = V/N_est
  return(c(V[1,1],V[2,2],V[3,3]))
}


mat_var_ps_est=matrix(0,nrow=rep,ncol=3)
mat_var_cal_est=matrix(0,nrow=rep,ncol=3)

for(i in 1:rep){
   try({
    data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,var_z_w_p,theta,dw,dwz1,dwz2,sigma_3)
    intdata=simu_int(data,gamma_int)
    total<-totals(data)
    estweights_ps=weights_ps(data,intdata)}
   )
    estweights_cal=calibration(intdata,total)
    mat_init[i,]= unweighted(intdata)
  intdata1=intdata[which(estweights_ps!=0),]
  estweights_ps=estweights_ps[which(estweights_ps!=0)]
  theta_hat_ps=mat_ps[i,]= weighted(intdata1,estweights_ps)[[1]]
  theta_hat_cal=mat_cal[i,]= weighted(intdata,estweights_cal)[[1]]
  mat_var_cal_est[i,]=sqrt(variance_cali(intdata,estweights_cal,theta_hat_cal))
  mat_var_ps_est[i,]=weighted(intdata1,estweights_ps)[[2]]
  }

ps_l=mat_ps-1.96*mat_var_ps_est
ps_u=mat_ps+1.96*mat_var_ps_est

cal_l=mat_cal-1.96*mat_var_cal_est
cal_u=mat_cal+1.96*mat_var_cal_est


list=list(ps_l,ps_u,cal_l,cal_u,mat_ps,mat_cal)
saveRDS(list,"summ4.rds")
