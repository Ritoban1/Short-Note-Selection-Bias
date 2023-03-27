library(nloptr)
library(nleqslv)
library(MASS)
library(nnet)
library(simplexreg)
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
gamma_ext = c(-0.6,1.2,0.4,0.5) #0.5

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

simu_ext<-function(data,gamma_ext){
  U2e <- runif(N)
  # Generate Sampling Status
  SELECT <- 0.75*expit(gamma_ext[1] +  gamma_ext[2]* data$D + gamma_ext[3] * data$W_p +  gamma_ext[4] * data$Z2)
  S_e  <- ifelse(SELECT > U2e, T, F)
  # Observed Data
  data_e <- data[which(S_e==1),]
  data_e$Select_Weights = 0.75*expit(gamma_ext[1] +  gamma_ext[2]* data_e$D + gamma_ext[3] * data_e$W_p +  gamma_ext[4] * data_e$Z2)
  return(data_e)
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

chen<-function(extdata,intdata){
  prop<-function(gamma){
    y=c(0,0,0,0)
    for(i in 1:nrow(extdata)){
      vec=c(1, extdata$D[i], extdata$W_p[i], extdata$Z2[i])
      y = y + as.vector(expit(gamma %*% vec)) * vec *
        1/(extdata$Select_Weights[i]) 
    }
    vec1=c(nrow(intdata),sum(intdata$D),sum(intdata$W_p),sum(intdata$Z2))
    y= vec1 - y
    y
  }
  
  start=c(0,0,0,0)
  z = nleqslv(x=start, fn=prop,method="Newton",global = "dbldog",control=list(trace=1,allowSingular=TRUE))
  estweights = 1/expit(cbind(1,intdata$D,intdata$W_p,intdata$Z2) %*% z$x)
  return(list(estweights,z$x))
}

lauren<-function(extdata,intdata){
  modsimp=simplexreg(data=extdata,formula=Select_Weights ~ (D+W_p +Z2) ,link="logit")
  wtintsimp=predict(modsimp,newdata=intdata,type="response")
  inter=intersect(intdata$id,extdata$id)
  intju=setdiff(intdata$id,inter)
  extju=setdiff(extdata$id,inter)
  justext=extdata[which(extdata$id %in% extju),]
  justint=intdata[which(intdata$id %in% intju),]
  bothie=extdata[which(extdata$id %in% inter),] 
  bothie$co = 1
  justint$co = 2
  justext$co = 3
  combdata=rbind(justint,bothie[,-6],justext[,-6])
  combdata$group=relevel(as.factor(combdata$co),ref=1)
  model1=multinom(group ~  D +  W_p + Z2, data= combdata)
  Pmult=model1$fitted.values
  Pmult1=Pmult[-c(1:nrow(justext)),]
  prob=rep(0,times=nrow(intdata))
  for(i in 1:nrow(intdata)){
    prob[i]= wtintsimp[i] * (Pmult1[i,1] + Pmult1[i,2])/(Pmult1[i,1] + Pmult1[i,3])
  }
  
  estweights=1/prob
  return(estweights)
}




weighted<-function(intdata,estweights){
  modelinit=glm(D ~ Z1 +Z2, family = stats::quasibinomial(),data=intdata)
  start <- coef(modelinit)
  design <- survey::svydesign(data = intdata,ids = 1:length(intdata$D), strata = NULL,
                              weights = estweights)
  mod=survey::svyglm(data=intdata,stats::formula(D ~ Z1 + Z2), design = design,
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

### variance estimation 
variance_chen<-function(intdata,extdata,estweights,theta_hat,gamma_hat){
  common_id=intersect(intdata$id,extdata$id)
  data_common=extdata[which(extdata$id %in% common_id),]
  N_est=sum(1/extdata$Select_Weights)
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
  
  ## M
  x=c(1, extdata$D[1], extdata$W_p[1], extdata$Z2[1])
  M = matrix(0,nrow=(length(x)),ncol=(length(x)))
  
  
  for(i in 1:nrow(extdata)){
    x=c(1, extdata$D[i], extdata$W_p[i], extdata$Z2[i])
    pi_a= exp(gamma_hat%*%x)/(1+ exp(gamma_hat%*%x))
    M = M + as.numeric(pi_a*(1-pi_a) * (1/extdata$Select_Weights[i]))* x%*%t(x)
  }
  M = -M/N_est
  
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
  E2 = matrix(0,nrow=(length(z)),ncol=(length(z)))
  sum1 = matrix(0,nrow=(length(x)),ncol=(length(z)))
  sum2 = matrix(0,nrow=(length(x)),ncol=(length(z)))
  
  for(i in 1:nrow(intdata)){
    x=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
    z=c(1,intdata$Z1[i],intdata$Z2[i])
    g = as.numeric(1/intprob[i] * (intdata$D[i] - exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * z
    sum1 = sum1 + x %*% t(g)
  }
  
  for(i in 1:nrow(data_common)){
    x=c(1, data_common$D[i], data_common$W_p[i], data_common$Z2[i])
    z=c(1,data_common$Z1[i],data_common$Z2[i])
    pi_a = exp(gamma_hat%*%x)/(1+ exp(gamma_hat%*%x))
    pi_e = data_common$Select_Weights[i]
    g = as.numeric(1/pi_a * (data_common$D[i] - exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * z
    sum2 = sum2 + as.numeric(pi_a/pi_e) * x%*%t(g)
  }
  E2= (g_alpha) %*%  (solve(M)) %*% (sum1-sum2)
  E2= E2/N_est
  
  ## E3
  E3 = t(E2)
  
  # E4
  x=c(1, intdata$D[1], intdata$W_p[1], intdata$Z2[1])
  z=c(1,intdata$Z1[1],intdata$Z2[1])
  E4 = matrix(0,nrow=(length(z)),ncol=(length(z)))
  sum1 = matrix(0,nrow=(length(x)),ncol=(length(x)))
  sum2 = matrix(0,nrow=(length(x)),ncol=(length(x)))
  sum3 = matrix(0,nrow=(length(x)),ncol=(length(x)))
  
  for(i in 1:nrow(intdata)){
    x=c(1, intdata$D[i], intdata$W_p[i], intdata$Z2[i])
    sum1 = sum1 + x %*% t(x)
  }
  
  for(i in 1:nrow(data_common)){
    x=c(1, data_common$D[i], data_common$W_p[i], data_common$Z2[i])
    pi_a = exp(gamma_hat%*%x)/(1+ exp(gamma_hat%*%x))
    pi_e = data_common$Select_Weights[i]
    sum2 = sum2 + as.numeric(pi_a/pi_e) * x%*%t(x)
  }
  
  for(i in 1:nrow(extdata)){
    x=c(1, extdata$D[i], extdata$W_p[i], extdata$Z2[i])
    pi_a = exp(gamma_hat%*%x)/(1+ exp(gamma_hat%*%x))
    pi_e = extdata$Select_Weights[i]
    r = as.numeric(pi_a/pi_e)
    sum3=sum3 + r^2 * x%*%t(x)
  }
  
  E4 =  (g_alpha) %*%  (solve(M)) %*% (sum1- 2*sum2 + sum3) %*% t((solve(M))) %*% t(g_alpha) 
  
  E4= E4/N_est
  
  E= E1-E2-E3+E4
  
  V = solve(g_theta) %*% E %*% t(solve(g_theta))
  V = V/N_est
  return(c(V[1,1],V[2,2],V[3,3]))
}





rep=100

mat_init=matrix(0,nrow=rep,ncol=3)
mat_chen=matrix(0,nrow=rep,ncol=3)
mat_lauren=matrix(0,nrow=rep,ncol=3)
mat_var_chen=matrix(0,nrow=rep,ncol=3)
mat_var_lauren=matrix(0,nrow=rep,ncol=3)



for(i in 1:rep){
  data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,var_z_w_p,theta,dw,dwz1,dwz2,sigma_3)
  extdata=simu_ext(data,gamma_ext)
  intdata=simu_int(data,gamma_int)
  chen_i = chen(extdata,intdata)
  estweights_chen=chen_i[[1]]
  estweights_lauren=lauren(extdata,intdata)
  mat_init[i,]= unweighted(intdata)
  theta_hat=mat_chen[i,]= weighted(intdata,estweights_chen)[[1]]
  lauren1=weighted(intdata,estweights_lauren)
  theta_hat=mat_lauren[i,]=lauren1[[1]]
  gamma_hat = chen_i[[2]]
  mat_var_chen[i,]= sqrt(variance_chen(intdata,extdata,estweights_chen,theta_hat,gamma_hat))
  mat_var_lauren[i,]= lauren1[[2]]
}

matlist=list(mat_init,mat_chen,mat_lauren)
saveRDS(matlist,"matcs4indi.rds")

