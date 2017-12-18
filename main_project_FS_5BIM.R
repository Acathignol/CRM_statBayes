rm(list=ls())

## PROJET STATS BAYESIENNES

library(rjags)
setwd(dir="C:/Users/Annie/Desktop/Bayes/")


 modelproject <-
   "model
   {
     # Definition des liens
     for (i in 1:n)
     {
        ptox[i] = exp(a0 + theta*dose[i])/(1+exp(a0 + theta*dose[i]))
        ntox[i] ~ dbinom(ptox[i],nb_ind)
     }
     # Definition des distributions a priori
    theta ~ dexp(1)
        
}"


################################################

### 1
 pprior_dose=c(0.05,0.15,0.333,0.333,0.333,0.333)
 a0 = 1
 
 data4jags <- list(
 ntox = c(0,1,1,0,1,2) ,
 a0 = 1,
 dose = log(pprior_dose/(1 - pprior_dose)) - a0,
 nb_ind = 3,
 n = 6
)
 
 
 inits1 <- 
   list(
     list(theta = 0.5),
     list(theta = 1.5),
     list(theta = 2)
   )
 
 
 
 # Simulation des chaines de Markov
 model <- jags.model(file =textConnection(modelproject),data = data4jags, n.chains=3, inits = inits1)
 
 # Phase de chauffe
 update(model,3000)
 
 # Simulations monitorees
 # Tirage de 2000 valeurs de theta
 mcmc <- coda.samples(model,c("theta"),n.iter = 2000)
 
 
 # Estimation des parametres a posteriori
 mean(mcmc[[1]][,"theta"])  #theta
 summary(mcmc)
 
 plot(mcmc,trace = T, density = F)
 gelman.diag(mcmc)
 autocorr.plot(mcmc[[1]])

 
##### Questions 1.a

dose_to_give <- c()
pprior_dose=c(0.05,0.15,0.333,0.333,0.333,0.333)
a0 = 1
ntox = c(0,1,1,0,1,2)
dose_tot = c(0.5,1,3,5,6)
pprior_dose_tot = c(0.05,0.1,0.15,0.333,0.5)
dose_star_tot = log(pprior_dose_tot/(1 - pprior_dose_tot)) - a0

nb_cohorte = 6

for(i in 1:(nb_cohorte-1)){
  
  data4jags <- list(
    ntox = ntox[1:i] ,
    a0 = 1,
    dose = log(pprior_dose[1:i]/(1 - pprior_dose[1:i])) - a0,
    nb_ind = 3,
    n = i
  )
  
  
  inits1 <- 
    list(
      list(theta = 0.5),
      list(theta = 1.5),
      list(theta = 2)
    )
  
  
  
  # Simulation des chaines de Markov
  model <- jags.model(file =textConnection(modelproject),data = data4jags, n.chains=3, inits = inits1)
  
  # Phase de chauffe
  update(model,3000)
  
  # Simulations monitorees
  # Tirage de 2000 valeurs de theta
  mcmc <- coda.samples(model,c("theta"),n.iter = 10000, thin = 3)
 
  # Estimation des parametres a posteriori
  # summary(mcmc)
  plot(mcmc,trace = T, density = F, main = as.character(i))
  # gelman.diag(mcmc)
  autocorr.plot(mcmc[[1]])
  
  mcmctot = as.data.frame(as.matrix(mcmc))
  
  ptox_tot = lapply(X = dose_star_tot, FUN = function(X) {median(exp(a0 + mcmctot$theta*X)/(1 +  exp(a0 + mcmctot$theta*X)))})
  dose_to_give[i] <- dose_tot[sum(ptox_tot<=0.33)] 
  
}

dose_to_give

##### Question 1.b

d0 <- list( a0 = 1,
            dose = log(pprior_dose/(1 - pprior_dose)) - a0,
            nb_ind = 3,
            n = 6)

model0 <- jags.model(file = textConnection(modelproject), data= d0,
                     inits = inits1, n.chains = 3)


update(model0, 3000)
mcmc0 <- coda.samples(model0, c("theta"), n.iter = 2000 )
mcmctot0 <- as.data.frame(as.matrix(mcmc0))
plot(mcmc0,trace=F)

ptox0 <- lapply(X = dose_star_tot, FUN = function(X) {exp(a0 + mcmctot0$theta*X)/(1 +  exp(a0 + mcmctot0$theta*X))})
ptox0 <- as.data.frame(ptox0)

for(p in 1:5){
  hist(ptox0[,p],main=paste("Dose :",as.character(dose_tot[p])),breaks=30)
  abline(v=pprior_dose_tot[p],col="red")
}


##### Question 1.c

modelproject2 <-
  "model
{
  # Definition des liens
  for (i in 1:n)
  {
    ptox[i] = exp(a0 + theta*dose[i])/(1+exp(a0 + theta*dose[i]))
    ntox[i] ~ dbinom(ptox[i],nb_ind)
  }
  # Definition des distributions a priori
  p0 ~ dunif(0,pm)
  gamma ~ dunif(dmin,dmax)
  theta <- log((p0*(1-pm))/(pm*(1-p0)))/(dmin-gamma)
  a0 <- log(pm/(1-pm))-log((p0*(1-pm))/(pm*(1-p0)))*gamma/(dmin-gamma)
}"

  dose_to_give <- c()
  dose_init=c(0.5,3,5,5,5,5)
  ntox = c(0,1,1,0,1,2)
  dose_tot = c(0.5,1,3,5,6)
  nb_cohorte = 6
  
  
  for(i in 1:(nb_cohorte-5)){
    
    data4jags <- list(
      ntox = ntox[1:i] ,
      dmin=0.55, dmax= 6 , pm=0.33,
      dose = dose_init[1:i],
      nb_ind = 3,
      n = i
    )
    
    
    inits1 <- 
      list(
        list(p0=0.1,gamma=4),
        list(p0=0.3,gamma=3.5),
        list(p0=0.2,gamma=1.5)
      )
    
    
    
    # Simulation des chaines de Markov
    model <- jags.model(file =textConnection(modelproject2),data = data4jags, n.chains=3,init=inits1)
    
    # Phase de chauffe
    update(model,3000)
    
    # Simulations monitorees
    # Tirage de 2000 valeurs de theta
    mcmc <- coda.samples(model,c("p0","gamma"),n.iter = 10000, thin = 3)
    
    # Estimation des parametres a posteriori
    S=summary(mcmc)
    plot(mcmc,trace = T, density = T)
    # gelman.diag(mcmc)
    #autocorr.plot(mcmc[[1]])
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    
    calc_p_tox <- function(mcmc,X){
      theta <- log((mcmc$p0*(1-0.33))/(0.33*(1-mcmc$p0)))/(0.5-mcmc$gamma)
      a0 <- log(0.33/(1-0.33))-log((mcmc$p0*(1-0.33))/(0.33*(1-mcmc$p0)))*mcmc$gamma/(0.5-mcmc$gamma)
      return (exp(a0 + theta*X)/(1 +  exp(a0 + theta*X)))
      
    }
    
    ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
    hist(ptox_tot[[4]])
    Snan=sum(is.nan(ptox_tot[[4]]))
    
    med_ptox=c()
    for(j in 1:5){
      med_ptox[j]<-median(ptox_tot[[j]])
    }
    
    dose_to_give[i] <- dose_tot[sum(med_ptox<=0.33)] 
    
  }  
  
  

