
rm(list=ls())

## Modele 2

library(rjags)


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
  regle1 <- c()
  regle2 <- c()
  
  dose_init=c(0.5,3,5,5,5,5)
  pprior_dose_tot = c(0.05,0.1,0.15,0.333,0.5)
  ntox = c(0,1,1,0,1,2)
  dose_tot = c(0.5,1,3,5,6)
  nb_cohorte = 6
  
  inits1 <- 
    list(
      list(p0=0.1,gamma=4),
      list(p0=0.3,gamma=3.5),
      list(p0=0.2,gamma=1.5)
    )
  
  calc_p_tox <- function(mcmc,X){
    theta <- log((mcmc$p0*(1-0.33))/(0.33*(1-mcmc$p0)))/(0.5-mcmc$gamma)
    a0 <- log(0.33/(1-0.33))-log((mcmc$p0*(1-0.33))/(0.33*(1-mcmc$p0)))*mcmc$gamma/(0.5-mcmc$gamma)
    return (exp(a0 + theta*X)/(1 +  exp(a0 + theta*X)))
    
  }
  
  par(mfrow=c(3,2))
  for(i in 1:(nb_cohorte)){
    
    data4jags <- list(
      ntox = ntox[1:i] ,
      dmin=0.55, dmax= 6 , pm=0.33,
      dose = dose_init[1:i],
      nb_ind = 3,
      n = i
    )
    
    # Simulation des chaines de Markov
    model <- jags.model(file =textConnection(modelproject2),data = data4jags, n.chains=3,init=inits1)
    
    # Phase de chauffe
    update(model,3000)
    
    # Simulations monitorees
    # Tirage de 2000 valeurs de theta
    mcmc <- coda.samples(model,c("p0","gamma"),n.iter = 10000, thin = 3)
    
    # Estimation des parametres a posteriori
    summary(mcmc)
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
    
    med_ptox=c()
    for(j in 1:5){
      med_ptox[j]<-median(ptox_tot[[j]])
    }
    
    dose_to_give[i] <- dose_tot[sum(med_ptox<=0.33)] 

    
    # Regle d'arret 1
    regle1[i] =  sum(mcmctot$p0>0.33)/length(mcmctot$p0)
    hist(mcmctot$p0,main= paste("Dose minimale, Cohorte ",as.character(i)),xlab="Probabilites de toxicite predites",col = "grey")
    abline(v=0.33,col="red", lwd=2)  
    
    # Regle d'arret 2
    ptox_dose_max = calc_p_tox(mcmctot,dose_tot[5])
    regle2[i] =  sum(ptox_dose_max<0.33)/length(ptox_dose_max)
    hist(ptox_dose_max,main=paste("Dose maximale, Cohorte ",as.character(i)),xlab="Probabilites de toxicite predites",col = "grey")
    abline(v=0.33,col="red", lwd=2) 
  }  
  
  dose_to_give
  
  
  ## Sans données
  
  d0 <- list( dmin=0.55, dmax= 6 , pm=0.33,
              dose = dose_init,
              nb_ind = 3,
              n = 6)
  
  model0 <- jags.model(file = textConnection(modelproject2), data= d0,
                       inits = inits1, n.chains = 3)
  
  
  update(model0, 3000)
  mcmc0 <- coda.samples(model0,c("p0","gamma"),n.iter = 10000, thin = 3)
  mcmctot0 <- as.data.frame(as.matrix(mcmc0))
  #plot(mcmc0,trace=F)
  
  ptox0 = sapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot0)
  ptox0 <- as.data.frame(ptox0)
  
  par(mfrow=c(3,2))  
  for(p in 1:5){  
    hist(ptox0[,p],main=paste("Dose :",as.character(dose_tot[p])), col = "grey",  
         xlab = "Probabilites de toxicite a priori") 
    abline(v=pprior_dose_tot[p],col="red", lwd=2)  
  }  
  
  ##### REGLE 1 ######
  # On calcule la ptox mediane de la dose minimale ainsi que son intervalle de confiance à 90%, 
  # il faut que l'intervalle soit superieur a 0.33 ie quantile(5%) >0.33
  if (sum(regle1>=0.9)==0) {
    "Regle 1 : Pas d'arrêt, la probabilite de toxicite severe est toujours inferieure a 33% pour la dose minimale"}else  {
      paste(" Regle 1 : On s'arrete a la cohorte",toString(min(which(regle1 >= 0.9))))}
  
  ##### REGLE 2 ######
  # On calcule la ptox mediane de la dose maximale ainsi que son intervalle de confiance à 90%, 
  # il faut que l'intervalle soit inferieur a 0.33 ie quantile(95%) <0.33
  if (sum(regle2>=0.9)==0) {
    "Regle 2 : Pas d'arrêt, la probabilite de toxicite severe est toujours superieur a 33% pour la dose maximale"}else  {
      paste(" Regle 2 : On s'arrete a la cohorte",toString(min(which(regle1 >= 0.9))))}
  
