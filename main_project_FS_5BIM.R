#### doesn't work
#### changes here
#### weird


rm(list=ls())

## PROJET STATS BAYESIENNES

library(rjags)

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
  
  inits1 <- 
    list(
      list(theta = 0.5),
      list(theta = 1.5),
      list(theta = 2)
    )
  
  for(i in 1:(nb_cohorte-1)){
    
    data4jags <- list(
      ntox = ntox[1:i] ,
      a0 = 1,
      dose = log(pprior_dose[1:i]/(1 - pprior_dose[1:i])) - a0,
      nb_ind = 3,
      n = i
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
  
  par(mfrow=c(3,2))  #### changes here
  for(p in 1:5){  #### changes here
    hist(ptox0[,p],main=paste("Dose :",as.character(dose_tot[p])), col = "grey",  #### changes here
         xlab = "probabilités de toxicité a priori")  #### changes here
    abline(v=pprior_dose_tot[p],col="red", lwd=2)  #### changes here
  }  #### changes here
  par(mfrow=c(1,1))  #### changes here
  
  
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
  
  
  for(i in 1:(nb_cohorte-1)){
    
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
    plot(mcmc,trace = F, density = T)
    # gelman.diag(mcmc)
    #autocorr.plot(mcmc[[1]])
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
    hist(ptox_tot[[3]])
    hist(ptox_tot[[4]])
    
    med_ptox=c()
    for(j in 1:5){
      med_ptox[j]<-median(ptox_tot[[j]])
    }
    
    dose_to_give[i] <- dose_tot[sum(med_ptox<=0.33)] 
    
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
  plot(mcmc0,trace=F)
  
  ptox2_0 = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot0)
  ptox2_0 <- as.data.frame(ptox2_0)
  
  par(mfrow=c(5,2)) #### doesn't work
  for(p in 1:5){
    hist(ptox0[,p],main=paste("Modele 1, Dose :",as.character(dose_tot[p])),breaks=30)
    abline(v=pprior_dose_tot[p],col="red")
    hist(ptox2_0[,p],main=paste("Modele 2, Dose :",as.character(dose_tot[p])),breaks=30)
    abline(v=pprior_dose_tot[p],col="red")
  }
  
  ################ CHANGEMENT REGLE D'ARRET = Question 2
  ##### REGLE 1 ######
  # On calcule la ptox mediane de la dose minimale ainsi que son intervalle de confiance à 90%, #### weird : crédibilité
  # il faut que l'intervalle soit superieur a 0.33 ie quantile(5%) >0.33
  
  # Modele 1 
  
  dose_to_give <- c()
  pprior_dose=c(0.05,0.15,0.333,0.333,0.333,0.333)
  a0 = 1
  ntox = c(0,1,1,0,1,2)
  dose_tot = c(0.5,1,3,5,6)
  pprior_dose_tot = c(0.05,0.1,0.15,0.333,0.5)
  dose_star_tot = log(pprior_dose_tot/(1 - pprior_dose_tot)) - a0
  
  nb_cohorte = 6
  
  for(i in 1:(nb_cohorte-1)){
    #for(i in 1){
    
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
    update(model,3000)
    mcmc <- coda.samples(model,c("theta"),n.iter = 10000, thin = 3)
    
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    ptox_tot = lapply(X = dose_star_tot, FUN = function(X) {median(exp(a0 + mcmctot$theta*X)/(1 +  exp(a0 + mcmctot$theta*X)))})
    dose_to_give[i] <- dose_tot[sum(ptox_tot<=0.33)] 
    
    ptox_dose_min=exp(a0 + mcmctot$theta*dose_star_tot[1])/(1 +  exp(a0 + mcmctot$theta*dose_star_tot[1]))
    hist(ptox_dose_min) #### doesn't work
    ptox_inf =  quantile( ptox_dose_min,0.05)
  }
  
  dose_to_give
  
  
  # Modele 2
  
  
  dose_to_give <- c()
  dose_init=c(0.5,3,5,5,5,5)
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
  
  
  for(i in 1:(nb_cohorte-1)){
    
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
    mcmc <- coda.samples(model,c("p0","gamma"),n.iter = 10000, thin = 3)
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_inf <- quantile(mcmctot$p0,0.05)
    
    dose_to_give[i] <- quantile( ptox_dose_min,0.05)
    
  }  
  
  dose_to_give  #### weird      1.874597e-06 1.874597e-06 1.874597e-06 1.874597e-06 1.874597e-06
  
  
  
  
  
  
  ##### REGLE 2 ######
  # On calcule la ptox mediane de la dose maximale ainsi que son intervalle de confiance à 90%, 
  # il faut que l'intervalle soit inferieur a 0.33 ie quantile(95%) <0.33
  
  # Modele 1 
  
  dose_to_give <- c()
  pprior_dose=c(0.05,0.15,0.333,0.333,0.333,0.333)
  a0 = 1
  ntox = c(0,1,1,0,1,2)
  dose_tot = c(0.5,1,3,5,6)
  pprior_dose_tot = c(0.05,0.1,0.15,0.333,0.5)
  dose_star_tot = log(pprior_dose_tot/(1 - pprior_dose_tot)) - a0
  
  nb_cohorte = 6
  
  for(i in 1:(nb_cohorte-1)){
    #for(i in 1){
    
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
    update(model,3000)
    mcmc <- coda.samples(model,c("theta"),n.iter = 10000, thin = 3)
    
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    ptox_tot = lapply(X = dose_star_tot, FUN = function(X) {median(exp(a0 + mcmctot$theta*X)/(1 +  exp(a0 + mcmctot$theta*X)))})
    dose_to_give[i] <- dose_tot[sum(ptox_tot<=0.33)] 
    
    ptox_dose_max=exp(a0 + mcmctot$theta*dose_star_tot[5])/(1 +  exp(a0 + mcmctot$theta*dose_star_tot[5]))
    hist(ptox_dose_max)  #### doesn't work
    ptox_sup =  quantile( ptox_dose_max,0.90)
  }
  
  dose_to_give
  
  
  # Modele 2
  
  
  dose_to_give <- c()
  dose_init=c(0.5,3,5,5,5,5)
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
  
  
  for(i in 1:(nb_cohorte-1)){
    
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
    mcmc <- coda.samples(model,c("p0","gamma"),n.iter = 10000, thin = 3)
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_dose_max = calc_p_tox(mcmctot,dose_tot[5])
    ptox_sup <- quantile( ptox_dose_max,0.05)
    
  }   #### doesn't work
  # Error in quantile.default(ptox_dose_max, 0.05) : 
  #  missing values and NaN's not allowed if 'na.rm' is FALSE
  
  dose_to_give #### doesn't work
  
  
  
  
  
  ##### REGLE 3 ######
  # La dose retenue restera incangée pour toutes les cohortes qui suivent.
  
  
  # Modele 2
  
  #################################################################
  
  
  dose_to_give <- c()
  # dose_init=c(0.5,3,5,5,5,5)
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
  
  
  
  data4jags <- list(
    ntox = ntox[1] ,
    dmin=0.55, dmax= 6 , pm=0.33,
    dose = c(0.5,rep(dose_inti_found,5))[1],
    nb_ind = 3,
    n = 1
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
  plot(mcmc,trace = F, density = T)
  # gelman.diag(mcmc)
  #autocorr.plot(mcmc[[1]])
  
  mcmctot = as.data.frame(as.matrix(mcmc))
  
  ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
  hist(ptox_tot[[3]]) #### doesn't work
  hist(ptox_tot[[4]]) #### doesn't work
  
  med_ptox=c()
  for(j in 1:5){
    med_ptox[j]<-median(ptox_tot[[j]])
  }
  
  dose_to_give[1] <- dose_tot[sum(med_ptox<=0.33)] 
  
  dose_init_found = dose_tot[sum(med_ptox<=0.33)] 
  
  
  
  
  for(i in 2:(nb_cohorte-1)){
    
    data4jags <- list(
      ntox = ntox[1:i] ,
      dmin=0.55, dmax= 6 , pm=0.33,
      dose = c(0.5,rep(dose_init_found,5))[1:i],
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
    plot(mcmc,trace = F, density = T)
    # gelman.diag(mcmc)
    #autocorr.plot(mcmc[[1]])
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
    # hist(ptox_tot[[3]]) #### doesn't work
    # hist(ptox_tot[[4]]) #### doesn't work
    
    med_ptox=c()
    for(j in 1:5){
      med_ptox[j]<-median(ptox_tot[[j]])
    }
    
    dose_to_give[i] <- dose_tot[sum(med_ptox<=0.33)] 
    
  }  
  
  dose_to_give #### weird
  
  
  
  
  
  
  
  
  ##### REGLE 4 ######
  # La variation relative de l'estimation de la probabilité de toxicité sévère entre 
  # la dernière cohorte et la cohorte en cours n'exèdera pas 5%
  
  
  # Modele 2
  
  #################################################################
  
  
  dose_to_give <- c()
  dose_init=c(0.5,3,5,5,5,5)
  ntox = c(0,1,1,0,1,2)
  dose_tot = c(0.5,1,3,5,6)
  nb_cohorte = 6
  variation = c()
  
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
  
  ptox_tab <- c()
  
  for(i in 1:(nb_cohorte-1)){
    
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
    plot(mcmc,trace = F, density = T)
    # gelman.diag(mcmc)
    #autocorr.plot(mcmc[[1]])
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
    # hist(ptox_tot[[3]])
    # hist(ptox_tot[[4]])
    
    med_ptox=c()
    for(j in 1:5){
      med_ptox[j]<-median(ptox_tot[[j]])
    }
    
    ptox_tab[i] <- list(med_ptox)
    
    if(i!=1){
      variation[i-1] = abs(ptox_tab[[i]][i-1]-ptox_tab[[i-1]][i-1])/(ptox_tab[[i]][i-1]+ptox_tab[[i-1]][i-1])
      # if (variation < 0.05)
      
    }
    
    dose_to_give[i] <- dose_tot[sum(med_ptox<=0.33)] 
    
  }  
  
  dose_to_give
  variation
  
  #################################################################
  

  
  
  
  
  
  
  ##### REGLE 4 ######
  # La variation relative de la largeur l'intervalle de crédibilité de 
  # la probabilité de toxicité sévère entre la dernière cohorte et 
  # la cohorte en cours n'exèdera pas 5%
  
  
  # Modele 2
  
  #################################################################
  
  
  dose_to_give <- c()
  dose_init=c(0.5,3,5,5,5,5)
  ntox = c(0,1,1,0,1,2)
  dose_tot = c(0.5,1,3,5,6)
  nb_cohorte = 6
  variation_cred = c()
  width_credibility = c()
  
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
  
  ptox_tab <- c()
  
  for(i in 1:(nb_cohorte-1)){
    
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
    plot(mcmc,trace = F, density = T)
    # gelman.diag(mcmc)
    #autocorr.plot(mcmc[[1]])
    
    mcmctot = as.data.frame(as.matrix(mcmc))
    
    ptox_tot = lapply(X = dose_tot, FUN = calc_p_tox, mcmc=mcmctot)
    # hist(ptox_tot[[3]])
    # hist(ptox_tot[[4]])
    
    med_ptox=c()
    for(j in 1:5){
      med_ptox[j]<-median(ptox_tot[[j]])
    }
    
    width_credibility[i] <- abs(quantile( med_ptox,0.025) - quantile( med_ptox,0.975))
    
    if(i!=1){
      variation_cred[i-1] = abs(width_credibility[i-1]-width_credibility[i])/(width_credibility[i-1] + width_credibility[i])
      # if (variation < 0.05)
      
    }
    
    dose_to_give[i] <- dose_tot[sum(med_ptox<=0.33)] 
    
  }  
  
  dose_to_give
  width_credibility
  variation_cred
  
  #################################################################
  
  
  
