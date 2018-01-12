
rm(list=ls())

## Modele 1

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


dose_to_give <- c()
regle1<- c()
regle2<- c()

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

for(i in 1:(nb_cohorte)){
  
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
  
  mcmc <- coda.samples(model,c("theta"),n.iter = 10000, thin = 3)
  
  mcmctot = as.data.frame(as.matrix(mcmc))
  
  ptox_tot = lapply(X = dose_star_tot, FUN = function(X) {median(exp(a0 + mcmctot$theta*X)/(1 +  exp(a0 + mcmctot$theta*X)))})
  dose_to_give[i] <- dose_tot[sum(ptox_tot<=0.33)] 
  
  
  #Regle arret 1
  ptox_dose_min=exp(a0 + mcmctot$theta*dose_star_tot[1])/(1 +  exp(a0 + mcmctot$theta*dose_star_tot[1]))
  regle1[i] =  sum(ptox_dose_min>0.33)/length(ptox_dose_min)
  
  #Regle arret 2
  ptox_dose_max=exp(a0 + mcmctot$theta*dose_star_tot[5])/(1 +  exp(a0 + mcmctot$theta*dose_star_tot[5]))
  regle2[i] =  sum(ptox_dose_max<0.33)/length(ptox_dose_max)
  hist(ptox_dose_max)
  
}

dose_to_give

#Sans donnees

d0 <- list( a0 = 1,
            dose = log(pprior_dose/(1 - pprior_dose)) - a0,
            nb_ind = 3,
            n = 6)

model0 <- jags.model(file = textConnection(modelproject), data= d0,
                     inits = inits1, n.chains = 3)


update(model0, 3000)
mcmc0 <- coda.samples(model0, c("theta"), n.iter = 2000 )
mcmctot0 <- as.data.frame(as.matrix(mcmc0))
#plot(mcmc0,trace=F)

ptox0 <- lapply(X = dose_star_tot, FUN = function(X) {exp(a0 + mcmctot0$theta*X)/(1 +  exp(a0 + mcmctot0$theta*X))})
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

