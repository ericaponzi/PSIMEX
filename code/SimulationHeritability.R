##############################################################
## The simulation extrapolation technique meets ecology and evolution: 
## A general and intuitive method to account for measurement error  
##
## R-code for simulation of pedigrees and P-SIMEX algorithm on 
## in simulated data on heritability
##############################################################


# load libraries
library("synbreed")
library("pedigree")
library("MCMCglmm")
library("mvtnorm")
library("GeneticsPed") 
library("pedigreemm")
library("AnimalINLA")
library("INLA")
library("Matrix")

# function to extract variances from INLA output (precisions)
get.variance<-function(prec)
{
  expect=inla.emarginal (function(x) 1/x, prec)
  stdev = sqrt((inla.emarginal(function(x) 1/x^2, prec)) 
               -(inla.emarginal(function(x) 1/x^1,prec))^2)
  q = inla.qmarginal(c(0.025,0.5,0.975),inla.tmarginal(function(x) 1/x, prec) )
  summarytab<-cbind(expect,stdev,q[1],q[2],q[3])
  colnames(summarytab)=c("mean","sd","0.025quant","0.5quant","0.975quant")
  return(summarytab)
  
}

#SIMEX FUNCTION
fun.simex <- function (lambda, nsim, known, unknown, ped0, founders_ped,  pedigree, trait){
  
  #setting initial values
  h <- matrix(data=NA, nrow=length(lambda), ncol=nsim)
  sd_h <- matrix(data=NA, nrow=length(lambda), ncol=nsim)
  VA <- matrix(data=NA, nrow=length(lambda), ncol=nsim)
  VE <- matrix(data=NA, nrow=length(lambda), ncol=nsim)
  sd_VA <- matrix(data=NA, nrow=length(lambda), ncol=nsim)
  sd_VE <- matrix(data=NA, nrow=length(lambda), ncol=nsim)
  
  for (k in 1:length(lambda)){
    for (j in 1:nsim) {
      animal <- known$animal
      dam <- known$dam
      sire <- known$sire  
      #add lambda% error
      #choose lambda% individuals and set their fathers as random
      #random fathers are chosen from the male individuals of the previous generation
      #excluding the real one
      m <- sample(1:length(known$animal), lambda[k]*length(known$animal))
      
      for (i in m) {
        generation <- ped0[which(ped0$id == known$animal[i]), 4]
        others <- ped0[ped0$generation == generation-1 & ped0$sex == 1 ,  ]$id
        others <- others[-which(others == known[i, 2])]
        sire[i] <- others[sample(1:length(others), 1)]
        
      }
      
      
      ped <- data.frame(animal = animal, dam = dam, sire = sire)
    
      #re calculate pedigree and inbreeding coefficient
      pedigree_sim <- rbind(unknown, ped)
      #prepare pedigree data
      pedigree <- data.frame(as.numeric(pedigree_sim$animal), as.numeric(pedigree_sim$sire), as.numeric(pedigree_sim$dam))
      for (i in 1:length(pedigree[,1])) {
        if (is.na(pedigree[i,1])) pedigree[i,1]<-0
        if (is.na(pedigree[i,2])) pedigree[i,2]<-0
        if (is.na(pedigree[i,3])) pedigree[i,3]<-0
      }
      names(pedigree) <- c('Individual', 'Parent1', 'Parent2')
      # calculate relatedness matrix from pedigree
      xx = compute.Ainverse(pedigree)
      Ainv = xx$Ainverse
      map  = xx$map
      Cmatrix = sparseMatrix(i=Ainv[,1],j=Ainv[,2],x=Ainv[,3])
      
      #check ID correspondence
      Ndata = dim(trait)[1]
      trait$IndexA = rep(0,Ndata)
      for(i in 1:Ndata)    trait$IndexA[i] = which(map[,1]==trait$ID[i])
      
      p.var <- var(trait$Y,na.rm = TRUE)
      
      #fit the inla model
      formula = Y ~ sex + 
        f(IndexA, model = "generic0", Cmatrix = Cmatrix, 
          hyper=list(prec=list(param = c(1/2,p.var/4))),
          constr=T) 
     
      model0 = inla(formula = formula, family = "gaussian",
                    data = trait,
                    control.family = list(hyper = list(prec =
                                                         list(param = c(p.var/2, p.var/2), fixed = FALSE))),
                    only.hyperparam =TRUE, control.compute = list(dic = T), num.threads = 2)
      
      #compute variances instead of precisions
      VA[k, j] <- get.variance(model0$marginals.hyperpar$"Precision for IndexA")[1]
      VE[k, j] <- get.variance(model0$marginals.hyperpar$"Precision for the Gaussian observations")[1]
      
      # extract samples 
      posterior.sample.var <- inla.hyperpar.sample(100000,model0)
      # compute heritability
      h2 <- (1/posterior.sample.var[,2]) /( 1/posterior.sample.var[,1] + 1/posterior.sample.var[,2] )
      h2 <- as.numeric(h2)
      h[k,j] <- posterior.mode(as.mcmc(h2))
      
      # standard errors
      sd_h[k, j] <- sd(h2)
      sd_VA[k, j] <- get.variance(model0$marginals.hyperpar$"Precision for IndexA")[2]
      sd_VE[k, j] <- get.variance(model0$marginals.hyperpar$"Precision for the Gaussian observations")[2]
      
      
      
      
    }
    
    
  }
  
  
  #store results
  res = list(h, sd_h, VA, 
             VE, sd_VA, sd_VE)
  names(res) <- c('h', 'sd_h', 'VA', 'VE', 'sd_VA', 'sd_VE')
  return(res)
}

###############################################################################
###############################SIMULATION#################################
###############################################################################

nped <- 50
h_true <- c()
VA_true <- c()
VE_true <- c()
sd_h_true <- c()
sd_VA_true <- c()
sd_VE_true <- c()
sd_h0 <- c()
h0 <- c()
VA0 <- c()
VE0 <- c()
sd_VA0 <- c()
sd_VE0 <- c()
h_pred <- c()
h_pred1 <- c()
h_pred2 <- c()
var_pred1 <- c()
var_pred2 <- c()
var_pred3 <- c()
h_ave  <- matrix(data  =  NA, ncol  =  10, nrow  =  nped)
sd_h <- matrix(data  =  NA, ncol  =  10, nrow  =  nped)
VA <- matrix(data  =  NA, ncol  =  10, nrow  =  nped)
VE <- matrix(data  =  NA, ncol  =  10, nrow  =  nped)
sd_VA <- matrix(data  =  NA, ncol  =  10, nrow  =  nped)
sd_VE <- matrix(data  =  NA, ncol  =  10, nrow  =  nped)


for (ii in 1: nped) {
  
  #simulate a pedigree
  ped0 <- generatePedigree(nId = 100, nGeneration = 30, nFather = 50, nMother = 50)
  
  #N individuals
  Nid <- length(ped0$id)
  animal.id <- 1 : Nid
  
  #fixed effects are mu, sex and inbreeding coefficient
  mu <- 10
  X <- rep(1, Nid)
  sex <- ped0$sex
  
  
  pedigree <- ped0[ , c(1,3,2)]
  names(pedigree) <- c("id", "dam", "sire")
  write.table(pedigree, paste("pedigree_", ii , ".txt", sep = '')) 
  
  #generate breeding values
  Z <- diag(Nid)
  sigmaA <- 0.3
  pedigree <- ped0[ ,1:3]
  u <- rbv(pedigree, sigmaA)
  
  #simulate environment
  sigmaE <- 0.1 # replace to obtain different values of h2
  R <- sigmaE*diag(Nid)
  e<- matrix(rmvnorm(n = 1, mean = rep(0,Nid), sigma = R))
  
  #simulate y
  Y <- X*mu+2*X*sex+Z%*%u+e
  
  #generate dataset
  trait <- data.frame(Y, sex, animal.id)
  trait$ID <- trait$animal.id
  trait$animal <- trait$animal.id
  trait$sex <- as.factor(sex)
  
  write.table(trait, paste('trait_', ii, ".txt", sep = ""))
  
  
  #INLA to compute heritability 
  #prepare pedigree data
  ped <- pedigree
  pedigree <- data.frame(as.numeric(pedigree$id), as.numeric(pedigree$father), as.numeric(pedigree$mother))
  for (i in 1:length(pedigree[,1])) {
    if (is.na(pedigree[i,1])) pedigree[i,1]<-0
    if (is.na(pedigree[i,2])) pedigree[i,2]<-0
    if (is.na(pedigree[i,3])) pedigree[i,3]<-0
  }
  names(pedigree) <- c('Individual', 'Parent1', 'Parent2')
  # compute relatedness matrix
  xx = compute.Ainverse(pedigree)
  Ainv = xx$Ainverse
  map  = xx$map
  Cmatrix = sparseMatrix(i=Ainv[,1],j=Ainv[,2],x=Ainv[,3])
  
  #ID
  Ndata = dim(trait)[1]
  trait$IndexA = rep(0,Ndata)
  for(i in 1:Ndata)    trait$IndexA[i] = which(map[,1]==trait$ID[i])
  
  p.var <- var(trait$Y,na.rm = TRUE)
  
  #fit the inla model
  formula = Y ~ sex + 
    f(IndexA, model = "generic0", Cmatrix = Cmatrix, 
      hyper=list(prec=list(param = c(1/2,p.var/4))),
      constr=T) 
  model = inla(formula = formula, family = "gaussian",
               data = trait,
               control.family = list(hyper = list(prec =
                                                  list(param = c(p.var/2, p.var/2), fixed = FALSE))),
               only.hyperparam =TRUE, control.compute = list(dic = T))
 
  #compute variances instead of precisions
  VA_true[ii] <- get.variance(model$marginals.hyperpar$"Precision for IndexA")[1]
  VE_true[ii] <- get.variance(model$marginals.hyperpar$"Precision for the Gaussian observations")[1]
  
  # extract samples 
  posterior.sample.var <- inla.hyperpar.sample(100000,model)
  # compute heritability
  h2 <- (1/posterior.sample.var[,2]) /( 1/posterior.sample.var[,1] + 1/posterior.sample.var[,2])
  h2 <- as.numeric(h2)
  h_true[ii] <- posterior.mode(as.mcmc(h2))
  # standard errors
  sd.true[ii] <- sd(h2)
  sd_VA_true[ii] <- get.variance(model$marginals.hyperpar$"Precision for IndexA")[2]
  sd_VE_true[ii] <- get.variance(model$marginals.hyperpar$"Precision for the Gaussian observations")[2]
 
  
  
   #generate some errors in the initial pedigree 
  #exclude founders from this
  founders <- ped0[which(is.na(ped0$father)), 1]
  founders_ped <- ped0[founders, 1:3 ]
  pedigree <- ped0[-founders, 1:3 ]
  names(founders_ped) <- c("animal", "sire", "dam")
  names(pedigree) <- c("animal","sire", "dam")
  
  
  #add 10% error
  #choose 10% individuals and set their fathers as random
  #random fathers are chosen from the male individuals of the previous generation
  #excluding the real one
  m <- sample(1:length(pedigree$animal), 0.1*length(pedigree$animal))
  
  for (i in m) {
    generation <- ped0[which(ped0$id == pedigree$animal[i]), 4]
    others <- ped0[ped0$generation == generation-1 & ped0$sex == 1 ,  ]$id
    others <- others[-which(others == pedigree[i, 2])]
    pedigree$sire[i] <- others[sample(1:length(others), 1)]
    
  }
  
  
  #re calculate pedigree and inbreeding coefficient
  pedigree_sim <- rbind(founders_ped, pedigree)
  #prepare pedigree data
  pedigree <- data.frame(as.numeric(pedigree_sim$animal), as.numeric(pedigree_sim$sire), as.numeric(pedigree_sim$dam))
  for (i in 1:length(pedigree[,1])) {
    if (is.na(pedigree[i,1])) pedigree[i,1]<-0
    if (is.na(pedigree[i,2])) pedigree[i,2]<-0
    if (is.na(pedigree[i,3])) pedigree[i,3]<-0
  }
  names(pedigree) <- c('Individual', 'Parent1', 'Parent2')
  # compute relatedness matrix
  xx = compute.Ainverse(pedigree)
  Ainv = xx$Ainverse
  map  = xx$map
  Cmatrix = sparseMatrix(i=Ainv[,1],j=Ainv[,2],x=Ainv[,3])
  
  #ID
  Ndata = dim(trait)[1]
  trait$IndexA = rep(0,Ndata)
  for(i in 1:Ndata)    trait$IndexA[i] = which(map[,1]==trait$ID[i])
  p.var <- var(trait$Y,na.rm = TRUE)
  
  #fit the inla model
  formula = Y ~ sex + 
    f(IndexA, model = "generic0", Cmatrix = Cmatrix, 
      hyper=list(prec=list(param = c(1/2,p.var/4))),
      constr=T) 
  
  model0 = inla(formula = formula, family = "gaussian",
               data = trait,
               control.family = list(hyper = list(prec =
                                                    list(param = c(p.var/2, p.var/2), fixed = FALSE))),
               only.hyperparam =TRUE, control.compute = list(dic = T))
  
  #compute variances instead of precisions
  VA0[ii] <- get.variance(model0$marginals.hyperpar$"Precision for IndexA")[1]
  VE0[ii] <- get.variance(model0$marginals.hyperpar$"Precision for the Gaussian observations")[1]
 
  # extract samples
  posterior.sample.var <- inla.hyperpar.sample(100000,model0)
  # compute heritability
  h2 <- (1/posterior.sample.var[,2]) /( 1/posterior.sample.var[,1] + 1/posterior.sample.var[,2])
  h0[ii] <- posterior.mode(as.mcmc(h2))
  
  # standard errors
  sd_h0[ii] <- sd(h2)
  sd_VA0[ii] <- get.variance(model0$marginals.hyperpar$"Precision for IndexA")[2]
  sd_VE0[ii] <- get.variance(model0$marginals.hyperpar$"Precision for the Gaussian observations")[2]
  
  #implement SIMEX on this data
  #library(parallel)
  nsim <- 100
  
  #where to generate more errors
  known <- pedigree_sim[which(!is.na(pedigree_sim$dam)), ]
  #already missing
  unknown <- pedigree_sim[which(is.na(pedigree_sim$dam)), ]
  
  #increasing percentages of errors
  lambda0 <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  #subtract the ones already missing
  lambda <-  1- (1-lambda0)/(1-0.1)
  
  
  results <- 
   fun.simex(lambda= lambda, nsim = nsim, known = known,
              unknown = unknown, pedigree = pedigree_sim, ped0 = ped0,  founders_ped = founders_ped,  
              trait = trait)

  #compute means
  h_ave_sim <- rowMeans(results$h, na.rm = TRUE)
  sd_h_ave <- rowMeans(results$sd_h, na.rm = TRUE)
  VA_ave <- rowMeans(results$VA, na.rm = TRUE)
  VE_ave <- rowMeans(results$VE, na.rm = TRUE)
  sd_VA_ave <- rowMeans(results$sd_VA, na.rm = TRUE)
  sd_VE_ave <- rowMeans(results$sd_VE, na.rm = TRUE)
  
  
  #merge all values
  lambda <- c(0.1, lambda0)
  h_ave[ii,] <- c(h0[ii], as.numeric(h_ave_sim))
  sd_h[ii,] <- c(sd_h0[ii], sd_h_ave)
  VA[ii,] <- c(VA0[ii],VA_ave)
  VE[ii,] <- c(VE0[ii], VE_ave)
  sd_VA[ii,] <- c(sd_VA0[ii], sd_VA_ave)
  sd_VE[ii,] <- c(sd_VE0[ii], sd_VE_ave)
  
  
  
  ############EXTRAPOLATION##########
  

  
  p.names <- c("h", "sd_h")
  estimates <- data.frame(h_ave[ii, ], sd_h[ii, ])
  colnames(estimates) <- p.names
  
  #linear case
  
  extrapolation_h <- lm(estimates[ ,1] ~ lambda)
  h_pred[ii] <- predict(extrapolation_h, newdata = data.frame(lambda = 0))
 
  #quadratic case
  
  extrapolation_h1 <- lm(estimates[ ,1] ~ lambda+ I(lambda^2))
  h_pred1[ii] <- predict(extrapolation_h1, newdata = data.frame(lambda = 0))
 
  #non linear case
  
  extrapolation_h2 <- lm(estimates[ ,1] ~ lambda+ I(lambda^2)+ I(lambda^3))
  h_pred2[ii] <- predict(extrapolation_h2, newdata = data.frame(lambda = 0))
 
  
  
  
  
  # calculate SIMEX SE
  
  # first component 
  # sampling variability
  
  S <- c()
  
  #one per lambda 
  for ( i in 1:9) {
    diff <- c()
    #calculate the differences per each simulation
    for ( j in 1: nsim ) {
      diff[j] <- results$h[i,j]-as.numeric(h_ave_sim[i])
    }
    
    diff <- na.omit(diff)
    S[i] <- 1/(nsim-1)*sum(diff^2) 
    
  }
 
  
  S <- c(0, S)
  
  S1 <- sd_h[ii, ]
  S2 <- S
  Stot <- S1 - S2
  
  # extrapolated value
  # linear case
  extrapolation_var<- lm(Stot ~ lambda)
  var_pred[ii] <- predict(extrapolation_var, newdata = data.frame(lambda = 0))
  
  #quadratic case
  extrapolation_var1<- lm(Stot ~ lambda+ I(lambda^2))
  var_pred1[ii] <- predict(extrapolation_var1, newdata = data.frame(lambda = 0))
  
  #cubic case
  extrapolation_var2<- lm(Stot ~ lambda+ I(lambda^2)+ I(lambda^3))
  var_pred2[ii] <- predict(extrapolation_var2, newdata = data.frame(lambda = 0))
  
  library(AICcmodavg)
  AICl[ii] <- AICc(extrapolation_h)
  
  AICq[ii] <- AICc(extrapolation_h1)
  
  AICc[ii] <- AICc(extrapolation_h2)
  
 
  
}


extrapolation_data <- data.frame(h_true, h_pred, h_pred1, h_pred2, sqrt(var_ext1), sqrt(var_ext2), sqrt(var_ext3), 
                                 AICl, AICq, AICc)

# bias 
bias_linear <- mean(extrapolation_data$h_pred - extrapolation_data$h_true)
bias_quad  <- mean(extrapolation_data$h_pred1 - extrapolation_data$h_true)
bias_cubic  <- mean(extrapolation_data$h_pred2 - extrapolation_data$h_true)


# MSE
MSE_linear <- mean((extrapolation_data$h_pred-extrapolation_data$h_true)^2)
MSE_quad <- mean((extrapolation_data$h_pred1-extrapolation_data$h_true)^2)
MSE_cubic <- mean((extrapolation_data$h_pred2-extrapolation_data$h_true)^2)

