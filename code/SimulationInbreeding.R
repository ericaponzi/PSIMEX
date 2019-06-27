##############################################################
## The simulation extrapolation technique meets ecology and evolution: 
## A general and intuitive method to account for measurement error  
##
## R-code for simulation of pedigrees and P-SIMEX algorithm on 
## in simulated data on inbreeding depression
##############################################################



library("synbreed")
library("pedigree")
library("MCMCglmm")
library("mvtnorm")
library("GeneticsPed") 
library("pedigreemm")


#SIMULATION FUNCTION
fun.simex <- function (lambda, nsim, known, unknown, ped0, founders_ped,  pedigree, trait){
  
  #setting initial values
  inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol  =  nsim)
  se_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol  =  nsim)
  mean_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)
  median_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)
  var_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)
  pval <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)
  
  for (k in 1:length(lambda)){
    for (j in 1:nsim) {
      animal <- known$animal
      dam <- known$dam
      sire <- known$sire  
      
      #initializing inbreeding coefficient
      trait$f_sim <- c()
      
      #choosing lambda over the total number of individuals 
      #and setting their parents as random
      m <- sample(1:length(known$animal), lambda[k]*length(known$animal))
      for (i in m) {
        
        generation <- ped0[which(ped0$id == known$animal[i]), 4]
        others <- ped0[ped0$generation == generation-1 & ped0$sex == 1,  ]$id
        others <- others[-which(others == known[i, 2])]
        sire[i] <- others[sample(1:length(others), 1)]
        
      }
      
      
      #recalculate pedigree
      ped <- data.frame(animal = animal, dam = dam, sire = sire)
      ped <- rbind(ped, unknown, founders_ped)
      ord <- orderPed(ped)
      ped <- ped[order(ord),]
      
      #re calculate inbreeding coefficient
      f_sim <- data.frame(ped$animal, calcInbreeding(ped))
      names(f_sim) <- c("animal", "f_sim")
      trait <- merge(trait, f_sim, by = "animal")
      trait$f_sim <- as.numeric(trait$f_sim)
      
      if(sum(trait$f_sim)>0){
        
        #linear model 
        #inbreeding depression is the regression coefficient of f_sim
        model <- lm(Y ~ f_sim +sex,
                    data  =  trait )
        
        #inbreeding depression
        inb[k, j] <- summary(model)$coefficients[2,1]
        se_inb[k,j] <- summary(model)$coefficients[2,2]
        pval [k, j] <- summary(model)$coefficients[2,4]
        
        #info on inbreeding coefficient
        mean_inb[k, j] <- mean(trait$f_sim)
        median_inb[k, j] <- median(trait$f_sim)
        var_inb[k, j] <- var(trait$f_sim)
        
      }
      
      #this must be checked
      if(sum(trait$f_sim) == 0) {
        inb[k, j] <- NA
        se_inb[k,j] <- NA
        mean_inb[k, j] <- NA
        median_inb[k, j] <- NA
        var_inb[k, j] <- NA
        pval [k, j] <- NA
        
      }
      
    }
    
  }
  
  #store results
  res = list(inb, se_inb, pval, 
             mean_inb, median_inb, var_inb)
  names(res) <- c('inb', 'se_inb', 'pval', 'mean_inb', 'median_inb', 'var_inb')
  return(res)
}

###############################################################################
###############################SIMULATION#################################
###############################################################################


nped <- 100
inb_dep1 <- c()
se_inb1 <- c()
p_val1 <- c()
mean_inb1 <- c()
median_inb1 <- c()
var_inb1 <-c()
inb_pred <- c()
inb_pred_up <- c()
inb_pred_low <- c()
inb_pred1 <- c()
inb_pred1_up <- c()
inb_pred1_low <- c()
inb_pred2 <- c()
inb_pred2_up <- c()
inb_pred2_low <- c()

inb_ave <- matrix(data  =  NA, ncol  =  8, nrow  =  nped)
se_inb <- matrix(data  =  NA, ncol  =  8, nrow  =  nped)
mean_inb <- matrix(data  =  NA, ncol  =  8, nrow = nped)
median_inb <- matrix(data  =  NA, ncol  =  8, nrow = nped)
var_inb <- matrix(data  =  NA, ncol  =  8, nrow = nped)
p_val <- matrix(data  =  NA, ncol  =  8, nrow = nped)




for (ii in 1: nped) {
  
  #simulate a pedigree
  ped0 <- generatePedigree(nId = 100, nGeneration = 30, nFather = 50, nMother = 50)
  # change nFather and nMother to obtain different Ne/Nc values
  
  #N individuals
  Nid <- length(ped0$id)
  animal.id <- 1 : Nid
  
  #fixed effects are mu, sex and inbreeding coefficient
  mu <- 10
  X <- rep(1, Nid)
  sex <- ped0$sex
  
  #calculate the inbreeding coefficient from pedigree
  pedigree <- ped0[ , c(1,3,2)]
  names(pedigree) <- c("id", "dam", "sire")
  f_inb <- calcInbreeding(pedigree)
  
  #generate breeding values
  Z <- diag(Nid)
  sigmaA <- 0.33
  pedigree <- ped0[ ,1:3]
  u <- rbv(pedigree, sigmaA)
  
  #simulate environment
  sigmaE <- 0.1
  R <- sigmaE*diag(Nid)
  e<- matrix(rmvnorm(n = 1, mean = rep(0,Nid), sigma = R))
  
  #simulate y
  Y <- X*mu+2*X*sex-7*X*f_inb+Z%*%u+e
  
  #generate dataset
  trait <- data.frame(Y, sex, f_inb, animal.id)
  trait$ID <- trait$animal.id
  trait$animal <- trait$animal.id
  trait$sex <- as.factor(sex)
  
  #linear model to compute inb depression
  #it is the beta of f_inb
  model <- lm(Y ~ f_inb +sex,
              data  =  trait )
  
  #inbreeding depression
  inb_dep1[ii] <- summary(model)$coefficients[2,1] 
  se_inb1[ii] <- summary(model)$coefficients[2,2]
  p_val1[ii] <- summary(model)$coefficients[2,4]
  
  #inbreeding coefficient
  mean_inb1[ii] <- mean(trait$f_inb)
  median_inb1[ii] <- median(trait$f_inb)
  var_inb1[ii] <- var(trait$f_inb)
  
  
  
  #generate some errors in the initial pedigree 
  #exclude founders from this
  founders <- ped0[which(is.na(ped0$father)), 1]
  founders_ped <- ped0[founders, 1:3 ]
  pedigree <- ped0[-founders, 1:3 ]
  names(founders_ped) <- c("animal", "dam", "sire")
  names(pedigree) <- c("animal", "dam", "sire")
  
  
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
  f_inb1 <- data.frame(pedigree_sim$animal, calcInbreeding(pedigree_sim))
  names(f_inb1) <- c("animal", "f_inb1")
  trait <- merge(trait, f_inb1, by = "animal")
  trait$f_sim <- as.numeric(trait$f_inb1)
  
  #compute initial values (error prone values)
  model0 <-  lm(Y ~ f_sim + sex, 
                data  =  trait )
  
  #inbreeding depression
  inb_dep0 <- summary(model0)$coefficients[2,1]
  se_inb0 <- summary(model0)$coefficients[2,2]
  pval0 <- summary(model0)$coefficients[2,4]
  
  #inbreeding coefficient
  mean_inb0 <- mean(trait$f_sim)
  median_inb0 <- median(trait$f_sim)
  var_inb0 <- var(trait$f_sim)
  
  
  #implement SIMEX 
  nsim <- 100
  
  #where to generate more errors
  known <- pedigree[which(!is.na(pedigree$dam)), ]
  #already missing
  unknown <- pedigree[which(is.na(pedigree$dam)), ]
  
  #increasing percentages of errors
  lambda0 <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  #subtract the ones already missing
  lambda <- lambda0-0.1
  
  #parallel computing
 results <- fun.simex(lambda= lambda, nsim = nsim, known = known,
              unknown = unknown, pedigree = pedigree_sim, ped0 = ped0,  founders_ped = founders_ped,  
              trait = trait)
 
  
  
  #compute means
  inb_ave_sim <- rowMeans(results$inb, na.rm = TRUE)
  se_inb_ave <- rowMeans(results$se_inb, na.rm = TRUE)
  mean_inb_ave <- rowMeans(results$mean_inb, na.rm = TRUE)
  median_inb_ave <- rowMeans(results$median_inb, na.rm = TRUE)
  var_inb_ave <- rowMeans(results$var_inb, na.rm = TRUE)

  
  #merge all values
  lambda <- c(0.1, lambda0)
  inb_ave[ii ,] <- c(inb_dep0, as.numeric(inb_ave_sim))
  se_inb[ii ,] <- c(se_inb0, se_inb_ave)
  mean_inb[ii ,] <- c(mean_inb0,mean_inb_ave)
  median_inb[ii ,] <- c(median_inb0, median_inb_ave)
  var_inb[ii ,] <- c(var_inb0, var_inb_ave)

  
  
  ############EXTRAPOLATION##########
  
  #inbreeding depression
  p.names <- c("inb", "up_inb", "low_inb")
  low <- inb_ave-1.96*se_inb
  up <- inb_ave+1.96*se_inb
  estimates <- data.frame(inb_ave[ii, ], up[ii, ], low[ii, ])
  colnames(estimates) <- p.names
  
  #linear case
  
  extrapolation_inb <- lm(estimates[ ,1] ~ lambda)
  inb_pred[ii] <- predict(extrapolation_inb, newdata = data.frame(lambda = 0))
  
  #quadratic case
  
  extrapolation_inb1 <- lm(estimates[ ,1] ~ lambda+ I(lambda^2))
  inb_pred1[ii] <- predict(extrapolation_inb1, newdata = data.frame(lambda = 0))
  
  #cubic case
  
  extrapolation_inb2 <- lm(estimates[ ,1] ~ lambda+ I(lambda^2)+ I(lambda^3))
  inb_pred2[ii] <- predict(extrapolation_inb2, newdata = data.frame(lambda = 0))
  
  
}

extrapolation_data <- data.frame(inb_dep1, inb_pred, inb_pred1, inb_pred2)




# bias 
bias_linear <- mean(extrapolation_data$inb_pred - extrapolation_data$inb_dep1)
bias_quad  <- mean(extrapolation_data$inb_pred1 - extrapolation_data$inb_dep1)
bias_cubic  <- mean(extrapolation_data$inb_pred2 - extrapolation_data$inb_dep1)


# MSE
MSE_linear <- mean((extrapolation_data$inb_pred-extrapolation_data$inb_dep1)^2)
MSE_quad <- mean((extrapolation_data$inb_pred1-extrapolation_data$inb_dep1)^2)
MSE_cubic <- mean((extrapolation_data$inb_pred2-extrapolation_data$inb_dep1)^2)


