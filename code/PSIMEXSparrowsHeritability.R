library("pedigree")
library("MCMCglmm")
library("mvtnorm")
library("GeneticsPed") 
library("pedigreemm")
library("AnimalINLA")
library("INLA")
library("Matrix")


set.seed(123452)

# control function
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

# load data
allbirds <- read.csv("SongSparrowData_noembr.csv", sep="\t")

# calculate generation as birth year
allbirds$generation_dad <- c()
for (jj in 1:length(allbirds[ ,1])){
  
  mm <- which(allbirds$id == allbirds$social.father[jj])
  if (is.integer0(mm)) allbirds$generation_dad[jj] <- NA
  if (!is.integer0(mm)) allbirds$generation_dad[jj] <- allbirds[mm, ]$birth.year
  
  
}



morpho <- read.csv("MandarteBiosBasicData2007.csv")





w <- intersect(allbirds$id, morpho$ninecode)
k <- which(is.element(morpho$ninecode, w))
morpho <- morpho[k, ]
morpho$sex<-c()
for ( i in 1:length(morpho[ ,1])) {
  w<-which(allbirds$id==morpho$ninecode[i])
  morpho$sex[i]<-allbirds$sex[w]
  
  
}

dati <- morpho[!is.na(morpho$sex), ]

# extract variables of interest
trait <-data.frame(dati$ninecode, dati$mintars, dati$sex)
names(trait)<-c("animal", "Y", "sex")

# calculate true value
# from genetic pedigree
all.ids <- c(allbirds$id, allbirds$gen.father, allbirds$gen.mother)
all.ids <- all.ids[!duplicated(all.ids)]
all <- which(!is.na(all.ids))
gen.dad <- c()
gen.mom <- c()
idg <- c()

for (i in all) {
  w <- which(allbirds$id == all.ids[i])
  idg[i] <- allbirds$id[w]
  gen.dad[i] <- allbirds$father[w]
  gen.mom[i] <- allbirds$mother[w]
  
}


# store genetic pedigree
gen <- data.frame(id = idg, dam = gen.mom, sire = gen.dad)

###########################################################################
### New: I use the MasterBayes package and a much more efficient version of inversA():
###########################################################################
# use orderPed from MasterBayes:
library(MasterBayes)
gen <- orderPed(gen)

# Replace the ringnumber by an id (1:nanimals)
d.map <- data.frame(id=gen$id,id_new=1:(nrow(gen)))

gen$id_new <- 1:(nrow(gen))
gen$dam_new <- d.map[match(gen$dam,gen$id),"id_new"]
gen$sire_new <- d.map[match(gen$sire,gen$id),"id_new"]

ainvOut <- inverseA(gen[, 4:6])
Cmatrix <- ainvOut$Ainv

# inbreeding was also directly calculated with inversA()
f <- data.frame(gen$id_new ,ainvOut$inbreeding)
names(f) <- c("id_new","f")


# Generate a new id, running from 1:no_of_animals (directly stor id_new and IndexA and IndexA.2 for later use in INLA)
trait$IndexA <- trait$IndexA.2 <- trait$id_new <- d.map[match(trait$animal,d.map$id),"id_new"]

trait$f <- f[match(trait$id_new,f$id_new),"f"] #merge(trait,f,by="id_new")
###########################################################################

p.var <- var(trait$Y, na.rm = TRUE)


#fit the inla model
formula = Y ~ sex + f +
  f(IndexA, model = "generic0", Cmatrix = Cmatrix, 
    hyper=list(prec=list(param = c(1/2,p.var/6))),
    constr=T) +
  f(IndexA.2, model = "iid",
    hyper=list(prec=list(param = c(1/2,p.var/6))), 
    constr = TRUE)

model = inla(formula=formula, family="gaussian",
             data =  trait, 
             control.family=list(hyper = list(prec =
                                                list(param = c(1/2, p.var/6), fixed = FALSE))),
             only.hyperparam =FALSE,control.compute=list(dic=T))
summary(model)
 


posterior.sample.var <- inla.hyperpar.sample(100000,model)
h2 <- (1/posterior.sample.var[,2]) /( 1/posterior.sample.var[,1] + 1/posterior.sample.var[,2] + 1/posterior.sample.var[,3])

h2 <- as.numeric(h2)
mean(h2)
(h.true <- posterior.mode(as.mcmc(h2)))
sd.true <- sd(h2)
low.true <- HPDinterval(as.mcmc(h2)) 
up.true <- HPDinterval(as.mcmc(h2))


# compute naive heritability with INLA


#calculate social inbreeding coefficient
all.ids1<-c(allbirds$id, allbirds$social.father, allbirds$social.mother)
all.ids1<-all.ids1[!duplicated(all.ids1)]
all1<-which(!is.na(all.ids1))
soc.dad<-c()
soc.mom<-c()
idg<-c()

for (i in all1) {
  w<-which(allbirds$id==all.ids1[i])
  idg[i]<-allbirds$id[w]
  soc.dad[i]<-allbirds$social.father[w]
  soc.mom[i]<-allbirds$social.mother[w]
  
}

soc<-data.frame(id=idg, dam=soc.mom, sire=soc.dad)
###########################################################################
### New: I use the MasterBayes package and a much more efficient version of inversA():
###########################################################################
# use orderPed from MasterBayes:
soc <- orderPed(soc)

# Replace the ringnumber by an id (1:nanimals)
d.map <- data.frame(id=soc$id,id_new=1:(nrow(soc)))

soc$id_new <- 1:(nrow(soc))
soc$dam_new <- d.map[match(soc$dam,soc$id),"id_new"]
soc$sire_new <- d.map[match(soc$sire,soc$id),"id_new"]

ainvOut <- inverseA(soc[4:6])
Cmatrix <- ainvOut$Ainv

# inbreeding was also directly calculated with inversA()
f <- data.frame(soc$id_new ,ainvOut$inbreeding)
names(f) <- c("id_new","f")

# Generate a new id, running from 1:no_of_animals (directly stor id_new and IndexA and IndexA.2 for later use in INLA)
trait$IndexA <- trait$IndexA.2 <- trait$id_new <- d.map[match(trait$animal,d.map$id),"id_new"]

trait$f <- f[match(trait$id_new,f$id_new),"f"] #merge(trait,f,by="id_new")
###########################################################################


p.var <- var(trait$Y, na.rm = TRUE)


#fit the inla model
formula.naive = Y ~ sex + f +
  f(IndexA, model = "generic0", Cmatrix = Cmatrix, 
    hyper=list(prec=list(param = c(1/2,p.var/6))),
    constr=T) +
  f(IndexA.2, model = "iid",
    hyper=list(prec=list(param = c(1/2,p.var/6))), 
    constr = TRUE)

model.naive = inla(formula=formula.naive, family="gaussian",
             data =  trait, 
             control.family=list(hyper = list(prec =
                                                list(param = c(1/2, p.var/6), fixed = FALSE))),
             only.hyperparam =FALSE,control.compute=list(dic=T))



posterior.sample.var.naive <- inla.hyperpar.sample(100000,model.naive)
h2 <- (1/posterior.sample.var.naive[,2]) /( 1/posterior.sample.var.naive[,1] + 1/posterior.sample.var.naive[,2] + 1/posterior.sample.var.naive[,3])


h2 <- as.numeric(h2)
mean(h2)
(h.naive <- posterior.mode(as.mcmc(h2)))
sd.naive <- sd(h2)
low.naive <- HPDinterval(as.mcmc(h2)) 
up.naive <- HPDinterval(as.mcmc(h2))

# apply PSIMEX
# increasing error proportions 
lambda0 <- seq(0.2, 0.8, 0.1)
# actual error proportion needed 
lambda <- 1- (1-lambda0)/(1-0.17)
# total number of PSIMEX iteration

nsim <- 100


  h<-matrix(data=NA, nrow=length(lambda), ncol=nsim)
  sd<-matrix(data=NA, nrow=length(lambda), ncol=nsim)
  low<-matrix(data=NA, nrow=length(lambda), ncol=nsim)
  up<-matrix(data=NA, nrow=length(lambda), ncol=nsim)
  
  
  for (k in 1:length(lambda)){
    for (j in 1:nsim) {
      all.ids <- c(allbirds$id, allbirds$social.father, allbirds$social.mother)
      all.ids <- all.ids[!duplicated(all.ids)]
      all <- which(!is.na(all.ids))
      
      
      dad<-c()
      mom<-c()
      id<-c()
      
      for (i in all) {
        w<-which(allbirds$id==all.ids[i])
        id[i] <- allbirds$id[w]
        dad[i] <- allbirds$social.father[w]
        mom[i] <- allbirds$social.mother[w]
        
      }
      
      m <- sample(1:length(allbirds$id), lambda[k]*length(allbirds$id))
      for (i in m) {
        generation <- allbirds[i, ]$generation_dad
        others <- allbirds[which(allbirds$birth.year == generation & allbirds$sex == "Male"),  ]$id
        if (!is.integer0(others)){
          if (is.element(dad[i], others))  others <- others[-which(others == dad[i])]
          dad[i] <- others[sample(1:length(others), 1)]
          
        }
      }
      pedsim <- data.frame(id = id, dam = mom, sire = dad)
      
      # ord <- orderPed(pedsim)
      # pedsim <- pedsim[order(ord),]
      # 
      # f <- data.frame(pedsim$id, calcInbreeding(pedsim))
      # names(f) <- c("Individual", "f")
      # 
      # trait <-data.frame(dati$ninecode, dati$mintars, dati$sex)
      # names(trait)<-c("animal", "Y", "sex")
      # 
      # 
      # trait$Individual <- trait$animal
      # trait <- merge(trait, f, by = "Individual")
      # 
      # 
      # 
      # 
      # pedigree <- data.frame(as.numeric(pedsim$id), as.numeric(pedsim$sire), as.numeric(pedsim$dam))
      # for (i in 1:length(pedigree[,1])) {
      #   if (is.na(pedigree[i,1])) pedigree[i,1]<-0
      #   if (is.na(pedigree[i,2])) pedigree[i,2]<-0
      #   if (is.na(pedigree[i,3])) pedigree[i,3]<-0
      # }
      # names(pedigree) <- c('Individual', 'Parent1', 'Parent2')
      # xx = compute.Ainverse(pedigree)
      # Ainv = xx$Ainverse
      # map  = xx$map
      # Cmatrix = sparseMatrix(i=Ainv[,1],j=Ainv[,2],x=Ainv[,3])
      # 
      # #ID
      # Ndata = dim(trait)[1]
      # trait$ID <- trait$Individual
      # trait$IndexA = rep(0,Ndata)
      # trait$IndexA.2 = rep(0,Ndata)
      # for(i in 1:Ndata)    trait$IndexA[i] = which(map[,1]==trait$ID[i])
      # for(i in 1:Ndata)    trait$IndexA.2[i] = which(map[,1]==trait$ID[i])
      ###########################################################################
      ### New: I use the MasterBayes package and a much more efficient version of inversA():
      ###########################################################################
      # use orderPed from MasterBayes:
      pedsim <- orderPed(pedsim)
      
      # Replace the ringnumber by an id (1:nanimals)
      d.map <- data.frame(id=pedsim$id,id_new=1:(nrow(pedsim)))
      
      pedsim$id_new <- 1:(nrow(pedsim))
      pedsim$dam_new <- d.map[match(pedsim$dam,pedsim$id),"id_new"]
      pedsim$sire_new <- d.map[match(pedsim$sire,pedsim$id),"id_new"]
      
      ainvOut <- inverseA(pedsim[4:6])
      Cmatrix <- ainvOut$Ainv
      
      # inbreeding was also directly calculated with inversA()
      f <- data.frame(pedsim$id_new ,ainvOut$inbreeding)
      names(f) <- c("id_new","f")
      
      # Generate a new id, running from 1:no_of_animals (directly stor id_new and IndexA and IndexA.2 for later use in INLA)
      trait$IndexA <- trait$IndexA.2 <- trait$id_new <- d.map[match(trait$animal,d.map$id),"id_new"]
      
      trait$f <- f[match(trait$id_new,f$id_new),"f"] #merge(trait,f,by="id_new")
      ###########################################################################
      
      p.var <- var(trait$Y, na.rm = TRUE)
      
      
      #fit the inla model
      formula.sim = Y ~ sex + f +
        f(IndexA, model = "generic0", Cmatrix = Cmatrix, 
          hyper=list(prec=list(param = c(1/2,p.var/6))),
          constr=T) +
        f(IndexA.2, model = "iid",
          hyper=list(prec=list(param = c(1/2,p.var/6))), 
          constr = TRUE)
      
      model.sim = inla(formula=formula.sim, family="gaussian",
                         data =  trait, 
                         control.family=list(hyper = list(prec =
                                                            list(param = c(1/2, p.var/6), fixed = FALSE))),
                         only.hyperparam =FALSE,control.compute=list(dic=T))
      
      
      
      posterior.sample.var.sim <- inla.hyperpar.sample(100000,model.sim)
      h2 <- (1/posterior.sample.var.sim[,2]) /( 1/posterior.sample.var.sim[,1] + 1/posterior.sample.var.sim[,2] + 1/posterior.sample.var.sim[,3])
      
      
      h2 <- as.numeric(h2)
      mean(h2)
      h[k,j] <- posterior.mode(as.mcmc(h2))
      sd[k,j] <- as.numeric(sd(h2))
      low[k,j] <- HPDinterval(as.mcmc(h2))[1]
      up[k,j] <- HPDinterval(as.mcmc(h2))[2]
      
      
    }
  }


write.table(h, 'h_inla2005.txt')
write.table(sd, 'sd_inla2005.txt')
write.table(low, 'low_inla2005.txt')
write.table(up, 'up_inla2005.txt')


# compute means across simulations
h_ave <- rowMeans(h)
sd_h <- rowMeans(sd)

# merge all values with naive ones
lambda <- c(0.17, lambda0)
h_ave <- c(h.naive, h_ave)
se_h <- c(sd.naive, sd_h)

# Extrapolation phase of PSIMEX
p.names <- c("h",  "se")
estimates <- data.frame(h_ave,  se_h)
colnames(estimates) <- p.names

#linear case

extrapolation_h<- lm(estimates[ ,1] ~ lambda)
(h_pred<-predict(extrapolation_h, newdata = data.frame(lambda = 0)))


#quadratic case

extrapolation_h1<- lm(estimates[ ,1] ~ lambda+ I(lambda^2))
(h_pred1<-predict(extrapolation_h1, newdata = data.frame(lambda = 0)))


#cubic case

extrapolation_h2<- lm(estimates[ ,1] ~ lambda+ I(lambda^2)+ I(lambda^3))
(h_pred2<-predict(extrapolation_h2, newdata = data.frame(lambda = 0)))


# calculate SIMEX SE

# first component 
# sampling variability
S <- c()
#one per lambda 
for ( i in 1:length(lambda0)) {
  diff <- c()
  #calculate the differences per each simulation
  for ( j in 1: nsim ) {
    diff[j] <- h[i,j]-h_ave[i]
  }
  diff <- na.omit(diff)
  S[i] <- 1/(nsim-1)*sum(diff^2) 
  
}


S <- c(0, S)

S1 <- se_h
S2 <- S
Stot <- S1 - S2

# extrapolated value
# linear case
extrapolation_var<- lm(Stot ~ lambda)
(var_pred<-predict(extrapolation_var, newdata = data.frame(lambda = 0)))

#quadratic case
extrapolation_var1<- lm(Stot ~ lambda+ I(lambda^2))
(var_pred1<-predict(extrapolation_var1, newdata = data.frame(lambda = 0)))

#cubic case
extrapolation_var2<- lm(Stot ~ lambda+ I(lambda^2)+ I(lambda^3))
(var_pred2<-predict(extrapolation_var2, newdata = data.frame(lambda = 0)))

library(AICcmodavg)
AICc(extrapolation_h)

AICc(extrapolation_h2)

AICc(extrapolation_h1)




# bias
bias_naive <- h.naive - h.true
bias_linear <- h_pred - h.true
bias_quad <- h_pred1 - h.true
bias_cubic <- h_pred2 - h.true

# MSE
MSE_naive <- bias_naive^2 + sd.naive^2
MSE_linear <- bias_linear^2 + var_pred^2
MSE_quad <- bias_quad^2 + var_pred1^2
MSE_linear <- bias_cubic^2 + var_pred2^2




