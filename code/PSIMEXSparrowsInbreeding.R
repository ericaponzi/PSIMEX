##############################################################
## The simulation extrapolation technique meets ecology and evolution: 
## A general and intuitive method to account for measurement error  
##
## R-code for P-SIMEX algorithm on inbreeding depression
## of juvenile survival in song sparrows data
##############################################################



set.seed(123452)


# control function
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}


# load libraries
library(pedigree)
library(MCMCglmm)
library(parallel)
library(lme4)
library(AICcmodavg)

#Load dataset
allbirds <- read.csv("Sparrows.csv")


# calculate generation as birth year
# used later to replace fathers with random individuals from the same generation
allbirds$generation_dad <- c()
for (jj in 1:length(allbirds[ ,1])){
  
  mm <- which(allbirds$id == allbirds$social.father[jj])
  if (!is.integer0(mm)) allbirds$generation_dad[jj] <- allbirds[mm, ]$birth.year
  
  
}

allbirds <- allbirds[which(!is.na(allbirds$generation_dad)), ]
allbirds <- allbirds[which(!duplicated(allbirds$id)), ]


# true value of inbreeding depression
# from genetic pedigree
# use f.corr from data set

# store covariates of interest
trait1 <- data.frame(allbirds$id, allbirds$js, allbirds$f.corr,  allbirds$sex, allbirds$Year)
names(trait1) <- c("ID", "h", "f", "sex", "Year")


# calculate true estimate of inbreeding depression
# glm model with cloglog link
# sex and year as additional covariates

model1 <- glm(h ~ f+sex+Year,
              data  =  trait1, 
              family = binomial(link = cloglog))

# extract inbreeding depression and its SE
true <- summary(model1)$coefficients[2,1]
se_true <- summary(model1)$coefficients[2,2]


# naive value of inbreeding depression
# from social pedigree 
pedigree <- data.frame(allbirds$id, allbirds$social.mother, allbirds$social.father) 
names(pedigree) <- c("id", "dam", "sire")

# order pedigree
ord <- pedigree::orderPed(pedigree)
pedigree <- pedigree[order(ord), ]
# calculate inbreeding coefficient from pedigree
f_inb <- data.frame(pedigree$id, calcInbreeding(pedigree))
names(f_inb) <- c("ID", "f_sim")

# store covariates of interest
trait <- data.frame(allbirds$id, allbirds$js)
names(trait) <- c("ID", "h")
trait <- merge(trait, f_inb, by = "ID")
trait$f <- as.numeric(trait$f_sim)
trait <- data.frame(allbirds$id, allbirds$js, trait$f, allbirds$sex, allbirds$Year)
names(trait) <- c("ID", "h", "f", "sex", "Year")

# calculate naive estimate of inbreeding depression
# glm model with cloglog link
# sex and year as additional covariates

model <- glm(h ~ f+sex+Year, 
            data  =  trait, 
            family = binomial(link = cloglog))

# extract inbreeding depression and its SE
inb0 <- summary(model)$coefficients[2,1]
se_inb0 <-summary(model)$coefficients[2,2]
pval0 <- summary(model)$coefficients[2,4]

# store inbreeding coefficient characteristics
mean_inb0 <- mean(trait$f)
median_inb0 <- median(trait$f)
var_inb0 <- var(trait$f)


# apply PSIMEX
# increasing error proportions 
lambda0 <- seq(0.2, 0.8, 0.1)
# actual error proportion needed 
lambda <- 1- (1-lambda0)/(1-0.17) 
# total number of PSIMEX iterations
nsim <- 100

# check if parentages are all known
known <- pedigree[which(!is.na(pedigree$dam)), ]
unknown <- pedigree[which(is.na(pedigree$dam)), ]


#setting initial values
inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol  =  nsim)
se_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol  =  nsim)
mean_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)
median_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)
var_inb <- matrix(data  =  NA, nrow  =  length(lambda), ncol = nsim)

# repeat for each error proportions
for (k in 1:length(lambda)){
  # repeat for nsim simulations
  for (j in 1:nsim) {
    animal <- known$id
    dam <- known$dam
    sire <- known$sire  
    #initializing inbreeding coefficient
    trait$f_sim <- c()
    
    # choose a proportion of lambda over the total number of individuals 
    # and set their parents as unknown
    m <- sample(1:length(known$id), lambda[k]*length(known$id))
    for (i in m) {
      # sample possible replacements from the sam generation
      generation <- allbirds[which(allbirds$id == known$id[i]), ]$generation_dad
      # exclude the real father
      others <- allbirds[which(allbirds$birth.year == generation & allbirds$sex == "Male"),  ]$id
      # check if there are other possible individuals and pick one to replace 
      if (!is.integer0(others)){
        if (is.element(known[i, 3], others))  others <- others[-which(others == known[i, 3])]
        sire[i] <- others[sample(1:length(others), 1)]
      } 
    }
    
    #recalculate pedigree
    ped <- data.frame(id = animal, dam = dam, sire = sire)
    ped <- rbind(ped, unknown)
    ord <- pedigree::orderPed(ped)
    ped <- ped[order(ord),]
    
    #re calculate inbreeding coefficient
    f_sim <- data.frame(ped$id, calcInbreeding(ped))
    names(f_sim) <- c("ID", "f_sim")
    # include it in the data
    trait <- merge(trait, f_sim, by = "ID")
    trait$f_sim <- as.numeric(trait$f_sim)
    
    # calculate inbreeding depression 
    if(sum(trait$f_sim)>0) {
      
      # glm for survival 
      # cloglog link function 
      # sex and year as covariates
      model <- glm(h ~ f_sim+sex+Year , 
                   data  =  trait, 
                   family = binomial(link = cloglog))
      
      # extract inbreeding depression and its SE
      inb[k, j] <- summary(model)$coefficients[2,1]
      se_inb[k,j] <- summary(model)$coefficients[2,2]
      
      # store inbreeding coefficient characteristics
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
    }
    
    
  }
  
}


#compute means across simulations
inb_ave <- rowMeans(inb, na.rm = TRUE)
se_inb_ave <- rowMeans(se_inb, na.rm = TRUE)
mean_inb_ave <- rowMeans(mean_inb, na.rm = TRUE)
median_inb_ave <- rowMeans(median_inb, na.rm = TRUE)
var_inb_ave <- rowMeans(var_inb, na.rm = TRUE)

# merge all values with naive ones
lambda <- c(0.17,lambda0)
inb_ave <- c(inb0, inb_ave)
se_inb <- c(se_inb0, se_inb_ave)
mean_inb <- c(mean_inb0,mean_inb_ave)
median_inb <- c(median_inb0, median_inb_ave)
var_inb <- c(var_inb0, var_inb_ave)



# Extrapolation phase of PSIMEX
p.names <- c("inb", "se_inb")
estimates <- data.frame(inb_ave, se_inb)

colnames(estimates) <- p.names

#linear case
extrapolation_inb <- lm(estimates[ ,1] ~ lambda)
inb_pred1 <- predict(extrapolation_inb, newdata = data.frame(lambda = 0))

#quadratic case
extrapolation_inb1 <- lm(estimates[ ,1] ~ lambda+ I(lambda^2))
inb_pred2 <- predict(extrapolation_inb1, newdata = data.frame(lambda = 0))

#cubic case
extrapolation_inb2 <- lm(estimates[ ,1] ~ lambda+ I(lambda^2)+ I(lambda^3))
inb_pred3 <- predict(extrapolation_inb2, newdata = data.frame(lambda = 0))

AICc(extrapolation_inb)

AICc(extrapolation_inb2)

AICc(extrapolation_inb1)


# calculate SIMEX SE


# first component is the sampling variability

S <- c()
#one per lambda 
for ( i in 1:length(lambda0)) {
  diff <- c()
  #calculate the differences per each simulation
  for ( j in 1: nsim ) {
    diff[j] <- inb[i,j]-inb_ave[i]
  }
  
  diff <- na.omit(diff)
  S[i] <- 1/(nsim-1)*sum(diff^2) 
  
}
S <- c(0, S)

# calculate the total standard error
# second component is the se from regressions
S1 <- se_inb
S2 <- S
Stot <- S1 - S2

# extrapolated value
# linear case
extrapolation_var<- lm(Stot ~ lambda)
var_pred<-predict(extrapolation_var, newdata = data.frame(lambda = 0))

#quadratic case
extrapolation_var1<- lm(Stot ~ lambda+ I(lambda^2))
var_pred1<-predict(extrapolation_var1, newdata = data.frame(lambda = 0))

#cubic case
extrapolation_var2<- lm(Stot ~ lambda+ I(lambda^2)+ I(lambda^3))
var_pred2<-predict(extrapolation_var2, newdata = data.frame(lambda = 0))


# bias
bias_naive <- inb0 - true
bias_linear <- inb_pred1 - true
bias_quad <- inb_pred2 - true
bias_cubic <- inb_pred3 - true

# MSE
MSE_naive <- bias_naive^2 + se_inb0^2
MSE_linear <- bias_linear^2 + var_pred^2
MSE_quad <- bias_quad^2 + var_pred1^2
MSE_linear <- bias_cubic^2 + var_pred2^2


