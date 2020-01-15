###########################
## test smooth.FEM.basis ##
###########################

library(fdaPDE)
setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

GCVFLAG=FALSE  # do not compute GCV (default)
#GCVFLAG=TRUE

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

#GCVMETHODFLAG="Exact" #for exact GCV
GCVMETHODFLAG="Stochastic" #for stochastic GCV (default)

# load the C-shaped mesh
load("Data/2D/meshC.RData")
FEMbasis = create.FEM.basis(mesh)

# exact solution f: Ramsay's test function

fs.test <- function (x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = TRUE) 
{
  
  q <- pi * r/2
  a <- d <- x * 0
  
  ind <- x >= 0 & y > 0
  a[ind] <- q + x[ind]
  d[ind] <- y[ind] - r
  
  ind <- x >= 0 & y <= 0
  a[ind] <- (-q - x[ind])
  d[ind] <- -r - y[ind]
  
  ind <- x < 0
  a[ind] <- -atan(y[ind]/x[ind]) * r 
  d[ind] <-( sqrt(x[ind]^2 + y[ind]^2) - r )* (y[ind]/r0*(as.numeric(abs(y[ind])<=r0 & x[ind]>-0.5))+(as.numeric(abs(y[ind])>r0 || x[ind]<(-0.5))))
  
  ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - r0)^2)
  
  f <- a * b + d^2
  
  if (exclude) 
    f[ind] <- NA
  
  attr(f, "exclude") <- ind
  f
}

### locations at mesh nodes

# generate exact data from the test function
dati_esatti = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE)

# add gaussian noise to data

sd = 5*abs(max(dati_esatti) - min(dati_esatti))/100

set.seed(5847947) # MEMO: do not forget the seed to replicate the simulation!

data = dati_esatti + rnorm(length(dati_esatti), sd = sd)

# choose value of lambda parameter
lambda = 10^-2

# solution
FEM_CPP = smooth.FEM.basis(observations = data, FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
#FEM_CPP = smooth.FEM.basis(observations = data, FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG)

plot(FEM_CPP$fit.FEM)

points=eval.FEM(FEM_CPP$fit.FEM, locations = mesh$nodes)
#write.table(points,file="smoothFEMbasis2D_nod_nocov.txt")
write.table(points,file="smoothFEMbasis2D_nod_nocov_exactGCV.txt")

#exact GCV: 0.2084202
#      edf: 58.4193
#     stderr:0.393723

#stoch GCV: 0.2103883
#      edf: 59.21434
#     stderr:0.3946491

### locations different from nodes

set.seed(5847947) 

loc1=cbind(runif(50,min=0.5,max=2.5),runif(50,min=0.1,max=0.9)) #braccio superiore
loc2=cbind(runif(50,min=0.5,max=2.5),runif(50,min=-0.9,max=-0.1)) # braccio inferiore
loc3=cbind(runif(50,min=-0.7,max=-0.1),runif(50,min=-0.5,max=0.5)) #circonferenza grande
noditest=mesh$nodes[1:50,]# alcune oss coincidenti con nodi

oss=rbind(loc1,loc2,loc3,noditest)

dati_esatti2 = fs.test(oss[,1], oss[,2], exclude = FALSE)

sd2 = 5*abs(max(dati_esatti2) - min(dati_esatti2))/100

data2 = dati_esatti2 + rnorm(length(dati_esatti2), sd = sd2)

lambda = 10^-2

#FEM_CPP2 = smooth.FEM.basis(locations=oss,observations = data2, FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG)
FEM_CPP2 = smooth.FEM.basis(locations=oss,observations = data2, FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(FEM_CPP2$fit.FEM)

points2=eval.FEM(FEM_CPP2$fit.FEM,locations=mesh$nodes)
#write.table(points2,"smoothFEMbasis2D_nonod_nocov.txt")
write.table(points2,"smoothFEMbasis2D_nonod_nocov_exactGCV.txt")


#exact GCV: 0.1284621
#      edf: 38.88559
#     stderr:0.3216916

#stoch GCV: 0.1280843
#      edf: 38.64815
#     stderr:0.3214548

### locations at nodes, with covariates

cov=cbind(rnorm(dim(mesh$nodes)[1],mean=0.5,sd=0.01), rnorm(dim(mesh$nodes)[1], mean=2.5,sd=0.05))
#FEM_CPP3= smooth.FEM.basis(observations = data, covariates=cov,FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG)
FEM_CPP3= smooth.FEM.basis(observations = data, covariates=cov,FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
plot(FEM_CPP3$fit.FEM)
FEM_CPP3$beta
#[1,] -0.7866522  coincidono nei due gcv e senza
#[2,] 0.1198352

points3=eval.FEM(FEM_CPP3$fit.FEM,locations=mesh$nodes)
#write.table(points3,"smoothFEMbasis2D_nod_cov.txt")
write.table(points3,"smoothFEMbasis2D_nod_cov_exactGCV.txt")
#write.table(solution$beta)

#exact GCV: 0.2126453
#      edf: 60.15952
#     stderr:0.3956479

#stoch GCV: 0.2061425
#      edf: 57.53278
#     stderr:0.3925878

### locations different from nodes, with covariates

cov=cbind(rnorm(dim(oss)[1],mean=0.5,sd=0.01), rnorm(dim(oss)[1], mean=2.5,sd=0.05))
#FEM_CPP4 = smooth.FEM.basis(locations=oss,observations = data2, covariates=cov,FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG)
FEM_CPP4 = smooth.FEM.basis(locations=oss,observations = data2, covariates=cov,FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
plot(FEM_CPP4$fit.FEM)

FEM_CPP4$beta
# [1,] -3.6750712 coincidono nei due GCV e senza
# [2,] 0.2618758

points4=eval.FEM(FEM_CPP4$fit.FEM,locations=oss)
#write.table(points4,"smoothFEMbasis2D_nonod_cov.txt")
write.table(points4,"smoothFEMbasis2D_nonod_cov_exactGCV.txt")

#exact GCV: 0.1295854
#      edf: 40.74686
#     stderr:0.3212233

#stoch GCV: 0.1187033
#      edf: 33.60719
#     stderr:0.3142561