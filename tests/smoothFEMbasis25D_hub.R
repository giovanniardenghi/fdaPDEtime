################################
## test smooth.FEM.basis 2.5D ##
################################

library(fdaPDE)

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

GCVFLAG=FALSE  # =FALSE do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

#GCVMETHODFLAG=1 #for exact GCV
GCVMETHODFLAG=2 #for stochastic GCV (default)

#set.seed(5847947)

data(hub)

cat('Plotting the mesh \n')
plot(hub)

### no covariates, observations at node locations ###

# Generate some random data

nnodes = hub$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*hub$nodes[3*i+1]) +  a2* sin(2*pi*hub$nodes[3*i+2]) +  a3*sin(2*pi*hub$nodes[3*i+3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
FEMbasis <- create.FEM.basis(hub)

lambda=c(0.00375)
# output_CPP =smooth.FEM.basis(observations = data,
#                              FEMbasis = FEMbasis, lambda = lambda,
#                              CPP_CODE = TRUE)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda,
                             CPP_CODE = TRUE, GCV = GCVFLAG)

cat("Showing result")
plot(output_CPP$fit.FEM)

# exact GCV 0.3584893
#       stderr 0.5348416
#       edf 68.69792
#       

# stoch GCV 0.3580569
#       stderr 0.5346803
#       edf 68.53417
#  

points=eval.FEM(output_CPP$fit.FEM, locations=hub$nodes) # non va bene!

nodesLocations=matrix(data=hub$nodes, ncol=3, nrow=nnodes, byrow=T) #sarà giusto? o devo distribuire per colonne?

points=eval.FEM(output_CPP$fit.FEM, locations=nodesLocations) # il primo è un NA
write.table(points, file="smoothFEMbasis25D_nod_nocov_stochGCV.txt")

### covariates, observations at node locations ###

cov1=3*nodesLocations[,1]+2*nodesLocations[,2]+5*nodesLocations[,3]+rnorm(nnodes,mean=0,sd=0.1)
cov2=rnorm(nnodes, mean=3, sd=1)

# output_CPP2=smooth.FEM.basis(observations = data, covariates = cbind(cov1,cov2), FEMbasis = FEMbasis, lambda =lambda)

output_CPP2=smooth.FEM.basis(observations = data, covariates = cbind(cov1,cov2), 
                             FEMbasis = FEMbasis, lambda =lambda, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)

plot(output_CPP2$fit.FEM)

# exact GCV 0.3558301
#       stderr 0.5314574
#       edf  70.11852

# exact GCV 0.3581721
#       stderr 0.5323297
#       edf  71.00231

points2=eval.FEM(output_CPP2$fit.FEM, locations=nodesLocations) # il primo è un NA
write.table(points2, file="smoothFEMbasis25D_nod_cov_stochGCV.txt")
