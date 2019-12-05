#################################
## smooth.FEM.basis meshC 2.5D ##
#################################

library(fdaPDE)

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

GCVFLAG=TRUE # =FALSE do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

GCVMETHODFLAG="Exact" #for exact GCV
#GCVMETHODFLAG="Stochastic" #for stochastic GCV (default)

set.seed(5847947)

# build the mesh using create.MESH.2.5D

nome_mesh="C2_5D"

vertici <- read.table(paste0("Data/2.5D/",nome_mesh,"_vertici.txt"), quote="\"", comment.char="")
triangoli <- read.table(paste0("Data/2.5D/",nome_mesh,"_triangoli.txt"), quote="\"", comment.char="")

mesh <- fdaPDE::create.MESH.2.5D(nodes = vertici[,1:3],triangles = triangoli[,1:3])
#mesh=second.order.MESH.2.5D(mesh)
FEMbasis <- fdaPDE::create.FEM.basis(mesh)

plot(mesh)

nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +  a2* sin(2*pi*mesh$nodes[i+1,2]) +  a3*sin(2*pi*mesh$nodes[i+1,3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)

lambda=c(0.00375)
# output_CPP =smooth.FEM.basis(observations = data,
#                              FEMbasis = FEMbasis, lambda = lambda,
#                              CPP_CODE = TRUE)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda, GCV = GCVFLAG, 
                             GCVmethod = GCVMETHODFLAG)

cat("Showing result")
plot(output_CPP$fit.FEM)

# exact GCV  0.8277302
#       stderr 0.6194387
#       edf  179.7067
#       

# stoch GCV  0.8242321
#       stderr 0.6187832 
#       edf 179.3775 
#  

points=eval.FEM(output_CPP$fit.FEM, locations=mesh$nodes) # il primo è un NA
write.table(points, file="smoothFEMbasis25D_meshC_nod_nocov.txt")

### covariates, observations at node locations ###
nodesLocations=mesh$nodes
cov1=3*nodesLocations[,1]+2*nodesLocations[,2]+5*nodesLocations[,3]+rnorm(nnodes,mean=0,sd=0.1)
cov2=rnorm(nnodes, mean=3, sd=1)

# output_CPP2=smooth.FEM.basis(observations = data, covariates = cbind(cov1,cov2), FEMbasis = FEMbasis, lambda =lambda)

output_CPP2=smooth.FEM.basis(observations = data, covariates = cbind(cov1,cov2), 
                             FEMbasis = FEMbasis, lambda =lambda, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)

plot(output_CPP2$fit.FEM)

# exact GCV  0.8333681
#       stderr  0.6198072
#       edf  180.5737 

# exact GCV 0.8599396
#       stderr 0.6246897 
#       edf 182.9783 

# output_CPP2$beta
# [1,] -0.18535120
# [2,]  0.03038352

points2=eval.FEM(output_CPP2$fit.FEM, locations=nodesLocations) # il primo è un NA
write.table(points2, file="smoothFEMbasis25D_meshC_nod_cov.txt")
