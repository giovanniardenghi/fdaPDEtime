####################################
## test smooth.FEM.basis 3D meshC ##
####################################

library(fdaPDE)

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

GCVFLAG=FALSE  # =FALSE do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

GCVMETHODFLAG=1 #for exact GCV
#GCVMETHODFLAG=2 #for stochastic GCV (default)

set.seed(5847947)

# create the mesh using create.mesh.3D

nome_mesh = "meshCcicciona"

vertici <- read.table(paste0("Data/3D/",nome_mesh,"_vertici.txt"), quote="\"", comment.char="")
tetraedri <- read.table(paste0("Data/3D/",nome_mesh,"_tetraedri.txt"), quote="\"", comment.char="")

mesh <- fdaPDE::create.MESH.3D(nodes = vertici[,1:3],tetrahedrons = tetraedri[,1:4])
FEMbasis <- fdaPDE::create.FEM.basis(mesh)

plot(mesh)

### no covariates, observations at node locations ###

nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +  a2* sin(2*pi*mesh$nodes[i+1,2]) +  a3*sin(2*pi*mesh$nodes[i+1,3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
FEMbasis <- create.FEM.basis(mesh)

lambda=c(10^-2)
# output_CPP =smooth.FEM.basis(observations = data,
#                              FEMbasis = FEMbasis, lambda = lambda,
#                              CPP_CODE = TRUE)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda,
                             CPP_CODE = TRUE,GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
plot(output_CPP$fit.FEM)

# exact GCV  3.848495
#       stderr  1.693281
#       edf  608.1299 

#stoch GCV  3.854978 
#      stderr  1.693994
#      edf  609.6247 

nodesLocations=matrix(data=mesh$nodes, ncol=3, nrow=nnodes, byrow=T) #sarà giusto? o devo distribuire per colonne?

points=eval.FEM(output_CPP$fit.FEM, locations=nodesLocations)

write.table(points, file="smoothFEMbasis3D_meshC_nod_nocov.txt")

### no covariates, observations different from node locations ###

loc1=cbind(runif(50,min=1,max=5), runif(50, min=-0.5,max=0.5),runif(50,min=1,max=1.5)) #braccio superiore
loc2=cbind(runif(50,min=1,max=5),runif(50,min=-0.5,max=0.5),runif(50,min=-1.5,max=-1)) # braccio inferiore
loc3=cbind(runif(50,min=-0.5,max=0.5),runif(50,min=-0.5,max=0.5),runif(50,min=-0.5,max=0.5)) #circonferenza
loc=rbind(loc1,loc2,loc3)

func_evaluation2=numeric(150)

for (i in 0:(150-1)){
  func_evaluation2[i+1] = a1* sin(2*pi*loc[3*i+1]) +  a2* sin(2*pi*loc[3*i+2]) +  a3*sin(2*pi*loc[3*i+3]) +1
}

data2=func_evaluation2+rnorm(150, mean=0, sd=0.5)

# output_CPP2 =smooth.FEM.basis(observations = data2,locations=loc,
#                              FEMbasis = FEMbasis, lambda = lambda,
#                              CPP_CODE = TRUE)
output_CPP2 =smooth.FEM.basis(observations = data2,locations=loc,
                              FEMbasis = FEMbasis, lambda = lambda,
                              CPP_CODE = TRUE,GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
plot(output_CPP2$fit.FEM)

#exact GCV 7.284317
#      stderr 2.363542 
#      edf  34.96529 

#stoch GCV  7.229292 
#      stderr 2.359065 
#      edf 34.52833  

points2=eval.FEM(output_CPP2$fit.FEM, locations=nodesLocations)
write.table(points2, file="smoothFEMbasis3D_meshC_nonod_nocov.txt")

### covariates, observations at node locations ###

cov1=3*nodesLocations[,1]+2*nodesLocations[,2]+5*nodesLocations[,3]+rnorm(nnodes,mean=0,sd=0.1)
cov2=rnorm(nnodes, mean=3, sd=1)

# output_CPP3 =smooth.FEM.basis(observations = data,covariates=cbind(cov1,cov2),
#                               FEMbasis = FEMbasis, lambda = lambda,
#                               CPP_CODE = TRUE)
output_CPP3 =smooth.FEM.basis(observations = data,covariates=cbind(cov1,cov2),
                              FEMbasis = FEMbasis, lambda = lambda,
                              CPP_CODE = TRUE, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
plot(output_CPP3$fit.FEM)

output_CPP3$beta
# [1,]  -0.25399166
# [2,]   0.05041702


#exact GCV   3.854487
#      stderr 1.693951 
#      edf  609.489 

#stoch GCV  3.867213 
#      stderr  1.695347 
#      edf  612.4129 

points3=eval.FEM(output_CPP3$fit.FEM, locations=nodesLocations)
write.table(points3, file="smoothFEMbasis3D_meshC_nod_cov.txt")

### covariates, observations different from node locations ###

cov1=3*loc[,1]+2*loc[,2]+5*loc[,3]+rnorm(150,mean=0,sd=0.1)
cov2=rnorm(150, mean=3, sd=1)

# output_CPP4 =smooth.FEM.basis(observations = data2,locations=loc,covariates=cbind(cov1,cov2),
#                               FEMbasis = FEMbasis, lambda = lambda,
#                               CPP_CODE = TRUE)
output_CPP4 =smooth.FEM.basis(observations = data2,locations=loc,covariates=cbind(cov1,cov2),
                              FEMbasis = FEMbasis, lambda = lambda,
                              CPP_CODE = TRUE,GCV = GCVFLAG,GCVmethod = GCVMETHODFLAG)
plot(output_CPP4$fit.FEM)

output_CPP4$beta
# [1,]  0.2298134
# [2,] 0.4508698 

#exact GCV 7.212591  
#      stderr 2.33946
#      edf  36.17667 

#stoch GCV  7.25451 
#      stderr 2.342852 
#      edf  36.506 

points4=eval.FEM(output_CPP4$fit.FEM, locations=nodesLocations)
write.table(points4, file="smoothFEMbasis3D_meshC_nonod_cov.txt")
