##############################
## test smooth.FEM.basis 3D ##
##############################

library(fdaPDE)

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

GCVFLAG=TRUE  # =FALSE do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

#GCVMETHODFLAG=1 #for exact GCV
GCVMETHODFLAG=2 #for stochastic GCV (default)

set.seed(5847947)

data("sphere3D")

cat('Plotting the mesh \n')
plot(sphere3D)

### no covariates, observations at node locations ###

nnodes = sphere3D$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*sphere3D$nodes[3*i+1]) +  a2* sin(2*pi*sphere3D$nodes[3*i+2]) +  a3*sin(2*pi*sphere3D$nodes[3*i+3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
FEMbasis <- create.FEM.basis(sphere3D)

lambda=c(10^-2)
# output_CPP =smooth.FEM.basis(observations = data,
#                              FEMbasis = FEMbasis, lambda = lambda,
#                              CPP_CODE = TRUE)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda,
                             CPP_CODE = TRUE,GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
plot(output_CPP$fit.FEM)

# exact GCV 1.089359
#       stderr 0.8726719
#       edf  176.6362

#stoch GCV  1.086962
#      stderr 0.8721914
#      edf  176.1839

nodesLocations=matrix(data=sphere3D$nodes, ncol=3, nrow=nnodes, byrow=T) #sarà giusto? o devo distribuire per colonne?

points=eval.FEM(output_CPP$fit.FEM, locations=nodesLocations)

write.table(points, file="smoothFEMbasis3D_nod_nocov_stochGCV.txt")

### no covariates, observations different from node locations ###

loc=matrix(data=runif(3*150, min=-0.5,max=0.5),nrow=150,ncol=3,byrow=T) # 150 punti casuali all'interno della sfera

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

#exact GCV  8.16833
#      stderr 2.590291
#      edf  26.78744

#stoch GCV  8.166831
#      stderr 2.590172
#      edf  26.77613

points2=eval.FEM(output_CPP2$fit.FEM, locations=nodesLocations)
write.table(points2, file="smoothFEMbasis3D_nonod_nocov_stochGCV.txt")

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
# [1,] -0.25748060
# [2,] -0.01327809


#exact GCV  1.094607
#      stderr 0.8734922
#      edf  177.8355

#stoch GCV  1.087523
#      stderr  0.8720756
#      edf  176.5052

points3=eval.FEM(output_CPP3$fit.FEM, locations=nodesLocations)
write.table(points3, file="smoothFEMbasis3D_nod_cov_stochGCV.txt")

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
# [1,] -0.4678160
# [2,] -0.2023703

#exact GCV  8.239225
#      stderr 2.589103
#      edf  27.95962

#stoch GCV  8.241929
#      stderr 2.589315
#      edf  27.97964

points4=eval.FEM(output_CPP4$fit.FEM, locations=nodesLocations)
write.table(points4, file="smoothFEMbasis3D_nonod_cov_stochGCV.txt")
