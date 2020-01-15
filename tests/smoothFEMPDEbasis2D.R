###############################
## test smooth.FEM.PDE.basis ##
###############################

library(fdaPDE)

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

GCVFLAG=TRUE  # =FALSE do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

#GCVMETHODFLAG="Exact" #for exact GCV
GCVMETHODFLAG="Stochastic" #for stochastic GCV (default)

set.seed(5847947)

# Load the mesh and plot it
data(mesh.2D.simple)
x11()
plot(mesh.2D.simple)


### observations at node locations

truefield=sin(pi*mesh.2D.simple$nodes[,1])
observations = sin(pi*mesh.2D.simple$nodes[,1]) + rnorm(n = nrow(mesh.2D.simple$nodes), sd = 0.1)

# Create the FEM basis object
FEMbasis = create.FEM.basis(mesh.2D.simple)
#image(FEM(truefield,FEMbasis = FEMbasis))

# Set a vector of smoothing coefficients
lambda = c(10^-4, 1, 10^4)

# Set the anysotropic smoothing matrix K
PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)

# Estimate one field for each smoothing parameter and plot these
#FEM_CPP_PDE = smooth.FEM.PDE.basis(observations = observations, 
#                                   FEMbasis = FEMbasis, lambda = lambda, 
#                                  PDE_parameters = PDE_parameters_anys, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)

FEM_CPP_PDE = smooth.FEM.basis(observations = observations, 
                                   FEMbasis = FEMbasis, lambda = lambda, 
                                   PDE_parameters = PDE_parameters_anys, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

plot(FEM_CPP_PDE$fit.FEM)

points=eval.FEM(FEM_CPP_PDE$fit.FEM, locations = mesh.2D.simple$nodes)
write.table(points, file="smoothFEMPDEbasis2D_nod_nocov.txt")

# exact GCV:0.03729336 0.01430638 0.12801092
#       edf:100.621313  15.451292   1.874006
#       stderr:0.05047699 0.11072310 0.35466841

# stoch GCV:0.03628253 0.01419894 0.12731236
#       edf:100.519234  15.101795   1.583249
#       stderr:0.05013141 0.11051462 0.35418356

### locations different from nodes

xobs=runif(min=-0.5,max=0.5,n=80)
yobs=runif(min=-0.5,max=0.5,n=80)
points(xobs,yobs,col='blue')
observations2 = sin(pi*xobs) + rnorm(n = length(xobs), sd = 0.1)
#FEM_CPP_PDE2 = smooth.FEM.basis(locations=cbind(xobs,yobs),observations = observations2, 
#                                    FEMbasis = FEMbasis, lambda = lambda, 
#                                    PDE_parameters = PDE_parameters_anys,GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
FEM_CPP_PDE2 = smooth.FEM.basis(locations=cbind(xobs,yobs),observations = observations2, 
                                     FEMbasis = FEMbasis, lambda = lambda, 
                                     PDE_parameters = PDE_parameters_anys,GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
                                    

plot(FEM_CPP_PDE2$fit.FEM)

points2=eval.FEM(FEM_CPP_PDE2$fit.FEM,locations=mesh.2D.simple$nodes)
write.table(points2, file="smoothFEMPDEbasis2D_nonod_nocov.txt")

# exact GCV:0.010726484 0.008111932 0.135202953
#       edf:43.844443  8.992484  1.549041
#       stderr:0.06962596 0.08485341 0.36412226

# stoch GCV:0.01008789 0.00807306 0.13434940
#       edf:42.717631  8.821738  1.300227
#       stderr:0.06856571 0.08475158 0.36354620

### observations at node locations with covariates

cov1=sin(pi*mesh.2D.simple$nodes[,1])
cov2=rnorm(mean=0, sd=0.5,n=length(mesh.2D.simple$nodes[,1]))
FEM_CPP_PDE3 = smooth.FEM.basis(observations = observations, covariates=cbind(cov1,cov2),
                                    FEMbasis = FEMbasis, lambda = lambda, 
                                    PDE_parameters = PDE_parameters_anys, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

FEM_CPP_PDE3 = smooth.FEM.basis(observations = observations, covariates=cbind(cov1,cov2),
                                    FEMbasis = FEMbasis, lambda = lambda, 
                                    PDE_parameters = PDE_parameters_anys, GCV=GCVFLAG)


points3=eval.FEM(FEM_CPP_PDE3$fit.FEM,locations=mesh.2D.simple$nodes)
write.table(points3, file="smoothFEMPDEbasis2D_nod_cov.txt")
plot(FEM_CPP_PDE3$fit.FEM)
FEM_CPP_PDE3$beta
#                     [,1]                 [,2]                 [,3]
# [1,]  0.22379931632932593  0.75498735993142074  1.00669466219865256
# [2,] -0.07480413617039104 -0.02786067961962656 -0.01943171430804428

# exact GCV:0.03470415 0.01448802 0.01340897
#       edf:100.918910  16.518433   3.392779
#       stderr:0.04770116 0.11077950 0.11396373

# stoch GCV:0.0006380062 0.0146805773 0.0133384873
#       edf:160.224995  17.120381   3.116763
#       stderr:NaN 0.1111458 0.1138137

### observations different from node locations, with covariates

cov1=sin(pi*xobs)
cov2=rnorm(mean=0, sd=0.5,n=length(xobs))
FEM_CPP_PDE4 = smooth.FEM.basis(observations = observations2, locations=cbind(xobs,yobs), covariates=cbind(cov1,cov2),
                                    FEMbasis = FEMbasis, lambda = lambda, 
                                    PDE_parameters = PDE_parameters_anys, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

FEM_CPP_PDE4 = smooth.FEM.basis(observations = observations2, locations=cbind(xobs,yobs), covariates=cbind(cov1,cov2),
                                    FEMbasis = FEMbasis, lambda = lambda, 
                                    PDE_parameters = PDE_parameters_anys, GCV=GCVFLAG)

points4=eval.FEM(FEM_CPP_PDE4$fit.FEM,locations=mesh.2D.simple$nodes)
write.table(points4, file="smoothFEMPDEbasis2D_nonod_cov.txt")
plot(FEM_CPP_PDE4$fit.FEM)
FEM_CPP_PDE4$beta

# [,1]        [,2]       [,3]
# [1,]  0.01343384431183954  0.37958127642305461  0.9839547016485253
# [2,] -0.04918467667551139 -0.04997457636598323 -0.0526896964475254

# exact GCV:0.010343063 0.007832468 0.007009382
#       edf:45.594815 10.202435  3.065789
#       stderr: 0.06669473 0.08266553 0.08210217

# stoch GCV:0.132739404 0.008928010 0.006973934
#       edf:70.396080 14.624919  2.870508
#       stderr:0.12623483 0.08541585 0.08199817