##################################
## test smooth.FEM.PDE.sv.basis ##
##################################

library(fdaPDE)

setwd("C:/Users/Gianmaria/Desktop/progetto pacs/script x test")

GCVFLAG=TRUE  # =FALSE do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

#GCVMETHODFLAG="Exact" #for exact GCV
GCVMETHODFLAG="Stochastic" #for stochastic GCV (default)

set.seed(5847947)

# Loading the mesh
data(mesh.2D.rectangular)
# Create the FEM basis object
FEMbasis = create.FEM.basis(mesh.2D.rectangular)

### observations at node locations

truefield=sin(0.2*pi*mesh.2D.rectangular$nodes[,1])

observations = sin(0.2*pi*mesh.2D.rectangular$nodes[,1]) + 
  rnorm(n = nrow(mesh.2D.rectangular$nodes), sd = 0.1)

# Set the smoothing coefficient
lambda = c(10^-2)

#Set the space variant coefficients of the penalizying PDE
K_func<-function(points)
{
  mat<-c(0.01,0,0,1)
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 0.5*mat %*% t(points[i,1]^2)
  output
}

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 0
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}

u_func<-function(points)
{
  rep(c(0), nrow(points))
}

# Assemble the parameters in one object
PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)

# Estimate the underlying spatial field and plot these
#FEM_CPP_PDE = smooth.FEM.PDE.sv.basis(observations = observations, 
#                                     FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV = GCVFLAG, GCVmethod = GCVMETHODFLAG)
FEM_CPP_PDE = smooth.FEM.basis(observations = observations, FEMbasis = FEMbasis, lambda = lambda, 
                               PDE_parameters = PDE_parameters, GCV = GCVFLAG, GCVmethod = GCVMETHODFLAG)

plot(FEM_CPP_PDE$fit.FEM)


points=eval.FEM(FEM_CPP_PDE$fit.FEM, locations = mesh.2D.rectangular$nodes)
write.table(points, file="smoothFEMPDEsvbasis2D_nod_nocov.txt")
image(FEM(points,FEMbasis = FEMbasis))

# exact GCV: 0.01076758
#       edf: 200.921
#       stderr: 0.09944809

# stoch GCV: 0.01078943
#       edf: 203.2143
#       stderr: 0.0994985

### observations different from node locations
xobs=runif(min=1,max=20,n=1000)
yobs=runif(min=1,max=5,n=1000)
observations2=sin(0.2*pi*xobs) + 
  rnorm(n = length(xobs), sd = 0.1)
# FEM_CPP_PDE2 = smooth.FEM.PDE.sv.basis(observations = observations2, locations=cbind(xobs,yobs),
#                                        FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)
FEM_CPP_PDE2 = smooth.FEM.basis(observations = observations2, locations=cbind(xobs,yobs),
                                       FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

plot(FEM_CPP_PDE2$fit.FEM)

points2=eval.FEM(FEM_CPP_PDE2$fit.FEM, locations = mesh.2D.rectangular$nodes)
write.table(points2, file="smoothFEMPDEsvbasis2D_nonod_nocov.txt")
image(FEM(points2,FEMbasis = FEMbasis))

# exact GCV: 0.01060132
#       edf: 115.2342
#       stderr:0.09684878

# stoch GCV: 0.01058419
#       edf: 114.5185
#       stderr: 0.09680962

### observations at node locations, with covariates

cov1=sin(pi*mesh.2D.rectangular$nodes[,1])
cov2=rnorm(mean=1,sd=0.1,n=length(mesh.2D.rectangular$nodes[,1]))

# FEM_CPP_PDE3 = smooth.FEM.PDE.sv.basis(observations = observations, covariates=cbind(cov1,cov2),
#                                        FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
FEM_CPP_PDE3 = smooth.FEM.basis(observations = observations, covariates=cbind(cov1,cov2),
                                       FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)

plot(FEM_CPP_PDE3$fit.FEM)

points3=eval.FEM(FEM_CPP_PDE3$fit.FEM, locations = mesh.2D.rectangular$nodes)
write.table(points3, file="smoothFEMPDEsvbasis2D_nod_cov.txt")
image(FEM(points3,FEMbasis = FEMbasis))
FEM_CPP_PDE3$beta
#[1,] 0.03085058
#[2,] 0.00321524

# exact GCV:0.0107662
#       edf: 202.0619
#       stderr: 0.09941666

# stoch GCV: 0.01080935
#       edf: 206.5831
#       stderr: 0.09951612

### observations different from node locations, with covariates

cov1=sin(pi*xobs)
cov2=rnorm(mean=1,sd=0.1,n=length(xobs))

# FEM_CPP_PDE4 = smooth.FEM.PDE.sv.basis(observations = observations2, locations=cbind(xobs,yobs),
#                                        covariates=cbind(cov1,cov2),FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)
FEM_CPP_PDE4 = smooth.FEM.basis(observations = observations2, locations=cbind(xobs,yobs),
                                       covariates=cbind(cov1,cov2),FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)

plot(FEM_CPP_PDE4$fit.FEM)

points4=eval.FEM(FEM_CPP_PDE4$fit.FEM, locations = mesh.2D.rectangular$nodes)
write.table(points4, file="smoothFEMPDEsvbasis2D_nonod_cov.txt")
image(FEM(points4,FEMbasis = FEMbasis))
FEM_CPP_PDE4$beta
#[1,]  0.01729489
#[2,] -0.05339664

# exact GCV: 0.0105862
#       edf: 116.4147
#       stderr: 0.0967151

# stoch GCV: 0.01052669
#       edf: 113.9208
#       stderr: 0.0965789