################
## test BC 2D ##
################

library(fdaPDE)
setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")
source('pointsFEM.R')

GCVFLAG=FALSE  # do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

# GCVMETHODFLAG=1 for exact GCV
# GCVMETHODFLAG=2 for stochastic GCV (default)

load("Data/2D/DatiEsempioSV.RData")
x11()
plot(mesh)

x11()
points.2D.data(SpacePoints = SpacePoints, Data = DatiEsatti)

FEMbasis = create.FEM.basis(mesh)

# Set BC 
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1)
BC$BC_values = rep(0,length(BC$BC_indices))

# Set PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

lambda = 10^-2

set.seed(5839745)

Data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))

### smooth.FEM.basis ###

Sol = smooth.FEM.basis(locations = SpacePoints,
                       observations = Data, 
                       FEMbasis = FEMbasis, 
                       lambda = lambda, 
                       BC = BC, GCV=GCVFLAG)
plot(Sol$fit.FEM)                                                

points=eval.FEM(Sol$fit.FEM, locations = mesh$nodes)
write.table(points, file="NEWsmoothFEMbasis2D_BC_nonod_nocov.txt")


### smooth.FEM.PDE.basis ###

PDE_parameters = list(K = matrix(10*c(K1,0,0,K2), nrow = 2), b = c(0,0), c = 0)

Sol = smooth.FEM.basis(locations = SpacePoints,
                           observations = Data, 
                           FEMbasis = FEMbasis, 
                           lambda = lambda, PDE_parameters = PDE_parameters,
                           BC = BC, GCV=GCVFLAG)
plot(Sol$fit.FEM)                                                

points=eval.FEM(Sol$fit.FEM, locations = mesh$nodes)
write.table(points, file="NEWsmoothFEMPDEbasis2D_BC_nonod_nocov.txt")


### smooth.FEM.PDE.sv.basis ###

K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
                             (K1-1)*points[i,1]*points[i,2]),
                           c((K1-1)*points[i,1]*points[i,2],
                             points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
  output
}

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
  output[,i] = 10*beta*c(points[i,1],points[i,2])  
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

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)


Data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))

Sol = smooth.FEM.basis(locations = SpacePoints,
                              observations = Data, 
                              FEMbasis = FEMbasis, 
                              lambda = lambda, 
                              BC = BC,
                              PDE_parameters = PDE_parameters,GCV=GCVFLAG) 
plot(Sol$fit.FEM)

points=eval.FEM(Sol$fit.FEM, locations = mesh$nodes)
write.table(points, file="NEWsmoothFEMPDEsvbasis2D_BC_nonod_nocov.txt")


