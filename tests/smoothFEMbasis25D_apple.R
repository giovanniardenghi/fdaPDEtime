################################
## test smooth.FEM.basis 2.5D ##
################################

library(fdaPDE)

GCVFLAG=TRUE  # =FALSE do not compute GCV (default)
#GCVFLAG=FALSE

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

GCVMETHODFLAG=1 #for exact GCV
#GCVMETHODFLAG=2 #for stochastic GCV (default)

set.seed(5847947)

load("Data/2.5D/apple.RData")
hub<-apple

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
write(func_evaluation,file="data_apple_25D.csv", sep=",",ncolumns=1)  #called F in matlab script

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
write(data,file="func_apple_25D.csv", sep=",",ncolumns=1)  # called z in matlab script

FEMbasis <- create.FEM.basis(hub)

lambda=c(0.00375)
# output_CPP =smooth.FEM.basis(observations = data,
#                              FEMbasis = FEMbasis, lambda = lambda,
#                              CPP_CODE = TRUE)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda,
                             CPP_CODE = TRUE, GCV = GCVFLAG,GCVmethod = GCVMETHODFLAG)

cat("Showing result")
plot(output_CPP$fit.FEM)

# exact GCV  0.254867
#       stderr 0.5016631 
#       edf  10.88973
#       

# stoch GCV 0.2542807
#       stderr  0.5013743
#       edf  9.903201
#  

points=eval.FEM(output_CPP$fit.FEM, locations=hub$nodes) # non va bene!

nodesLocations=matrix(data=hub$nodes, ncol=3, nrow=nnodes, byrow=T) #sarà giusto? o devo distribuire per colonne?

points=eval.FEM(output_CPP$fit.FEM, locations=nodesLocations) # il primo è un NA
write.table(points, file="smoothFEMbasis25D_apple_nod_nocov_exactGCV.txt")

### covariates, observations at node locations ###

cov1=3*nodesLocations[,1]+2*nodesLocations[,2]+5*nodesLocations[,3]+rnorm(nnodes,mean=0,sd=0.1)
cov2=rnorm(nnodes, mean=3, sd=1)

# output_CPP2=smooth.FEM.basis(observations = data, covariates = cbind(cov1,cov2), FEMbasis = FEMbasis, lambda =lambda)

output_CPP2=smooth.FEM.basis(observations = data, covariates = cbind(cov1,cov2), 
                             FEMbasis = FEMbasis, lambda =lambda, GCV=GCVFLAG,GCVmethod = GCVMETHODFLAG)

plot(output_CPP2$fit.FEM)

# exact GCV  0.256168
#       stderr 0.5024047 
#       edf 12.71736  

# stochastic GCV  0.255355
#       stderr 0.5020056
#       edf  11.35848 

points2=eval.FEM(output_CPP2$fit.FEM, locations=nodesLocations) # il primo è un NA
write.table(points2, file="smoothFEMbasis25D_apple_nod_cov_exactGCV.txt")


##### differenze
# # salvo i coefficienti
# write.table(output_CPP$fit.FEM$coeff, file=".txt") #coeff_25D_nod_cov_Hugo_mod/orig_apple
# 
# 
# coeff_Hugo_orig_matlab=read.csv("file:///C:/Users/unico/Desktop/PROGETTO PACS/TEST/2.5D/Apple/coeff25D_nod_nocov_matlab_Hugo_orig_apple.csv", header = FALSE)
# coeff_Hugo_mod=read.table("file:///C:/Users/unico/Desktop/PROGETTO PACS/TEST/2.5D/Apple/coeff_25D_nod_cov_Hugo_mod_apple.txt")
# coeff_Hugo_orig=read.table("file:///C:/Users/unico/Desktop/PROGETTO PACS/TEST/2.5D/Apple/coeff_25D_nod_cov_Hugo_orig_apple.txt")
# 
# # coefficienti R vs coefficienti matlab (da versione originale)
# diff=abs(coeff_Hugo_orig_matlab-coeff_Hugo_mod)
# maxdif=max(diff)
# maxdif
# write.table(diff, file="diff_coeff_25D_nod_cov_apple_Hugo_mod_vs_Matlab_Lila.txt")
# 
# # coefficienti R da versione originale vs coefficienti R da versione modificata
# diff2=abs(coeff_Hugo_orig-output_CPP$fit.FEM$coeff)
# maxdif2=max(diff2)
# maxdif2
# write.table(diff2, file="diff_coeff_25D_nod_cov_apple_Hugo_mod_vs_Hugo_orig.txt")
