#########################
## test fPCA 3D sphere ##
#########################

library("fdaPDE")

setwd("C:/Users/Gianmaria/Desktop/progetto pacs/script x test")

nome_mesh = "Sfera"
nome_eigenf1 = "Sfera_eigenf_1.RData"
nome_eigenf2 = "Sfera_eigenf_2.RData"
nome_eigenf3 = "Sfera_eigenf_3.RData"

seed = 5728969
sd_score1 = 0.1
sd_score2 = 0.05
sd_score3 = 0.0125
sd_errore = 0.075

# we create the 3D mesh with create.FEM.3D

vertici <- read.table(paste0("Data/3D/",nome_mesh,"_vertici.txt"), quote="\"", comment.char="")
tetraedri <- read.table(paste0("Data/3D/",nome_mesh,"_tetraedri.txt"), quote="\"", comment.char="")

mesh <- fdaPDE::create.MESH.3D(nodes = vertici[,1:3],tetrahedrons = tetraedri[,1:4])
FEMbasis <- fdaPDE::create.FEM.basis(mesh)

# Load data for exact solution

load(paste0("Data/3D/",nome_eigenf1))
eigenfunc1 = FEM(coeff = sol_esatta, FEMbasis = FEMbasis)

load(paste0("Data/3D/",nome_eigenf2))
eigenfunc2 = FEM(coeff = sol_esatta, FEMbasis = FEMbasis)

load(paste0("Data/3D/",nome_eigenf3))
eigenfunc3 = FEM(coeff = sol_esatta, FEMbasis = FEMbasis)

truedatarange = max(c(eigenfunc1$coeff, eigenfunc2$coeff, eigenfunc3$coeff)) - min(c(eigenfunc1$coeff, eigenfunc2$coeff, eigenfunc3$coeff))
# MEMO: forse va cambiato a seconda dei dati e della mesh:
# lambda_fpca = 10^seq(-4,0, by = 0.1)  # errore 0.08 Sfera
lambda_fpca = 10^-2

# 
nlocations = nrow(vertici)
nsamples = 50
norm_vec <- function(x) sqrt(sum(x^2))

set.seed(seed) 

# Generate random data

score1 = rnorm(n = nsamples, sd = sd_score1*truedatarange)
score2 = rnorm(n = nsamples, sd = sd_score2*truedatarange)
score3 = rnorm(n = nsamples, sd = sd_score3*truedatarange)
errore = rnorm(n = nsamples*nlocations, sd = sd_errore*truedatarange)

datamatrix_esatta = matrix(score1) %*% t(matrix(eigenfunc1$coeff)) + matrix(score2) %*% t(matrix(eigenfunc2$coeff)) + matrix(score3) %*% t(matrix(eigenfunc3$coeff))
datamatrix = datamatrix_esatta + errore

# Center data
data_bar = colMeans(datamatrix, na.rm = TRUE)
data_bar = matrix(rep(data_bar, nrow(datamatrix)), nrow = nrow(datamatrix), byrow = TRUE)
datamatrix_no_mean = datamatrix - data_bar

# write(datamatrix,file="datamatrix_hub.csv", sep=",",ncolumns=50)
# write(datamatrix_no_mean,file="datamatrix_meshC3D_demeaned.csv", sep=",",ncolumns=50)

### No validation ###

output_CPP1 = smooth.FEM.FPCA(datamatrix = datamatrix_no_mean, FEMbasis = FEMbasis, lambda = lambda_fpca, 
                              nPC = 3)

write(output_CPP1$loadings.FEM$coeff[,1],file="loadings1_sphere3D.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,1],file="scores1_sphere3D.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,2],file="loadings2_sphere3D.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,2],file="scores2_sphere3D.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,3],file="loadings3_sphere3D.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,3],file="scores3_sphere3D.csv", sep=",",ncolumns=50)

plot(output_CPP1$loadings.FEM)

x11()
plot(output_CPP1$variance_explained,type="o")

x11()
plot(output_CPP1$cumsum_percentage,type="o")


### K-fold cross validation ###

lambda_fpca = 10^seq(-4,0, by = 0.1)

output_CPP2 = smooth.FEM.FPCA(datamatrix = datamatrix_no_mean, FEMbasis = FEMbasis, lambda = lambda_fpca, 
                              nPC = 3, validation = "KFold")

write(output_CPP2$loadings.FEM$coeff[,1],file="loadings1_sphere3D_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,1],file="scores1_sphere3D_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,2],file="loadings2_sphere3D_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,2],file="scores2_sphere3D_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,3],file="loadings3_sphere3D_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,3],file="scores3_sphere3D_KF.csv", sep=",",ncolumns=50)

plot(output_CPP2$loadings.FEM)

x11()
plot(output_CPP2$variance_explained,type="o")

x11()
plot(output_CPP2$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP2$lambda)))


### exact GCV validation ###

lambda_fpca = 10^seq(-4,0, by = 0.1)

output_CPP3 = smooth.FEM.FPCA(datamatrix = datamatrix_no_mean, FEMbasis = FEMbasis, lambda = lambda_fpca, 
                              nPC = 3, validation = "GCV",GCVmethod = 1)

write(output_CPP3$loadings.FEM$coeff[,1],file="loadings1_sphere3D_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,1],file="scores1_sphere3D_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,2],file="loadings2_sphere3D_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,2],file="scores2_sphere3D_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,3],file="loadings3_sphere3D_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,3],file="scores3_sphere3D_GCVexact.csv", sep=",",ncolumns=50)

plot(output_CPP3$loadings.FEM)

x11()
plot(output_CPP3$variance_explained,type="o")

x11()
plot(output_CPP3$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP3$lambda)))


### stochastic GCV validation ###

lambda_fpca = 10^seq(-4,0, by = 0.1)

output_CPP4 = smooth.FEM.FPCA(datamatrix = datamatrix_no_mean, FEMbasis = FEMbasis, lambda = lambda_fpca, 
                              nPC = 3, validation = "GCV",GCVmethod = 2)

write(output_CPP4$loadings.FEM$coeff[,1],file="loadings1_sphere3D_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,1],file="scores1_sphere3D_GCVstoch.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,2],file="loadings2_sphere3D_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,2],file="scores2_sphere3D_GCVstoch.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,3],file="loadings3_sphere3D_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,3],file="scores3_sphere3D_GCVstoch.csv", sep=",",ncolumns=50)

plot(output_CPP4$loadings.FEM)

x11()
plot(output_CPP4$variance_explained,type="o")

x11()
plot(output_CPP4$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP4$lambda)))

