######################
## test fPCA hub #####
######################

library(fdaPDE)

data(hub)

## NOTA: i write.table producono file .csv che possono essere importati in Matlab per un confronto

#cat('Plotting the mesh \n')
#plot.MESH2.5D(hub)

## Generate some random data ##

nnodes = hub$nnodes
cat("Nnodes: \n")
nnodes
cat("\n")

set.seed(5847947)

datamatrix<-NULL;
for(ii in 1:50){
  a1 = rnorm(1,mean = 1, sd = 1)
  a2 = rnorm(1,mean = 1, sd = 1)
  a3 = rnorm(1,mean = 1, sd = 1)
  
  func_evaluation = numeric(nnodes)
  
  for (i in 0:(nnodes-1)){
    func_evaluation[i+1] = a1* sin(2*pi*hub$nodes[3*i+1]) +  a2* sin(2*pi*hub$nodes[3*i+2]) +  a3*sin(2*pi*hub$nodes[3*i+3]) +1
  }
  
  data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
  datamatrix<-rbind(datamatrix,data)
}

#datamatrix

data_bar=colMeans(datamatrix)
data_demean=matrix(rep(data_bar,50),nrow=50,byrow=TRUE)

#data_demean

datamatrix_demeaned=datamatrix-data_demean

# write(datamatrix,file="datamatrix_hub.csv", sep=",",ncolumns=50)
# write(datamatrix_demeaned,file="datamatrix_hub_demeaned.csv", sep=",",ncolumns=50)

FEMbasis <- create.FEM.basis(hub)

### No validation ###

lambda=c(0.00375)
output_CPP1 =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda,
                             nPC=5)

str(output_CPP1)

#cat("Loadings\n")
#output_CPP$loadings
#cat("Scores\n")
#output_CPP$scores
write(output_CPP1$loadings.FEM$coeff[,1],file="loadings1_hub.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,1],file="scores1_hub.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,2],file="loadings2_hub.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,2],file="scores2_hub.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,3],file="loadings3_hub.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,3],file="scores3_hub.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,4],file="loadings4_hub.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,4],file="scores4_hub.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,5],file="loadings5_hub.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,5],file="scores5_hub.csv", sep=",",ncolumns=50)


plot(output_CPP1$loadings.FEM)

x11()
plot(output_CPP1$variance_explained,type="o")

x11()
plot(output_CPP1$cumsum_percentage,type="o")


### K-fold cross validation ###

lambda1=c(-4,-3,-2,-1,0,1,2,3,4)
lambda1=10^lambda1
output_CPP2 =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda1,
                             nPC=5,validation="KFold",NFolds=5)

write(output_CPP2$loadings.FEM$coeff[,1],file="loadings1_hub_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,1],file="scores1_hub_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,2],file="loadings2_hub_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,2],file="scores2_hub_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,3],file="loadings3_hub_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,3],file="scores3_hub_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,4],file="loadings4_hub_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,4],file="scores4_hub_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,5],file="loadings5_hub_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,5],file="scores5_hub_KF.csv", sep=",",ncolumns=50)

plot(output_CPP2$loadings.FEM)

x11()
plot(output_CPP2$variance_explained,type="o")

x11()
plot(output_CPP2$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP2$lambda)))

### exact GCV validation ###

output_CPP3 =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda1,
                             nPC=5,validation="GCV",GCVmethod = 1)

write(output_CPP3$loadings.FEM$coeff[,1],file="loadings1_hub_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,1],file="scores1_hub_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,2],file="loadings2_hub_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,2],file="scores2_hub_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,3],file="loadings3_hub_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,3],file="scores3_hub_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,4],file="loadings4_hub_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,4],file="scores4_hub_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,5],file="loadings5_hub_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,5],file="scores5_hub_GCVexact.csv", sep=",",ncolumns=50)

plot(output_CPP3$loadings.FEM)

x11()
plot(output_CPP3$variance_explained,type="o")

x11()
plot(output_CPP3$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP3$lambda)))

### stochastic GCV validation ###

output_CPP4 =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda1,
                             nPC=5,validation="GCV",GCVmethod = 2)

write(output_CPP4$loadings.FEM$coeff[,1],file="loadings1_hub_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,1],file="scores1_hub_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,2],file="loadings2_hub_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,2],file="scores2_hub_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,3],file="loadings3_hub_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,3],file="scores3_hub_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,4],file="loadings4_hub_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,4],file="scores4_hub_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,5],file="loadings5_hub_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,5],file="scores5_hub_GCV.csv", sep=",",ncolumns=50)

plot(output_CPP4$loadings.FEM)

x11()
plot(output_CPP4$variance_explained,type="o")

x11()
plot(output_CPP4$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP4$lambda)))