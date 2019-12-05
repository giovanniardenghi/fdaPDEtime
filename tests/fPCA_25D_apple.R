########################
## test fPCA apple #####
########################

library(fdaPDE)

data(apple)

nodes=matrix(apple$nodes,ncol=3,byrow=T)
triangles=matrix(apple$triangles,ncol=3,byrow=T)
hub<-create.MESH.2.5D(nodes,triangles)


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

# write(datamatrix,file="datamatrix_apple.csv", sep=",",ncolumns=50)
# write(datamatrix_demeaned,file="datamatrix_apple_demeaned.csv", sep=",",ncolumns=50)

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
write(output_CPP1$loadings.FEM$coeff[,1],file="loadings1_apple.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,1],file="scores1_apple.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,2],file="loadings2_apple.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,2],file="scores2_apple.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,3],file="loadings3_apple.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,3],file="scores3_apple.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,4],file="loadings4_apple.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,4],file="scores4_apple.csv", sep=",",ncolumns=50)
write(output_CPP1$loadings.FEM$coeff[,5],file="loadings5_apple.csv", sep=",",ncolumns=nnodes)
write(output_CPP1$scores[,5],file="scores5_apple.csv", sep=",",ncolumns=50)


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

write(output_CPP2$loadings.FEM$coeff[,1],file="loadings1_apple_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,1],file="scores1_apple_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,2],file="loadings2_apple_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,2],file="scores2_apple_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,3],file="loadings3_apple_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,3],file="scores3_apple_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,4],file="loadings4_apple_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,4],file="scores4_apple_KF.csv", sep=",",ncolumns=50)
write(output_CPP2$loadings.FEM$coeff[,5],file="loadings5_apple_KF.csv", sep=",",ncolumns=nnodes)
write(output_CPP2$scores[,5],file="scores5_apple_KF.csv", sep=",",ncolumns=50)

plot(output_CPP2$loadings.FEM)

x11()
plot(output_CPP2$variance_explained,type="o")

x11()
plot(output_CPP2$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP2$lambda)))

### exact GCV validation ###

output_CPP3 =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda1,
                             nPC=5,validation="GCV",GCVmethod = "Exact")

write(output_CPP3$loadings.FEM$coeff[,1],file="loadings1_apple_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,1],file="scores1_apple_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,2],file="loadings2_apple_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,2],file="scores2_apple_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,3],file="loadings3_apple_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,3],file="scores3_apple_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,4],file="loadings4_apple_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,4],file="scores4_apple_GCVexact.csv", sep=",",ncolumns=50)
write(output_CPP3$loadings.FEM$coeff[,5],file="loadings5_apple_GCVexact.csv", sep=",",ncolumns=nnodes)
write(output_CPP3$scores[,5],file="scores5_apple_GCVexact.csv", sep=",",ncolumns=50)

plot(output_CPP3$loadings.FEM)

x11()
plot(output_CPP3$variance_explained,type="o")

x11()
plot(output_CPP3$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP3$lambda)))


### stochastic GCV validation ###

output_CPP4 =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda1,
                             nPC=5,validation="GCV",GCVmethod = "Stochastic")

write(output_CPP4$loadings.FEM$coeff[,1],file="loadings1_apple_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,1],file="scores1_apple_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,2],file="loadings2_apple_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,2],file="scores2_apple_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,3],file="loadings3_apple_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,3],file="scores3_apple_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,4],file="loadings4_apple_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,4],file="scores4_apple_GCV.csv", sep=",",ncolumns=50)
write(output_CPP4$loadings.FEM$coeff[,5],file="loadings5_apple_GCVstoch.csv", sep=",",ncolumns=nnodes)
write(output_CPP4$scores[,5],file="scores5_apple_GCV.csv", sep=",",ncolumns=50)

plot(output_CPP4$loadings.FEM)

x11()
plot(output_CPP4$variance_explained,type="o")

x11()
plot(output_CPP4$cumsum_percentage,type="o")

print(paste("Best lambda:", log10(output_CPP4$lambda)))
