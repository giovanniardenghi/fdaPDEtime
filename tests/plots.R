##################
### plot tests ###
##################

#### plot 3D tests ####

data(sphere3Ddata)
nodes=sphere3Ddata$nodes
tetrahedrons=sphere3Ddata$tetrahedrons

#Create the triangulated mesh from the connectivity matrix and nodes locations
sphere3D=create.MESH.3D(nodes,tetrahedrons)

##Define the number of nodes
nnodes=sphere3D$nnodes
## Create a linear Finite Element basis
FEMbasis = create.FEM.basis(sphere3D)
## Define a function as expansion of this basis
coeff=numeric(nnodes)
for(i in 0:nnodes-1){
  coeff[i+1] <- 0.5*sqrt(15/pi)*sphere3D$nodes[i+1,1]*sphere3D$nodes[i+1,2]
}
FEM_object<- FEM(coeff, FEMbasis)
##Create the grid of evaluated points
FEMevaluated=evaluate.grid.FEM(FEM_object)

##Plot the slices
slices.FEMevaluated(FEMevaluated,xs=c(-0.5,0,0.5),ys=c(-0.5,0,0.5),NAcol="transparent")

slicecontours.FEMevaluated(FEMevaluated,xs=c(-0.5,0,0.5),ys=0,NAcol="blue")

##Plot the isosurfaces
isosurfaces.FEMevaluated(FEMevaluated,levels=seq(-0.45,0.45,by=0.15))


#### plots mesh2D ####

# order 1

data(meshC)
plot.MESH2D(mesh)
#image.FEM(mesh)

## Create a linear Finite Element basis
FEMbasis = create.FEM.basis(mesh)
## Define a function as expansion of this basis
coeff=numeric(dim(mesh$nodes)[1])
for(i in 0:dim(mesh$nodes)[1]-1){
  coeff[i+1] <- 0.5*sqrt(15/pi)*mesh$nodes[i+1,1]*mesh$nodes[i+1,2]
}
FEM_object<- FEM(coeff, FEMbasis)
plot.FEM(FEM_object)

# order 2

data("mesh.2D.simple")
plot.MESH2D(mesh.2D.simple)
#image.FEM(mesh)

FEMbasis = create.FEM.basis(mesh.2D.simple)
## Define a function as expansion of this basis
coeff=numeric(108)
for(i in 0:(108-1)){
  coeff[i+1] <- 0.5*sqrt(15/pi)*mesh.2D.simple$nodes[i+1,1]*mesh.2D.simple$nodes[i+1,2]
}
FEM_object<- FEM(coeff, FEMbasis)
plot.FEM(FEM_object)

#### plots mesh.2.5D ####

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

nome_mesh="C2_5D"

vertici <- read.table(paste0("Data/2.5D/",nome_mesh,"_vertici.txt"), quote="\"", comment.char="")
triangoli <- read.table(paste0("Data/2.5D/",nome_mesh,"_triangoli.txt"), quote="\"", comment.char="")

#order 1

mesh <- fdaPDE::create.MESH.2.5D(nodes = vertici[,1:3],triangles = triangoli[,1:3])
plot.MESH.2.5D(mesh)

FEMbasis = create.FEM.basis(mesh)
## Define a function as expansion of this basis
coeff=numeric(dim(mesh$nodes)[1])
for(i in 0:dim(mesh$nodes)[1]-1){
  coeff[i+1] <- 0.5*sqrt(15/pi)*mesh$nodes[i+1,1]*mesh$nodes[i+1,2]
}
FEM_object<- FEM(coeff, FEMbasis)
plot.FEM(FEM_object)

#order 2

mesh=second.order.MESH.2.5D(mesh)
plot.MESH.2.5D(mesh)

FEMbasis = create.FEM.basis(mesh)
## Define a function as expansion of this basis
coeff=numeric(dim(mesh$nodes)[1])
for(i in 0:dim(mesh$nodes)[1]-1){
  coeff[i+1] <- 0.5*sqrt(15/pi)*mesh$nodes[i+1,1]*mesh$nodes[i+1,2]
}
FEM_object<- FEM(coeff, FEMbasis)
plot.FEM(FEM_object)

#### plots mesh.3D ####

setwd("C:/Users/Gianmaria/Desktop/progetto_pacs/script x test")

nome_mesh = "meshCcicciona"

vertici <- read.table(paste0("Data/3D/",nome_mesh,"_vertici.txt"), quote="\"", comment.char="")
tetraedri <- read.table(paste0("Data/3D/",nome_mesh,"_tetraedri.txt"), quote="\"", comment.char="")

mesh <- fdaPDE::create.MESH.3D(nodes = vertici[,1:3],tetrahedrons = tetraedri[,1:4])

plot.MESH.3D(mesh)

FEMbasis = create.FEM.basis(mesh)
## Define a function as expansion of this basis
coeff=numeric(dim(mesh$nodes)[1])
for(i in 0:dim(mesh$nodes)[1]-1){
  coeff[i+1] <- 0.5*sqrt(15/pi)*mesh$nodes[i+1,1]*mesh$nodes[i+1,2]
}
FEM_object<- FEM(coeff, FEMbasis)
plot.FEM(FEM_object)



tetrahedrons = FEM_object$FEMbasis$mesh$tetrahedrons

nodes=FEM_object$FEMbasis$mesh$nodes

ntetrahedrons = FEM_object$FEMbasis$mesh$ntetrahedrons

tet=t(rbind(tetrahedrons[,-1],tetrahedrons[,-2],tetrahedrons[,-3],tetrahedrons[,-4]))

coeff = FEM_object$coeff

nsurf = dim(coeff)[[2]]
cols=NULL
p=jet.col(n=128,alpha=0.8)
#p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(128)
ncolors=length(p)
for (isurf in 1:nsurf)
{	col=rep(0,ntetrahedrons)
for(j in 1:ntetrahedrons)
  col[j]=mean(c(coeff[tetrahedrons[j,1],isurf],coeff[tetrahedrons[j,2],isurf],
                coeff[tetrahedrons[j,3],isurf],coeff[tetrahedrons[j,4],isurf]))


col=rep(1,3) %o% col


diffrange = max(col)-min(col)

col= (col - min(col))/diffrange*(ncolors-1)+1
col=p[col]
cols=c(cols,col)
open3d()
axes3d()
rgl.pop("lights") 
light3d(specular="black") 

rgl.triangles(nodes[tet,1],nodes[tet,2],nodes[tet,3],col=col,alpha=0.7)

#tetramesh(tetrahedrons,nodes,col=col,alpha=0.7)
}
legend3d("topleft", legend=levels(cut(c(-7,7), 10)), col=rainbow(10), pch=20)
