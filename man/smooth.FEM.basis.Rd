% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothing.R
\name{smooth.FEM.time}
\alias{smooth.FEM.time}
\title{Spatial regression with differential regularization}
\usage{
smooth.FEM.time(locations = NULL, observations, FEMbasis, lambda, 
       covariates = NULL, BC = NULL, GCV = FALSE)
}
\arguments{
\item{locations}{A #observations-by-ndim matrix where each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} if ndim=3) of the corresponding observations in the vector \code{observations}.
This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.}

\item{observations}{A vector of length #observations with the observed data values over the domain. 
The locations of the observations can be specified with the \code{locations} argument. 
Otherwise if only the vector of observations is given, these are considered to be located in the corresponding node in the table
\code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
 node.}

\item{FEMbasis}{A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.}

\item{lambda}{A scalar or vector of smoothing parameters.}

\item{covariates}{A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.}

\item{PDE_parameters}{A list specifying the parameters of the PDE in the regularizing term. Default is NULL, i.e. regularization is by means of the Laplacian (stationary, isotropic case). 
 If the PDE is elliptic it must contain: \code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K; \code{b}, a vector of length 2 of advection coefficients. This induces a 
smoothing only in the direction specified by the vector \code{b}; \code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero
If the PDE is space-varying it must contain: \code{K}, a function that for each spatial location in the spatial domain 
(indicated by the vector of the 2 spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
smoothing with a local preferential direction that corresponds to the first eigenvector of the diffusion matrix K.The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
an array with dimensions 2-by-2-by-#points.\code{b}, a function that for each spatial location in the spatial domain returns 
a vector of length 2 of transport coefficients. This induces a local smoothing only in the direction specified by the vector \code{b}. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
a matrix with dimensions 2-by-#points; \code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
\code{c} induces a shrinkage of the surface to zero. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
a vector with length #points; \code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
\code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
a vector with length #points.
For 2.5D and 3D only the Laplacian is available (\code{PDE_parameters=NULL})}

\item{incidence_matrix}{A #regions-by-#triangles/tetrahedrons matrix where the element (i,j) equals 1 if the j-th triangle/tetrahedron is in the i-th region and 0 otherwise.
This is only for areal data. In case of pointwise data, this parameter is set to \code{NULL}.}

\item{BC}{A list with two vectors: 
\code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
\code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.}

\item{GCV}{Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.}

\item{GCVmethod}{either "Exact" or "Stochastic". If set to "Exact" perform an exact (but possibly slow) computation of the GCV index. If set to "Stochastic" approximate the GCV with a stochastic algorithm.
This parameter is considered only when \code{GCV=TRUE}}
}
\value{
A list with the following variables:
\item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
\item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
\item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
the smoothing parameter is equal to \code{lambda[j]}.}
\item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
\item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
\item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
}
\description{
This function implements a spatial regression model with differential regularization; isotropic and stationary case. In particular, the regularizing term involves the Laplacian of the spatial field. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
}
\examples{
library(fdaPDE)
## Load the Meuse data and a domain boundary for these data
data(MeuseData)
data(MeuseBorder)
## Create a triangular mesh for these data with the provided boundary and plot it
order=1
mesh <- create.mesh.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
plot(mesh)
## Create the Finite Element basis 
FEMbasis = create.FEM.basis(mesh)
## Estimate zync field without using covariates, setting the smoothing parameter to 10^3.5
data = log(MeuseData[,"zinc"])
lambda = 10^3.5
ZincMeuse = smooth.FEM.time(observations = data, 
                             FEMbasis = FEMbasis, lambda = lambda)
## Plot the estimated spatial field 
plot(ZincMeuse$fit.FEM)
# Now repeat the analysis using as covariates the square root of the log-distance 
# from river \\code{sqrt(dist.log(m))} and the altitude \\code{elev}
desmat = matrix(1,nrow=nrow(MeuseData),ncol=2)
desmat[,1] = sqrt(MeuseData[,"dist.log(m)"])
desmat[,2] = MeuseData[,"elev"]
ZincMeuseCovar = smooth.FEM.time(observations = data, covariates = desmat, 
                                   FEMbasis = FEMbasis, lambda = lambda)
# Plot of the non parametric part (f) of the regression model y_i = beta_1 x_i1 + beta_2 x_i2 + f
plot(ZincMeuseCovar$fit.FEM)
# Print covariates' regression coefficients
print(ZincMeuseCovar$beta)
}
\references{
Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
}
