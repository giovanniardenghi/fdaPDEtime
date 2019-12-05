#' Spatial regression with differential regularization
#' 
#' @param observations A vector of length #observations with the observed data values over the domain. 
#' The locations of the observations can be specified with the \code{locations} argument. 
#' Otherwise if only the vector of observations is given, these are considered to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, an \code{NA} value in the \code{observations} vector indicates that there is no observation associated to the corresponding
#'  node.
#' @param locations A #observations-by-ndim matrix where each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} if ndim=3) of the corresponding observations in the vector \code{observations}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{observations}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param covariates A #observations-by-#covariates matrix where each row represents the covariates associated with the corresponding observed data value in \code{observations}.
#' @param PDE_parameters A list specifying the parameters of the PDE in the regularizing term. Default is NULL, i.e. regularization is by means of the Laplacian (stationary, isotropic case). 
#'  If the PDE is elliptic it must contain: \code{K}, a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a preferential direction that corresponds to the first eigenvector of the diffusion matrix K; \code{b}, a vector of length 2 of advection coefficients. This induces a 
#' smoothing only in the direction specified by the vector \code{b}; \code{c}, a scalar reaction coefficient. \code{c} induces a shrinkage of the surface to zero
#' If the PDE is space-varying it must contain: \code{K}, a function that for each spatial location in the spatial domain 
#' (indicated by the vector of the 2 spatial coordinates) returns a 2-by-2 matrix of diffusion coefficients. This induces an anisotropic 
#' smoothing with a local preferential direction that corresponds to the first eigenvector of the diffusion matrix K.The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' an array with dimensions 2-by-2-by-#points.\code{b}, a function that for each spatial location in the spatial domain returns 
#' a vector of length 2 of transport coefficients. This induces a local smoothing only in the direction specified by the vector \code{b}. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a matrix with dimensions 2-by-#points; \code{c}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{c} induces a shrinkage of the surface to zero. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points; \code{u}, a function that for each spatial location in the spatial domain  returns a scalar reaction coefficient.
#' \code{u} induces a reaction effect. The function must support recycling for efficiency reasons, thus if the input parameter is a #point-by-2 matrix, the output should be
#' a vector with length #points.
#' For 2.5D and 3D only the Laplacian is available (\code{PDE_parameters=NULL})

#' @param incidence_matrix A #regions-by-#triangles/tetrahedrons matrix where the element (i,j) equals 1 if the j-th triangle/tetrahedron is in the i-th region and 0 otherwise.
#' This is only for areal data. In case of pointwise data, this parameter is set to \code{NULL}.
#' @param BC A list with two vectors: 
#'  \code{BC_indices}, a vector with the indices in \code{nodes} of boundary nodes where a Dirichlet Boundary Condition should be applied;
#'  \code{BC_values}, a vector with the values that the spatial field must take at the nodes indicated in \code{BC_indices}.
#' @param GCV Boolean. If \code{TRUE} the following quantities are computed: the trace of the smoothing matrix, the estimated error standard deviation,  and 
#'        the Generalized Cross Validation criterion, for each value of the smoothing parameter specified in \code{lambda}.
#' @param GCVmethod either "Exact" or "Stochastic". If set to "Exact" perform an exact (but possibly slow) computation of the GCV index. If set to "Stochastic" approximate the GCV with a stochastic algorithm.
#'        This parameter is considered only when \code{GCV=TRUE}

#' @return A list with the following variables:
#' \item{\code{fit.FEM}}{A \code{FEM} object that represents the fitted spatial field.}
#' \item{\code{PDEmisfit.FEM}}{A \code{FEM} object that represents the Laplacian of the estimated spatial field.}
#' \item{\code{beta}}{If covariates is not \code{NULL}, a matrix with number of rows equal to the number of covariates and numer of columns equal to length of lambda.  The \code{j}th column represents the vector of regression coefficients when 
#' the smoothing parameter is equal to \code{lambda[j]}.}
#' \item{\code{edf}}{If GCV is \code{TRUE}, a scalar or vector with the trace of the smoothing matrix for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{stderr}}{If GCV is \code{TRUE}, a scalar or vector with the estimate of the standard deviation of the error for each value of the smoothing parameter specified in \code{lambda}.}
#' \item{\code{GCV}}{If GCV is \code{TRUE}, a  scalar or vector with the value of the GCV criterion for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a spatial regression model with differential regularization; isotropic and stationary case. In particular, the regularizing term involves the Laplacian of the spatial field. Space-varying covariates can be included in the model. The technique accurately handle data distributed over irregularly shaped domains. Moreover, various conditions can be imposed at the domain boundaries.
#' @usage smooth.FEM.basis(locations = NULL, observations, FEMbasis, lambda, 
#'        covariates = NULL, BC = NULL, GCV = FALSE)

#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
#' @examples
#' library(fdaPDE)
#' ## Load the Meuse data and a domain boundary for these data
#' data(MeuseData)
#' data(MeuseBorder)
#' ## Create a triangular mesh for these data with the provided boundary and plot it
#' order=1
#' mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = order)
#' plot(mesh)
#' ## Create the Finite Element basis 
#' FEMbasis = create.FEM.basis(mesh)
#' ## Estimate zync field without using covariates, setting the smoothing parameter to 10^3.5
#' data = log(MeuseData[,"zinc"])
#' lambda = 10^3.5
#' ZincMeuse = smooth.FEM.basis(observations = data, 
#'                              FEMbasis = FEMbasis, lambda = lambda)
#' ## Plot the estimated spatial field 
#' plot(ZincMeuse$fit.FEM)
#' # Now repeat the analysis using as covariates the square root of the log-distance 
#' # from river \code{sqrt(dist.log(m))} and the altitude \code{elev}
#' desmat = matrix(1,nrow=nrow(MeuseData),ncol=2)
#' desmat[,1] = sqrt(MeuseData[,"dist.log(m)"])
#' desmat[,2] = MeuseData[,"elev"]
#' ZincMeuseCovar = smooth.FEM.basis(observations = data, covariates = desmat, 
#'                                    FEMbasis = FEMbasis, lambda = lambda)
#' # Plot of the non parametric part (f) of the regression model y_i = beta_1 x_i1 + beta_2 x_i2 + f
#' plot(ZincMeuseCovar$fit.FEM)
#' # Print covariates' regression coefficients
#' print(ZincMeuseCovar$beta)


smooth.FEM.basis<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, PDE_parameters=NULL, incidence_matrix = NULL, BC = NULL, GCV = FALSE, GCVmethod = "Stochastic", nrealizations = 100)
{
  if(class(FEMbasis$mesh) == "MESH.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "MESH.2.5D"){
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "MESH.3D"){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }
  ##################### Checking parameters, sizes and conversion ################################
  
  if(GCVmethod=="Stochastic")
    GCVMETHOD=2
  else if(GCVmethod=="Exact")
    GCVMETHOD=1
  else{
    stop("GCVmethod must be either Stochastic or Exact")
  }
  
  space_varying=checkSmoothingParameters(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, BC, GCV, PDE_parameters, GCVMETHOD , nrealizations)
  
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
    covariates = as.matrix(covariates)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }
  
  # if I have PDE non-sv case I need (constant) matrices as parameters
  
  if(!is.null(PDE_parameters) & space_varying==FALSE) 
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }
  
  
  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, BC, GCV, space_varying, PDE_parameters,ndim, mydim)
  
  ################## End checking parameters, sizes and conversion #############################
  
  if(class(FEMbasis$mesh) == 'MESH.2D' & is.null(PDE_parameters)){	
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
                                  covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                  BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations)
  
    numnodes = nrow(FEMbasis$mesh$nodes)
    
  } else if(class(FEMbasis$mesh) == 'MESH.2D' & !is.null(PDE_parameters) & space_varying==FALSE){
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
                                      PDE_parameters = PDE_parameters,
                                      covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                      BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations)
    
    numnodes = nrow(FEMbasis$mesh$nodes)

  } else if(class(FEMbasis$mesh) == 'MESH.2D' & !is.null(PDE_parameters) & space_varying==TRUE){	
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.PDE.sv.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
                                         PDE_parameters = PDE_parameters,
                                         covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
                                         BC=BC, GCV=GCV,GCVMETHOD=GCVMETHOD, nrealizations=nrealizations)
  
    numnodes = nrow(FEMbasis$mesh$nodes)
  
  }else if(class(FEMbasis$mesh) == 'MESH.2.5D'){
    
    bigsol = NULL  
    print('C++ Code Execution')
    bigsol = CPP_smooth.manifold.FEM.basis(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, ndim, mydim, BC, GCV,GCVMETHOD, nrealizations)
    
    numnodes = FEMbasis$mesh$nnodes
    
  }else if(class(FEMbasis$mesh) == 'MESH.3D'){
    
    bigsol = NULL  
    print('C++ Code Execution')
    bigsol = CPP_smooth.volume.FEM.basis(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, ndim, mydim, BC, GCV,GCVMETHOD, nrealizations)
    
    numnodes = FEMbasis$mesh$nnodes
  }
  
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)
  
  reslist = NULL
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates, incidence_matrix, ndim, mydim)
  if(GCV == TRUE)
  {
    seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, incidence_matrix = incidence_matrix, edf = bigsol[[2]], ndim, mydim)
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
  }else{
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta)
  }
  
  return(reslist)
}

getBetaCoefficients<-function(locations, observations, fit.FEM, covariates, incidence_matrix = NULL, ndim, mydim)
{
  loc_nodes = NULL
  fnhat = NULL
  betahat = NULL
  
  if(!is.null(covariates))
  {
    if(is.null(locations))
    {
      loc_nodes = (1:length(observations))[!is.na(observations)]
      fnhat = as.matrix(fit.FEM$coeff[loc_nodes,])
    }else{
      loc_nodes = 1:length(observations)
      fnhat = eval.FEM(FEM = fit.FEM, locations = locations, incidence_matrix = incidence_matrix)
    }
    ## #row number of covariates, #col number of functions
    betahat = matrix(0, nrow = ncol(covariates), ncol = ncol(fnhat))
    for(i in 1:ncol(fnhat))
      betahat[,i] = as.vector(lm.fit(covariates,as.vector(observations-fnhat[,i]))$coefficients)
  }
  
  return(betahat)
}


getGCV<-function(locations, observations, fit.FEM, covariates = NULL, incidence_matrix = NULL, edf, ndim, mydim)
{
  loc_nodes = NULL
  fnhat = NULL
  
  edf = as.vector(edf)
  
  if(is.null(locations) && is.null(incidence_matrix))
  {
    loc_nodes = (1:length(observations))[!is.na(observations)]
    fnhat = as.matrix(fit.FEM$coeff[loc_nodes,])
  }else{
    loc_nodes = 1:length(observations)
    fnhat = eval.FEM(FEM = fit.FEM, locations = locations, incidence_matrix = incidence_matrix)
  }
  
  zhat = NULL
  zhat = matrix(nrow = length(loc_nodes), ncol = length(edf))
  if(!is.null(covariates))
  {
    desmatprod = ( solve( t(covariates) %*% covariates ) ) %*% t(covariates)
    for ( i in 1:length(edf))
    {
      betahat  = desmatprod %*% (observations-fnhat[,i])
      zhat[,i] = covariates %*% betahat + fnhat[,i]
    }
  }else{
    zhat = fnhat
  }
  
  np = length(loc_nodes)
  
  stderr2 = numeric(length(edf))
  GCV = numeric(length(edf))
  
  zhat <- as.matrix(zhat)
  
  if(any(np - edf <= 0))
  {
    warning("Some values of 'edf' are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.")  
  }
  
  for (i in 1:length(edf))
  {
    stderr2[i] = t(observations[loc_nodes] - zhat[,i]) %*% (observations[loc_nodes] - zhat[,i]) / ( np - edf[i] )
    GCV[i] = ( np / ( np - edf[i] )) * stderr2[i]
  }
  
  # NA if stderr2 is negative
  stderr = vector('numeric', length(stderr2));
  stderr[stderr2>=0] = sqrt(stderr2[stderr2>=0]);
  stderr[stderr2<0] = NaN;
  
  return(list(stderr = stderr, GCV = GCV))
}
