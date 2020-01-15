CPP_smooth.manifold.FEM.basis<-function(locations, time_locations, observations, FEMbasis,time_mesh, lambdaS, lambdaT,
                                    covariates, incidence_matrix, ndim, mydim, BC, FLAG_MASS, FLAG_PARABOLIC, IC,GCV, GCVMETHOD, nrealizations)
{

  # C++ function for manifold works with vectors not with matrices

  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1


  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(IC))
  {
    IC<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(time_locations))
  {
    if(FLAG_PARABOLIC==FALSE)
    {
      time_locations<-time_mesh
    }
    else
    {
      NumTimeInstants = length(time_mesh)
      time_locations<-time_mesh[1:NumTimeInstants-1]
    }
  }

  if(is.null(time_mesh))
  {
    if(FLAG_PARABOLIC==FALSE)
    {
      time_mesh<-time_locations
    }
    else
    {
      NumTimeInstants = length(time_locations)
      lastTimeDiff = time_locations[NumTimeInstants]-time_locations[NumTimeInstants-1]
      time_mesh<-c(time_locations,time_locations[NumTimeInstants]+lastTimeDiff)
    }
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1

  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }

  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  data <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntriangles) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  GCV <- as.integer(GCV)
  storage.mode(GCV) <- "integer"

  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <-"integer"

  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"

  ## Call C++ function


  bigsol <- .Call("regression_Laplace", locations, time_locations, data, FEMbasis$mesh, time_mesh, FEMbasis$order,
                  mydim, ndim, lambdaS, lambdaT, covariates, incidence_matrix, BC$BC_indices, BC$BC_values, FLAG_MASS, FLAG_PARABOLIC,
                  IC, GCV, GCVMETHOD, nrealizations, PACKAGE = "fdaPDEtime")

  return(bigsol)
}

CPP_eval.manifold.FEM = function(FEM_time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim)
{
  FEMbasis = FEM_time$FEMbasis
  # C++ function for manifold works with vectors not with matrices

  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  #NEW: riporto shift indici in R
  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(FEM_time$mesh_time) <- "double"
  coeff <- as.matrix(FEM_time$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  storage.mode(FLAG_PARABOLIC) <- "integer"

  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_time", FEMbasis$mesh, FEM_time$mesh_time, locations, time_locations, incidence_matrix, coeff[,i],
                         FEMbasis$order, redundancy, FLAG_PARABOLIC, mydim, ndim, package = "fdaPDEtime")
  }

  #Returning the evaluation matrix
  evalmat
}
