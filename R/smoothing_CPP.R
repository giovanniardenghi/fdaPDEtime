CPP_smooth.FEM.basis<-function(locations, time_locations, observations, FEMbasis, time_mesh, lambdaS, lambdaT, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, FLAG_MASS, FLAG_PARABOLIC, IC, GCV ,GCVMETHOD = 2, nrealizations = 100, DOF=TRUE,DOF_matrix=NULL)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF_matrix))
  {
    DOF_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
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
    if(FLAG_PARABOLIC)
      time_locations <- time_mesh[2:length(time_mesh)]
    else
      time_locations <- time_mesh
  }

  if(is.null(time_mesh))
  {
    time_mesh<-time_locations
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

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  DOF_matrix <- as.matrix(DOF_matrix)
  storage.mode(DOF_matrix) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"

  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"

  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <-"integer"

  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
                  mydim, ndim, lambdaS, lambdaT, covariates, incidence_matrix, BC$BC_indices, BC$BC_values, FLAG_MASS, FLAG_PARABOLIC,
                  IC, GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, PACKAGE = "fdaPDEtime")
  return(bigsol)
}

CPP_smooth.FEM.PDE.basis<-function(locations, time_locations, observations, FEMbasis, time_mesh, lambdaS, lambdaT, PDE_parameters, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, FLAG_MASS, FLAG_PARABOLIC, IC, GCV,GCVMETHOD = 2, nrealizations = 100, DOF=TRUE,DOF_matrix=NULL)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }


  if(is.null(DOF_matrix))
  {
    DOF_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
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
    if(FLAG_PARABOLIC)
      time_locations <- time_mesh[2:length(time_mesh)]
    else
      time_locations <- time_mesh
  }

  if(is.null(time_mesh))
  {
    time_mesh<-time_locations
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

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  DOF_matrix <- as.matrix(DOF_matrix)
  storage.mode(DOF_matrix) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"

  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"

  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <-"integer"

  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <-"integer"

  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_PDE", locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
                  mydim, ndim, lambdaS, lambdaT, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c,
                  covariates, incidence_matrix, BC$BC_indices, BC$BC_values, FLAG_MASS, FLAG_PARABOLIC,
                  IC, GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, PACKAGE = "fdaPDEtime")
  return(bigsol)
}

CPP_smooth.FEM.PDE.sv.basis<-function(locations, time_locations, observations, FEMbasis, time_mesh, lambdaS, lambdaT, PDE_parameters, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, FLAG_MASS, FLAG_PARABOLIC, IC, GCV,GCVMETHOD = 2, nrealizations = 100, DOF=TRUE,DOF_matrix=NULL)
{

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF_matrix))
  {
    DOF_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
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
    if(FLAG_PARABOLIC)
      time_locations <- time_mesh[2:length(time_mesh)]
    else
      time_locations <- time_mesh
  }

  if(is.null(time_mesh))
  {
    if(FLAG_PARABOLIC)
      time_mesh <- rbind(0,time_locations)
    else
      time_mesh<-time_locations
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

  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = FEMbasis$mesh, order = FEMbasis$order),ncol = 2)
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$b = (PDE_parameters$b)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  time_locations <- as.matrix(time_locations)
  storage.mode(time_locations) <- "double"
  time_mesh <- as.matrix(time_mesh)
  storage.mode(time_mesh) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  DOF_matrix <- as.matrix(DOF_matrix)
  storage.mode(DOF_matrix) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambdaS) <- "double"
  storage.mode(lambdaT) <- "double"
  IC <- as.matrix(IC)
  storage.mode(IC) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"

  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"

  FLAG_MASS <- as.integer(FLAG_MASS)
  storage.mode(FLAG_MASS) <-"integer"

  FLAG_PARABOLIC <- as.integer(FLAG_PARABOLIC)
  storage.mode(FLAG_PARABOLIC) <-"integer"

  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"

  ## Call C++ function
  bigsol <- .Call("regression_PDE_space_varying", locations, time_locations, observations, FEMbasis$mesh, time_mesh, FEMbasis$order,
                  mydim, ndim, lambdaS, lambdaT, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u,
                  covariates, incidence_matrix, BC$BC_indices, BC$BC_values, FLAG_MASS, FLAG_PARABOLIC,
                  IC, GCV, GCVMETHOD, nrealizations, DOF, DOF_matrix, PACKAGE = "fdaPDEtime")
  return(bigsol)
}

CPP_eval.FEM = function(FEM, locations, incidence_matrix, redundancy, ndim, mydim)
{

  # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
  #
  #        arguments:
  # X         an array of x-coordinates.
  # Y         an array of y-coordinates.
  # FELSPLOBJ a FELspline object
  # FAST      a boolean indicating if the walking algorithm should be apply
  #        output:
  # EVALMAT   an array of the same size as X and Y containing the value of
  #           FELSPLOBJ at (X,Y).

  FEMbasis = FEM$FEMbasis
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff <- as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"

  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations, incidence_matrix, coeff[,i],
                         FEMbasis$order, redundancy, mydim, ndim, package = "fdaPDE")
  }

  #Returning the evaluation matrix
  evalmat
}

CPP_eval.FEM_time <- function(FEM_time, locations, time_locations, incidence_matrix, FLAG_PARABOLIC, redundancy, ndim, mydim)
{

  # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
  #
  #        arguments:
  # X         an array of x-coordinates.
  # Y         an array of y-coordinates.
  # FELSPLOBJ a FELspline object
  # FAST      a boolean indicating if the walking algorithm should be apply
  #        output:
  # EVALMAT   an array of the same size as X and Y containing the value of
  #           FELSPLOBJ at (X,Y).

  FEMbasis = FEM_time$FEMbasis
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

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
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
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


CPP_get_evaluations_points = function(mesh, order)
{
  #here we do not shift indices since this function is called inside CPP_smooth.FEM.PDE.sv.basis

  # Imposing types, this is necessary for correct reading from C++
  if(class(mesh)=="MESH.2D"){
    ndim=2
    mydim=2
  }else if(class(mesh) == "MESH.2.5D" || class(mesh) == "MESH.3D"){
    stop('Function not yet implemented for this mesh class')
  }else{
    stop('Unknown mesh class')
  }

  storage.mode(ndim)<-"integer"
  storage.mode(mydim)<-"integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  storage.mode(mesh$edges) <- "integer"
  storage.mode(mesh$neighbors) <- "integer"
  storage.mode(order) <- "integer"


  points <- .Call("get_integration_points",mesh, order,mydim, ndim,
                  PACKAGE = "fdaPDEtime")

  #Returning the evaluation matrix
  points
}

CPP_get.FEM.Mass.Matrix<-function(FEMbasis)
{
  if(class(FEMbasis$mesh) == "MESH.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(mesh) == "MESH.3D"){
    stop('Function not yet implemented for this mesh class')
  }else{
    stop('Unknown mesh class')
  }


  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  ## Set propr type for correct C++ reading
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(ndim)<-"integer"
  storage.mode(mydim)<-"integer"

  ## Call C++ function
  triplets <- .Call("get_FEM_mass_matrix", FEMbasis$mesh,
                    FEMbasis$order,mydim, ndim,
                    PACKAGE = "fdaPDE")

  A = sparseMatrix(i = triplets[[1]][,1], j=triplets[[1]][,2], x = triplets[[2]], dims = c(nrow(FEMbasis$mesh$nodes),nrow(FEMbasis$mesh$nodes)))
  return(A)
}

CPP_get.FEM.Stiff.Matrix<-function(FEMbasis)
{
  if(class(FEMbasis$mesh) == "MESH.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(mesh) == "MESH.3D"){
    stop('Function not yet implemented for this mesh class')
  }else{
    stop('Unknown mesh class')
  }

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  ## Set propr type for correct C++ reading
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(ndim)<-"integer"
  storage.mode(mydim)<-"integer"

  ## Call C++ function
  triplets <- .Call("get_FEM_stiff_matrix", FEMbasis$mesh,
                    FEMbasis$order, mydim, ndim,
                    PACKAGE = "fdaPDE")

  A = sparseMatrix(i = triplets[[1]][,1], j=triplets[[1]][,2], x = triplets[[2]], dims = c(nrow(FEMbasis$mesh$nodes),nrow(FEMbasis$mesh$nodes)))
  return(A)
}

CPP_get.FEM.PDE.Matrix<-function(FEMbasis, PDE_parameters)
{
  if(class(FEMbasis$mesh) == "MESH.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(mesh) == "MESH.3D"){
    stop('Function not yet implemented for this mesh class')
  }else{
    stop('Unknown mesh class')
  }
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  covariates<-matrix(nrow = 0, ncol = 1)
  locations<-matrix(nrow = 0, ncol = 2)
  incidence_matrix<-matrix(nrow = 0, ncol = 1)
  BC$BC_indices<-vector(length=0)
  BC$BC_values<-vector(length=0)
  lambda = 0
  GCV = 0
  GCVmethod = 0
  nrealizations = 0

  ## Set proper type for correct C++ reading

  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  storage.mode(GCV) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"

  storage.mode(PDE_parameters$K) <- "double"
  storage.mode(PDE_parameters$b) <- "double"
  storage.mode(PDE_parameters$c) <- "double"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVmethod) <- "integer"

  ## Call C++ function
  triplets <- .Call("get_FEM_PDE_matrix", locations, observations, FEMbasis$mesh,
                    FEMbasis$order,mydim, ndim, lambda, PDE_parameters$K, PDE_parameters$b, PDE_parameters$c, covariates,
                    incidence_matrix, BC$BC_indices, BC$BC_values, GCV,GCVmethod, nrealizations,
                    PACKAGE = "fdaPDE")

  A = sparseMatrix(i = triplets[[1]][,1], j=triplets[[1]][,2], x = triplets[[2]], dims = c(nrow(FEMbasis$mesh$nodes),nrow(FEMbasis$mesh$nodes)))
  return(A)
}


CPP_get.FEM.PDE.sv.Matrix<-function(FEMbasis, PDE_parameters)
{

  if(class(FEMbasis$mesh) == "MESH.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(mesh) == "MESH.3D"){
    stop('Function not yet implemented for this mesh class')
  }else{
    stop('Unknown mesh class')
  }

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  covariates<-matrix(nrow = 0, ncol = 1)
  locations<-matrix(nrow = 0, ncol = 2)
  incidence_matrix<-matrix(nrow = 0, ncol = 1)
  BC$BC_indices<-vector(length=0)
  BC$BC_values<-vector(length=0)
  lambda = 0
  GCV = 0
  GCVmethod = 0
  nrealizations = 0

  PDE_param_eval = NULL
  points_eval = matrix(CPP_get_evaluations_points(mesh = FEMbasis$mesh, order = FEMbasis$order),ncol = 2)
  PDE_param_eval$K = (PDE_parameters$K)(points_eval)
  PDE_param_eval$b = (PDE_parameters$b)(points_eval)
  PDE_param_eval$c = (PDE_parameters$c)(points_eval)
  PDE_param_eval$u = (PDE_parameters$u)(points_eval)

  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  storage.mode(GCV) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"

  storage.mode(PDE_param_eval$K) <- "double"
  storage.mode(PDE_param_eval$b) <- "double"
  storage.mode(PDE_param_eval$c) <- "double"
  storage.mode(PDE_param_eval$u) <- "double"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVmethod) <- "integer"

  ## Call C++ function
  triplets <- .Call("get_FEM_PDE_space_varying_matrix", locations, observations, FEMbasis$mesh,
                    FEMbasis$order,mydim, ndim, lambda, PDE_param_eval$K, PDE_param_eval$b, PDE_param_eval$c, PDE_param_eval$u, covariates,
                    incidence_matrix, BC$BC_indices, BC$BC_values, GCV,GCVmethod, nrealizations,
                    PACKAGE = "fdaPDE")

  A = sparseMatrix(i = triplets[[1]][,1], j=triplets[[1]][,2], x = triplets[[2]], dims = c(nrow(FEMbasis$mesh$nodes),nrow(FEMbasis$mesh$nodes)))
  return(A)
}
