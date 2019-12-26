/*
 * FEMeval.cpp
 *
 *  Created on: Aug 16, 2015
 *      Author: eardi
 */


#define R_VERSION_

#include "fdaPDE.h"
//#include "IO_handler.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "evaluator.h"
#include "spline.h"



SEXP CPP_eval_FEM_fd(SEXP Rmesh, double* X,  double* Y,  double* Z, UInt n_X, UInt** incidenceMatrix, UInt nRegions, UInt nElements, double* coef, UInt order, UInt fast, UInt mydim, UInt ndim)
{


    SEXP result;

	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocVector(REALSXP, n_X));
		std::vector<bool> isinside(n_X);
		//Set the mesh
		//std::cout<<"Length "<<n_X<<"--X0 "<<X[0]<<"--Y0 "<<Y[0];
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			//std::cout<<"Starting evaluation from FEMeval \n";
			evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh);
			//mesh.printTriangles(std::cout);
			//mesh.printPoints(std::cout);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh);
			//mesh.printTriangles(std::cout);
			//mesh.printPoints(std::cout);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
		}

		for (int i=0; i<n_X; ++i)
		{
			if(!(isinside[i]))
			{
				REAL(result)[i]=NA_REAL;
			}
		}
	}
	else //areal data
	{
		PROTECT(result = Rf_allocVector(REALSXP, nRegions));
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}
	}

	UNPROTECT(1);
    // result list
  return(result);
}

extern "C" {
//! This function manages the various option for the solution evaluation.
/*!
	This function is then called from R code.
	Calls the walking algoritm for efficient point location inside the mesh in 2D.

	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rlocations an R-matrix (seen as an array) containing the xyz coordinates of the points where the solution has to be evaluated
	\param RincidenceMatrix an R-matrix for the incidence matrix defining the regions in the case of areal data
	\param Rcoef an R-vector the coefficients of the solution
	\param Rorder an R integer containg the order of the solution
	\param Rfast an R integer 0 for Naive location algorithm, 1 for Walking Algorithm (can miss location for non convex meshes)
*/


SEXP eval_FEM_fd(SEXP Rmesh, SEXP Rlocations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rmydim, SEXP Rndim)
{
	int n_X = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	int nRegions = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	int nElements = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1]; //number of triangles/tetrahedron if areal data
	//Declare pointer to access data from C++
	double *X, *Y, *Z;
	UInt **incidenceMatrix;
	double *coef;
	int order, mydim, ndim;
	bool fast;

	coef  = REAL(Rcoef);
  order = INTEGER(Rorder)[0];
  fast  = INTEGER(Rfast)[0];
  mydim = INTEGER(Rmydim)[0];
  ndim  = INTEGER(Rndim)[0];

	X = (double*) malloc(sizeof(double)*n_X);
	Y = (double*) malloc(sizeof(double)*n_X);
	Z = (double*) malloc(sizeof(double)*n_X);
	incidenceMatrix = (UInt**) malloc(sizeof(UInt*)*nRegions);

    // Cast all computation parameters
	if (ndim==3)
	{
		for (int i=0; i<n_X; i++)
		{
			X[i] = REAL(Rlocations)[i + n_X*0];
			//Rprintf("X[%i]= %d", i, X[i]);
			Y[i] = REAL(Rlocations)[i + n_X*1];
			Z[i] = REAL(Rlocations)[i + n_X*2];
		}
	}
	else //ndim==2
	{
		for (int i=0; i<n_X; i++)
		{
			X[i] = REAL(Rlocations)[i + n_X*0];
			Y[i] = REAL(Rlocations)[i + n_X*1];
			Z[i] = 0;
		}
	}
	for (int i=0; i<nRegions; i++)
	{
		incidenceMatrix[i] = (UInt*) malloc(sizeof(UInt)*nElements);
		for (int j=0; j<nElements; j++)
		{
			incidenceMatrix[i][j] = INTEGER(RincidenceMatrix)[i+nRegions*j];
		}
	}

    SEXP result;

	if (n_X>0) //pointwise data
	{
		PROTECT(result = Rf_allocVector(REALSXP, n_X));
		std::vector<bool> isinside(n_X);
		//Set the mesh
		//std::cout<<"Length "<<n_X<<"--X0 "<<X[0]<<"--Y0 "<<Y[0];
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			//std::cout<<"Starting evaluation from FEMeval \n";
			evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.eval(X, Y, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh);
			//mesh.printTriangles(std::cout);
			//mesh.printPoints(std::cout);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh);
			//mesh.printTriangles(std::cout);
			//mesh.printPoints(std::cout);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.eval(X, Y, Z, n_X, coef, fast, REAL(result), isinside);
		}

		for (int i=0; i<n_X; ++i)
		{
			if(!(isinside[i]))
			{
				REAL(result)[i]=NA_REAL;
			}
		}
	}
	else //areal data
	{
		PROTECT(result = Rf_allocVector(REALSXP, nRegions));
		if(order==1 && mydim==2 && ndim==2)
		{
			MeshHandler<1,2,2> mesh(Rmesh);
			Evaluator<1,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==2)
		{
			MeshHandler<2,2,2> mesh(Rmesh);
			Evaluator<2,2,2> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==2 && ndim==3)
		{
			MeshHandler<1,2,3> mesh(Rmesh);
			Evaluator<1,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==2 && mydim==2 && ndim==3)
		{
			MeshHandler<2,2,3> mesh(Rmesh);
			Evaluator<2,2,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));
		}
		else if(order==1 && mydim==3 && ndim==3)
		{
			MeshHandler<1,3,3> mesh(Rmesh);
			Evaluator<1,3,3> evaluator(mesh);
			evaluator.integrate(incidenceMatrix, nRegions, nElements, coef, REAL(result));

		}
	}

	free(X); free(Y); free(Z);
	for (int i=0; i<nRegions; i++)
	{
		free(incidenceMatrix[i]);
	}
	free(incidenceMatrix);


	UNPROTECT(1);
    // result list
    return(result);
}

SEXP eval_FEM_time(SEXP Rmesh, SEXP Rmesh_time, SEXP Rlocations, SEXP Rtime_locations, SEXP RincidenceMatrix, SEXP Rcoef, SEXP Rorder, SEXP Rfast, SEXP Rflag_parabolic, SEXP Rmydim, SEXP Rndim)
{
	UInt n = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	UInt ns = INTEGER(Rf_getAttrib(VECTOR_ELT(Rmesh, 0), R_DimSymbol))[0];
	UInt nt = Rf_length(Rmesh_time);
	UInt nRegions = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt nElements = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1]; //number of triangles/tetrahedron if areal data


	//Declare pointer to access data from C++
	Real *X, *Y, *Z, *mesh_time, *t;
	UInt **incidenceMatrix;
	double *coef;
	int order, mydim, ndim;
	bool fast,flag_par;

	coef  = REAL(Rcoef);
  order = INTEGER(Rorder)[0];
  fast  = INTEGER(Rfast)[0];
  mydim = INTEGER(Rmydim)[0];
  ndim  = INTEGER(Rndim)[0];
	flag_par = INTEGER(Rflag_parabolic)[0];
	mesh_time = REAL(Rmesh_time);
	t = REAL(Rtime_locations);

	X = (double*) malloc(sizeof(double)*n);
	Y = (double*) malloc(sizeof(double)*n);
	Z = (double*) malloc(sizeof(double)*n);
	incidenceMatrix = (UInt**) malloc(sizeof(UInt*)*nRegions);

    // Cast all computation parameters
	if (ndim==3)
	{
		for (int i=0; i<n; i++)
		{
			X[i] = REAL(Rlocations)[i + n*0];
			//Rprintf("X[%i]= %d", i, X[i]);
			Y[i] = REAL(Rlocations)[i + n*1];
			Z[i] = REAL(Rlocations)[i + n*2];
		}
	}
	else //ndim==2
	{
		for (int i=0; i<n; i++)
		{
			X[i] = REAL(Rlocations)[i + n*0];
			Y[i] = REAL(Rlocations)[i + n*1];
			Z[i] = 0;
		}
	}
	for (int i=0; i<nRegions; i++)
	{
		incidenceMatrix[i] = (UInt*) malloc(sizeof(UInt)*nElements);
		for (int j=0; j<nElements; j++)
		{
			incidenceMatrix[i][j] = INTEGER(RincidenceMatrix)[i+nRegions*j];
		}
	}


	UInt DEGREE = flag_par ? 1 : 3;
	UInt M = nt + DEGREE - 1;
	SpMat phi(n,M);

	if(flag_par)
	{
		Spline<IntegratorGaussP5,1,0>spline(mesh_time,nt);
		Real value;
		for (UInt i = 0; i < n; ++i)
		{
			for (UInt j = 0; j < M; ++j)
			{
				value = spline.BasisFunction(DEGREE, j, t[i]);
				if (value!=0)
				{
					phi.coeffRef(i,j) = value;
				}
			}
		}
	}
	else
	{
		Spline<IntegratorGaussP5,3,2>spline(mesh_time,nt);
		Real value;
		for (UInt i = 0; i < n; ++i)
		{
			for (UInt j = 0; j < M; ++j)
			{
				value = spline.BasisFunction(DEGREE, j, t[i]);
				if (value!=0)
				{
					phi.coeffRef(i,j) = value;
				}
			}
		}
	}
	phi.makeCompressed();

	SEXP result;

	PROTECT(result=Rf_allocVector(REALSXP, n));
	Real* COEFF;
	COEFF = (double*) malloc(sizeof(double)*ns);
	std::vector<Real> XX,YY,ZZ;
	std::vector<UInt> indices;

	for(UInt j=0; j<ns; ++j)
	{
		COEFF[j] = coef[j];
	}
	// for(UInt i=0; i<n; i++)
	// {
	// 	if(phi.coeff(i,0)!=0)
	// 	{
	// 		if (ndim==3)
	// 		{
	// 			XX.push_back(X[i]);
	// 			YY.push_back(Y[i]);
	// 			ZZ.push_back(Z[i]);
	// 			indices.push_back(i);
	// 		}
	// 		else //ndim==2
	// 		{
	// 			XX.push_back(X[i]);
	// 			YY.push_back(Y[i]);
	// 			indices.push_back(i);
	// 		}
	// 	}
	// }
	SEXP temp = CPP_eval_FEM_fd(Rmesh, X, Y, Z, n, incidenceMatrix, nRegions, nElements, COEFF, order, fast, mydim, ndim);
	for(UInt k=0; k < n; k++)
	{
		REAL(result)[k] = REAL(temp)[k];
		if(!ISNA(REAL(result)[k]))
			REAL(result)[k] = REAL(result)[k]*phi.coeff(k,0);
	}

	for(UInt i=1; i<M; ++i)
	{
		for(UInt j=0; j<ns; ++j)
		{
			COEFF[j] = coef[i*ns+j];
		}
		for(UInt k=0; k<n; k++)
		{
			if(phi.coeff(k,i)!=0)
			{
				if (ndim==3)
				{
					XX.push_back(X[k]);
					YY.push_back(Y[k]);
					ZZ.push_back(Z[k]);
					indices.push_back(k);
				}
				else //ndim==2
				{
					XX.push_back(X[k]);
					YY.push_back(Y[k]);
					indices.push_back(k);
				}
			}
		}
		temp = CPP_eval_FEM_fd(Rmesh, XX.data(), YY.data(), ZZ.data(), XX.size(), incidenceMatrix, nRegions, nElements, COEFF, order, fast, mydim, ndim);
		for(UInt k=0; k<indices.size(); ++k)
		{
			if(!ISNA(REAL(result)[indices[k]]))
				REAL(result)[indices[k]] = REAL(result)[indices[k]] + REAL(temp)[k]*phi.coeff(indices[k],i);
		}
		XX.clear();YY.clear();ZZ.clear();indices.clear();
	}

	free(X); free(Y); free(Z);
	for (int i=0; i<nRegions; i++)
	{
		free(incidenceMatrix[i]);
	}
	free(incidenceMatrix);

	UNPROTECT(1);
	return(result);
}

SEXP eval_FEM_time_nodes(SEXP Rns, SEXP Rmesh_time, SEXP Rtime, SEXP Rcoef, SEXP Rflag_parabolic)
{
	UInt ns = INTEGER(Rns)[0];
	UInt nt = Rf_length(Rmesh_time);

	Real *mesh_time = REAL(Rmesh_time);
	Real t = REAL(Rtime)[0];
	bool flag_par = INTEGER(Rflag_parabolic)[0];

	UInt DEGREE = flag_par ? 1 : 3;
	UInt M = nt + DEGREE - 1;
	VectorXr phi(M);

	if(flag_par)
	{
		Spline<IntegratorGaussP5,1,0>spline(mesh_time,nt);
		for (UInt j = 0; j < M; ++j)
			{
				phi(j) = spline.BasisFunction(DEGREE, j, t);
			}
	}
	else
	{
		Spline<IntegratorGaussP5,3,2>spline(mesh_time,nt);
		for (UInt j = 0; j < M; ++j)
			{
				phi(j) = spline.BasisFunction(DEGREE, j, t);
			}
	}

	SEXP result;

	PROTECT(result=Rf_allocVector(REALSXP, ns));

	for(UInt k=0; k<ns; ++k)
	{
		REAL(result)[k] = REAL(Rcoef)[k]*phi(0);
	}
	for(UInt i=1; i < M; i++)
	{
		for(UInt k=0; k<ns; ++k)
		{
			REAL(result)[k] = REAL(result)[k] + REAL(Rcoef)[k+ns*i]*phi(i);
		}
	}

	UNPROTECT(1);
	return(result);
}

}
