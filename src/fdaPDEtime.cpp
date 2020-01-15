
#define R_VERSION_

#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "solverdefinitions.h"
//#include <chrono>

#include "mixedFERegression.h"

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
SEXP regression_skeleton(InputHandler &regressionData, SEXP Rmesh, SEXP Rmesh_time)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	UInt n_time = Rf_length(Rmesh_time);
	std::vector<Real> mesh_time(n_time);
	for(UInt i=0; i<n_time; ++i)
	{
		mesh_time[i] = REAL(Rmesh_time)[i];
	}
	SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim> regression(mesh, mesh_time,regressionData);

	regression.apply();

	MatrixXv const & solution = regression.getSolution();
	MatrixXr const & dof = regression.getDOF();
	MatrixXr const & GCV = regression.getGCV();
	UInt bestLambdaS = regression.getBestLambdaS();
	UInt bestLambdaT = regression.getBestLambdaT();
	MatrixXv beta;
	if(regressionData.getCovariates().rows()==0)
	{
		beta.resize(1,1);
		beta(0,0).resize(1);
		beta(0,0)(0) = 10e20;
	}
	else
		 beta = regression.getBeta();

	//Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution(0,0).size(), solution.rows()*solution.cols()));
	SET_VECTOR_ELT(result, 1, Rf_allocMatrix(REALSXP, dof.rows(), dof.cols()));
	SET_VECTOR_ELT(result, 2, Rf_allocMatrix(REALSXP, GCV.rows(), GCV.cols()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(INTSXP, 2));
	SET_VECTOR_ELT(result, 4, Rf_allocMatrix(REALSXP, beta(0,0).size(), beta.rows()*beta.cols()));

	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt i = 0; i < solution.rows(); i++)
	{
		for(UInt j = 0; j < solution.cols(); j++)
		{
			for(UInt k = 0; k < solution(0,0).size(); k++)
				rans[k + solution(0,0).size()*i + solution(0,0).size()*solution.rows()*j] = solution.coeff(i,j)(k);
		}
	}

	Real *rans2 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < dof.rows(); i++)
	{
		for(UInt j = 0; j < dof.cols(); j++)
		{
		rans2[i + dof.rows()*j] = dof.coeff(i,j);
		}
	}

	Real *rans3 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < GCV.rows(); i++)
	{
		for(UInt j = 0; j < GCV.cols(); j++)
		{
		rans3[i + GCV.rows()*j] = GCV.coeff(i,j);
		}
	}

	UInt *rans4 = INTEGER(VECTOR_ELT(result, 3));
	rans4[0] = bestLambdaS;
	rans4[1] = bestLambdaT;

	Real *rans5 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < beta.rows(); i++)
	{
		for(UInt j = 0; j < beta.cols(); j++)
		{
			for(UInt k = 0; k < beta(0,0).size(); k++)
				rans5[k + beta(0,0).size()*i + beta(0,0).size()*beta.rows()*j] = beta.coeff(i,j)(k);
		}
	}
	UNPROTECT(1);
	return(result);
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
SEXP get_integration_points_skeleton(SEXP Rmesh)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);
	FiniteElement<Integrator,ORDER, mydim, ndim> fe;

	SEXP result;
	PROTECT(result=Rf_allocVector(REALSXP, 2*Integrator::NNODES*mesh.num_elements()));
	for(UInt i=0; i<mesh.num_elements(); i++)
	{
		fe.updateElement(mesh.getElement(i));
		for(UInt l = 0;l < Integrator::NNODES; l++)
		{
			Point p = fe.coorQuadPt(l);
			REAL(result)[i*Integrator::NNODES + l] = p[0];
			REAL(result)[mesh.num_elements()*Integrator::NNODES + i*Integrator::NNODES + l] = p[1];
		}
	}

	UNPROTECT(1);
	return(result);
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim, typename A>
SEXP get_FEM_Matrix_skeleton(SEXP Rmesh, EOExpr<A> oper)
{
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	SpMat AMat;
	Assembler::operKernel(oper, mesh, fe, AMat);

	//Copy result in R memory
	SEXP result;
	result = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(INTSXP, AMat.nonZeros() , 2));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, AMat.nonZeros()));

	int *rans = INTEGER(VECTOR_ELT(result, 0));
	Real  *rans2 = REAL(VECTOR_ELT(result, 1));
	UInt i = 0;
	for (UInt k=0; k < AMat.outerSize(); ++k)
		{
			for (SpMat::InnerIterator it(AMat,k); it; ++it)
			{
				//std::cout << "(" << it.row() <<","<< it.col() <<","<< it.value() <<")\n";
				rans[i] = 1+it.row();
				rans[i + AMat.nonZeros()] = 1+it.col();
				rans2[i] = it.value();
				i++;
			}
		}
	UNPROTECT(1);
	return(result);
}

extern "C" {

//! This function manages the various options for Spatial Regression, Sangalli et al version
/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_Laplace(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,
					SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations)
{
    //Set input data
	RegressionDataTime regressionData(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(regression_skeleton<RegressionDataTime,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(regression_skeleton<RegressionDataTime,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
    else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
		return(regression_skeleton<RegressionDataTime,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
   else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
		return(regression_skeleton<RegressionDataTime,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
		return(regression_skeleton<RegressionDataTime,IntegratorTetrahedronP2, 1, IntegratorGaussP5, 3, 2, 3, 3>(regressionData, Rmesh, Rmesh_time));
    return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R-matrix representing the diffusivity matrix of the model
	\param Rbeta an R-vector representing the advection term of the model
	\param Rc an R-double representing the reaction term of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\return R-vector containg the coefficients of the solution
*/

SEXP regression_PDE(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations)
{
	RegressionDataTimeElliptic regressionData(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton<RegressionDataTimeElliptic,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton<RegressionDataTimeElliptic,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton<RegressionDataTimeElliptic,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton<RegressionDataTimeElliptic,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	return(NILSXP);
}

/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param RK an R object representing the diffusivity tensor of the model
	\param Rbeta an R object representing the advection function of the model
	\param Rc an R object representing the reaction function of the model
	\param Ru an R object representing the forcing function of the model
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\return R-vector containg the coefficients of the solution
*/


SEXP regression_PDE_space_varying(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations)
{
    //Set data
	RegressionDataTimeEllipticSpaceVarying regressionData(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations);

	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	if(regressionData.getOrder() == 1 && ndim==2)
		return(regression_skeleton<RegressionDataTimeEllipticSpaceVarying,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==2)
		return(regression_skeleton<RegressionDataTimeEllipticSpaceVarying,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 2>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 1 && ndim==3)
		return(regression_skeleton<RegressionDataTimeEllipticSpaceVarying,IntegratorTriangleP2, 1, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	else if(regressionData.getOrder() == 2 && ndim==3)
		return(regression_skeleton<RegressionDataTimeEllipticSpaceVarying,IntegratorTriangleP4, 2, IntegratorGaussP5, 3, 2, 2, 3>(regressionData, Rmesh, Rmesh_time));
	return(NILSXP);
}

//! A function required for anysotropic and nonstationary regression (only 2D)
/*!
    \return points where the PDE space-varying params are evaluated in the R code
*/
SEXP get_integration_points(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	//Declare pointer to access data from C++
	int order = INTEGER(Rorder)[0];

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];
//Not implemented for ndim==3
    if(order == 1 && ndim ==2)
    	return(get_integration_points_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh));
    else if(order == 2 && ndim==2)
    	return(get_integration_points_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh));
    return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging

SEXP get_FEM_mass_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	int order = INTEGER(Rorder)[0];

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);

    if(order==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, mass));
	if(order==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, mass));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_stiff_matrix(SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim)
{
	int order = INTEGER(Rorder)[0];

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Stiff> ETMass;   Stiff EStiff;   ETMass stiff(EStiff);

    if(order==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, stiff));
	if(order==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, stiff));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_PDE_matrix(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations)
{
	RegressionDataTimeElliptic regressionData(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, RK, Rbeta, Rc, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations);

	//Get mydim and ndim
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Real& c = regressionData.getC();
	const Eigen::Matrix<Real,2,2>& K = regressionData.getK();
	const Eigen::Matrix<Real,2,1>& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

//! A utility, not used for system solution, may be used for debugging
SEXP get_FEM_PDE_space_varying_matrix(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rmesh, SEXP Rmesh_time, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
					SEXP RlambdaS, SEXP RlambdaT, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix,
					SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations)
{
	RegressionDataTimeEllipticSpaceVarying regressionData(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations);

	//Get mydim and ndim
	//UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

	typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	const Reaction& c = regressionData.getC();
	const Diffusivity& K = regressionData.getK();
	const Advection& beta = regressionData.getBeta();

    if(regressionData.getOrder()==1 && ndim==2)
    	return(get_FEM_Matrix_skeleton<IntegratorTriangleP2, 1,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	if(regressionData.getOrder()==2 && ndim==2)
		return(get_FEM_Matrix_skeleton<IntegratorTriangleP4, 2,2,2>(Rmesh, c*mass+stiff[K]+dot(beta,grad)));
	return(NILSXP);
}

}
