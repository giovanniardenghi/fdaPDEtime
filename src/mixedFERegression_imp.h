#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

//#include <libseq/mpi.h>

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::addDirichletBC()
{
	UInt id1,id3;

	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bc_indices[i];
			id3=id1+N*M;

			matrixNoCov_.coeffRef(id1,id1)=pen;
			matrixNoCov_.coeffRef(id3,id3)=pen;


			_rightHandSide(id1)=bc_values[i]*pen;
			_rightHandSide(id3)=0;
	 }

	matrixNoCov_.makeCompressed();
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::addNA()
{

	UInt id;

	UInt nnodes = mesh_.num_nodes();
	//std::cout << "nnodes: " << nnodes << std::endl;

	const std::vector<UInt>& observations_na= regressionData_.getObservationsNA();

	UInt n_NA = observations_na.size();

	for(UInt i = 0; i < n_NA; ++i)
	{
		id = observations_na[i];
		matrixNoCov_.coeffRef(id, id) = 0;
	}
	//std::cout << _solution << std::endl;
	matrixNoCov_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setPsi(){

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofSpaceObservations();

	psi_.resize(nlocations, nnodes);
	if (regressionData_.isLocationsByNodes()) //pointwise data
	{
		psi_.resize(nnodes, nnodes);
		// std::vector<coeff> tripletAll;
		// auto k = regressionData_.getObservationsIndices();
		// tripletAll.reserve(k.size());
		// for (int i = 0; i< k.size(); ++i){
		// 	tripletAll.push_back(coeff(i,k[i],1.0));
		// }
		// psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		// psi_.makeCompressed();
		psi_.resize(nnodes,nnodes);
		psi_.setIdentity();
	}
	else if (regressionData_.getNumberOfRegions() == 0)
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;
		Element<Nodes, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,Nodes,1> coefficients;

		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else
			{
				for(UInt node = 0; node < Nodes ; ++node)
				{
					coefficients = Eigen::Matrix<Real,Nodes,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<Nodes,mydim,ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}
		psi_.makeCompressed();
	}
	else //areal data
	{
		constexpr UInt Nodes = mydim==2 ? 3*ORDER : 6*ORDER-2;

		Real *tab; //Psi_i
		tab = (Real*) malloc(sizeof(Real)*nnodes);
		for(UInt i=0; i<nlocations;i++) //nlocations = number of regions
		{
			for (UInt k=0; k<nnodes; k++) {tab[k]=0;}
			for (UInt j=0; j<mesh_.num_elements(); j++)
			{
				if (regressionData_.getIncidenceMatrix()(i,j) == 1) //element j is in region i
				{
					Element<Nodes, mydim, ndim> tri = mesh_.getElement(j);
					for (UInt k=0; k<Nodes; k++)
					{
						tab[tri[k].getId()] += integratePsi(tri,k); // integral over tri of psi_k
					}
				}
			}
			for (int k=0; k<nnodes; k++)
			{
				if (tab[k] != 0)
				{
					psi_.insert(i,k) = tab[k]/A_.coeffRef(i,i); //divide by |D_i|
				}
			}
		}
		free(tab);
		psi_.makeCompressed();
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
MatrixXr SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::LeftMultiplybyQ(const MatrixXr& u)
{
	if (regressionData_.getCovariates().rows() == 0){
		return u;
	}
	else{
		MatrixXr W(this->regressionData_.getCovariates());
		if (isWTWfactorized_ == false ){
			WTW_.compute(W.transpose()*W);
			isWTWfactorized_=true;
		}
		MatrixXr Pu= W*WTW_.solve(W.transpose()*u);
		return u-Pu;
	}

}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::buildMatrixNoCov(const SpMat& DMat,  const SpMat& SWblock,  const SpMat& SEblock)
{
	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;
	UInt nnodes = M*N;

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*SWblock.nonZeros() + SEblock.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for (int k=0; k<SEblock.outerSize(); ++k)
		for (SpMat::InnerIterator it(SEblock,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for (int k=0; k<SWblock.outerSize(); ++k)
	  for (SpMat::InnerIterator it(SWblock,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<SWblock.outerSize(); ++k)
	  for (SpMat::InnerIterator it(SWblock,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	matrixNoCov_.setZero();
	matrixNoCov_.resize(2*nnodes,2*nnodes);
	matrixNoCov_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	matrixNoCov_.makeCompressed();
	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
}

// template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
// void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::buildMatrixOnlyCov(const SpMat& B,  const MatrixXr& H)
// {
//
// 	UInt N = mesh_.num_nodes();
// 	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;
// 	UInt nnodes = M*N;
//
// 	MatrixXr NWblock= MatrixXr::Zero(2*nnodes,2*nnodes);
//
// 	if(regressionData_.getNumberOfRegions()==0)
//     	NWblock.topLeftCorner(nnodes,nnodes)=B.transpose()*(-H)*B;
//   else
// 	    NWblock.topLeftCorner(nnodes,nnodes)=B.transpose()*Ak_*(-H)*B;
//
//   matrixOnlyCov_=NWblock.sparseView();
// 	matrixOnlyCov_.makeCompressed();
// }

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::system_factorize()
{
	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;
	UInt nnodes = M*N;

	// First phase: Factorization of matrixNoCov
	// matrixNoCovdec_.compute(matrixNoCov_);

	if (regressionData_.getCovariates().rows() != 0)
	{
		// Second phase: factorization of matrix  G =  C + [V * matrixNoCov^-1 * U]= C + D

		// Definition of matrix U = [ psi^T * A * W | 0 ]^T and V= [ W^T*psi| 0]

		MatrixXr W(this->regressionData_.getCovariates());

		U_ = MatrixXr::Zero(2*nnodes, W.cols());

		V_ = MatrixXr::Zero(W.cols(),2*nnodes);
		V_.leftCols(nnodes)=W.transpose()*B_;

		if(regressionData_.getNumberOfRegions()==0)
		{ // pointwise data
		  U_.topRows(nnodes) = B_.transpose()*W;
		}
		else
		{                                          //areal data
		  U_.topRows(nnodes) = B_.transpose()*Ak_*W;
    }

		// MatrixXr D = V_*matrixNoCovdec_.solve(U_);
		MatrixXr temp(U_.rows(),U_.cols());
		Mumps::template solve(matrixNoCov_,U_,temp);
		MatrixXr D = V_*temp;

		// G = C + D
		MatrixXr G = -W.transpose()*W + D;
		Gdec_.compute(G);

	}
}

template<typename InputHandler,typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> &b) {

	// Resolution of the system matrixNoCov * x1 = b
	// MatrixXr x1 = matrixNoCovdec_.solve(b);
	MatrixXr x1(b.rows(),b.cols());
	Mumps::template solve(matrixNoCov_,b,x1);

	if (regressionData_.getCovariates().rows() != 0)
	{
		// Resolution of G * x2 = V * x1

		MatrixXr x2 = Gdec_.solve(V_*x1);

		// Resolution of the system matrixNoCov * x3 = U * x2
		MatrixXr xtemp(b.rows(),b.cols());
		Mumps::template solve(matrixNoCov_,U_*x2,xtemp);
		x1 -= xtemp;
		// x1 -= matrixNoCovdec_.solve(U_*x2);
	}

	return x1;
}

// template<typename InputHandler,typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
// void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::setQ()
// {
//  	//std::cout<<"Computing Orthogonal Space Projection Matrix"<<std::endl;
//  	Q_.resize(H_.rows(),H_.cols());
//  	Q_ = -H_;
//  	for (int i=0; i<H_.rows();++i)
//  	{
//  		Q_(i,i) += 1;
//  	}
//  }
//
// template<typename InputHandler,typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
// void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::setH()
// {
// 	MatrixXr W(this->regressionData_.getCovariates());
	// UInt nlocations = regressionData_.getNumberofObservations();
 	// if(regressionData_.isLocationsByNodes())
 	// {
 	// 	MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
 	// 	for (auto i=0; i<nlocations;++i)
 	// 	{
 	// 		auto index_i = regressionData_.getObservationsIndices()[i];
 	// 		for (auto j=0; j<W.cols();++j)
 	// 		{
 	// 			W_reduced(i,j) = W(index_i,j);
 	// 		}
 	// 	}
 	// 	W = W_reduced;
 	// }
 // 	MatrixXr WTW(W.transpose()*W);
 // 	H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix
 // }

// template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>addDirichletBC()
// void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, UInt mydim, UInt ndim>::setAk()
// {
// 	UInt nRegions = regressionData_.getNumberOfRegions();
// 	UInt M = mesh_time_.size();
//
// 	Ak_.resize(nRegions*M,1);
// 	for (int i=0; i<nRegions; i++)
// 	{
// 		Ak_(i)=0;
// 		for (int j=0; j<regressionData_.getIncidenceMatrix().cols(); j++)
// 		{
// 			if (regressionData_.getIncidenceMatrix()(i,j) == 1)
// 			{
// 				for (int k=0; k<M; k++)
// 				{
// 					Ak_(i+k*nRegions)+=mesh_.elementMeasure(j);
// 				}
// 			}
// 		}
// 	}
// }
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim, ndim>::setA()
{
	UInt nRegions = regressionData_.getNumberOfRegions();
	A_.resize(nRegions,nRegions);
	for (int i=0; i<nRegions; i++)
	{
		A_.coeffRef(i,i)=0;
		for (int j=0; j<regressionData_.getIncidenceMatrix().cols(); j++)
		{
			if (regressionData_.getIncidenceMatrix()(i,j) == 1)
			{
				A_.coeffRef(i,i)+=mesh_.elementMeasure(j);
			}
		}
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::getRightHandData(VectorXr& rightHandData)
{
	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;

	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(N*M);

	if (regressionData_.getCovariates().rows() == 0) //no covariate
	{
		if (regressionData_.isLocationsByNodes() && regressionData_.getFlagParabolic())
		{
			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else if (regressionData_.getNumberOfRegions() == 0) //pointwise data
		{
			rightHandData=B_.transpose()*regressionData_.getObservations();
		}
		else //areal data
		{
			rightHandData=B_.transpose()*Ak_*regressionData_.getObservations();
		}
	}
	else if (regressionData_.getNumberOfRegions() == 0) //with covariates, pointwise data
	{
		rightHandData=B_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
	else //with covariates, areal data
	{
		rightHandData=B_.transpose()*Ak_*LeftMultiplybyQ(regressionData_.getObservations());
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{
	int GCVmethod = regressionData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_indexS, output_indexT, lambdaS, lambdaT);
			break;
	}
	VectorXr dataHat;
	VectorXr z = regressionData_.getObservations();
	if(regressionData_.getCovariates().rows()==0)
		dataHat = B_*_solution(output_indexS,output_indexT).topRows(B_.cols());
	else
		dataHat = z - LeftMultiplybyQ(z) + LeftMultiplybyQ(B_*_solution(output_indexS,output_indexT).topRows(B_.cols()));
	UInt n = dataHat.rows();

	_GCV(output_indexS,output_indexT) = (n / ((n-_dof(output_indexS,output_indexT)) * (n-_dof(output_indexS,output_indexT)))) * (z-dataHat).dot(z-dataHat);
	if (_GCV(output_indexS,output_indexT) < _bestGCV)
	{
		bestLambdaS_ = output_indexS;
		bestLambdaT_ = output_indexT;
		_bestGCV = _GCV(output_indexS,output_indexT);
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{

	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;

	// UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees=0;
	SpMat X;
	SpMat BBsmall(M*N,M*N);
	if(regressionData_.getCovariates().rows() == 0)
	{
		BBsmall = B_.transpose()*B_;
	}
	else
	{
		BBsmall = (SpMat(B_.transpose())*LeftMultiplybyQ(B_)).sparseView();
	}

	SpMat BB(2*M*N,2*M*N);

	std::vector<coeff> tripletAll;
	tripletAll.reserve(BBsmall.nonZeros());

	for (int k=0; k<BBsmall.outerSize(); ++k)
		for (SpMat::InnerIterator it(BBsmall,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}

	BB.setFromTriplets(tripletAll.begin(),tripletAll.end());
	BB.makeCompressed();

	Real *values = matrixNoCov_.valuePtr();
	UInt *inner = matrixNoCov_.innerIndexPtr();
	UInt *outer = matrixNoCov_.outerIndexPtr();

	UInt nz = matrixNoCov_.nonZeros();

	UInt nzj;
	UInt counter = 0;
	UInt *jcn = new UInt[nz];

	for(UInt i = 1; i < 2*M*N+1; ++i)
	{
		nzj = outer[i]-outer[i-1];

		for(UInt j=0;j<nzj;++j)
		{
			jcn[counter+j] = i;
			//cout << jcn[counter+j] << " ";
		}
		counter+=nzj;
	}

	UInt *irn = new UInt[nz];

	for(UInt i=0; i<nz; ++i)
	{
		irn[i]=inner[i]+1;
		//cout << irn[i] << " ";
	}

	DMUMPS_STRUC_C id;


	UInt nz_rhs = BB.nonZeros();
	UInt *innerBB = BB.innerIndexPtr();
	UInt *outerBB = BB.outerIndexPtr();

	UInt irhs_sparse[nz_rhs];

	for(UInt i=0; i<nz_rhs; ++i)
	{
		irhs_sparse[i]=innerBB[i]+1;
		//cout << irn[i] << " ";
	}

	UInt irhs_ptr[BB.cols()+1];

	for(UInt i=0; i<BB.cols()+1; ++i)
	{
		irhs_ptr[i]=outerBB[i]+1;
		//cout << irn[i] << " ";
	}

	Real rhs_sparse[nz_rhs];

	// Initialize a MUMPS instance. Use MPI_COMM_WORLD
	id.job=JOB_INIT; id.par=1; id.sym=0;id.comm_fortran=USE_COMM_WORLD;
	dmumps_c(&id);

	//Define the problem on the host
	id.n = BB.cols(); id.nz = nz; id.irn=irn; id.jcn=jcn;
	id.a = values;
	id.nz_rhs = nz_rhs; id.nrhs = BB.cols();
	id.rhs_sparse = rhs_sparse;
	id.irhs_sparse = irhs_sparse;
	id.irhs_ptr = irhs_ptr;

	#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
	/* No outputs */
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
	id.ICNTL(14)=200;
	id.ICNTL(20)=1; id.ICNTL(30)=1;

	/* Call the MUMPS package. */
	id.job=6;
	dmumps_c(&id);

	/* Terminate instance */
	id.job=JOB_END; dmumps_c(&id);

	delete[] irn;
	delete[] jcn;
	//std::cout<<"delete ok"<<std::endl;

	SpMat BBPinv = BB;

	Real *valueBBPinv = BBPinv.valuePtr();

	for(UInt i = 0; i < nz_rhs; ++i)
		valueBBPinv[i] = rhs_sparse[i];

	X = BBPinv*BB;

	// MatrixXr X1;
	// if (regressionData_.getNumberOfRegions() == 0)
	// { //pointwise data
	// 	X1 = B_.transpose() * LeftMultiplybyQ(B_);
	// }
	// else
	// { //areal data
	// 	X1 = B_.transpose() * Ak_ * LeftMultiplybyQ(B_);
	// }
	//
	// if(isRcomputed_==false && regressionData_.getFlagParabolic())
	// {
	// 	isRcomputed_ = true;
	// 	R_.compute(R0k_);
	// }
	//
	// SpMat P;
	// if(regressionData_.getFlagParabolic())
	// {
	// 	SpMat X2 = R1k_+lambdaT*LR0k_;
	// 	P = lambdaS * X2.transpose() * R_.solve(X2);
	// }
	// else
	// {
	// 	P = lambdaS*Psk_ + lambdaT*Ptk_;
	// }
	// MatrixXr X3 = X1 + P;
	// Eigen::LDLT<MatrixXr> Dsolver(X3);
	//
	// MatrixXr X;
	// X = Dsolver.solve(MatrixXr(X1));
	//
	if (regressionData_.getCovariates().rows() != 0) {
		MatrixXr x1 = Gdec_.solve(V_);
		SpMat x2 = (BBPinv * U_ * x1).sparseView();
		X -= x2 * BBPinv * BB;
		degrees += regressionData_.getCovariates().cols();
	}
	for (int i = 0; i<M*N; ++i) {
		degrees += X.coeff(i,i);
	}
	_dof(output_indexS,output_indexT) = degrees;
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT)
{

	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;
	UInt nnodes = M*N;
	UInt nlocations = B_.rows();

	std::default_random_engine generator;
	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = regressionData_.getNrealizations();
	MatrixXr u(nlocations, nrealizations);
	for (int j=0; j<nrealizations; ++j) {
		for (int i=0; i<nlocations; ++i) {
			if (distribution(generator)) {
				u(i,j) = 1.0;
			}
			else {
				u(i,j) = -1.0;
			}
		}
	}

	// Define the first right hand side : | I  0 |^T * B^T * A * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	if (regressionData_.getNumberOfRegions() == 0){
		b.topRows(nnodes) = B_.transpose() * LeftMultiplybyQ(u);
	}else{
		b.topRows(nnodes) = B_.transpose() * Ak_ * LeftMultiplybyQ(u);
	}

	// Resolution of the system
	// system_factorize();

	// MatrixXr x = system_solve(b);
	Real *values = matrixNoCov_.valuePtr();
	UInt *inner = matrixNoCov_.innerIndexPtr();
	UInt *outer = matrixNoCov_.outerIndexPtr();

	UInt nz = matrixNoCov_.nonZeros();

	UInt nzj;
	UInt counter = 0;
	UInt *jcn = new UInt[nz];

	for(UInt i = 1; i < 2*nnodes+1; ++i)
	{
		nzj = outer[i]-outer[i-1];

		for(UInt j=0;j<nzj;++j)
		{
			jcn[counter+j] = i;
			//cout << jcn[counter+j] << " ";
		}
		counter+=nzj;
	}

	UInt *irn = new UInt[nz];

	for(UInt i=0; i<nz; ++i)
	{
		irn[i]=inner[i]+1;
		//cout << irn[i] << " ";
	}

	DMUMPS_STRUC_C id;


	UInt nrhs = b.cols();
	UInt lrhs = b.rows();

	Real* rhs = b.data();

	// Initialize a MUMPS instance. Use MPI_COMM_WORLD
	id.job=JOB_INIT; id.par=1; id.sym=0;id.comm_fortran=USE_COMM_WORLD;
	dmumps_c(&id);

	//Define the problem on the host
	id.n = b.rows(); id.nz = nz; id.irn=irn; id.jcn=jcn;
	id.a = values;
	id.lrhs = lrhs; id.nrhs = nrhs;
	id.rhs = rhs;

	#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
	/* No outputs */
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
	id.ICNTL(14)=1000;
	id.ICNTL(20)=0;

	/* Call the MUMPS package. */
	id.job=6;
	dmumps_c(&id);

	/* Terminate instance */
	id.job=JOB_END; dmumps_c(&id);

	delete[] irn;
	delete[] jcn;
	//std::cout<<"delete ok"<<std::endl;

	MatrixXr x(b.rows(),b.cols());

	for(UInt j = 0; j < nrhs; ++j)
		for(UInt i = 0; i < lrhs; ++i)
			x(i,j) = rhs[i+j*lrhs];

	MatrixXr uTB = u.transpose()*B_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * B * | I  0 |* x ]
	if (regressionData_.getCovariates().rows() != 0) {
		q = regressionData_.getCovariates().cols();
	}
	// For any realization we compute the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {

		edf_vect(i) = uTB.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof(output_indexS,output_indexT) = mean;
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::buildMatrices()
{
	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;

	MixedFERegression<InputHandler, IntegratorSpace, ORDER, mydim, ndim> RegressionSpace(mesh_,regressionData_);

	RegressionSpace.buildSpaceMatrices();
	SpMat psi = RegressionSpace.getPsi();
	SpMat R0 = RegressionSpace.getR0();
	SpMat R1 = RegressionSpace.getR1();

	MixedSplineRegression <InputHandler, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE> Spline(mesh_time_,regressionData_);
	MixedFDRegression <InputHandler> FiniteDifference(mesh_time_,regressionData_);

	SpMat IM(M,M);
	SpMat phi;

	if(regressionData_.getFlagParabolic())
	{
		FiniteDifference.setDerOperator();
		SpMat L = FiniteDifference.getDerOpL();
		IM.setIdentity();
		LR0k_ = kroneckerProduct(L,R0);
		Ptk_.resize(N*M,N*M);
		phi = IM;
		//! right hand side correction for the initial condition:
		rhs_ic_correction_ = (1/(mesh_time_[1]-mesh_time_[0]))*(R0*regressionData_.getInitialValues());
	}
	else
	{
		SpMat IN(N,N);
		Spline.setPhi();
		Spline.setTimeMass();
		Spline.smoothSecondDerivative();
		if(regressionData_.getFlagMass())
		{
			IM = Spline.getTimeMass();
			IN = R0;
		}
		else
		{
			IM.setIdentity();
			IN.setIdentity();
		}
		phi = Spline.getPhi();
		SpMat Pt = Spline.getPt();
		Ptk_ = kroneckerProduct(Pt,IN);
		// if(regressionData_.computeDOF())
		// {
		// 	MatrixXr Ps(N,N);
		// 	Eigen::SparseLU<SpMat> solver;
		// 	solver.compute(R0);
		// 	auto X1 = solver.solve(R1);
		// 	Ps = R1.transpose() * X1;
		// 	Psk_ = kroneckerProduct(IM,Ps.sparseView());
		// }
		LR0k_.resize(N*M,N*M);
	}

	// if(regressionData_.getFlagParabolic() && regressionData_.isLocationsByNodes())
	// 	{
	// 		B.resize(N*M,N*M);
	// 		B.setIdentity();
	// 	}
	B_ = kroneckerProduct(phi,psi);
	B_.makeCompressed();
	R1k_ = kroneckerProduct(IM,R1);
	R1k_.makeCompressed();
	R0k_ = kroneckerProduct(IM,R0);
	R0k_.makeCompressed();
	Ak_ = kroneckerProduct(IM,RegressionSpace.getA());
	Ak_.makeCompressed();

	//! right hand side correction for the forcing term:
	rhs_ft_correction_.resize(M*N);

	if(regressionData_.isSV())
	{
		VectorXr ft = RegressionSpace.getForcingTerm();
		for(UInt i=0; i<N; i++)
		{
			for(UInt j=0; j<M; j++)
			{
				rhs_ft_correction_(i+j*N) = ft(i);
			}
		}
	}
}

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
void SpaceTimeRegression<InputHandler, IntegratorSpace, ORDER, IntegratorTime, SPLINE_DEGREE, ORDER_DERIVATIVE, mydim, ndim>::apply()
{
	UInt N = mesh_.num_nodes();
	UInt M = regressionData_.getFlagParabolic() ? mesh_time_.size()-1 : mesh_time_.size()+SPLINE_DEGREE-1;

	buildMatrices();
	//
	// if(!regressionData_.getCovariates().rows() == 0)
	// {
	// 	setH();
	// 	setQ();
	// }

	VectorXr rightHandData;
	getRightHandData(rightHandData); //updated
	this->_rightHandSide = VectorXr::Zero(2*M*N);
	this->_rightHandSide.topRows(M*N)=rightHandData;

	this->_solution.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	this->_dof.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	this->_GCV.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	if(regressionData_.getCovariates().rows()!=0)
	{
		this->_beta.resize(regressionData_.getLambdaS().size(),regressionData_.getLambdaT().size());
	}

	VectorXr rhs=_rightHandSide;

	for(UInt s = 0; s<regressionData_.getLambdaS().size(); ++s)
	{
		for(UInt t = 0; t<regressionData_.getLambdaT().size(); ++t)
		{
			Real lambdaS = regressionData_.getLambdaS()[s];
			Real lambdaT = regressionData_.getLambdaT()[t];
			_rightHandSide=rhs;
			SpMat R1k_lambda = (-lambdaS)*(R1k_+lambdaT*LR0k_);
			SpMat R0k_lambda = (-lambdaS)*R0k_;
			SpMat BTB_lambda ;
			if(regressionData_.getNumberOfRegions()==0) // pointwise data
			    BTB_lambda=B_.transpose()*B_ + lambdaT*Ptk_;
			else                                        // areal data: need to add the diag(|A_1|,...,|A_N|)
			    BTB_lambda=B_.transpose()*Ak_*B_ + lambdaT*Ptk_;

			this->buildMatrixNoCov(BTB_lambda, R1k_lambda, R0k_lambda);

			if(regressionData_.isSV())
			{
				_rightHandSide.bottomRows(M*N) = -lambdaS*rhs_ft_correction_;
			}
			if(regressionData_.getFlagParabolic())
			{
				for(UInt i = 0; i<regressionData_.getInitialValues().rows(); i++)
				{
					_rightHandSide(M*N+i) -= lambdaS*rhs_ic_correction_(i);
				}
			}

			addNA();
			//Applying boundary conditions if necessary
			if(regressionData_.getDirichletIndices().size() != 0)  // if areal data NO BOUNDARY CONDITIONS
				addDirichletBC();

			system_factorize();
			//
		  _solution(s,t) = this->template system_solve(this->_rightHandSide);
			// Mumps::template solve(matrixNoCov_,_rightHandSide,_solution(s,t));
			//
			if(regressionData_.computeDOF())
			{
				computeDegreesOfFreedom(s,t,lambdaS,lambdaT);
			}
			else
			{
				_dof(s,t) = -1;
				_GCV(s,t) = -1;
			}

			if(regressionData_.getCovariates().rows()!=0)
			{
				MatrixXr W(this->regressionData_.getCovariates());
				VectorXr beta_rhs = W.transpose()*(regressionData_.getObservations() - B_*_solution(s,t).topRows(B_.cols()));
				_beta(s,t) = WTW_.solve(beta_rhs);
			}
			// Rprintf("s:%d		t:%d		dof:%f\n",s,t,_dof(s,t));
		}
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildSpaceMatrices(EOExpr<A> oper, const ForcingTerm & u)
{
	UInt nnodes=mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	setA();
	setPsi();

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);


	if(regressionData_.isSV())
	{
		Assembler::forcingTerm(mesh_, fe, u, this->forcingTerm_);
	}
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataTime, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataTime, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataTime& regressionData):MixedFERegressionBase<RegressionDataTime, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void buildSpaceMatrices()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionDataTime, Integrator, ORDER, mydim, ndim>::buildSpaceMatrices(stiff, ForcingTerm(std::vector<Real>(1)));
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataTimeElliptic, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataTimeElliptic, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataTimeElliptic& regressionData):MixedFERegressionBase<RegressionDataTimeElliptic, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void buildSpaceMatrices()
	{
		if(mydim!=2 || ndim !=2)
		{

		#ifdef R_VERSION_
			Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
		#else
			std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
		#endif

		}
		else
		{
			typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
			typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
			typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = this->regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& b = this->regressionData_.getBeta();

	    MixedFERegressionBase<RegressionDataTimeElliptic, Integrator, ORDER, mydim, ndim>::buildSpaceMatrices(c*mass+stiff[K]+dot(b,grad), ForcingTerm(std::vector<Real>(1)));
		}
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataTimeEllipticSpaceVarying, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataTimeEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>
{

	public:
		MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataTimeEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataTimeEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

		void buildSpaceMatrices()
		{
			if(mydim!=2 || ndim !=2)
			{
				#ifdef R_VERSION_
					Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
				#else
					std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
				#endif
			}
			else
			{
				typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
				typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
				typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

				const Reaction& c = this->regressionData_.getC();
				const Diffusivity& K = this->regressionData_.getK();
				const Advection& b = this->regressionData_.getBeta();
				const ForcingTerm& u= this->regressionData_.getU();

				MixedFERegressionBase<RegressionDataTimeEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>::buildSpaceMatrices(c*mass+stiff[K]+dot(b,grad), u);
			}
		}
};

// Parte temporale:
template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::setPhi(){

		Spline<Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);
		UInt M = spline.num_knots()-SPLINE_DEGREE-1;
		UInt m = regressionData_.getNumberofTimeObservations();

		phi_.resize(m, M);
		Real value;

    for (UInt i = 0; i < m; ++i)
        for (UInt j = 0; j < M; ++j)
        {
					//value = spline.BasisFunction(SPLINE_DEGREE, j, mesh_time_[i]);
					value = spline.BasisFunction(SPLINE_DEGREE, j, this->regressionData_.getTimeLocations()[i]);
					if (value!=0)
					{
						phi_.coeffRef(i,j) = value;
					}
				}
    phi_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::setTimeMass(){

    using ETTimeMass = EOExpr<TimeMass>;

    Spline<Integrator, SPLINE_DEGREE, 0> spline(mesh_time_);

    TimeMass ETimeMass;
    ETTimeMass timeMass(ETimeMass);

    Assembler::operKernel(timeMass, spline, timeMass_);
}

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
void MixedSplineRegression<InputHandler, Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE>::smoothSecondDerivative(){

    using ETTimeMass = EOExpr<TimeMass>;

    Spline<Integrator, SPLINE_DEGREE, ORDER_DERIVATIVE> spline(mesh_time_);

    TimeMass ETimeMass;
    ETTimeMass timeMass(ETimeMass);

    Assembler::operKernel(timeMass, spline, Pt_);
}

// Parabolic
template<typename InputHandler>
void MixedFDRegression<InputHandler>::setDerOperator(){

	UInt M = mesh_time_.size()-1;
	derOpL_.resize(M, M);

	// set the first and the last rows
	Real delta = mesh_time_[1] - mesh_time_[0];
	derOpL_.coeffRef(0,0) = 1/delta;

	delta = mesh_time_[M-1] - mesh_time_[M-2];
	derOpL_.coeffRef(M-1,M-1) = 1/delta;
	derOpL_.coeffRef(M-1,M-2) = -1/delta;

	for (UInt i = 1; i < M-1; ++i)
	{
		delta = mesh_time_[i] - mesh_time_[i-1];
		derOpL_.coeffRef(i,i-1) = -1/delta;
		derOpL_.coeffRef(i,i) 	= 1/delta;
	}

	derOpL_.makeCompressed();

}


#endif
