#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include "integratePsi.h"
#include "kronecker_product.h"
#include <memory>
#include "../inst/include/dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

//! A base class for the smooth regression.
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegressionBase
{
	protected:

		const MeshHandler<ORDER, mydim, ndim> &mesh_;
		const InputHandler& regressionData_;

		SpMat R1_;	//! Finite Element matrix associated to the bilinear form a
		SpMat R0_;	//! Mass matrix in space
		SpMat psi_; //! Matrix of the evaluations of the spatial basis functions in the space locations
		SpMat A_;		//! A_.asDiagonal() = diag(|A_1|,...,|A_N|) areal matrix

		VectorXr forcingTerm_;

		void setPsi(); 	//! A method computing psi_
		void setA();		//! A method computing A_

	public:
		MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& regressionData): mesh_(mesh), regressionData_(regressionData) {};

		SpMat const & getPsi() const {return psi_;}
		SpMat const & getR0() const {return R0_;}
		SpMat const & getR1() const {return R1_;}
		SpMat const & getA() const {return A_;}
		VectorXr const & getForcingTerm() const {return forcingTerm_;}

		template<typename A>
		void buildSpaceMatrices(EOExpr<A> oper, const ForcingTerm & u);
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>
{
	public:
		MixedFERegression(const MeshHandler<ORDER, ndim, mydim>& mesh, const InputHandler& regressionData):MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

		void buildSpaceMatrices()
		{
			std::cout << "Option not implemented! \n";
		}
};

//! A class for the construction of the temporal matrices needed for the parabolic case
template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
class MixedSplineRegression
{
	private:
		const std::vector<Real>& mesh_time_;
		const InputHandler& regressionData_;

		SpMat phi_;   //! Matrix of the evaluations of the spline basis functions in the time locations
		SpMat Pt_;
		SpMat timeMass_; //! Mass matrix in time

	public:
		MixedSplineRegression(const std::vector<Real>& mesh_time, const InputHandler& regressionData):mesh_time_(mesh_time), regressionData_(regressionData){};

    void setPhi();
		void setTimeMass();
    void smoothSecondDerivative();

		inline SpMat const & getPt() const { return Pt_; }
		inline SpMat const & getPhi() const { return phi_; }
		inline SpMat const & getTimeMass() const { return timeMass_; }

};

//! A class for the construction of the temporal matrices needed for the separable case
template<typename InputHandler>
class MixedFDRegression
{
	private:
		const std::vector<Real>& mesh_time_;
		const InputHandler& regressionData_;

		SpMat derOpL_; //!matrix associated with derivation in time

	public:
		MixedFDRegression(const std::vector<Real>& mesh_time, const InputHandler& regressionData):mesh_time_(mesh_time), regressionData_(regressionData){};

    void setDerOperator(); //! sets derOpL_
		inline SpMat const & getDerOpL() const { return derOpL_; }

};

//! A class that stores all the needed matrices and contains the methods necessary to assemble and solve the final system
template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class SpaceTimeRegression
{
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const std::vector<Real>& mesh_time_;
	const UInt N_; //! Number of spatial basis functions. The final system ha dimension 2N_*M_
	const UInt M_;
	/*!
		M_= number of temporal basis functions for parabolic problem.
		M_= number of time istants -1 for parabolic problem.
		M_= 1. This happens only when initial condition for parabolic problem are not provided and has to be estimated solving a "only-space" problem
	*/
  const InputHandler& regressionData_;
	// Separable case:
	//  system matrix= 	| B^T * Ak *B + lambdaT*Ptk |  -lambdaS*R1k^T  |   +  |B^T * Ak * (-H) * B |  O |   =  matrixNoCov + matrixOnlyCov
	//	                |      -lambdaS*R1k^T       |  -lambdaS*R0k	   |      |         O          |  O |

	// Parabolic case:
	//  system matrix= 	|          B^T * Ak *B           | -lambdaS*(R1k^T+lambdaT*LR0k)  |   +  |B^T * Ak * (-H) * B |  O |   =  matrixNoCov + matrixOnlyCov
	//	                | -lambdaS*(R1k^T+lambdaT*LR0k)  |        -lambdaS*R0k	          |      |         O          |  O |


	SpMat matrixNoCov_;	//! System matrix with psi^T*psi or psi^T*A*psi in north-west block  (is the full system matrix if no covariates)

	//! kron(IM,Ps) (separable version)
	SpMat Psk_;
	//! kron(Pt,IN) (separable version)
	SpMat Ptk_;
	//! kron(IM,R1)
	SpMat R1k_;
	//! kron(L,R0) (parabolic version)
	SpMat LR0k_;
	//! kron(IM,R0)
	SpMat R0k_;
	/*! kron(phi,psi) for separable case
	 		kron(IM,psi) 	for parabolic case
	*/
	SpMat B_;
	//! Kronecker product of the matrix W (1/domainArea) and identity
	SpMat Ak_; //Ak_.asDiagonal() = kron(IM,diag(|A_1|,...,|A_N|)) areal matrix

	MatrixXr U_;	//! psi^T * W or psi^T * A * W padded with zeros, needed for Woodbury decomposition
	MatrixXr V_;   //! W^T*psi, if pointwise data is U^T, needed for Woodbury decomposition
	VectorXr z_; //! Observations

	Eigen::SparseLU<SpMat> matrixNoCovdec_; //! Stores the factorization of matrixNoCov_
	Eigen::PartialPivLU<MatrixXr> Gdec_;	//! Stores factorization of G =  C + [V * matrixNoCov^-1 * U]
	Eigen::PartialPivLU<MatrixXr> WTW_;	//! Stores the factorization of W^T * W
	bool isWTWfactorized_=false; //! true if the matrix W^T*W has been factorized, false otherwise
	bool isRcomputed_=false; //! true if the matrix R0k has been factorized, false otherwise
	Eigen::SparseLU<SpMat> R_; //! Stores the factorization of R0k_

	VectorXr rhs_ft_correction_;	//! right hand side correction for the forcing term:
	VectorXr rhs_ic_correction_;	//!Initial condition correction (parabolic case)
	VectorXr _rightHandSide;         //!A Eigen::VectorXr: Stores the system right hand side.
	MatrixXv _solution; //!A Eigen::VectorXv: Stores the system solution.
	MatrixXr _dof;          //! A Eigen::MatrixXr storing the computed dofs
	MatrixXr _GCV;	 //! A Eigen::MatrixXr storing the computed GCV
	UInt bestLambdaS_=0;	//!Stores the index of the best lambdaS according to GCV
	UInt bestLambdaT_=0;	//!Stores the index of the best lambdaT according to GCV
	Real _bestGCV=10e20;	//!Stores the value of the best GCV
	MatrixXv _beta;		//! A Eigen::MatrixXv storing the computed beta coefficients

	//! A method computing the no-covariates version of the system matrix
	void buildMatrixNoCov(const SpMat& B,  const SpMat& SWblock,  const SpMat& SEblock);
	//! A function that given a vector u, performs Q*u efficiently
	MatrixXr LeftMultiplybyQ(const MatrixXr& u);
	//! A method which adds Dirichlet boundary conditions before solving the system ( Remark: BC for areal data are not implemented!)
	void addDirichletBC();
	//! A method which takes care of missing values setting to 0 the corresponding rows of B_
	void addNA();
	//! A method returning the system right hand data
	void getRightHandData(VectorXr& rightHandData);
	//! A method which builds all the matrices needed for assembling matrixNoCov_
	void buildMatrices();
	//! A method computing the dofs
	void computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method computing dofs in case of exact GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method computing dofs in case of stochastic GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method computing GCV from the dofs
	void computeGeneralizedCrossValidation(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A method to factorize the system, using Woodbury decomposition useful when there are covariates
	void system_factorize();
	//! A function which solves the factorized system
	template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);

public:
	SpaceTimeRegression(const MeshHandler<ORDER,mydim,ndim>& mesh, const std::vector<Real>& mesh_time, const InputHandler& regressionData):
		    mesh_(mesh), mesh_time_(mesh_time), N_(mesh_.num_nodes()), M_(mesh_time.size()==1? 1: regressionData.getFlagParabolic() ? mesh_time.size()-1 : mesh_time.size()+SPLINE_DEGREE-1),
				regressionData_(regressionData),_dof(regressionData_.getDOF_matrix()){}
	/*!
		The function that builds the needed matrices, for each combination of lambdaS
		and lambdaT assembles the system matrix and right hand side and solves the
		corresponding system
	*/
	void apply();
	//! A  method returning the computed solution
	MatrixXv const & getSolution() const{return _solution;}
	//! A method returning the computed dofs of the model
	MatrixXr const & getDOF() const{return _dof;}
	//! A method returning the computed GCV of the model
	MatrixXr const & getGCV() const{return _GCV;}
	//! A method returning the computed beta coefficients of the model
	MatrixXv const & getBeta() const{return _beta;}
	//! A method returning the index of the best lambdaS according to GCV
	UInt getBestLambdaS(){return bestLambdaS_;}
	//! A method returning the index of the best lambdaT according to GCV
	UInt getBestLambdaT(){return bestLambdaT_;}


};


#include "mixedFERegression_imp.h"

#endif
