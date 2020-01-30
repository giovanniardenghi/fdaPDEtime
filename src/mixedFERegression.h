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

		SpMat R1_;	//! North-east block of system matrix matrixNoCov_
		SpMat R0_;	//! South-east block of system matrix matrixNoCov_
		SpMat psi_; //! Psi matrix of the model
		SpMat A_;	//A_.asDiagonal() = diag(|A_1|,...,|A_N|) areal matrix

		VectorXr forcingTerm_;

		void setPsi(); 	//! A member function computing the Psi matrix
		void setA();

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

template<typename InputHandler, typename Integrator, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE>
class MixedSplineRegression
{
	private:
		const std::vector<Real>& mesh_time_;
		const InputHandler& regressionData_;

		SpMat phi_;
		SpMat Pt_;
		SpMat timeMass_;

	public:
		MixedSplineRegression(const std::vector<Real>& mesh_time, const InputHandler& regressionData):mesh_time_(mesh_time), regressionData_(regressionData){};

    void setPhi();
		void setTimeMass();
    void smoothSecondDerivative();

		inline SpMat const & getPt() const { return Pt_; }
		inline SpMat const & getPhi() const { return phi_; }
		inline SpMat const & getTimeMass() const { return timeMass_; }

};

//! A LinearSystem class: A class for the construction of the temporal matrices for the separable case (FLAG_PARABOLIC = TRUE)
template<typename InputHandler>
class MixedFDRegression
{
	private:
		const std::vector<Real>& mesh_time_;
		const InputHandler& regressionData_;

		SpMat derOpL_;

	public:
		MixedFDRegression(const std::vector<Real>& mesh_time, const InputHandler& regressionData):mesh_time_(mesh_time), regressionData_(regressionData){};

    void setDerOperator();

		inline SpMat const & getDerOpL() const { return derOpL_; }

};

template<typename InputHandler, typename IntegratorSpace, UInt ORDER, typename IntegratorTime, UInt SPLINE_DEGREE, UInt ORDER_DERIVATIVE, UInt mydim, UInt ndim>
class SpaceTimeRegression
{
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const std::vector<Real>& mesh_time_;
	const UInt N_;
	UInt M_;
  const InputHandler& regressionData_;
	// Separable case:
	//  system matrix= 	| B^T * Ak *B + lambdaT*Ptk |  -lambdaS*R1k^T  |   +  |B^T * Ak * (-H) * B |  O |   =  matrixNoCov + matrixOnlyCov
	//	                |      -lambdaS*R1k^T       |  -lambdaS*R0k	   |      |         O          |  O |

	// Parabolic case:
	//  system matrix= 	|          B^T * Ak *B           | -lambdaS*(R1k^T+lambdaT*LR0k)  |   +  |B^T * Ak * (-H) * B |  O |   =  matrixNoCov + matrixOnlyCov
	//	                | -lambdaS*(R1k^T+lambdaT*LR0k)  |        -lambdaS*R0k	          |      |         O          |  O |


	SpMat matrixNoCov_;	//! System matrix with psi^T*psi or psi^T*A*psi in north-west block  (is the full system matrix if no covariates)
	// SpMat matrixOnlyCov_; //! coeffmatrix=matrixNoCov+matrixOnlyCov

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

	// MatrixXr Q_;  //! Identity - H, projects onto the orthogonal subspace
	// MatrixXr H_; //! The hat matrix of the regression

	VectorXr rhs_ft_correction_;	//! right hand side correction for the forcing term:
	VectorXr rhs_ic_correction_;	//!Initial condition correction (parabolic case)
	VectorXr _rightHandSide;         //!A Eigen::VectorXr: Stores the system right hand side.
	MatrixXv _solution; //!A Eigen::VectorXr: Stores the system solution.
	MatrixXr _dof;          //! A Eigen::MatrixXr storing the computed dofs
	MatrixXr _GCV;	 //! A Eigen::MatrixXr storing the computed GCV
	UInt bestLambdaS_=0;	//!Stores the index of the lambdaS of best GCV
	UInt bestLambdaT_=0;	//!Stores the index of the lambdaT of best GCV
	Real _bestGCV=10e20;	//!Stores the value of the bestGCV
	MatrixXv _beta;		//! A Eigen::MatrixXv storing the computed beta coefficients

	//! A member function computing the no-covariates version of the system matrix
	void buildMatrixNoCov(const SpMat& B,  const SpMat& SWblock,  const SpMat& SEblock);
	//! A member function computing the matrix to be added to matrixNoCov_ to obtain the full system matrix
	// void buildMatrixOnlyCov(const SpMat& B, const MatrixXr& H);
	//! A function that given a vector u, performs Q*u efficiently
	MatrixXr LeftMultiplybyQ(const MatrixXr& u);
	//! A function which adds Dirichlet boundary conditions before solving the system ( Remark: BC for areal data are not implemented!)
	void addDirichletBC();
	//! A function which takes care of missing values setting to 0 the corresponding rows of B_
	void addNA();
	//! A member function which builds the Q matrix
	// void setQ();
	// //! A member function which builds the H matrix
	// void setH();
	// //! A member function which builds the A vector containing the areas of the regions in case of areal data
	// void setAk();

	//! A member function returning the system right hand data
	void getRightHandData(VectorXr& rightHandData);
	//! A member function which builds all the matrices needed for assembling matrixNoCov_
	void buildMatrices();
	//! A member function computing the dofs
	void computeDegreesOfFreedom(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A function computing dofs in case of exact GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomExact(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A function computing dofs in case of stochastic GCV, it is called by computeDegreesOfFreedom
	void computeDegreesOfFreedomStochastic(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);
	//! A function computing GCV from the dofs
	void computeGeneralizedCrossValidation(UInt output_indexS, UInt output_indexT, Real lambdaS, Real lambdaT);

	//! A function to factorize the system, using Woodbury decomposition when there are covariates
	void system_factorize();
	//! A function which solves the factorized system
	template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);

public:
	SpaceTimeRegression(const MeshHandler<ORDER,mydim,ndim>& mesh, const std::vector<Real>& mesh_time, const InputHandler& regressionData):
		    mesh_(mesh), mesh_time_(mesh_time), N_(mesh_.num_nodes()), regressionData_(regressionData),_dof(regressionData_.getDOF_matrix())
				{
					if(mesh_time.size()==1)
						M_=1;
					else
						M_ = regressionData.getFlagParabolic() ? mesh_time.size()-1 : mesh_time.size()+SPLINE_DEGREE-1;
				};
	//! The function solving the system, used by the children classes. Saves the result in _solution
	/*!
	    \param oper an operator, which is the Stiffness operator in case of Laplacian regularization
	    \param u the forcing term, will be used only in case of anysotropic nonstationary regression
	*/
	void apply();

	//! A inline member that returns a VectorXr, returns the whole solution_.
	MatrixXv const & getSolution() const{return _solution;}
	//! A function returning the computed dofs of the model
	MatrixXr const & getDOF() const{return _dof;}
	//! A function returning the computed GCV of the model
	MatrixXr const & getGCV() const{return _GCV;}
	//! A function returning the computed beta coefficients of the model
	MatrixXv const & getBeta() const{return _beta;}
	//! A function returning the index of the lambdaS of the best GCV
	UInt getBestLambdaS(){return bestLambdaS_;}
	//! A function returning the index of the lambdaT of the best GCV
	UInt getBestLambdaT(){return bestLambdaT_;}


};





#include "mixedFERegression_imp.h"
// #include "timeRegression_imp.h"s

#endif
