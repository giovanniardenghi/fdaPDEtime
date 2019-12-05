#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

#include "fdaPDE.h"
#include "../inst/include/dmumps_c.h"

//!  A Linear System QR solver class
/*!
 * This class gives offers a standard interface to the QR resolutor for dense matrices
*/
class QR{
	public:
	static void solve(const MatrixXr & A, const VectorXr & b,VectorXr &x){x=A.householderQr().solve(b);};
};

//!  A Linear System LU Partial Pivoting solver class
/*!
 * This class gives offers a standard interface to the LU Partial Pivoting resolutor for dense matrices.
 * OBS: The matrix should be invertible.
*/
class LUPV{
	public:
	static void solve(MatrixXr const & A, VectorXr const & b,VectorXr &x){x=A.partialPivLu().solve(b);};
};

//!  A Linear System LDLT solver class
/*!
 * This class gives offers a standard interface to the LDLT resolutor for dense matrices.
 * OBS: The matrix should be symmetric and SDP.
*/
class Symmetric{
	public:
	static void solve(MatrixXr const & A, VectorXr const & b,VectorXr &x){x=A.ldlt().solve(b);};
};

//!  A Linear System Cholesky solver class
/*!
 * This class gives offers a standard interface to the Cholesky resolutor for dense matrices.
 * OBS: The matrix should be symmetric and SDP, faster and more stable than others.
*/
class Cholesky{
	public:
	static void solve(MatrixXr const & A, VectorXr const & b,VectorXr &x){x=A.ldlt().solve(b);};
};

//!  A Linear System LU sparse solver class
/*!
 * This class gives offers a standard interface to the LU resolutor for sparse matrices.
*/
class SpLU{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::SparseLU<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System QR sparse solver class
/*!
 * This class gives offers a standard interface to the QR resolutor for sparse matrices.
*/
class SpQR{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::SparseQR<SpMat,Eigen::COLAMDOrdering<int> > solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System Cholesky sparse solver class
/*!
 * This class gives offers a standard interface to the Cholesky resolutor for sparse matrices.
*/
class SpCholesky{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::SimplicialLDLT<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System Conjugate Gradient sparse solver class
/*!
 * This class gives offers a standard interface to the Conjugate Gradient resolutor for sparse matrices.
*/
class SpConjGrad{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::ConjugateGradient<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System BiConjugate Gradient stabilized sparse solver class
/*!
 * This class gives offers a standard interface to the BiConjugate Gradient stabilized resolutor for sparse matrices.
*/

class BiCGSTAB{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::BiCGSTAB<SpMat> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};

//!  A Linear System BiConjugate Gradient stabilized with Incomplete LUT preconditioner sparse solver class
/*!
 * This class gives offers a standard interface to the BiConjugate Gradient stabilized BiConjugate Gradient stabilized with Incomplete LUT preconditioner resolutor for sparse matrices.
*/

class BiCGSTABILUT{
	public:
	static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{
		Eigen::BiCGSTAB<SpMat,Eigen::IncompleteLUT<Real>> solver;
		solver.compute(A);
		if(solver.info()!=Eigen::Success){
		//std::cerr<<"Decomposition failed!"<<std::endl;
		}
		x=solver.solve(b);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
	};
};


class Mumps{
	public:
    static void solve(SpMat const & A, VectorXr const & b, VectorXr &x )
	{

		const Real *values = A.valuePtr();
		const UInt *inner = A.innerIndexPtr();
		const UInt *outer = A.outerIndexPtr();

		UInt nz = A.nonZeros();
		UInt n = A.cols();

		UInt nzj;
		UInt counter = 0;
		UInt *jcn = new UInt[nz];

		for(UInt i = 1; i <n+1; ++i)
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

		Real *a = new Real[nz];

		for(UInt i=0; i<nz; ++i)
		{
	    	a[i]=values[i];
	    	//cout << irn[i] << " ";
		}

	    DMUMPS_STRUC_C id;

		//Real *rhs = b.array();
        Real *rhs = new Real[n];
        for(UInt i = 0; i < n; ++i) rhs[i] = b(i);

	    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
	    id.job=JOB_INIT; id.par=1; id.sym=0; id.comm_fortran=USE_COMM_WORLD;
	    dmumps_c(&id);
	    /* Define the problem on the host */

		id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
		id.a = a; id.rhs = rhs;

	    #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
	    /* No outputs */
	    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
	    /* Call the MUMPS package. */
	    id.job=6;
	    dmumps_c(&id);
	    id.job=JOB_END; dmumps_c(&id); /* Terminate instance */

		for(UInt i=0; i<n; ++i)
		{
	    	x(i)=rhs[i];
		}

     };
};


#endif
