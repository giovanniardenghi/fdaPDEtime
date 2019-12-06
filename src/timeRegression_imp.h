#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

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
