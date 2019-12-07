#ifndef __REGRESSIONDATA_IMP_HPP__
#define __REGRESSIONDATA_IMP_HPP__

RegressionDataTime::RegressionDataTime(std::vector<Point>& locations, std::vector<Real>& time_locations, VectorXr& observations, UInt order, std::vector<Real>& lambdaS, std::vector<Real>& lambdaT, MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF):
					locations_(locations), time_locations_(time_locations), observations_(observations), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
					order_(order), lambdaS_(lambdaS), lambdaT_(lambdaT), bc_values_(bc_values), bc_indices_(bc_indices), ic_(ic), flag_mass_(flag_mass), flag_parabolic_(flag_parabolic), DOF_(DOF)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(int i = 0; i<observations_.size();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

RegressionDataTimeElliptic::RegressionDataTimeElliptic(std::vector<Point>& locations, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
												std::vector<Real>& lambdaS, std::vector<Real>& lambdaT, Eigen::Matrix<Real,2,2>& K,
												Eigen::Matrix<Real,2,1>& beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF):
		 RegressionDataTime(locations, time_locations, observations, order, lambdaS, lambdaT, covariates, incidenceMatrix, bc_indices, bc_values, ic, flag_mass, flag_parabolic, DOF),
		 										K_(K), beta_(beta), c_(c)
{;}

RegressionDataTimeEllipticSpaceVarying::RegressionDataTimeEllipticSpaceVarying(std::vector<Point>& locations, std::vector<Real>& time_locations, VectorXr& observations, UInt order,
												std::vector<Real>& lambdaS, std::vector<Real>& lambdaT,
									const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
									const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
									const std::vector<Real>& c, const std::vector<Real>& u,
									MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, VectorXr& ic, bool flag_mass, bool flag_parabolic, bool DOF):
		 RegressionDataTime(locations, time_locations, observations, order, lambdaS, lambdaT, covariates, incidenceMatrix, bc_indices, bc_values, ic, flag_mass, flag_parabolic, DOF),
											K_(K), beta_(beta), c_(c), u_(u)
{;}


#ifdef R_VERSION_
RegressionDataTime::RegressionDataTime(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP Rcovariates,
						SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP DOF, SEXP RGCVmethod,
						SEXP Rnrealizations)
{
	setLocations(Rlocations);
	setTimeLocations(Rtime_locations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);
	setNrealizations(Rnrealizations);

	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	flag_mass_ = LOGICAL(Rflag_mass)[0];
	flag_parabolic_ = LOGICAL(Rflag_parabolic)[0];

	DOF_ = INTEGER(DOF)[0];
	UInt length_indexes = Rf_length(RBCIndices);

	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);

	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

	UInt length_lambdaS = Rf_length(RlambdaS);
	for (UInt i = 0; i<length_lambdaS; ++i)  lambdaS_.push_back(REAL(RlambdaS)[i]);

	UInt length_lambdaT = Rf_length(RlambdaT);
	for (UInt i = 0; i<length_lambdaT; ++i)  lambdaT_.push_back(REAL(RlambdaT)[i]);

	UInt length_ic = Rf_length(Ric);
	ic_.resize(length_ic);
	for (UInt i = 0; i<length_ic; ++i)  ic_[i]=REAL(Ric)[i];

}

RegressionDataTimeElliptic::RegressionDataTimeElliptic(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP RK,
		SEXP Rbeta, SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric,	SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations):
	RegressionDataTime(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations)
{
	K_.resize(2, 2);
	for(auto i=0; i<2; ++i)
	{
		for(auto j=0; j<2 ; ++j)
		{
			K_(i,j)=REAL(RK)[i+ 2*j];
		}
	}

	beta_.resize(2);
	for(auto i=0; i<2 ; ++i)
	{
		beta_(i)=REAL(Rbeta)[i];
	}

	c_ =  REAL(Rc)[0];
}

RegressionDataTimeEllipticSpaceVarying::RegressionDataTimeEllipticSpaceVarying(SEXP Rlocations, SEXP Rtime_locations, SEXP Robservations, SEXP Rorder, SEXP RlambdaS, SEXP RlambdaT, SEXP RK,
		SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP Rflag_mass, SEXP Rflag_parabolic, SEXP Ric, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations):
					 RegressionDataTime(Rlocations, Rtime_locations, Robservations, Rorder, RlambdaS, RlambdaT, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, Rflag_mass, Rflag_parabolic, Ric, DOF, RGCVmethod, Rnrealizations),
					 K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{this->isSpaceVarying=TRUE;}

void RegressionDataTime::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(locations_.size() == 0 && nRegions_ == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				observations_[i] = REAL(Robservations)[i];
				observations_indices_.push_back(i);
			}
			else
			{
				observations_[i] = 0.0;
				observations_na_.push_back(i);
			}
		}
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		locations_by_nodes_ = false;
		for(auto i=0;i<n_obs_;++i)
		{
			observations_[i] = REAL(Robservations)[i];
		}
	}

	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

void RegressionDataTime::setCovariates(SEXP Rcovariates)
{
	n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];

	covariates_.resize(n_, p_);

	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			covariates_(i,j)=REAL(Rcovariates)[i+ n_*j];
		}
	}
}

void RegressionDataTime::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0){
		int ndim = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[1];

	  if (ndim == 2){
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
			}
		}else{ //ndim == 3
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1],REAL(Rlocations)[i+ n_*2]);
			}
		}
	}
}

void RegressionDataTime::setTimeLocations(SEXP Rtime_locations)
{
	UInt n_time_loc_ = Rf_length(Rtime_locations);
  time_locations_.resize(n_time_loc_);

	for(auto i=0;i<n_time_loc_;++i)
	{
		time_locations_[i] = REAL(Rtime_locations)[i];
	}
}

void RegressionDataTime::setNrealizations(SEXP Rnrealizations) {
	nrealizations_ = INTEGER(Rnrealizations)[0];
}

void RegressionDataTime::setIncidenceMatrix(SEXP RincidenceMatrix)
{
	nRegions_ = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt p = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1];

	incidenceMatrix_.resize(nRegions_, p);

	for(auto i=0; i<nRegions_; ++i)
	{
		for(auto j=0; j<p; ++j)
		{
			incidenceMatrix_(i,j) = INTEGER(RincidenceMatrix)[i+nRegions_*j];
		}
	}
}

#endif

void RegressionDataTime::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.size(); i++)
	{
		out<<i<<"\t"<<observations_(i)<<std::endl;
	}
}

void RegressionDataTime::printCovariates(std::ostream & out) const
{

	for(auto i=0;i<covariates_.rows(); i++)
	{
		for(auto j=0; j<covariates_.cols(); j++)
		{
			out<<covariates_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

void RegressionDataTime::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

void RegressionDataTime::printIncidenceMatrix(std::ostream & out) const
{
	for (auto i=0; i<incidenceMatrix_.rows(); i++)
	{
		for (auto j=0; j<incidenceMatrix_.cols(); j++)
		{
			out << incidenceMatrix_(i,j) << "\t";
		}
		out << std::endl;
	}
}

#endif
