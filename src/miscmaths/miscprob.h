/*  miscprob.h

    Christian Beckmann & Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

// Miscellaneous maths functions that rely on libprob build ontop of miscmaths


#if !defined(__miscprob_h)
#define __miscprob_h

#include "miscmaths.h"
#include "libprob.h"
#include "stdlib.h"

using namespace NEWMAT;

namespace MISCMATHS {

//   ReturnMatrix betarnd(const int dim1, const int dim2, 
// 		       const float a, const float b); 

  ReturnMatrix betapdf(const RowVector& vals, 
		       const float a, const float b); 

  ReturnMatrix unifrnd(const int dim1 = 1, const int dim2 = -1, 
		       const float start = 0, const float end = 1);
  
  ReturnMatrix normrnd(const int dim1 = 1, const int dim2 = -1, 
		       const float mu = 0, const float sigma = 1);

  // returns nsamps*nparams matrix:
  ReturnMatrix mvnrnd(const RowVector& mu, const SymmetricMatrix& covar, int nsamp = 1);
  
  float mvnpdf(const RowVector& vals, const RowVector& mu, const SymmetricMatrix& covar);

  float bvnpdf(const RowVector& vals, const RowVector& mu, const SymmetricMatrix& covar);

  float normpdf(const float val, const float mu = 0, const float var = 1);
  float lognormpdf(const float val, const float mu = 0, const float var = 1);

  ReturnMatrix normpdf(const RowVector& vals, const float mu = 0, const float var = 1);

  ReturnMatrix normpdf(const RowVector& vals, const RowVector& mus, 
		       const RowVector& vars);

  ReturnMatrix normcdf(const RowVector& vals, const float mu = 0, const float var = 1);

  ReturnMatrix gammapdf(const RowVector& vals, const float mu = 0, const float var = 1);

  ReturnMatrix gammacdf(const RowVector& vals, const float mu = 0, const float var = 1);

//   ReturnMatrix gammarnd(const int dim1, const int dim2, 
// 			const float a, const float b);

  // returns n! * n matrix of all possible permutations
  ReturnMatrix perms(const int n);

  
  class Mvnormrandm
    {
    public:
      Mvnormrandm(){}

      Mvnormrandm(const RowVector& pmu, const SymmetricMatrix& pcovar) :
	mu(pmu),
	covar(pcovar)
	{
	  Matrix eig_vec;
	  DiagonalMatrix eig_val;
	  EigenValues(covar,eig_val,eig_vec);

	  covarw = sqrt(eig_val)*eig_vec.t();
	}

      ReturnMatrix next(int nsamp = 1) const 
	{
	  Matrix ret = ones(nsamp, 1)*mu + normrnd(nsamp,mu.Ncols())*covarw;
	  ret.Release();
	  return ret;
	}

      ReturnMatrix next(const RowVector& pmu, int nsamp = 1)  
	{
	  mu=pmu;

	  Matrix ret = ones(nsamp, 1)*mu + normrnd(nsamp,mu.Ncols())*covarw;
	  ret.Release();
	  return ret;
	}

      void setcovar(const SymmetricMatrix& pcovar)
	{
	  covar=pcovar;

	  mu.ReSize(covar.Nrows());
	  mu=0;

	  Matrix eig_vec;
	  DiagonalMatrix eig_val;
	  EigenValues(covar,eig_val,eig_vec);

	  covarw = sqrt(eig_val)*eig_vec.t();
	}

    private:      

      RowVector mu;
      SymmetricMatrix covar;

      Matrix covarw;

    };
}
#endif






