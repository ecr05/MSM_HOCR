/*  optimise.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

// Mathematical optimisation functions


#if !defined(__optimise_h)
#define __optimise_h

#include <cmath>
#include "newmatap.h"
#include "string"

using namespace NEWMAT;

namespace MISCMATHS {

float optimise1d(ColumnVector &pt, const ColumnVector dir, 
		const ColumnVector tol, int &iterations_done, 
		float (*func)(const ColumnVector &), int max_iter,
		float &init_value, float boundguess);


 float optimise(ColumnVector &pt, int numopt, const ColumnVector &tol, 
		float (*func)(const ColumnVector &), int &iterations_done, 
		int max_iter, const ColumnVector& boundguess, 
		const std::string type="brent");

}

#endif
