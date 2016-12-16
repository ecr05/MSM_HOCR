/*  rungekutta.h

    Mark Woolrich - FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(rungekutta_h)
#define rungekutta_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "newmatap.h"
#include "newmatio.h"

using namespace NEWMAT;

namespace MISCMATHS {

class Derivative
{
public:
  Derivative(int pny) : ny(pny), dy(pny) {}

  // x is time point to evaluate at
  // y is state variables
  // paramvalues are "constants" in the diff eqn
  virtual const ColumnVector& evaluate(float x,const ColumnVector& y,const ColumnVector& paramvalues) const = 0;

  virtual ~Derivative(){};

protected:
  int ny;
  mutable ColumnVector dy;
};

void rk(ColumnVector& ret, const ColumnVector& y, const ColumnVector& dy, float x, float h, const Derivative& deriv,const ColumnVector& paramvalues);

void rkqc(ColumnVector& y, float& x, float& hnext, ColumnVector& dy, float htry, float eps, const Derivative& deriv,const ColumnVector& paramvalues);

void runge_kutta(Matrix& yp, ColumnVector& xp, ColumnVector& hp, const ColumnVector& ystart, float x1, float x2, float eps, float hmin, const Derivative& deriv,const ColumnVector& paramvalues);

}

#endif
