/*  cspline
    
    Cubic spline fitting and interpolation
    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(__cspline_h)
#define __cspline_h

#include <string>
#include <iostream>
#include <fstream>

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"

#define WANT_STREAM
#define WANT_MATH

using namespace NEWMAT;
using namespace std;
///////////////////////////////////////////////////////


namespace MISCMATHS { 
  class Cspline{
  public: 
    Cspline(){}
    Cspline(ColumnVector& pnodes,ColumnVector& pvals):
      nodes(pnodes),
      vals(pvals),
      n(nodes.Nrows())
    {
      fit();
      fitted=true;
    }

    Cspline(ColumnVector& pnodes, Matrix& pcoefs) : 
      nodes(pnodes),
      coefs(pcoefs),
      n(nodes.Nrows())
    { fitted=true;}

    ~Cspline(){
      fitted=false;
    };
    
    void set(ColumnVector& pnodes,ColumnVector& pvals);
    void set(ColumnVector& pnodes, Matrix& pcoefs);
    
    void fit();
    float interpolate(float xx) const;
    float interpolate(float xx,int ind) const;
    ColumnVector interpolate(const ColumnVector& x) const;
    ColumnVector interpolate(const ColumnVector& x, const ColumnVector& indvec) const;
    
  protected:

    bool fitted;
    ColumnVector nodes;
    ColumnVector vals;
    Matrix coefs;
    int n;
    void diff(const ColumnVector& x, ColumnVector& dx );

  };
}


#endif
