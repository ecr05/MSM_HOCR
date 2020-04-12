/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT  */

#ifndef _point
#define _point

#include <iostream>
#include <cmath>
#include "newmat.h"

using namespace std;
using namespace NEWMAT;

namespace NEWMESH{

  class Pt {

  public:

  double X;
  double Y; 
  double Z;
  
  Pt() : X(0), Y(0), Z(0){};
  Pt(double x, double y, double z) : X(x), Y(y), Z(z){};
  Pt (const Pt& p) : X(p.X), Y(p.Y), Z(p.Z){};


  inline Pt operator =(const Pt& p)
    {
	X = p.X;
	Y = p.Y;
	Z = p.Z;
	return *this;
    }
    
    
  
    inline void operator+=(const Pt p)      
    {
      X=X+p.X;
      Y=Y+p.Y;
      Z=Z+p.Z;
    }

    inline void operator*=(const double d)    
    {
      X*=d;
      Y*=d;
      Z*=d;
    }

  inline void operator/=(const double d) 
    {
      if (d!=0)
	{
	  X/=d;
	  Y/=d;
	  Z/=d;
	}
      else cerr << "division by zero" << endl;
    }

   const inline double norm() const /// transferred from Vec by emma
   {
      return (sqrt(X*X + Y*Y + Z*Z));
   }

  void normalize(){      /// transferred from Vec by emma
    double n = norm();
    if (n!=0){
    X/=n;
    Y/=n;
    Z/=n;}
  }

  inline bool operator==(const Pt &p) const
    {
      return((fabs(X-p.X)<1e-8) && (fabs(Y-p.Y)<1e-8) && (fabs(Z-p.Z)<1e-8));
    }

};

// Pt operators added by emma

  const double operator|(const Pt &p1, const Pt &p2);    // dot product
  const Pt operator*(const Pt &p1, const Pt &p2);        // cross product
  const Pt operator*(const Pt &v, const double &d);
  const Pt operator/(const Pt &v, const double &d);   //NB all these were const in original meshclass, but I don't think they need/should be
  const Pt operator -(const Pt&p1, const Pt &p2);
  const Pt operator +(const Pt&p1, const Pt &p2);
  const Pt operator*(const Matrix &M,const Pt &p1); 
  const Pt operator*(const Pt &p1, const Matrix &M); 
}

#endif
