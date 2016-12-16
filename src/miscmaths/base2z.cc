/*  base2z.cc

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <cmath>
#include "base2z.h"
#include "libprob.h"

namespace MISCMATHS {

  Base2z* Base2z::base2z = NULL;

  float Base2z::logbeta(float v, float w)
  {
    return MISCMATHS::lgam(v)+MISCMATHS::lgam(w)-MISCMATHS::lgam(v+w);
  }

  float Base2z::logp2largez(float logp) 
  {
    // Large Z extrapolation routine for converting log(p) to Z values
    //  written by Mark Jenkinson, March 2000
    //
    // Equations were derived by using integration by parts and give the
    //  following formulae:
    //   Z to log(p)
    //       log(p) = -1/2*z*z - 1/2*log(2*pi) - log(z) 
    //                + log(1 - 1/(z*z) + 3/(z*z*z*z))
    // this equation is then solved by the recursion:
    //   z_0 = sqrt(2*(-log(p) - 1/2*log(2*pi)))
    //   z_{n+1} = sqrt(2*(-log(p) - 1/2*log(2*pi) - log(z_n)  
    //             + log(1 - 1/(zn*zn) + 3/(zn*zn*zn*zn)) ))
    // In practice this recursion is quite accurate in 3 to 5 iterations
    // The equation is accurate to 1 part in 10^3 for Z>3.12  (3 iterations)


    static const float pi = 3.141592653590;
    static const float log2pi = log(2*pi);
    
    float z0, zn;
    // iteratively solve for z given log p
    float b = -2*logp - log2pi; 
    z0 = sqrt(b);
    zn = z0;
    for (int m=1; m<=3; m++) {
      // zn = sqrt(b + 2*log(1/zn - 1/(zn*zn*zn) + 3/(zn*zn*zn*zn*zn)));
      zn = sqrt(b + 2*log(((3/(zn*zn) - 1)/(zn*zn) + 1)/zn) );
    }
  
    return zn;
  }
  
  float Base2z::convertlogp2z(float logp) 
    {
      // logp must be the *natural* logarithm of p, not base 10
      float z = 0.0;
      
      if(!issmalllogp(logp)) {
	z = MISCMATHS::ndtri(exp(logp));
      }
      else {
	z = logp2largez(logp);
      }

      return z;
    }

}






























