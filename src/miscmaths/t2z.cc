/*  t2z.cc

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <cmath>
#include "t2z.h"
#include "newmat.h"
#include "utils/tracer_plus.h"
#include "libprob.h"

using namespace NEWMAT;
using namespace Utilities;

namespace MISCMATHS {

  T2z* T2z::t2z = NULL;
  Z2t* Z2t::z2t = NULL;
 
  float Z2t::convert(float z, int dof) 
    {
      float t = 0.0;

      if(z>8)
	throw Exception("z is too large to convert to t");

      double p = MISCMATHS::ndtr(z);
      cerr << "p = " << p << endl;
      t = MISCMATHS::stdtri(dof,p);

      return t;
    }

  float T2z::larget2logp(float t, int dof)
  {
    //    static float logbeta[] = {   1.144729885849, 0.693147180560,
    //			 0.451582705289, 0.287682072452, 
    //			 0.163900632838, 0.064538521138,
    //			 -0.018420923956, -0.089612158690, 
    //			 -0.151952316581, -0.207395194346 } ;

    //static const float pi = 3.141592653590;
    //    static const float log2pi = log(2*pi);

    // Large T extrapolation routine for converting T to Z values
    //  written by Mark Jenkinson, March 2000
    //
    // It does T to Z via log(p) rather than p, since p becomes very
    //  small and underflows the arithmetic
    // Equations were derived by using integration by parts and give the
    //  following formulae:
    //  (1) T to log(p)    NB: n = DOF
    //       log(p) = -1/2*log(n) - log(beta(n/2,1/2)) - (n-1)/2*log(1+t*t/n)
    //                + log(1 - (n/(n+2))/(t*t) + 3*n*n/((n+2)*(n+4)*t*t*t*t)) 
    //  (2) Z to log(p)
    //       log(p) = -1/2*z*z - 1/2*log(2*pi) - log(z) 
    //                + log(1 - 1/(z*z) + 3/(z*z*z*z))
    // equation (2) is then solved by the recursion:
    //   z_0 = sqrt(2*(-log(p) - 1/2*log(2*pi)))
    //   z_{n+1} = sqrt(2*(-log(p) - 1/2*log(2*pi) - log(z_n) 
    //             + log(1 - 1/(zn*zn) + 3/(zn*zn*zn*zn))
    // In practice this recursion is quite accurate in 3 to 5 iterations
    // Equation (1) is accurate to 1 part in 10^3 for T>7.5  (any n)
    // Equation (2) is accurate to 1 part in 10^3 for Z>3.12  (3 iterations)


    if (t<0) {
      return larget2logp(-t,dof);
    }
  
    float logp, lbeta;

    if (dof<=0) { 
      cerr << "DOF cannot be zero or negative!" << endl;
      return 0.0; 
    }

    float n = (float) dof;

    // complete Beta function
    lbeta = this->logbeta(1/2.0,n/2.0);
    //if (dof<=10) {
    //lbeta = logbeta[dof-1];
    //} else {
    //lbeta = log2pi/2 - log(n)/2 + 1/(4*n);
    //}
    
    // log p from t value
    // logp = log( (1 - n/((n+2)*t*t) + 3*n*n/((n+2)*(n+4)*t*t*t*t))/(sqrt(n)*t))
    //          - ((n-1)/2)*log(1 + t*t/n) - lbeta;
    logp = log(( (3*n*n/((n+2)*(n+4)*t*t) - n/(n+2))/(t*t) + 1)/(sqrt(n)*t))
      - ((n-1)/2)*log(1 + t*t/n) - lbeta;

    return logp;
  }

  bool T2z::islarget(float t, int dof, float &logp) {
    // aymptotic formalae are valid if 
    //   log(p) < -14.5  (derived from Z-statistic approximation error)
    // For dof>=15, can guarantee that log(p)>-33 (p > 1e-14) if T<7.5
    //  and so in this region use conventional means, not asymptotic
    if ((dof>=15) && (fabs(t)<7.5)) { return false; }
    logp=larget2logp(t,dof);
    if (dof>=15)  return true;  // force asymptotic calc for all T>=7.5, D>=15
    return issmalllogp(logp);
  }

  bool T2z::issmalllogp(float logp) {
	// aymptotic formula accurate to 1 in 10^3 for Z>4.9
	// which corresponds to log(p)=-14.5
	return (logp < -14.5);
      }

  float T2z::convert(float t, int dof) {

      float z = 0.0, logp=0.0;
      
      if(!islarget(t,dof,logp)) {
	//	cerr << "t = " << t << endl;
	double p = MISCMATHS::stdtr(dof, t);
	//cerr << "p = " << p << endl;
	z = MISCMATHS::ndtri(p);
      }
      else {
	
	z = logp2largez(logp);

	//	cerr<<endl<<"logp="<<logp<<endl;

	if (t<0) z=-z;
      }

      return z;

    }


  float T2z::converttologp(float t, int dof) 
    {
      float logp=0.0;
      
      if(!islarget(t,dof,logp)) {
	logp = log(1-MISCMATHS::stdtr(dof, t));
      }
      else if(t<0) {
	// t < 0 and abs(t) is large enough to require asymptotic approx.
	// but t to logp is not available for negative t 
	// so just hardcode it to be -1e-12
	logp=-1e-12;
      }

//       cerr << "logp = " << logp << endl;
//       cerr << "exp(logp) = " << std::exp(logp) << endl;

      return logp;
    }

  void T2z::ComputePs(const ColumnVector& p_vars, const ColumnVector& p_cbs, int p_dof, ColumnVector& p_ps)
    {
      Tracer ts("T2z::ComputePs");

      int numTS = p_vars.Nrows();

      T2z& t2z = T2z::getInstance();

      p_ps.ReSize(numTS);

      for(int i = 1; i <= numTS; i++)
	{
	  //cerr << "var = " << p_vars(i) << " at index "<< i << endl;
	  //cerr << "cb = " << p_cbs(i) << " at index "<< i << endl;
	  if ( (p_vars(i) != 0.0) && (p_cbs(i) != 0.0) )
	    {
	      if(p_vars(i) < 0.0)
		{
		  //cerr << "var = " << p_vars(i) << " at index "<< i << endl;
		  p_ps(i) = 0.0;
		}
	      else
		{
		  p_ps(i) = t2z.converttologp(p_cbs(i)/sqrt(p_vars(i)),p_dof);
		  
		  //if(p_zs(i) == 0.0)
		  //cerr << " at index " << i << endl;
		}
	    }
	  else
	    p_ps(i) = 0.0;
	}
    }

  void T2z::ComputeZStats(const ColumnVector& p_vars, const ColumnVector& p_cbs, int p_dof, ColumnVector& p_zs)
    {
      ColumnVector dof = p_vars;
      dof = p_dof;
      ComputeZStats(p_vars,p_cbs,dof,p_zs);
    }

  void T2z::ComputeZStats(const ColumnVector& p_vars, const ColumnVector& p_cbs, const ColumnVector& p_dof, ColumnVector& p_zs)
    {
      Tracer ts("T2z::ComputeStats");

      int numTS = p_vars.Nrows();

      T2z& t2z = T2z::getInstance();

      p_zs.ReSize(numTS);

      for(int i = 1; i <= numTS; i++)
	{
	  //cerr << "var = " << p_vars(i) << " at index "<< i << endl;
	  //cerr << "cb = " << p_cbs(i) << " at index "<< i << endl;
	  if ( (p_vars(i) != 0.0) && (p_cbs(i) != 0.0) )
	    {
	      if(p_vars(i) < 0.0)
		{
		  //cerr << "var = " << p_vars(i) << " at index "<< i << endl;
		  p_zs(i) = 0.0;
		}
	      else
		{
		  p_zs(i) = t2z.convert(p_cbs(i)/sqrt(p_vars(i)),int(p_dof(i)));
		  
		  //if(p_zs(i) == 0.0)
		  //cerr << " at index " << i << endl;
		}
	    }
	  else
	    p_zs(i) = 0.0;
	}
    }
}






























