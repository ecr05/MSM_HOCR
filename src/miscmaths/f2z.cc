/*  f2z.cc

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <cmath>
#include "f2z.h"
#include "utils/log.h"
#include "utils/tracer_plus.h"
#include <stdexcept>
#include "libprob.h"

using namespace NEWMAT;
using namespace Utilities;

namespace MISCMATHS {

  F2z* F2z::f2z = NULL;
 
  float F2z::largef2logp(float f, int d1, int d2)
  {
    Tracer_Plus ts("F2z::largef2logp");

    // no of iterations:
    int N = 20;

//     cout << f<< endl;
//     cout << d1<< endl;
//     cout << d2<< endl;

    if (f<=0.0) {
      cerr << "f cannot be zero or negative!" << endl;
      return 0.0;
    }
  
    if (d1<=0 || d2<=0) { 
      cerr << "DOFs cannot be zero or negative!" << endl;
      return 0.0; 
    }
 
    double alpha=d1/(double)d2;
    double m=(d1+d2)/2.0;
    double n=(1-d1/2.0);
    double loggam = (d1/2.0)*(::log(d1/(double)d2)-logbeta(d2/2.0,d1/2.0));

    //iter=f^(-n)/(alpha*(n+m-1)) + n*f^(-(n+1))/(alpha^2*(n+m-1)*(n+m)) + n*(n+1)*f^(-(n+2))/(alpha^3*(n+m-1)*(n+m)*(n+m+1));

    double top = 1.0;
    double bot = n+m-1;
    double iter = 0.0;
 
//     cerr << "logbeta(d2/2.0,d1/2.0)=" << logbeta(d2/2.0,d1/2.0) << endl;
//     cerr << "loggam = " << loggam << endl;
//     cerr << "n = " << n << endl;
//     cerr << "m = " << m << endl;

    for(int i = 1; i <= N; i++)
      {
	// cerr << "i=" << i;
		  iter = iter + top* ( std::pow( f,float(-(n+i-1)) ) / ( std::pow(alpha,double(i))*bot ) );	
	top = top*(n-1+i)*(-1);
	bot = bot*(n+m-1+i);
// 	cerr << "iter=" << iter;
      }


    if(iter <= 0) throw Exception("iter negative");

    float logp = loggam-(m-1)*(::log(1+alpha*f))+::log(iter);

//     cerr << "iter = " << iter << endl;
//     cerr << "logp = " << logp << endl;

    return logp;
  }

  bool F2z::islargef(float f, int d1, int d2, float &logp) {
   
    if(f > 2.0 && d1>1)
      {

	try
	  {
	    logp=largef2logp(f,d1,d2);	
	  }
	catch(Exception& p_excp) 
	  {
	    cerr << "Negative iter in F2z::largef2logp" << endl;
	    return false;
	  }

	return issmalllogp(logp);
      }
    else
      return false;
  }

  bool F2z::issmalllogp(float logp) {
    return (logp < -14.5);
  }

  float F2z::convert(float f, int d1, int d2) 
  {
    Tracer_Plus ts("F2z::convert");

    float z = 0.0, logp=0.0;

    if(!islargef(f,d1,d2,logp)) {

      double p = MISCMATHS::fdtr(d1, d2, f);

      z = MISCMATHS::ndtri(p);
    }
      else {

	z = logp2largez(logp);
      }

      return z;
    }

  void F2z::ComputeFStats(const ColumnVector& p_fs, int p_dof1, int p_dof2, ColumnVector& p_zs)
  {
    ColumnVector dof2 = p_fs;
    dof2 = p_dof2;
    ComputeFStats(p_fs,p_dof1,dof2,p_zs);
  }
  
  void F2z::ComputeFStats(const ColumnVector& p_fs, int p_dof1, const ColumnVector& p_dof2, ColumnVector& p_zs)
  {
    Tracer_Plus ts("F2z::ComputeFStats");
    
    int numTS = p_fs.Nrows();

    p_zs.ReSize(numTS);
    F2z& f2z = F2z::getInstance();
    
    for(int i = 1; i <= numTS; i++)
      {		  	
	if (p_fs(i) > 0.0)
	  {

// 	    cerr << "i=" << i;
// 	    cerr << ",p_fs(i)=" << p_fs(i);
// 	    cerr << ",p_dof1=" << p_dof1;
// 	    cerr << ",p_dof2=" << p_dof2(i) << endl;

	    p_zs(i) = f2z.convert(p_fs(i),int(p_dof1),int(p_dof2(i)));  
	  }
	else
	  {
	    p_zs(i) = 0.0;
	  }     
      }
  }
   void F2z::ComputeFStats(const ColumnVector& p_fs, const ColumnVector& p_dof1, const ColumnVector& p_dof2, ColumnVector& p_zs)
  {
    Tracer_Plus ts("F2z::ComputeFStats");
    
    int numTS = p_fs.Nrows();

    p_zs.ReSize(numTS);
    F2z& f2z = F2z::getInstance();
    
    for(int i = 1; i <= numTS; i++)
      {		  	
	if (p_fs(i) > 0.0)
	  {

// 	    cerr << "i=" << i;
// 	    cerr << ",p_fs(i)=" << p_fs(i);
// 	    cerr << ",p_dof1=" << p_dof1;
// 	    cerr << ",p_dof2=" << p_dof2(i) << endl;

	    p_zs(i) = f2z.convert(p_fs(i),int(p_dof1(i)),int(p_dof2(i)));  
	  }
	else
	  {
	    p_zs(i) = 0.0;
	  }     
      }
  }
  
}






























