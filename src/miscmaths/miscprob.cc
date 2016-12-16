/*  miscprob.cc

    Christian Beckmann & Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

// Miscellaneous maths functions that rely on libprob

#include "miscprob.h"
#include "stdlib.h"
#include "newmatio.h"
#include <iostream>
// #include "gam.h"

using namespace NEWMAT;

namespace MISCMATHS {

//   ReturnMatrix betarnd(const int dim1, const int dim2, const float a, const float b)
//   {
//     // Devroye, L. (1986) Non-Uniform Random Variate Generation, Springer-Verlag.

//     int tdim = dim2;

//     if(tdim<0){tdim=dim1;}

    
//     Matrix g1=gammarnd(dim1, tdim, a, 1);
//     Matrix g2=gammarnd(dim1, tdim, b, 1);

//     Matrix res(dim1,tdim);
//     for (int mc=1; mc<=res.Ncols(); mc++) {
//       for (int mr=1; mr<=res.Nrows(); mr++) {
// 	res(mr,mc)=g1(mr,mc)/(g1(mr,mc)+g2(mr,mc));
//       }
//     }
    
//     res.Release();
//     return res;
//   }

ReturnMatrix betapdf(const RowVector& vals, const float a, const float b)
{

  RowVector res(vals);

  if(a<0 || b<0)
    {
      throw Exception("Negative a or b in call to Miscprob::betapdf");
    }

  for (int mc=1; mc<=res.Ncols(); mc++)
    {
      float x=vals(mc);
      if(x<0)
	{
	  res(mc)=0;
	}
      else
	{
	  float logkerna=(a-1)*std::log(x);
	  float logkernb=(b-1)*std::log(1-x);
	  float betaln_ab=lgam(a)+lgam(b)-lgam(a+b);
	  res(mc)=std::exp(logkerna+logkernb-betaln_ab);	 
	}
    }
  
  res.Release();
  return res;
}

ReturnMatrix unifrnd(const int dim1, const int dim2, const float start, const float end)
{
  int tdim = dim2;
  double tmpD=1.0;
  if(tdim<0){tdim=dim1;}
  Matrix res(dim1,tdim);
  for (int mc=1; mc<=res.Ncols(); mc++) {
    for (int mr=1; mr<=res.Nrows(); mr++) {
      tmpD = (rand()+1)/double(RAND_MAX+2.0);
      res(mr,mc)=(tmpD)*(end-start)+start;
      //drand(&tmpD);
      //res(mr,mc)=(tmpD-1)*(end-start)+start;
    }
  }

  res.Release();
  return res;
}

ReturnMatrix normrnd(const int dim1, const int dim2, const float mu, const float sigma)
{
  int tdim = dim2;
  double tmpD=1.0;
  if(tdim<0){tdim=dim1;}
  Matrix res(dim1,tdim);
  for (int mc=1; mc<=res.Ncols(); mc++) {
    for (int mr=1; mr<=res.Nrows(); mr++) {
      tmpD = (rand()+1)/double(RAND_MAX+2.0);
      res(mr,mc)=ndtri(tmpD)*sigma+mu ;
      //drand(&tmpD);
      //res(mr,mc)=ndtri(tmpD-1)*sigma+mu ;
    }
  }

  res.Release();

  return res;
}

ReturnMatrix normpdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  for (int mc=1; mc<=res.Ncols(); mc++){
    res(mc) = std::exp(-0.5*(std::pow(vals(mc)-mu,2)/var))*std::pow(2*M_PI*var,-0.5);
  }

  res.Release();
  return res;
}


ReturnMatrix normcdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  RowVector tmp;
  tmp = (vals-mu)/std::sqrt(var);
  for (int mc=1; mc<=res.Ncols(); mc++){
    res(mc) = ndtr(tmp(mc));
  }

  res.Release();
  return res;
}

ReturnMatrix gammacdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  res=0;
  if((mu>0)&&(var>0)){
    float b = std::pow(mu,2)/var;
    float a = mu/var;  
    for (int mc=1; mc<=res.Ncols(); mc++){
      if(vals(mc)>0)
	res(mc) = gdtr(a,b,vals(mc));
    }
  }
  res.Release();
  return res;
}

ReturnMatrix gammapdf(const RowVector& vals, const float mu, const float var)
{
  RowVector res(vals);
  res=0;
  if((mu>0)&&(var>0.00001)){
    float a = std::pow(mu,2)/var;
    float b = mu/var;
    float c = lgam(a);
    if(std::abs(c) < 150){
      for (int mc=1; mc<=res.Ncols(); mc++){
	if(vals(mc)>0.000001){
	  res(mc) = std::exp(a*std::log(b) + 
			     (a-1) * std::log(vals(mc)) 
			     - b*vals(mc) - c);
	}
      }
    }
  }
  res.Release();
  return res;
}

  float normpdf(const float val, const float mu, const float var)
  {
    return std::exp(-0.5*(std::pow(val-mu,2)/var))*std::pow(2*M_PI*var,-0.5);
  }
  float lognormpdf(const float val, const float mu, const float var)
  {
    return -0.5*(std::pow(val-mu,2)/var+std::log(2*M_PI*var));
  }
    
ReturnMatrix normpdf(const RowVector& vals, const RowVector& mu, const RowVector& var)
{
  Matrix res(mu.Ncols(),vals.Ncols());
  for (int mc=1; mc<=res.Ncols(); mc++){
    for (int mr=1; mr<=res.Nrows(); mr++){
      res(mr,mc) = std::exp(-0.5*(std::pow(vals(mc)-mu(mr),2)/var(mr)))*std::pow(2*M_PI*var(mr),-0.5);
    }
  }

  res.Release();
  return res;
}


ReturnMatrix mvnrnd(const RowVector& mu, const SymmetricMatrix& covar, int nsamp) 
{     
//   Matrix eig_vec; 
//   DiagonalMatrix eig_val;
//   EigenValues(covar,eig_val,eig_vec);

//   Matrix ret = ones(nsamp, 1)*mu + dnormrandm(nsamp,mu.Ncols())*sqrt(eig_val)*eig_vec.t();
  Mvnormrandm mvn(mu, covar);

  return mvn.next(nsamp);
}

 float mvnpdf(const RowVector& vals, const RowVector& mu, const SymmetricMatrix& covar)
 {
   if(vals.Ncols()==2)
     return bvnpdf(vals,mu,covar);
   else
     return std::exp(-0.5*((vals-mu)*covar.i()*(vals-mu).t()).AsScalar())/(std::pow(covar.Determinant(),0.5)*std::pow(2*M_PI,vals.Ncols()/2.0));
 }

 float bvnpdf(const RowVector& vals, const RowVector& mu, const SymmetricMatrix& covar)
 {
   // bivariate normal pdf
   double det=covar(1,1)*covar(2,2)-Sqr(covar(1,2));
   float m1=vals(1)-mu(1);
   float m2=vals(2)-mu(2);
   float ss=(Sqr(m1)*covar(2,2)-2*m1*m2*covar(1,2)+Sqr(m2)*covar(1,1))/det;

   return std::exp(-0.5*ss)/(std::pow(det,0.5)*std::pow(2*M_PI,vals.Ncols()/2.0));
 }

// ReturnMatrix gammarnd(const int dim1, const int dim2, 
// 			const float a, const float b)
// {
//   // Marsaglia, G. and Tsang, W.W. (2000) "A Simple Method for Generating Gamma Variables", Acm Trans. Math. Soft. 26(3):363-372.

//   int tdim = dim2;
//   if(tdim<0){tdim=dim1;}
//   Matrix res(dim1,tdim);

//   Gam& gam=Gam::getInstance();
//   gam.setParams(a,b);

//   for (int mc=1; mc<=res.Ncols(); mc++) {
//     for (int mr=1; mr<=res.Nrows(); mr++) {
//       res(mr,mc)=gam.rnd();
//     }
//   }
//   res.Release();
//   return res;
// }

ReturnMatrix perms(const int n){
  if(n<=1){
    Matrix P(1,1);
    P << n;
    P.Release();
    return P;
  }
  Matrix Q = perms(n-1);  // recursive calls
  int m = Q.Nrows();
  Matrix P(n*m,n);
  for(int i=1;i<=m;i++){
    P(i,1)=n;
    for(int j=1;j<=Q.Ncols();j++)
      P(i,j+1)=Q(i,j);
  }
  for(int i=n-1;i>=1;i--){
    int jj=1;
    for(int j=(n-i)*m+1;j<=(n-i+1)*m;j++){
      P(j,1)=i;
      for(int k=1;k<=n-1;k++){
	P(j,k+1)= (Q(jj,k)==i) ? n : Q(jj,k);
      }
      jj++;
    } 
  }
  P.Release();
  return P;
}

}
