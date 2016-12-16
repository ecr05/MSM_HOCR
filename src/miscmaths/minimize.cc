/*  minimize
 
    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <string>
#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <vector>
#include <algorithm>
#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"
#include "minimize.h"

#define WANT_STREAM
#define WANT_MATH

using namespace NEWMAT;
using namespace std;
///////////////////////////////////////////////////////

namespace MISCMATHS {

float diff1(const ColumnVector& x, const EvalFunction& func, int i,float h,int errorord)
{
  //computes the first derivative of "eval" wrt the i^th parameter at point "x" with step size h
  ColumnVector xtmp=x;
  float deriv=0;
  if(errorord==1){
    xtmp(i)=xtmp(i)+h;
    float en_plus=func.evaluate(xtmp);
    float en=func.evaluate(x);
    deriv=(en_plus-en)/h;
  }
  else if(errorord==2){
    xtmp(i)=xtmp(i)+h;
    float en_plus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-(2*h);
    float en_minus=func.evaluate(xtmp);
    deriv=(en_plus-en_minus)/(2*h);
  }
  else{
    xtmp(i)=xtmp(i)+(2*h);
    float en_2plus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;
    float en_plus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-(2*h);
    float en_minus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;
    float en_2minus=func.evaluate(xtmp);
    deriv=(-en_2plus+8*en_plus-8*en_minus+en_2minus)/(12*h);
  }
  return deriv;
}

float diff2(const ColumnVector& x, const EvalFunction& func, int i,float h,int errorord)
{
  //computes the second derivative of "eval" wrt the i^th parameter at point "x" with step size h 
   ColumnVector xtmp=x;
   float deriv=0;
   if(errorord==1){
     xtmp(i)=xtmp(i)+(2*h);
     float en_2plus=func.evaluate(xtmp);
     xtmp(i)=xtmp(i)-h;
     float en_plus=func.evaluate(xtmp);
     float en=func.evaluate(x);
     deriv=(en_2plus-2*en_plus+en)/(h*h);
   }
   else if(errorord==2){
     xtmp(i)=xtmp(i)+h;
     float en_plus=func.evaluate(xtmp);
     xtmp(i)=xtmp(i)-(2*h);
     float en_minus=func.evaluate(xtmp);
     float en=func.evaluate(x);
     deriv=(en_plus-2*en+en_minus)/(h*h);
   }
   else{
     xtmp(i)=xtmp(i)+(2*h);
     float en_2plus=func.evaluate(xtmp);
     xtmp(i)=xtmp(i)-h;
     float en_plus=func.evaluate(xtmp);
     xtmp(i)=xtmp(i)-(2*h);
     float en_minus=func.evaluate(xtmp);
     xtmp(i)=xtmp(i)-h;
     float en_2minus=func.evaluate(xtmp);
     float en=func.evaluate(x);
     deriv=(-en_2plus+16*en_plus-30*en+16*en_minus-en_2minus)/(12*h*h);
   }
   return deriv;
}

float diff2(const ColumnVector& x, const EvalFunction& func, int i,int j,float h,int errorord)
{//computes the cross derivative of "eval" wrt the i^th and j^th  parameter at point "x" with step size h 
  ColumnVector xtmp=x;
  float deriv=0;
  if(errorord==1){ 
    xtmp(i)=xtmp(i)+h; xtmp(j)=xtmp(j)+h;
    float  en_iplus_jplus=func.evaluate(xtmp);
    xtmp(j)=xtmp(j)-h;
    float en_iplus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;xtmp(j)=xtmp(j)+h;
    float en_jplus=func.evaluate(xtmp);
    float en=func.evaluate(x);
    deriv=(en_iplus_jplus-en_iplus-en_jplus+en)/(h*h);}
  else if(errorord==2){
    xtmp(i)=xtmp(i)+h; xtmp(j)=xtmp(j)+h;
    float  en_iplus_jplus=func.evaluate(xtmp);
    xtmp(j)=xtmp(j)-2*h;
    float en_iplus_jminus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-2*h;xtmp(j)=xtmp(j)+2*h;
    float en_iminus_jplus=func.evaluate(xtmp);
    xtmp(j)=xtmp(j)-2*h;
    float en_iminus_jminus=func.evaluate(xtmp);
    deriv=(en_iplus_jplus-en_iplus_jminus-en_iminus_jplus+en_iminus_jminus)/(4*h*h); 
  }
  else{
    xtmp(i)=xtmp(i)+2*h;xtmp(j)=xtmp(j)+2*h;
    float en_i2plus_j2plus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;
    float en_iplus_j2plus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-2*h;
    float en_iminus_j2plus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;
    float en_i2minus_j2plus=func.evaluate(xtmp);
    xtmp(j)=xtmp(j)-h;
    float en_i2minus_jplus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)+h;
    float en_iminus_jplus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)+2*h;  
    float en_iplus_jplus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)+h;
    float en_i2plus_jplus=func.evaluate(xtmp);
    xtmp(j)=xtmp(j)-2*h;
    float en_i2plus_jminus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;
    float en_iplus_jminus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-2*h;
    float en_iminus_jminus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)-h;
    float en_i2minus_jminus=func.evaluate(xtmp);
    xtmp(j)=xtmp(j)-h;
    float en_i2minus_j2minus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)+h;
    float en_iminus_j2minus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)+2*h;
    float en_iplus_j2minus=func.evaluate(xtmp);
    xtmp(i)=xtmp(i)+h;
    float en_i2plus_j2minus=func.evaluate(xtmp);
    deriv=(en_i2plus_j2plus-8*en_iplus_j2plus+8*en_iminus_j2plus-en_i2minus_j2plus
	   -8*en_i2plus_jplus+64*en_iplus_jplus-64*en_iminus_jplus+8*en_i2minus_jplus
	   +8*en_i2plus_jminus-64*en_iplus_jminus+64*en_iminus_jminus-8*en_i2minus_jminus
	   -en_i2plus_j2minus+8*en_iplus_j2minus-8*en_iminus_j2minus+en_i2minus_j2minus)/(144*h*h);
  }
  return deriv;
}

ReturnMatrix gradient(const ColumnVector& x, const EvalFunction& func, float h,int errorord){
  ColumnVector deriv(x.Nrows());
  deriv = 0;

  for(int i=1;i<=x.Nrows();i++){
    deriv(i) = diff1(x,func,i,h,errorord);
  }

  deriv.Release();
  return deriv;
}

ReturnMatrix hessian(const ColumnVector& x, const EvalFunction& func, float h,int errorord)
{ //evaluates the hessian of function "eval" at x in parameter space
  //errorord=4 requires something like 8n^2-3n evaluations
  //errorord=2 requires something like 2n^2+n evaluations
  //errorord=1 requires same as errorord=2. no point really.

  // NB Hessian will compute _all_ derivatives even if non_varying parameters exist. The user must prune out rows/colums that are non needed.



  SymmetricMatrix hess(x.Nrows());
  for(int i=1;i<=x.Nrows();i++){
    for(int j=1;j<=i;j++){
      if(i!=j) hess(i,j)=diff2(x,func,i,j,h,errorord);
      else hess(i,j)=diff2(x,func,i,h,errorord);
    }
  }
  hess.Release();
  return hess;
  
}


void minsearch(ColumnVector& x, const EvalFunction& func, ColumnVector& paramstovary){


  //perform generic function minimization without gradient info


  // Number of nonvarying parameters
  int n_nonvary=0;
  for(int i=1;i<=paramstovary.Nrows();i++){
    if(paramstovary(i)>0){
      paramstovary(i)=1;
    }
    else{
      paramstovary(i)=0;
      n_nonvary++;
    }
  }
  
  
  //Number of parameters to estimate
  int n=x.Nrows()-n_nonvary, maxiter=200*n,iter=0;
  int ntot=x.Nrows();
  int func_evals=0;  
  // Some things we'll need.
  float rho=1,chi=2,psi=0.5,sigma=0.5;
  float tolx=1e-6,tolf=1e-6;
  ColumnVector onesn(ntot);
  onesn=1;
  ColumnVector one2n(ntot),two2np1(ntot);
  for(int i=1;i<=ntot;i++){
    one2n(i)=i;
    two2np1(i)=i+1;
  }
  
  // We want to store the best n+1 parameter estimates
  // I'm going to store them as a vector of pairs of floats and ColVecs 
  // so I can sort them based on energy
  
  vector<pair<float, ColumnVector> > v;
  float en=func.evaluate(x);
  func_evals++;
  pair<float, ColumnVector> tmppair;
  tmppair.first=en;
  tmppair.second=x;
  v.push_back(tmppair);
  
  float usual_delta=0.05,zero_term_delta=0.00025;
  //perturb each parameter by a bit, and store the cost.
  ColumnVector y=x;

  for(int i=1;i<=ntot;i++){
    // The values of nonvarying parameters should be the same in 
    // all of the optional param vectors and therefore in all
    // combinations of them in the remainder of the code.
    
    if(paramstovary(i)==1){
      if(y(i)!=0){y(i)=(1+usual_delta)*y(i);}
      else{y(i)=(1+zero_term_delta);}
      en=func.evaluate(y);
      
      func_evals++;
      tmppair.first=en;
      tmppair.second=y;
      v.push_back(tmppair);
    }
    
  }

  
  sort(v.begin(),v.end(),pair_comparer()); //wasn't that easy...
  string how="";
  ColumnVector xbar(ntot),xr(ntot),xe(ntot),xc(ntot),xcc(ntot),xtmp(ntot);
  //cerr<<"starting loop"<<endl;
  while(iter<=maxiter){
    iter++;
    if(v[n].first-v[0].first< tolf){
      ColumnVector tmpvec1,tmpvec2;
      bool stopsearch=true;
      for(int i=0;i<n;i++){//iterate over paramsets
	tmpvec1=v[i].second;
	tmpvec2=v[i+1].second;
	for(int j=1;j<=n;j++){//iterate over n params
	  if(fabs( tmpvec1(j)-tmpvec2(j) ) >tolx){stopsearch=false;}
	}
      }
      if(stopsearch){break;}
    } 
    //compute reflection point
  
  // xbar is average of best n paramsets.
    xbar=0;
    for(int i=0;i<n;i++){
      xbar=xbar+v[i].second;
    }
    xbar=xbar/n;
    xr=(1+rho)*xbar-rho*v[n].second; //reflection point
    float en_xr=func.evaluate(xr);func_evals++;
    if(en_xr < v[0].first){ //en_xr is better than our current best
 
      //compute expansion point
      xe=(1+rho*chi)*xbar-rho*chi*v[n].second;
      float en_xe=func.evaluate(xe);func_evals++;
    
      if(en_xe<en_xr){
	tmppair.first=en_xe;
	tmppair.second=xe;
	v[n]=tmppair;
	how="expand";
      }
      else{ //en_xr < en_xe
	tmppair.first=en_xr;
	tmppair.second=xr;
	v[n]=tmppair;
	how="reflect1";
      }
    }
    else{ //en_xr is worse than our current best
      if(en_xr<=v[n-1].first){ //en_xr is better than our current second worst
	tmppair.first=en_xr;
	tmppair.second=xr;
	v[n]=tmppair;
	how="reflect2";
      }
      else{//en_xr is worse than our current secind worst
	//perform contraction

	if(en_xr<v[n].first){//en_xr better than current worst
	  //perform outside contraction
	  xc=(1+rho*psi)*xbar-rho*psi*v[n].second;
	  float en_xc = func.evaluate(xc); func_evals++;

	  if(en_xc<=v[n].first){ //en_xr is better than our current worst
	    tmppair.first=en_xc;
	    tmppair.second=xc;
	    v[n]=tmppair;
	    how="contract outside";
	  }
	  else{ //xc no good
	    //perform a shrink
	    how="shrink";
	  }
	}
	else{//en_xr worse than currenst worst
	  //perform inside contraction
	  xcc = (1-psi)*xbar + psi*v[n].second;
	  float en_xcc=func.evaluate(xcc);func_evals++;
	  if(en_xcc<v[n].first){
	    tmppair.first=en_xcc;
	    tmppair.second=xcc;
	    v[n]=tmppair;
	    how="contract inside";
	  }
	  else{
	    //perform a shrink
	    how="shrink";
	  }
	}
	if(how=="shrink"){
	  for(int i=1;i<n+1;i++){
	    tmppair.second=v[0].second+sigma*(v[i].second-v[0].second);
	    tmppair.first=func.evaluate(xtmp);func_evals++;
	    v[i]=tmppair;
	  }
	}


      }
    }

    sort(v.begin(),v.end(),pair_comparer()); //double bracks constructs a temporary object


  } //closing while
 
  x=v[0].second;
    
}

void scg(ColumnVector& x,const gEvalFunction& func, ColumnVector& paramstovary, float tol, float eps, int niters){
  
  //number of nonvarying parameters.
  int n_nonvary=0;
  for(int i=1;i<=paramstovary.Nrows();i++){
    if(paramstovary(i)>0){
      paramstovary(i)=1;
    }
    else{
      paramstovary(i)=0;
      n_nonvary++;
    }
  }
  
  int fevals=0;
  int gevals=0;
  int nparams=x.Nrows();
  float sigma0 = 1.0e-4;
  float  fold=func.evaluate(x); fevals++;
  float fnow=0,fnew=0;
  //Comput gradient wrt all parameters which we wish to estimate
  ColumnVector gradold=func.g_evaluate(x);gevals++;
  gradold=SP(gradold,paramstovary);
  
  ColumnVector gradnew=gradold;
  ColumnVector d=-gradnew;       // search direction
  ColumnVector xplus,xnew,gplus;
  bool success=true;
  int nsuccess=0;
  float lambda=1.0;
  float lambdamin = 1.0e-15; 
  float lambdamax = 1.0e15;			
  int j = 1;	
  float mu=0,kappa=0,sigma=0,gamma=0,alpha=0,delta=0,Delta,beta=0;

  // main loop..
  while(j<niters){

    if(success){
      mu=(d.t()*gradnew).AsScalar();
      if(mu >= 0){
	d=-gradnew;
	mu=(d.t()*gradnew).AsScalar();
      }

      kappa=(d.t()*d).AsScalar();
      if(kappa<eps){
	break;
      }

      sigma=sigma0/std::sqrt(kappa);      
      xplus = x + sigma*d;

      gplus=func.g_evaluate(xplus);gevals++;
      gplus=SP(gplus,paramstovary);
	
      gamma = (d.t()*(gplus - gradnew)).AsScalar()/sigma;

    }

    delta = gamma + lambda*kappa;
    if (delta <= 0){ 
      delta = lambda*kappa;
      lambda = lambda - gamma/kappa;
    }
    alpha = - mu/delta;
    
    xnew = x + alpha*d;
    fnew=func.evaluate(xnew);fevals++;
    Delta = 2*(fnew - fold)/(alpha*mu);

    if (Delta  >= 0){
      success = true;
      nsuccess = nsuccess + 1;
       x = xnew;
       fnow = fnew;}
    else{
      success = false;
      fnow = fold;
    }    
    
    if (success == 1){
      
      //Test for termination...
      
      if ( (max(abs(d*alpha))).AsScalar() < tol && std::abs(fnew-fold) < tol){
	break;
      }
      else{
	fold = fnew;
	gradold = gradnew;
	
	gradnew=func.g_evaluate(x);gevals++;
	gradnew=SP(gradnew,paramstovary);
	
	if ((gradnew.t()*gradnew).AsScalar() == 0){
	  break;
	}
      }
      
    }
    if (Delta < 0.25){
      //   lambda = min(4.0*lambda, lambdamax);
      lambda=4.0*lambda<lambdamax ? 4.0*lambda : lambdamax;
    }
     if (Delta > 0.75){
       //lambda = max(0.5*lambda, lambdamin);
       lambda = 0.5*lambda > lambdamin ? 0.5*lambda : lambdamin;
     }

     if (nsuccess == nparams){
       d = -gradnew;
       nsuccess = 0;
     }
     else{
       if (success == 1){	 
	 beta = ((gradold - gradnew).t()*gradnew).AsScalar()/mu;
	 d = beta*d - gradnew;
       }
     }
     j++; 
  }

}




}








