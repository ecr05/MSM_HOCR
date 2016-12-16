/*  rungekutta.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 2002 University of Oxford  */

/*  CCOPYRIGHT  */

#include "rungekutta.h"

using namespace std;

namespace MISCMATHS {

void rk(ColumnVector& ret, const ColumnVector& y, const ColumnVector& dy, float x, float h, const Derivative& deriv,const ColumnVector& paramvalues)
{ 
  Tracer tr("rk"); 

  float hh=h*0.5;
  float xh=x+hh;
  
  //first step
  ColumnVector yt=y+hh*dy;
  
  //second step
  ColumnVector dyt = deriv.evaluate(xh,yt,paramvalues);
  yt=y+hh*dyt;
  
  //third step
  ColumnVector dym = deriv.evaluate(xh,yt,paramvalues);
  yt=y+h*dym;
  dym=dym+dyt;
  
  //fourth step
  dyt = deriv.evaluate(x+h,yt,paramvalues);
  
  //addup
  ret = y+h*(dy+dyt+2*dym)/6;
}

void rkqc(ColumnVector& y, float& x, float& hnext, ColumnVector& dy, float htry, float eps, const Derivative& deriv,const ColumnVector& paramvalues)
{
  Tracer tr("rkqc"); 

  float xsav = x;
  ColumnVector dysav = dy;
  ColumnVector ysav = y;
  float h = htry;
  float hdid;
  ColumnVector ytemp;

  while(true)
    {
      // take 2 1/2 step sizes
  
      // first 1/2 step
      float hh=h*0.5;

      rk(ytemp,ysav,dysav,xsav,hh,deriv,paramvalues);
  
      // second 1/2 step
      x=xsav+hh;
      dy = deriv.evaluate(x,ytemp,paramvalues);
      rk(y,ytemp,dysav,xsav,hh,deriv,paramvalues);
      x=xsav+h;
      if(x==xsav) cerr << "step size too small" << endl;

      // take large step size
      rk(ytemp,ysav,dysav,xsav,h,deriv,paramvalues);
   
      // eval accuracy
      float errmax = 0.0;
      for(int i=1; i<=y.Nrows(); i++)
	{
	  //errmax=max(abs((y-ytemp)./y));
	  
	  float tmp = fabs((y(i)-ytemp(i))/y(i));
	  if(tmp > errmax) errmax = tmp;
	}

      errmax=errmax/eps;
      
      if(errmax <=1.0) 
	{
	  // step OK, compute step size for next step
	  hdid=h;
	  
	  if(errmax>6e-4)
	    hnext=h*std::exp(-0.2*std::log(errmax));
	  else
	    hnext=4*h;
	  
	  break;
      }
      else 
	{
	  // step too large,
	  h=h*std::exp(-0.25*std::log(errmax));
	}
    }

  y = y+(y-ytemp)/15;
}

void runge_kutta(Matrix& yp, ColumnVector& xp, ColumnVector& hp, const ColumnVector& ystart, float x1, float x2, float eps, float hmin, const Derivative& deriv,const ColumnVector& paramvalues)
{
  Tracer tr("runge_kutta"); 

  int MAXSTEP=1000;

  ColumnVector y = ystart;

  float x=x1;
  xp.ReSize(MAXSTEP,1);
  xp = 0;
  xp(1) =x1;

  float h=hp(1);
  hp.ReSize(MAXSTEP,1);
  hp = 0;

  yp.ReSize(MAXSTEP,y.Nrows());
  yp = 0;

  int kout=1;
  
  ColumnVector dy;

  for(int k=1; k <= MAXSTEP; k++)
    { 
      dy = deriv.evaluate(x,y,paramvalues);

      // store results:
      xp(kout)=x;
      yp.Row(kout)=y;
      hp(kout)=h;
  
      kout=kout+1;
    
      // stop overshoot of step past x2:
      if((x+h-x2)*(x+h-x1)>0) h=x2-x;

      float hnext = 0.0;
      rkqc(y,x,hnext,dy,h,eps,deriv,paramvalues);

      if((x-x2)*(x2-x1) >= 0.0)
      {
	xp(kout)=x;
	yp.Row(kout)=y;
	hp(kout)=h;
	//kout=kout+1;
	
        xp = xp.Rows(1,kout);
	yp = yp.Rows(1,kout);

	return;
      }
      else
      {
	if(hnext<=hmin) cerr << "step size too small" << endl;
	h=hnext;
      } 
      
    }
  cerr << "too many steps" << endl;
}

}
