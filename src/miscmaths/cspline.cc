/*  cspline
   
    Cubic spline fitting and interpolation
    
    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <string>
#include <iostream>
#include <fstream>

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"
#include "cspline.h"

#define WANT_STREAM
#define WANT_MATH

using namespace NEWMAT;
using namespace std;

///////////////////////////////////////////////////////



namespace MISCMATHS{

  //  void Cspline::Cspline(){}
  
  void Cspline::set(ColumnVector& pnodes,ColumnVector& pvals){
    nodes=pnodes;vals=pvals;
    fitted=false;
    n=vals.Nrows();
  }
  void Cspline::set(ColumnVector& pnodes, Matrix& pcoefs){
    nodes=pnodes;coefs=pcoefs;
    fitted=false;
    n=vals.Nrows();
  }
  
  void Cspline::diff(const ColumnVector& x, ColumnVector& dx ){
    // dx should be of length length(x)-1
    dx.ReSize(x.Nrows()-1);
    for(int i=2;i<=x.Nrows();i++){
      dx(i-1)=x(i)-x(i-1);
    }						    
  }
  
  
  void Cspline::fit(){
    if(vals.Nrows()<4){
      cerr<<"Cspline::fit - You have less than 4 data pts for spline fitting."<<endl;
      exit(-1);
    }
    if(nodes.Nrows()!=vals.Nrows()){
      cerr<<"Nodes and VALS must be the same length in your spline"<<endl;
      exit(-1);
    }
    int n=vals.Nrows();
    ColumnVector s(n);
    ColumnVector dx,dy,dydx(n-1);
    diff(nodes,dx);
    diff(vals,dy);
    
    for(int i=1;i<=n-1;i++){
      dydx(i)=dy(i)/dx(i);
    }
    ColumnVector b(n);
    b=0;
    for(int i=2;i<=b.Nrows()-1;i++){
      b(i)=3*(dx(i)*dydx(i-1)+dx(i-1)*dydx(i));
    }
    
    float x31=nodes(3)-nodes(1),xn=nodes(n)-nodes(n-2);
    b(1)=((dx(1)+2*x31)*dx(2)*dydx(1)+dx(1)*dx(1)*dydx(2))/x31;
    b(n)=(dx(n-1)*dx(n-1)*dydx(n-2)+(2*xn+dx(n-1))*dx(n-2)*dydx(n-1))/xn;
    
    Matrix tridiag(n,n);
    tridiag=0;
    ColumnVector  y3(n);
    for(int j=2;j<=n-1;j++){
      tridiag(j,j-1)=dx(j);
      tridiag(j,j)=2*(dx(j)+dx(j-1));
      tridiag(j,j+1)=dx(j-1);
    }
    tridiag(1,1)=dx(2);tridiag(1,2)=x31;
    tridiag(n,n-1)=xn;tridiag(n,n)=dx(n-2);
    s=tridiag.i()*b; 
    
    
    ColumnVector d(n-1),c(n-1);
    for(int j=1;j<n;j++){
      d(j)=(s(j)+s(j+1)-2*dydx(j))/dx(j);
      c(j)=(dydx(j)-s(j))/dx(j)-d(j);
    }
    
    coefs.ReSize(n-1,4);
    for(int j=1;j<n;j++){
      coefs(j,1)=vals(j);
      coefs(j,2)=s(j);
      coefs(j,3)=c(j);
      coefs(j,4)=d(j)/dx(j);
    }
    fitted=true;
  }
  
 
  float Cspline::interpolate(float xx) const{
    // nodes must be monotonically increasing. I don't check this. 
    // On your head be it if you don't.
    if(nodes.Nrows()!=vals.Nrows()){
      cerr<<"Cspline:interpolate: Nodes and Vals should be the same length"<<endl;
      exit(-1);
    }

    float ret;
    if(!fitted){
      cerr<<"Cspline::interpolate - Cspline has not been fitted"<<endl;
      exit(-1);
    }
    else{
      
      bool stop=false;
      int ind=0;
      
      if(xx<nodes(1)){
	ind=1;
      }
      else if(xx>nodes(nodes.Nrows())){
	ind=nodes.Nrows()-1;
      }
      else{
	for(int i=1;i<nodes.Nrows();i++){
	  if(!stop){
	    if( (nodes(i)<=xx) && (nodes(i+1)>xx) ){
	      ind=i;
	      stop=true;
	    }
	  }
	}
      }
      
       float a=coefs(ind,1);
      float b=coefs(ind,2);
      float c=coefs(ind,3);
      float d=coefs(ind,4);
      float t=xx-nodes(ind);
      ret=a+b*t+c*t*t+d*t*t*t;
      
    }
    return ret;  
}

  float Cspline::interpolate(float xx, int ind) const{
    float ret;
    if(!fitted){
      cerr<<"Cspline::interpolate - Cspline has not been fitted"<<endl;
      exit(-1);
    }
    else{
      if(ind>n-1){
	cerr<<"Cspline::interpolate - segment index is greater than number of segments - exiting"<<endl;
	exit(-1);
      }
      else if(ind<1){
      	cerr<<"Cspline::interpolate - segment index is less than 1 - exiting"<<endl;
	exit(-1);
      }
      float a=coefs(ind,1);
      float b=coefs(ind,2);
      float c=coefs(ind,3);
      float d=coefs(ind,4);
      float t=xx-nodes(ind);
      ret=a+b*t+c*t*t+d*t*t*t;
    }
    return ret;
  }


  ColumnVector Cspline::interpolate(const ColumnVector& x) const{
    // nodes must be monotonically increasing. I don't check this. 
    // On your head be it if you don't.
    
    if(nodes.Nrows()!=vals.Nrows()){
      cerr<<"Cspline::interpolate -  Nodes and Vals should be the same length"<<endl;
      exit(-1);
    }
    
    ColumnVector ret(x.Nrows());  

    if(!fitted){
      cerr<<"Cspline::interpolate - Cspline has not been fitted"<<endl;
      exit(-1);
    }
    else{

      for(int xnum=1;xnum<=x.Nrows();xnum++){
	
	float xx=x(xnum);
	
	bool stop=false;
	int ind=0;
	
	if(xx<nodes(1)){
	  ind=1;
	}
	else if(xx>=nodes(nodes.Nrows())){
	  ind=nodes.Nrows()-1;
	}
	else{
	  for(int i=1;i<nodes.Nrows();i++){
	    if(!stop){
	      if( (nodes(i)<=xx) && (nodes(i+1)>xx) ){
		ind=i;
		stop=true;
	      }
	    }
	  }
	}

	float a=coefs(ind,1);
	float b=coefs(ind,2);
	float c=coefs(ind,3);
	float d=coefs(ind,4);
	float t=xx-nodes(ind);
	ret(xnum)=a+b*t+c*t*t+d*t*t*t;
      }
    }
    return ret;   
  }
  
 
  ColumnVector Cspline::interpolate(const ColumnVector& x,const ColumnVector& indvec) const{
    // nodes must be monotonically increasing. I don't check this. 
    // On your head be it if you don't.
    
    if(nodes.Nrows()!=vals.Nrows()){
      cerr<<"Cspline::interpolate - Nodes and Vals should be the same length"<<endl;
      exit(-1);
    }
    
    ColumnVector ret(x.Nrows());  

    if(!fitted){
      cerr<<"Cspline::interpolate - Cspline has not been fitted"<<endl;
      exit(-1);
    }
    else{
      for(int xnum=1;xnum<=x.Nrows();xnum++){
	
	float xx=x(xnum);
	
	int ind=int(indvec(xnum));

	float a=coefs(ind,1);
	float b=coefs(ind,2);
	float c=coefs(ind,3);
	float d=coefs(ind,4);
	float t=xx-nodes(ind);
	ret(xnum)=a+b*t+c*t*t+d*t*t*t;
      }
    }
    return ret;   
  }
  



}





























