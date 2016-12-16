/*  minimize
 
    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(minimize_h)
#define minimize_h

#include <string>
#include <iostream>
#include <fstream>
//#include <unistd.h>
#include <vector>
#include <algorithm>
#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"


#define WANT_STREAM
#define WANT_MATH

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace std;
///////////////////////////////////////////////////////

//fminsearch.m

namespace MISCMATHS {

class pair_comparer
{
public:
  bool operator()(const pair<float,ColumnVector>& p1,const pair<float,ColumnVector>& p2) const
  {
    return p1.first < p2.first;
  }    
};

  class EvalFunction;
  class gEvalFunction;

float diff1(const ColumnVector& x, const EvalFunction& func, int i,float h,int errorord=4);// finite diff derivative

float diff2(const ColumnVector& x, const EvalFunction& func, int i,float h,int errorord=4);// finite diff 2nd derivative

float diff2(const ColumnVector& x, const EvalFunction& func, int i,int j,float h,int errorord=4);// finite diff cross derivative

ReturnMatrix gradient(const ColumnVector& x, const EvalFunction& func,float h,int errorord=4);// finite diff derivative vector 

ReturnMatrix hessian(const ColumnVector& x, const EvalFunction& func,float h,int errorord=4);// finite diff hessian

void minsearch(ColumnVector& x, const EvalFunction& func, ColumnVector& paramstovary);

void scg(ColumnVector& x, const gEvalFunction& func, ColumnVector& paramstovary, float tol = 0.0000001, float eps=1e-16, int niters=500);

class EvalFunction
{//Function where gradient is not analytic (or you are too lazy to work it out) (required for fminsearch)
public:
  EvalFunction(){}
  virtual float evaluate(const ColumnVector& x) const = 0; //evaluate the function
  virtual ~EvalFunction(){};

  virtual void minimize(ColumnVector& x)
  {
    ColumnVector paramstovary(x.Nrows());
    paramstovary = 1;
    minsearch(x,*this,paramstovary);
  }

  virtual void minimize(ColumnVector& x, ColumnVector& paramstovary)
  {
    minsearch(x,*this,paramstovary);
  }

private:
  const EvalFunction& operator=(EvalFunction& par);
  EvalFunction(const EvalFunction&);
};

class gEvalFunction : public EvalFunction
{//Function where gradient is analytic (required for scg)
public:
  gEvalFunction() : EvalFunction(){}
  // evaluate is inherited from EvalFunction
  
  virtual ReturnMatrix g_evaluate(const ColumnVector& x) const = 0; //evaluate the gradient
  virtual ~gEvalFunction(){};

  virtual void minimize(ColumnVector& x)
  {
    ColumnVector paramstovary(x.Nrows());
    paramstovary = 1;
    scg(x,*this,paramstovary);
  }

  virtual void minimize(ColumnVector& x, ColumnVector& paramstovary)
  {
    scg(x,*this,paramstovary);
  }

private:
    
  const gEvalFunction& operator=(gEvalFunction& par);
  gEvalFunction(const gEvalFunction&);
};

}
   

#endif







