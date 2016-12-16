/*  sparsefn.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#define WANT_STREAM
#define WANT_MATH

#include "sparse_matrix.h"
#include "sparsefn.h"
#include "newmatio.h"
#include "newmat.h"
#include "miscmaths.h"
#include "utils/tracer_plus.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace Utilities;

namespace MISCMATHS {

  float quadratic(const ColumnVector& m, const SparseMatrix& C)
  {
    Tracer_Plus trace("sparsefns::quadratic");

    // computes m'*C*m
    // assumes that C is symmetric

    float sum = 0;

    for(int j = 1; j<=m.Nrows(); j++)
      {
	// do diagonal
	sum += C(j,j)*m(j)*m(j);

	// do off-diagonal
	const SparseMatrix::Row& row = C.row(j);

	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;

	    if(c>=j) break;

	    double val = (*it).second;
	    sum += 2*val*m(j)*m(c); 
	  }
      }
    return sum;
  }

  void speye(int n, SparseMatrix& ret)
  {
    ret.ReSize(n,n);
    for(int j = 1; j<=n; j++)
      {
	ret.insert(j,j,1);
      }
  }

  void addto(SparseMatrix::Row& A, const SparseMatrix::Row& B, float S)
  {
    // computes A = A+B*S
    if(S!=0)
      {
	for(SparseMatrix::Row::const_iterator it=B.begin();it!=B.end();it++)
	  {
	    int c = (*it).first;
	    double val = (*it).second;
	    A[c] += val*S;
	  }      
      }
  }

  void addto(SparseMatrix& A, const SparseMatrix& B, float S)
  {
    Tracer_Plus trace("sparsefns::addto");
    // computes A+B*S
    if(S!=0)
      {
	for(int j = 1; j<=B.Nrows(); j++)
	  {
	    const SparseMatrix::Row& row = B.row(j);
	    for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	      {
		int c = (*it).first+1;
		double val = (*it).second*S;
		A.addto(j,c,val);
	      }
	  }
      }
  }

  void symmetric_addto(SparseMatrix& A, const SparseMatrix& B, float S)
  {
    Tracer_Plus trace("sparsefns::symmetric_addto");
    // computes A+B*S
    if(S!=0)
      {
	for(int j = 1; j<=B.Nrows(); j++)
	  {
	    const SparseMatrix::Row& row = B.row(j);
	    
	    A.addto(j,j,B(j,j)*S);

	    for(SparseMatrix::Row::const_iterator it=row.lower_bound(j);it!=row.end();it++)
	      {		
		int c = (*it).first+1;
		double val = (*it).second*S;

		A.addto(j,c,val);
		A.addto(c,j,val);
	      }
	  }
      }
  }

  void addto(SparseMatrix& A, const Matrix& B)
  {
    Tracer_Plus trace("sparsefns::addto2");
    for(int r=1; r <= B.Nrows(); r++)
      for(int c=1; c <= B.Ncols(); c++)
	{
	  if(B(r,c)!=0)
	    A.addto(r,c,B(r,c));
	}
  }
  
  
  void chol(const SparseMatrix& A, SparseMatrix& U, SparseMatrix& L)
  {
    Tracer_Plus trace("sparsefns::chol");

    int length = A.Nrows();
    U.ReSize(length,length);

    for(int j = 1; j<=length; j++)
      {
      
	const SparseMatrix::Row& rowAj = A.row(j);
	SparseMatrix::Row& rowUj = U.row(j);

	for(SparseMatrix::Row::const_iterator it=rowAj.lower_bound(j-1);it!=rowAj.end();it++)
	  {
	    int c = (*it).first;
	    double val = (*it).second;
	    rowUj[c] = val;
	  }

	for(int k = 1; k<=j-1; k++)
	  {
	    SparseMatrix::Row& rowk = U.row(k);
	    double Ukj = U(k,j);
	    if(Ukj!=0)
	      for(SparseMatrix::Row::iterator it=rowk.lower_bound(j-1);it!=rowk.end();it++)
		{	  
		  int c = (*it).first+1;      
		  double val = (*it).second*Ukj;
		  U.addto(j,c,-val);
		}
	  }

	double sqrtUjj = std::sqrt(Max(U(j,j),1e-6));

	for(SparseMatrix::Row::iterator it=rowUj.lower_bound(j-1);it!=rowUj.end();it++)
	  {   	      	  
	    (*it).second /= sqrtUjj;
	  }
      }

    U.transpose(L);
  }

  void inv(const SparseMatrix& U, const SparseMatrix& L, SparseMatrix& ret)
  {
    Tracer_Plus trace("sparsefns::inv");

  // assumes A=LU is symmetric

    int length = U.Nrows();

    ret.ReSize(length,length);

    SparseMatrix b;
    speye(length,b);

    for(int bi=1;bi<=b.Ncols();bi++)
      {      
	// solve for y (L*y=b)
	ColumnVector y(length);
	y = 0;
	y(1) = b(1,bi)/L(1,1);

	bool compute = false;
	if(b(1,bi)!=0) compute = true;

	for(int r = 2; r<=length; r++)
	  {
	    if(!compute && b(r,bi)!=0) compute = true;
	  
	    if(compute)
	      {
		float sum = 0.0;
		const SparseMatrix::Row& row = L.row(r);
		for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
		  {
		    int c = (*it).first+1;
		  
		    if(c > r-1) break;
		  
		    double val = (*it).second;
		    sum += val*y(c);
		  }

		y(r) = (b(r,bi)-sum)/L(r,r);
	      }
	  }
      
	// solve for x(bi) (U*x=y)

	ret.set(length,bi,y(length)/U(length,length));
	compute = false;
	if(y(length)!=0) compute = true;

	// do not do rows which we already have from symmetry
	// therefore end at r=bi and not r=1
	for(int r = length; r>=bi; r--)
	  {
	    if(!compute && y(r)!=0) compute = true;
	  
	    if(compute)
	      {
		float sum = 0.0;
		const SparseMatrix::Row& row = U.row(r);
		for(SparseMatrix::Row::const_iterator it=row.lower_bound(r);it!=row.end();it++)
		  {
		    int c = (*it).first+1;
		  
		    double val = (*it).second;
		    sum += val*ret(c,bi);
		  }	
		ret.set(r,bi,(y(r)-sum)/U(r,r));
		ret.set(bi,r,(y(r)-sum)/U(r,r));     
	      }
	  }     
      }
  }

  void solvefortracex(const SparseMatrix& U, const SparseMatrix& L, const SparseMatrix& b1, const SparseMatrix& b2, float& tr1, float& tr2)
  {
    Tracer_Plus trace("sparsefns::solvefortracex");

    int length = U.Nrows();

    tr1 = 0.0;
    tr2 = 0.0;
    
    for(int bi=1;bi<=b1.Ncols();bi++)
      {
	// solve for y (L*y=b)
	ColumnVector y1(length);
	ColumnVector y2(length);
	y1 = 0;
	y2 = 0;
	y1(1) = b1(1,bi)/L(1,1);
	y2(1) = b2(1,bi)/L(1,1);
      
	bool compute1 = false;
	if(b1(1,bi)!=0) compute1 = true;

	bool compute2 = false;
	if(b2(1,bi)!=0) compute2 = true;

	for(int r = 2; r<=length; r++)
	  {
	    if(!compute1 && b1(r,bi)!=0) compute1 = true;
	    if(!compute2 && b2(r,bi)!=0) compute2 = true;
	  
	    if(compute1 || compute2)
	      {
		float sum1 = 0.0;
		float sum2 = 0.0;
		const SparseMatrix::Row& row = L.row(r);
		for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
		  {
		    int c = (*it).first+1;
		  
		    if(c > r-1) break;
		  
		    double val = (*it).second;
		    if(compute1) sum1 += val*y1(c);
		    if(compute2) sum2 += val*y2(c);
		  }

		if(compute1) y1(r) = (b1(r,bi)-sum1)/L(r,r);
		if(compute2) y2(r) = (b2(r,bi)-sum2)/L(r,r);
	      }
	  }
      
	// solve for x(bi) (U*x=y)
	ColumnVector x1(length);
	ColumnVector x2(length);
	x1 = 0;
	x2 = 0;

	x1(length) = y1(length)/U(length,length);
	x2(length) = y2(length)/U(length,length);

	compute1 = false;
	if(y1(length)!=0) compute1 = true;
	compute2 = false;
	if(y2(length)!=0) compute2 = true;

	for(int r = length; r>=bi; r--)
	  {
	    if(!compute1 && y1(r)!=0) compute1 = true;
	    if(!compute2 && y2(r)!=0) compute2 = true;
	  
	    if(compute1 || compute2)
	      {
		float sum1 = 0.0;
		float sum2 = 0.0;
		const SparseMatrix::Row& row = U.row(r);
		for(SparseMatrix::Row::const_iterator it=row.lower_bound(r);it!=row.end();it++)
		  {
		    int c = (*it).first+1;
		  
		    double val = (*it).second;
		    if(compute1) sum1 += val*x1(c);
		    if(compute2) sum2 += val*x2(c);
		  }
	      
		if(compute1) x1(r) = (y1(r)-sum1)/U(r,r);     
		if(compute2) x2(r) = (y2(r)-sum2)/U(r,r);
	      }
	  }
      
	tr1 += x1(bi);
	tr2 += x2(bi);
      }

  }

  float solvefortracex(const SparseMatrix& A, const SparseMatrix& b, SparseMatrix& x, int nsamps, float tol)
  {
    Tracer_Plus trace("sparsefns::solvefortracex");
  
    int every = Max(1,A.Ncols()/nsamps);
    //    int every = 1;
    //    OUT(every);

    float tr = 0.0;

    // assumes symmetric A and b
    for(int r = every; r<=A.Ncols();  r+=every)
      {
//  	cout << float(r)/A.Ncols() << "\r"; 
//  	cout.flush();	

	ColumnVector br = b.RowAsColumn(r);      
	ColumnVector xr = x.RowAsColumn(r);
      
	solveforx(A,br,xr,tol);
	
	for(int c = 1; c<=b.Ncols(); c++)
	  {
	    if(xr(c)!=0)
	      {
		x.set(r,c,xr(c));
	      }
	  }

	tr += xr(r);
      }

    cout << endl;

    tr *= every;

    return tr;
  }

  void solveforx(const SparseMatrix& A, const SparseMatrix& b, SparseMatrix& x)
  {
    Tracer_Plus trace("sparsefns::solveforx");
  
    // assumes symmetric A and b
    for(int r = 1; r<=A.Ncols();  r++)
      {
	cout << float(r)/A.Ncols() << "\r"; 
	cout.flush();	

	ColumnVector br = b.RowAsColumn(r);      
	ColumnVector xr = x.RowAsColumn(r);
      
	solveforx(A,br,xr);

	for(int c = 1; c<=b.Ncols(); c++)
	  {
	    if(xr(c)!=0)
	      {
		x.set(r,c,xr(c));
	      }
	  }      
      }
    cout << endl;
  }

  void solveforx(const SparseMatrix& A, const ColumnVector& b, ColumnVector& x, float tol, int kmax)
  {
    //
    // Algorithm based on Golub & van Loan, chapter 10, page 527.
    //
    Tracer_Plus trace("sparsefns::solveforx");

    if(norm2(b)==0)
      {
	x = 0;
      }
    else
      {
	int k = 2;
	kmax = Max(b.Nrows(),kmax);
  
	ColumnVector tmp;
	multiply(A,x,tmp);

	ColumnVector r = b-tmp;
	ColumnVector rho(kmax);
	rho = Sqr(norm2(r));
	ColumnVector w;
	ColumnVector p = r;

	while(std::sqrt(rho(k))>tol*norm2(b) && k < kmax)
	  {
	    k++;
	    //if(k>2)
	    p = r + p*rho(k-1)/rho(k-2);
	    //else
	    //	p = r;

	    // 	SparseMatrix::Row passparserow;
	    // 	colvectosparserow(p,passparserow);
	    // 	multiply(A,passparserow,w);

	    multiply(A,p,w);
   
	    float alpha = 0.0;
	    //if(k>1)
	    alpha = rho(k-1)/(p.t()*w).AsScalar();
	    //else
	    //alpha = 1;
   
	    x += alpha*p;
	    r -= alpha*w;
	    rho(k) = Sqr(norm2(r));            
	  }
    

	if(k>kmax/2.0)
	  {
  	    OUT(std::sqrt(rho(k-1)));
  	    OUT(norm2(b));
	    OUT(k);
	    cout.flush();
	  }

      }


    //    write_ascii_matrix("rho",rho);
  }

  void solveforx(const SparseMatrix& U, const SparseMatrix& L, const ColumnVector& b, ColumnVector& x)
  {
    Tracer_Plus trace("sparsefns::solveforx");

    int length = U.Nrows();

    x.ReSize(length);

    // solve for y (L*y=b)
    ColumnVector y(length);
    y = 0;
    y(1) = b(1)/L(1,1);
    bool compute = false;
    if(b(1)!=0) compute = true;

    for(int r = 2; r<=length; r++)
      {
	if(!compute && b(r)!=0) compute = true;
	  
	if(compute)
	  {
	    float sum = 0.0;
	    const SparseMatrix::Row& row = L.row(r);
	    for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	      {
		int c = (*it).first+1;
		  
		if(c > r-1) break;
		  
		double val = (*it).second;
		sum += val*y(c);

	      }

	    y(r) = (b(r)-sum)/L(r,r);
	  }
      }
      
    // solve for x (U*x=y)
    x(length) = y(length)/U(length,length);
    compute = false;
    if(y(length)!=0) compute = true;

    for(int r = length; r>=1; r--)
      {
	if(!compute && y(r)!=0) compute = true;
	  
	if(compute)
	  {
	    float sum = 0.0;
	    const SparseMatrix::Row& row = U.row(r);
	    for(SparseMatrix::Row::const_iterator it=row.lower_bound(r);it!=row.end();it++)
	      {
		int c = (*it).first+1;
		  
		double val = (*it).second;
		sum += val*x(c);
	      }
	      
	    x(r) = (y(r)-sum)/U(r,r);     
	  }
      }

  }

  void solve(const SparseMatrix& A, const Matrix& b, SparseMatrix& x)
  {
    Tracer_Plus trace("sparsefns::solve");
    
    int length = A.Nrows();
    
    SparseMatrix U;
    SparseMatrix L;
    chol(A,U,L);
    
    x.ReSize(length,b.Ncols());
  
    for(int bi=1;bi<=b.Ncols();bi++)
      {
	// solve for y (L*y=b)
	ColumnVector y(length);
	y = 0;
	y(1) = b(1,bi)/L(1,1);
	bool compute = false;
	if(b(1,bi)!=0) compute = true;

	for(int r = 2; r<=length; r++)
	  {
	    if(!compute && b(r,bi)!=0) compute = true;
	  
	    if(compute)
	      {
		float sum = 0.0;
		SparseMatrix::Row& row = L.row(r);
		for(SparseMatrix::Row::iterator it=row.begin();it!=row.end();it++)
		  {
		    int c = (*it).first+1;
		  
		    if(c > r-1) break;
		  
		    double val = (*it).second;
		    sum += val*y(c);

		  }

		y(r) = (b(r,bi)-sum)/L(r,r);
	      }
	  }
      
	// solve for x (U*x=y)
	x.set(length,bi,y(length)/U(length,length));
	compute = false;
	if(y(length)!=0) compute = true;

	for(int r = length; r>=1; r--)
	  {
	    if(!compute && y(r)!=0) compute = true;
	  
	    if(compute)
	      {
		float sum = 0.0;
		SparseMatrix::Row& row = U.row(r);
		for(SparseMatrix::Row::iterator it=row.lower_bound(r);it!=row.end();it++)
		  {
		    int c = (*it).first+1;
		  
		    double val = (*it).second;
		    sum += val*x(c,bi);
		  }
	      
		x.set(r,bi,(y(r)-sum)/U(r,r));     
	      }
	  }
      }
  }


  void cov(const ColumnVector& A, SparseMatrix& ret)
  {
    Tracer_Plus trace("sparsefns::cov");
  
    ret.ReSize(A.Nrows(),A.Nrows());

    for(int r=1; r <= A.Nrows(); r++)
      {
	// diagonal
	if(A(r) != 0)
	  {
	    ret.set(r,r,Sqr(A(r)));

	    // off-diagonal
	    for(int c=r+1; c <= A.Nrows(); c++)
	      {
		if(A(c) != 0)
		  {
		    ret.set(r,c,A(r)*A(c));
		    ret.set(c,r,A(r)*A(c));
		  }
	      }
	  }
      }
  }  
}
