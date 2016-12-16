/*  sparse_matrix.h

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
#include "newmatio.h"
#include "newmat.h"
#include "miscmaths.h"
#include "utils/tracer_plus.h"

using namespace std;
using namespace Utilities;
using namespace NEWMAT;
using namespace MISCMATHS;

namespace MISCMATHS {

  SparseMatrix::SparseMatrix(int pnrows, int pncols) : 
	nrows(pnrows),
	ncols(pncols),
	data(pnrows) 
  {
  }

  void SparseMatrix::ReSize(int pnrows, int pncols)
  {
    nrows = pnrows;
    ncols = pncols;
    
    data.clear();
    data.resize(nrows);
  }

  const SparseMatrix& SparseMatrix::operator=(const Matrix& pmatin) 
  {
    data.clear();
    data.resize(pmatin.Nrows());
    nrows = pmatin.Nrows();
    ncols = pmatin.Ncols();

    for(int r=1; r <= pmatin.Nrows(); r++)
      {	
	for(int c=1; c <= pmatin.Ncols(); c++)
	  {
	    if(pmatin(r,c)!=0)
	      insert(r,c,pmatin(r,c));
	  }
      }
    return *this;
  }     
 
  void SparseMatrix::transpose(SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::transpose");

    ret.ReSize(ncols,nrows);
    
    for(int r=1; r <= nrows; r++)
      for(map<int,double>::const_iterator it=data[r-1].begin(); it!=data[r-1].end(); it++)
	ret.insert((*it).first+1, r, (*it).second);
    
  }
  
  int SparseMatrix::maxnonzerosinrow() const
  {
    int mx = 0;
    for(int r=1; r <= nrows; r++)
      {
	int si = data[r-1].size();
	if(si > mx)
	  mx = si;
      }
    return mx;
  }
  
  void SparseMatrix::permute(const ColumnVector& p, SparseMatrix& pA)
  {
    Tracer_Plus tr("SparseMatrix::permute");

    pA.ReSize(nrows,ncols);
    
    ColumnVector ip(p.Nrows());
    for(int r=1; r <= nrows; r++)
      ip(int(p(r))) = r;
    
    for(int r=1; r <= nrows; r++)
      for(map<int,double>::const_iterator it=data[r-1].begin(); it!=data[r-1].end(); it++)
	{
	  pA.insert(int(ip(r)), int(ip((*it).first+1)), (*it).second); 
	}
  }
  
  ReturnMatrix SparseMatrix::AsMatrix() const
  {
    Matrix ret(nrows,ncols);
    ret = 0;	  
    
    for(int r=1; r <= nrows; r++)
      for(map<int,double>::const_iterator it=data[r-1].begin(); it!=data[r-1].end(); it++)
	ret(r,(*it).first+1) = (*it).second;
    
    ret.Release();
    return ret;
  }

  float SparseMatrix::trace() const
  {   
    float tr = 0.0;
    for(int k = 1; k<=Nrows(); k++)
      {
	tr += (*this)(k,k);
      }

    return tr;
  }

 void SparseMatrix::vertconcatbelowme(const SparseMatrix& B)
 {
   Tracer_Plus tr("SparseMatrix::vertconcatbelowme");

   if (Ncols() != B.Ncols()) {throw Exception("Cols don't match in SparseMatrix::vertconcatbelowme");}

   data.resize(Nrows()+B.Nrows());
   for (int i=1; i<=B.Nrows(); i++) {
     this->row(Nrows()+i) = B.row(i);
   }

   nrows += B.Nrows();
 }

 void SparseMatrix::vertconcataboveme(const SparseMatrix& A)
 {
   Tracer_Plus tr("SparseMatrix::vertconcataboveme");

   if (Ncols() != A.Ncols()) {throw Exception("Cols don't match in SparseMatrix::vertconcataboveme");}

   data.resize(Nrows()+A.Nrows());
   for (int i=Nrows(); i>=1; i--) {
     this->row(i+A.Nrows()) = this->row(i);
   }
   for (int i=1; i<=A.Nrows(); i++) {
     this->row(i) = A.row(i);
   }

   nrows += A.Nrows();
 }
  
 void SparseMatrix::horconcat2myright(const SparseMatrix& B)
 {
   Tracer_Plus tr("SparseMatrix::horconcat2myright");

   if (Nrows() != B.Nrows()) {throw Exception("Rows don't match in SparseMatrix::vertconcat2myright");}

   for (int i=1; i<=Nrows(); i++) {
     const SparseMatrix::Row& tmpRow = B.row(i);
     for (SparseMatrix::Row::const_iterator it=tmpRow.begin(); it!=tmpRow.end(); it++) {
       this->insert(i,Ncols()+int(it->first)+1,double(it->second));
     }
   }
   ncols += B.Ncols();
 }

 void SparseMatrix::horconcat2myleft(const SparseMatrix& A)
 {
   Tracer_Plus tr("SparseMatrix::horconcat2myright");

   if (Nrows() != A.Nrows()) {throw Exception("Rows don't match in SparseMatrix::vertconcat2myleft");}

   for (int i=1; i<=Nrows(); i++) {
     SparseMatrix::Row oldRow = this->row(i);
     this->row(i) = SparseMatrix::Row();  // Empty row.
     const SparseMatrix::Row& tmpRow = A.row(i);
     for (SparseMatrix::Row::const_iterator it=tmpRow.begin(); it!=tmpRow.end(); it++) {
       this->insert(i,int(it->first)+1,double(it->second));
     }
     for (SparseMatrix::Row::const_iterator it=oldRow.begin(); it!=oldRow.end(); it++) {
       this->insert(i,A.Ncols()+int(it->first)+1,double(it->second));
     }
   }     
   ncols += A.Ncols();
 }

 ReturnMatrix SparseMatrix::RowAsColumn(int r) const
  {
    Tracer_Plus tr("SparseMatrix::RowAsColumn");

    ColumnVector ret;
    ret.ReSize(ncols);
    ret = 0;
    
    const SparseMatrix::Row& rowtmp = row(r);
    for(SparseMatrix::Row::const_iterator it=rowtmp.begin();it!=rowtmp.end();it++)
      {
	int c = (*it).first+1;	     	      
	double val = (*it).second;
	ret(c) = val;
      }
    
    ret.Release();
    return ret;
  }
  
  void colvectosparserow(const ColumnVector& col, SparseMatrix::Row& row)
  {
    Tracer_Plus tr("SparseMatrix::colvectosparserow");
    for(int j = 1; j<=col.Nrows(); j++)
      {
	if(std::abs(col(j))>1e-4)
	  row[j-1] = col(j);
      }
  }

  void SparseMatrix::multiplyby(double S)
  {
    Tracer_Plus tr("SparseMatrix::multiplyby");

    for(int j = 1; j<=Nrows(); j++)
      {
	SparseMatrix::Row& row = (*this).row(j);
	for(SparseMatrix::Row::iterator it=row.begin();it!=row.end();it++)
	  {
	    (*it).second *= S;
	  }
      }
  }
  
  void multiply(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply");

    int nrows = lm.Nrows();
    int ncols = rm.Ncols();

    if(lm.Ncols() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::multiply");

    ret.ReSize(nrows,ncols);

    for(int j = 1; j<=nrows; j++)
      {
	const SparseMatrix::Row& row = lm.row(j);	
	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;
	    double val = (*it).second;
	    for(int k = 1; k<=ncols; k++)
	      {
		ret.addto(j,k,val*rm(c,k));
	      }
	  }
      }

  }
  
  void multiply(const SparseMatrix& lm, const ColumnVector& rm, ColumnVector& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply2");

    int nrows = lm.Nrows();   
    
    if(lm.Ncols() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::multiply");
    
    ret.ReSize(nrows);
    
    for(int j = 1; j<=nrows; j++)
      {
	float sum = 0.0;
	const SparseMatrix::Row& row = lm.row(j);	
	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;
	    double val = (*it).second;
	    sum += val*rm(c);
	  }

	ret(j) = sum;
      }
  }

  void multiply(const DiagonalMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply");

    int nrows = lm.Nrows();
    int ncols = rm.Ncols();

    if(lm.Ncols() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::multiply");

    ret.ReSize(nrows,ncols);

    for(int j = 1; j<=nrows; j++)
      {
	const SparseMatrix::Row& row = rm.row(j);	
	for(SparseMatrix::Row::const_iterator it=row.begin();it!=row.end();it++)
	  {
	    int c = (*it).first+1;
	    double val = (*it).second;
	    ret.insert(j,c,val*lm(j,j));
	  }
      }

  }

  void multiply(const SparseMatrix& lm, const SparseMatrix::Row& rm, ColumnVector& ret)
  {
    Tracer_Plus tr("SparseMatrix::multiply3");

    int nrows = lm.Nrows();   
       
    ret.ReSize(nrows);

    for(int j = 1; j<=nrows; j++)
      {
	float sum = 0.0;
	const SparseMatrix::Row& row = lm.row(j);

	SparseMatrix::Row::const_iterator it=row.begin();
	SparseMatrix::Row::const_iterator itrm=rm.begin();

	while(it!=row.end() && itrm!=rm.end())
	  {
	    int crm = (*itrm).first;
	    int c = (*it).first;
	    if(c==crm)
	      {
		sum += ((*itrm).second)*((*it).second);
		it++;
		itrm++;
	      }
	    else if(c < crm)
	      {
		it++;
	      }
	    else
	      {
		itrm++;
	      }
	  }

	ret(j) = sum;
      }
  }

  void add(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret)
  {
    Tracer_Plus tr("SparseMatrix::add");

    int nrows = lm.Nrows();
    int ncols = lm.Ncols();

    if(lm.Ncols() != rm.Ncols() || lm.Nrows() != rm.Nrows()) throw Exception("Rows and cols don't match in SparseMatrix::add"); 

    ret.ReSize(nrows,ncols);

    for(int j = 1; j<=nrows; j++)
      {
	const SparseMatrix::Row& lmrow = lm.row(j);
	const SparseMatrix::Row& rmrow = rm.row(j);

	SparseMatrix::Row::const_iterator lmit = lmrow.begin();
	SparseMatrix::Row::const_iterator rmit = rmrow.begin();
	int lmc = (*lmit).first+1;
	int rmc = (*rmit).first+1;

	while(lmit!=lmrow.end() || rmit!=rmrow.end())
	  {
	    if((lmc<rmc && lmit!=lmrow.end()) || rmit==rmrow.end())
	      {		
		ret.insert(j,lmc,(*lmit).second+rm(j,lmc));
		lmit++;
		lmc = (*lmit).first+1;
	      }
	    else if((rmc<lmc && rmit!=rmrow.end()) || lmit==lmrow.end())
	      {
		ret.insert(j,rmc,(*rmit).second+lm(j,rmc));
		rmit++;
		rmc = (*rmit).first+1;
	      }
	    else
	      {
		//lmc==rmc
		ret.insert(j,rmc,(*lmit).second+(*rmit).second);
		lmit++;
		lmc = (*lmit).first+1;
		rmit++;
		rmc = (*rmit).first+1;
	      }
	  }
      }
  }

  // Concatenation. Note that these work also for concatenating sparse and full (newmat) 
  // matrices by virtue of the Matrix->SparseMatrix constructor/converter.

  // ret = [A; B]; % Matlab lingo
  void vertconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret)
  {
    if (A.Ncols() != B.Ncols()) {throw Exception("Cols don't match in SparseMatrix::vertconcat");}     

    ret.ReSize(A.Nrows()+B.Nrows(),A.Ncols());

    for (int i=1; i<=A.Nrows(); i++) {ret.row(i) = A.row(i);}
    for (int i=1; i<=B.Nrows(); i++) {ret.row(i+A.Nrows()) = B.row(i);}
  }

  // ret = [A B]; % Matlab lingo
  void horconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret)
  {
    if (A.Nrows() != B.Nrows()) {throw Exception("Rows don't match in SparseMatrix::horconcat");}     

    ret.ReSize(A.Nrows(),A.Ncols()+B.Ncols());

    for (int i=1; i<=A.Nrows(); i++) {
      ret.row(i) = A.row(i);
      const SparseMatrix::Row& tmpRow = B.row(i);
      for (SparseMatrix::Row::const_iterator it=tmpRow.begin(); it!=tmpRow.end(); it++) {
        ret.insert(i,A.Ncols()+int(it->first)+1,double(it->second));
      }
    }         
  }

}










