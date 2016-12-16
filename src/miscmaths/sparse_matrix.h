/*  Sparse_Matrix.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(Sparse_Matrix_h)
#define Sparse_Matrix_h

#include <iostream>
#include "newmat.h"
#include <math.h>
#include <map>
#include <vector>
#include "newmatio.h"

using namespace NEWMAT;
using namespace std;

namespace MISCMATHS {  
  
  class SparseMatrix
    {
    public:

      typedef map<int,double> Row;
 
      SparseMatrix() : nrows(0), ncols(0) {}

      SparseMatrix(int pnrows, int pncols);

      SparseMatrix(const SparseMatrix& psm) 
	{
	  operator=(psm);
	}

      const SparseMatrix& operator=(const SparseMatrix& psm) 
	{
	  nrows = psm.nrows;
	  ncols = psm.ncols;
	  data = psm.data;

	  return *this;
	}

      SparseMatrix(const Matrix& pmatin) 
	{
	  operator=(pmatin);
	}

      const SparseMatrix& operator=(const Matrix& pmatin);

      //      void ReSize(int pnrows, int pncols)
      void ReSize(int pnrows, int pncols);

      void clear()
	{
	  ReSize(0,0);
	}
      
      void transpose(SparseMatrix& ret);
      
      ReturnMatrix RowAsColumn(int r) const;      

      int maxnonzerosinrow() const;

      void permute(const ColumnVector& p, SparseMatrix& pA);

      const double operator()(int x, int y) const
	{
	  double ret = 0.0;
	  map<int,double>::const_iterator it=data[x-1].find(y-1);
	  if(it != data[x-1].end())
	    ret = (*it).second;

	  return ret;
	}
      
      void set(int x, int y, double val) 
	{
	  data[x-1][y-1] = val;
	}

      void update(int x, int y, double val) 
	{
	  data[x-1][y-1] = val;
	}

      void insert(int x, int y, double val) 
	{
	  data[x-1].insert(Row::value_type(y-1,val));
	}

      void addto(int x, int y, double val) 
	{
	  if(val!=0)
	    data[x-1][y-1] += val;
	}

      void multiplyby(int x, int y, double val) 
	{
	  if((*this)(x,y)!=0)
	    data[x-1][y-1] *= val;
	}
            
      float trace() const;      

      Row& row(int r) { return data[r-1]; }

      const Row& row(int r) const { return data[r-1]; }

      ReturnMatrix AsMatrix() const;

      int Nrows() const { return nrows; }
      int Ncols() const { return ncols; }

      void multiplyby(double S);

      void vertconcatbelowme(const SparseMatrix& B); // me -> [me; B]
      void vertconcataboveme(const SparseMatrix& A); // me -> [A; me]
      void horconcat2myright(const SparseMatrix& B); // me -> [me B]
      void horconcat2myleft(const SparseMatrix& A);  // me -> [A me]

    private:
      
      int nrows;
      int ncols;

      vector<map<int,double> > data;

    };   

  void multiply(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret);
  void multiply(const DiagonalMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret);

  void multiply(const SparseMatrix& lm, const ColumnVector& rm, ColumnVector& ret);

  void multiply(const SparseMatrix& lm, const SparseMatrix::Row& rm, ColumnVector& ret);

  void add(const SparseMatrix& lm, const SparseMatrix& rm, SparseMatrix& ret);

  void colvectosparserow(const ColumnVector& col, SparseMatrix::Row& row);

  void vertconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret);

  void horconcat(const SparseMatrix& A, const SparseMatrix& B, SparseMatrix& ret);
}

#endif
