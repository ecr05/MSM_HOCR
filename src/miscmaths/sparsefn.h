/*  sparsefn.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

// Miscellaneous maths functions


#if !defined(__sparsefn_h)
#define __sparsefn_h

#define WANT_STREAM
#define WANT_MATH

#include "sparse_matrix.h"
#include "newmat.h"

using namespace NEWMAT;

namespace MISCMATHS {

  float quadratic(const ColumnVector& m, const SparseMatrix& C);
  void speye(int n, SparseMatrix& ret);
  void chol(const SparseMatrix& A, SparseMatrix& U, SparseMatrix& L);
  void inv(const SparseMatrix& U, const SparseMatrix& L, SparseMatrix& ret);
  void solvefortracex(const SparseMatrix& U, const SparseMatrix& L, const SparseMatrix& b1, const SparseMatrix& b2, float& tr1, float& tr2);
  void solveforx(const SparseMatrix& U, const SparseMatrix& L, const ColumnVector& b, ColumnVector& x);
  void solveforx(const SparseMatrix& A, const ColumnVector& b, ColumnVector& x, float tol = 0.001, int kmax = 500);
  void solveforx(const SparseMatrix& A, const ColumnVector& b, SparseMatrix& x);
  void solveforx(const SparseMatrix& A, const SparseMatrix& b, SparseMatrix& x);
  float solvefortracex(const SparseMatrix& A, const SparseMatrix& b, SparseMatrix& x, int nsamps = 50, float tol = 0.001);
  void solve(const SparseMatrix& A, const Matrix& b, SparseMatrix& x);
  void addto(SparseMatrix& A, const SparseMatrix& B, float S);
  void symmetric_addto(SparseMatrix& A, const SparseMatrix& B, float S);
  void addto(const SparseMatrix::Row& A, const SparseMatrix::Row& B, float S);
  void addto(SparseMatrix& A, const Matrix& B);
  void cov(const ColumnVector& A, SparseMatrix& ret);
}

#endif
