//*****************************************************************
// Iterative template routine -- BiCG
//
// BiCG solves the unsymmetric linear system Ax = b 
// using the Preconditioned BiConjugate Gradient method
//
// BiCG follows the algorithm described on p. 22 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************
//
// Slightly modified version of IML++ template. See ReadMe file.
//
// Jesper Andersson
//
/*    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#ifndef bicg_h
#define bicg_h

namespace MISCMATHS {


template < class Matrix, class Vector, class Preconditioner, class Real >
int 
BiCG(const Matrix &A, Vector &x, const Vector &b,
     const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1);
  Vector z, ztilde, p, ptilde, q, qtilde;

  Real normb = b.NormFrobenius();
  Vector r = b - A * x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = r.NormFrobenius() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    z = M.solve(r);
    ztilde = M.trans_solve(rtilde);
    rho_1(1) = DotProduct(z, rtilde);
    if (rho_1(1) == 0) { 
      tol = r.NormFrobenius() / normb;
      max_iter = i;
      return 2;
    }
    if (i == 1) {
      p = z;
      ptilde = ztilde;
    } else {
      beta(1) = rho_1(1) / rho_2(1);
      p = z + beta(1) * p;
      ptilde = ztilde + beta(1) * ptilde;
    }
    q = A * p;
    qtilde = A.trans_mult(ptilde);
    alpha(1) = rho_1(1) / DotProduct(ptilde, q);
    x += alpha(1) * p;
    r -= alpha(1) * q;
    rtilde -= alpha(1) * qtilde;

    rho_2(1) = rho_1(1);
    if ((resid = r.NormFrobenius() / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;
}
  
} // End namespace MISCMATHS

#endif // End #ifndef cg_h
