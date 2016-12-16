//*****************************************************************
// Iterative template routine -- CGS
//
// CGS solves the unsymmetric linear system Ax = b 
// using the Conjugate Gradient Squared method
//
// CGS follows the algorithm described on p. 26 of the 
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

#ifndef cgs_h
#define cgs_h

namespace MISCMATHS {

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
CGS(const Matrix &A, Vector &x, const Vector &b,
    const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1);
  Vector p, phat, q, qhat, vhat, u, uhat;

  Real normb = b.NormFrobenius();
  Vector r = b - A*x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = r.NormFrobenius() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    rho_1(1) = DotProduct(rtilde, r);
    if (rho_1(1) == 0) {
      tol = r.NormFrobenius() / normb;
      return 2;
    }
    if (i == 1) {
      u = r;
      p = u;
    } else {
      beta(1) = rho_1(1) / rho_2(1);
      u = r + beta(1) * q;
      p = u + beta(1) * (q + beta(1) * p);
    }
    phat = M.solve(p);
    vhat = A*phat;
    alpha(1) = rho_1(1) / DotProduct(rtilde, vhat);
    q = u - alpha(1) * vhat;
    uhat = M.solve(u + q);
    x += alpha(1) * uhat;
    qhat = A * uhat;
    r -= alpha(1) * qhat;
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

#endif // End #ifndef cgs_h
