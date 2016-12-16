//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the 
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

#ifndef cg_h
#define cg_h

namespace MISCMATHS {

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
CG(const Matrix &A, Vector &x, const Vector &b,
   const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector p, z, q;
  Vector alpha(1), beta(1), rho(1), rho_1(1);

  Real normb = b.NormFrobenius();
  Vector r = b - A*x;

  if (normb == 0.0) 
    normb = 1;
  
  if ((resid = r.NormFrobenius() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    z = M.solve(r);
    rho(1) = DotProduct(r, z);
    
    if (i == 1)
      p = z;
    else {
      beta(1) = rho(1) / rho_1(1);
      p = z + beta(1) * p;
    }
    
    q = A*p;
    alpha(1) = rho(1) / DotProduct(p, q);
    
    x += alpha(1) * p;
    r -= alpha(1) * q;
    
    if ((resid = r.NormFrobenius() / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;     
    }

    rho_1(1) = rho(1);
  }
  
  tol = resid;
  return 1;
}

} // End namespace MISCMATHS

#endif // End #ifndef cg_h
