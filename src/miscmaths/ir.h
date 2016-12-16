//*****************************************************************
// Iterative template routine -- Preconditioned Richardson
//
// IR solves the unsymmetric linear system Ax = b using 
// Iterative Refinement (preconditioned Richardson iteration).
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

#ifndef ir_h
#define ir_h

namespace MISCMATHS {

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
IR(const Matrix &A, Vector &x, const Vector &b,
   const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector z;

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
    x += z;
    r = b - A * x;
    
    if ((resid = r.NormFrobenius() / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;
}

} // End namespace MISCMATHS

#endif // End #ifndef ir_h


