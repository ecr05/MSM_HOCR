//*****************************************************************
// Iterative template routine -- QMR
//
// QMR.h solves the unsymmetric linear system Ax = b using the
// Quasi-Minimal Residual method following the algorithm as described
// on p. 24 in the SIAM Templates book.
//
//   -------------------------------------------------------------
//   return value     indicates
//   ------------     ---------------------
//        0           convergence within max_iter iterations
//        1           no convergence after max_iter iterations
//                    breakdown in:
//        2             rho
//        3             beta
//        4             gamma
//        5             delta
//        6             ep
//        7             xi
//   -------------------------------------------------------------
//   
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax=b
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

#ifndef qmr_h
#define qmr_h

#include <math.h>

namespace MISCMATHS {

template < class Matrix, class Vector, class Preconditioner1,
           class Preconditioner2, class Real >
int 
QMR(const Matrix &A, Vector &x, const Vector &b, const Preconditioner1 &M1, 
    const Preconditioner2 &M2, int &max_iter, Real &tol)
{
  Real resid;

  Vector rho(1), rho_1(1), xi(1), gamma(1), gamma_1(1), theta(1), theta_1(1);
  Vector eta(1), delta(1), ep(1), beta(1);

  Vector r, v_tld, y, w_tld, z;
  Vector v, w, y_tld, z_tld;
  Vector p, q, p_tld, d, s;

  Real normb = b.NormFrobenius();

  r = b - A * x;

  if (normb == 0.0)
    normb = 1;

  if ((resid = r.NormFrobenius() / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  v_tld = r;
  y = M1.solve(v_tld);
  rho(1) = y.NormFrobenius();

  w_tld = r;
  z = M2.trans_solve(w_tld);
  xi(1) = z.NormFrobenius();

  gamma(1) = 1.0;
  eta(1) = -1.0;
  theta(1) = 0.0;

  for (int i = 1; i <= max_iter; i++) {

    if (rho(1) == 0.0)
      return 2;                        // return on breakdown

    if (xi(1) == 0.0)
      return 7;                        // return on breakdown

    v = (1. / rho(1)) * v_tld;
    y = (1. / rho(1)) * y;

    w = (1. / xi(1)) * w_tld;
    z = (1. / xi(1)) * z;

    delta(1) = DotProduct(z, y);
    if (delta(1) == 0.0)
      return 5;                        // return on breakdown

    y_tld = M2.solve(y);               // apply preconditioners
    z_tld = M1.trans_solve(z);

    if (i > 1) {
      p = y_tld - (xi(1) * delta(1) / ep(1)) * p;
      q = z_tld - (rho(1) * delta(1) / ep(1)) * q;
    } else {
      p = y_tld;
      q = z_tld;
    }

    p_tld = A * p;
    ep(1) = DotProduct(q, p_tld);
    if (ep(1) == 0.0)
      return 6;                        // return on breakdown

    beta(1) = ep(1) / delta(1);
    if (beta(1) == 0.0)
      return 3;                        // return on breakdown

    v_tld = p_tld - beta(1) * v;
    y = M1.solve(v_tld);

    rho_1(1) = rho(1);
    rho(1) = y.NormFrobenius();
    w_tld = A.trans_mult(q) - beta(1) * w;
    z = M2.trans_solve(w_tld);

    xi(1) = z.NormFrobenius();

    gamma_1(1) = gamma(1);
    theta_1(1) = theta(1);

    theta(1) = rho(1) / (gamma_1(1) * beta(1));
    gamma(1) = 1.0 / sqrt(1.0 + theta(1) * theta(1));

    if (gamma(1) == 0.0)
      return 4;                        // return on breakdown

    eta(1) = -eta(1) * rho_1(1) * gamma(1) * gamma(1) / 
      (beta(1) * gamma_1(1) * gamma_1(1));

    if (i > 1) {
      d = eta(1) * p + (theta_1(1) * theta_1(1) * gamma(1) * gamma(1)) * d;
      s = eta(1) * p_tld + (theta_1(1) * theta_1(1) * gamma(1) * gamma(1)) * s;
    } else {
      d = eta(1) * p;
      s = eta(1) * p_tld;
    }

    x += d;                            // update approximation vector
    r -= s;                            // compute residual

    if ((resid = r.NormFrobenius() / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;                            // no convergence
}

} // End namespace MISCMATHS

#endif // End #ifndef qmr_h
