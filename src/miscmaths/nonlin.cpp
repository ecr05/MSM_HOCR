/*    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
// Definitions for module nonlin

#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "newmat.h"
#include "newmatio.h"
#include "bfmatrix.h"
#include "nonlin.h"
#include "Simplex.h"
#include "utils/fsl_isfinite.h"

using namespace std;
using namespace NEWMAT;

namespace MISCMATHS {

// Declarations of routines for use only in this module

// Main routine for Variable-Metric optimisation
NonlinOut varmet(const NonlinParam& p, const NonlinCF& cfo);

// Main routine for Gradient-descent optimisation
NonlinOut grades(const NonlinParam& p, const NonlinCF& cfo);

// Main routine for Conjugate-Gradient optimisation
NonlinOut congra(const NonlinParam& p, const NonlinCF& cfo);

// Main routine for scaled conjugate-gradient optimisation
NonlinOut sccngr(const NonlinParam& p, const NonlinCF& cfo);

// Main routine for Levenberg-Marquardt optimisation
NonlinOut levmar(const NonlinParam& p, const NonlinCF& cfo);

// Main routine for amoeba (Nelder-Mead) optimisation
NonlinOut amoeba(const NonlinParam& p, const NonlinCF& cfo);

LinOut linsrch(// Input
               const ColumnVector&  pdir,    // Search direction
               const ColumnVector&  p0,      // Current parameter values
               const ColumnVector&  grad,    // Gradient at p0
               const NonlinCF&      cfo,     // Cost-function object
               double               f0,      // Current cost-function value
               double               sf,      // Scale factor for cost-function
               double               maxiter, // Max # of iterations in line minimisation
               double               sm,      // Stepmax
               double               alpha,   // Alpha (sorry).
               double               ptol,    // Tolerance in parameter space
               // Output
               double               *lambda,// Resulting step length
               double               *of,    // Value of cost-function on output
               ColumnVector         *np);   // New parameters

double scale_factor(const ColumnVector&  p,       // Current parameter values
                    const ColumnVector&  pdir,    // Search direction
                    const NonlinCF&      cfo,     // Cost-function object
                    int                  maxiter, // Max # of iterations
                    double               sf);     // Scale factor.

LinOut linmin(// Input
              const ColumnVector&   p,      // Current parameter values
              const ColumnVector&   pdir,   // Search direction
              const NonlinCF&       cfo,    // Cost-function object
              double                sf,     // Scale factor for cost-function
              pair<double,double>   lp,     // Left point
              pair<double,double>   mp,     // Point somewhere in interval
              pair<double,double>   rp,     // Right point
              double                ftol,   // Fractional tolerance
              int                   maxiter,// Max # of iterations
              // Output
              pair<double,double>   *x);    // Best point

pair<double,double> bracket(// Input 
                            const ColumnVector& p,      // Current parameter values
                            const ColumnVector& pdir,   // Search direction
                            const NonlinCF&     cfo,    // Cost-function object
                            double              ptol,   // Relative tolerance for parameter values
                            double              sf,     // Scale factor of cost-function
                            // Output
                            pair<double,double> *p_0,   // Cost-function value at p
                            pair<double,double> *p_m);  // Point between p_0 and p_l

// Utility routine that checks for convergence based on "zero"-gradient

// Utility routines that checks for convergence based on various criteria
// Based on zero (neglible) gradient
bool zero_grad_conv(const ColumnVector&   par,
                    const ColumnVector&   grad,
                    double                cf,
                    double                gtol);
// Based on zero (neglible) decrease in cost-function
bool zero_cf_diff_conv(double cfo,
                       double cfn,
                       double cftol);
// Based on zero (neglible) step in parameter space
bool zero_par_step_conv(const ColumnVector&   par,
                        const ColumnVector&   step,
                        double                ptol);

void print_newmat(const NEWMAT::GeneralMatrix&  m,
                  std::string                   fname);


std::string NonlinParam::TextStatus() const
{
  switch (status) {
  case NL_UNDEFINED:
    return(std::string("Status is undefined. Object has been created but no minimisation has been performed"));
    break;
  case NL_MAXITER:
    return(std::string("The optimisation did not converge because the maximum number of iterations was exceeded"));
    break;
  case NL_LM_MAXITER:
    return(std::string("The optimisation did not converge because the maximum number of iterations for a single line minimisation was exceeded"));
    break;
  case NL_PARCONV:
    return(std::string("The optimisation converged. The convergence criterion was that the last step in parameter space was very short"));
    break;
  case NL_GRADCONV:
    return(std::string("The optimisation converged. The convergence criterion was that all the elements of the gradient were very small"));
    break;
  case NL_CFCONV:
    return(std::string("The optimisation converged. The convergence criterion was that the last step changed the cost-function by an insignificant amount"));
    break;
  case NL_LCONV:
    return(std::string("The optimisation converged. The convergence criterion was that lambda became too large"));
    break;
  default:
    return(std::string("Impossible status. This indicates there is a bug"));
    break;
  }
}

// If user choses not to overide the grad-method of the NonlinCF
// base class this routine will be used to calculate numerical derivatives.

ReturnMatrix NonlinCF::grad(const ColumnVector& p) const
{
  ColumnVector gradv(p.Nrows());
  ColumnVector tmpp = p;
  double tiny = 1e-8;
  double cf0 = cf(tmpp);
  for (int i=0; i<p.Nrows(); i++) {
    double step = tiny * std::max(tmpp.element(i),1.0);
    tmpp.element(i) += step;
    gradv.element(i) = (cf(tmpp) - cf0) / step;
    tmpp.element(i) -= step;
  }
  gradv.Release();
  return(gradv);    
}

// If user choses not to overide the hess-method of the NonlinCF
// base class this routine will calculate numerical 2nd derivatives.
// Note also that the hessian will only be used for the Levenberg-
// Marquardt minimisation.
// It is in general _not_ a good idea to use a method that explicitly
// uses the Hessian (i.e. LM) when calculating it numerically.
// Note that it returns a (safe) pointer to BFMatrix. BFMatrix has two
// subclasses FullBFMatrix, and SparseBFMatrix which can be used 
// used interchangeably depending on if the structure of the Hessian
// renders it sparse or not.

boost::shared_ptr<BFMatrix> NonlinCF::hess(const ColumnVector&         p,
                                           boost::shared_ptr<BFMatrix> iptr) const
{
  boost::shared_ptr<BFMatrix>    hessm;
  if (iptr && int(iptr->Nrows())==p.Nrows() && int(iptr->Ncols())==p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));
  ColumnVector tmpp = p;
  double tiny = 1e-4;
  double fx0y0 = cf(tmpp);
  ColumnVector fdx(p.Nrows());
  ColumnVector step(p.Nrows());

  // First calculate all f(x+dx_i) values

  for (int i=0; i<p.Nrows(); i++) {         
    step.element(i) = tiny * std::max(tmpp.element(i),1.0);
    tmpp.element(i) += step.element(i);
    fdx.element(i) = cf(tmpp);
    tmpp.element(i) -= step.element(i);
  }
  
  // Then values of matrix

  for (int i=0; i<p.Nrows(); i++) {
    for (int j=i; j<p.Nrows(); j++) {
      if (i==j) {   // If diagonal element
        tmpp.element(i) -= step.element(i);
        double tmp = cf(tmpp);
        hessm->Set(i+1,i+1,(fdx.element(i) + tmp - 2.0*fx0y0) / (step.element(i)*step.element(i)));
        tmpp.element(i) += step.element(i);
      }
      else {        // If off-diagonal element
        tmpp.element(i) += step.element(i);
        tmpp.element(j) += step.element(j);
        hessm->Set(i+1,j+1,(cf(tmpp)+fx0y0-fdx.element(i)-fdx.element(j)) / (step.element(i)*step.element(j)));
        hessm->Set(j+1,i+1,hessm->Peek(i+1,j+1));
        tmpp.element(i) -= step.element(i);
        tmpp.element(j) -= step.element(j);
      }
    }
  }

  return(hessm);
}


// Display (for debug purposes) matrix if it is small enough for that to make sense

void VarmetMatrix::print() const
{
  if (sz > 10) {
    cout << "Matrix too big to be meaningful to display" << endl;
    return;
  }
  else {
    if (mtp == VM_FULL) {
      cout << setw(10) << setprecision(5) << mat;
    }
    else if (mtp == VM_COL) {
      Matrix  tmp = IdentityMatrix(sz);
      for (unsigned int i=0; i<sf.size(); i++) {
        tmp += sf[i] * vec[i]*(vec[i]).t();
      }
      cout << setw(10) << setprecision(5) << tmp;
    }
  }
  return;
}

// Update estimate of inverse Hessian based on latest step

void VarmetMatrix::update(const NEWMAT::ColumnVector& pdiff,  // x_{i+1} - x_i
                          const NEWMAT::ColumnVector& gdiff)  // \nabla f_{i+1} - \nabla f_i
{
  // Sort out if this call defines size of problem
  if (pdiff.Nrows() != sz || gdiff.Nrows() != sz) {
    if (sf.size() == 0 && pdiff.Nrows()==gdiff.Nrows()) {
      sz = pdiff.Nrows();
      if (mtp == VM_OPT) {
        if (sz < 100) {mtp = VM_FULL;}
        else {mtp = VM_COL;}
      }
    }
    else {throw NonlinException("VarmetMatrix::update: mismatch between vector and matrix sizes");}
  }
  // Now do the actual update
  double sf1 = DotProduct(pdiff,gdiff);
  if ((sf1*sf1) > MISCMATHS::EPS*DotProduct(pdiff,pdiff)*DotProduct(gdiff,gdiff)) {
    sf1 = 1.0 / sf1;
    ColumnVector v2 = (*this) * gdiff;
    double sf2 = -1.0 / DotProduct(gdiff,v2);
    if (mtp == VM_FULL) {
      mat += sf1 * pdiff * pdiff.t();
      mat += sf2 * v2 * v2.t();
    }
    else {
      vec.push_back(pdiff);
      vec.push_back(v2);
      sf.push_back(sf1);
      sf.push_back(sf2);
    }
    if (utp == VM_BFGS) {
      if (mtp == VM_FULL) {
        ColumnVector u = sf1*pdiff + sf2*v2;
        mat -= (1.0/sf2) * u * u.t();
      }
      else {
        vec.push_back(sf1*pdiff + sf2*v2);
        sf.push_back(-1.0/sf2);
      }
    }
  }
}

// Multiply representation of matrix with vector

ColumnVector operator*(const VarmetMatrix& m, const ColumnVector& v)
{
  if (m.mtp == VM_FULL) {return(m.mat*v);}
  else {
    ColumnVector ov = v; // Multiplication with unity matrix
    if (m.sf.size() != 0) {
      std::vector<double>::const_iterator                sfp;
      std::vector<NEWMAT::ColumnVector>::const_iterator  vep;
      for (sfp=m.sf.begin(), vep=m.vec.begin(); sfp!=m.sf.end(); ++sfp, ++vep) {
        double tmp = (*sfp) * DotProduct((*vep),v);
        ov += tmp * (*vep);
      }
    }
    return(ov);
  }
}


// Gateway function to routines for non-linear optimisation

NonlinOut nonlin(const NonlinParam& p, const NonlinCF& cfo)
{
  NonlinOut status = NL_MAXITER;

  // Call functions that actually do the job

  switch (p.Method()) {
  case NL_VM:
    status = varmet(p,cfo);
    break;
  case NL_CG:
    status = congra(p,cfo);
    break;
  case NL_SCG:
    status = sccngr(p,cfo);
    break;
  case NL_LM:
    status = levmar(p,cfo);
    break;
  case NL_GD:
    status = grades(p,cfo);
    break;
  case NL_NM:
    status = amoeba(p,cfo);
    break;
  }

  return(status);
}

// Main routine for Levenberg-Marquardt optimisation

NonlinOut levmar(const NonlinParam& p, const NonlinCF& cfo)
{
  // Calculate initial values
  p.SetCF(cfo.cf(p.Par()));                                   // Cost-function evaluated at current parameters
  bool                            success = true;             // True if last step decreased CF
  double                          olambda = 0.0;              // How much the diagonal of H was nudged last time
  ColumnVector                    g;                          // Gradient
  boost::shared_ptr<BFMatrix>     H;                          // Hessian

  while (p.NextIter(success)) {
    if (success) {                                            // If last attempt decreased cost-function
      g = cfo.grad(p.Par());                                  // Gradient evaluated at current parameters
      H = cfo.hess(p.Par(),H);                                // Hessian evaluated at current parameters
    }
    for (int i=1; i<=p.NPar(); i++) {                         // Nudge it
      if (p.GaussNewtonType() == LM_LM) {                     // If Levenberg-Marquardt
        // H->AddTo(i,i,(p.Lambda()-olambda)*H->Peek(i,i));
        H->Set(i,i,((1.0+p.Lambda())/(1.0+olambda))*H->Peek(i,i));
      }
      else if (p.GaussNewtonType() == LM_L) {                // If Levenberg
        H->AddTo(i,i,p.Lambda()-olambda);             
      }
    }
    ColumnVector step;
    double ncf = 0.0;
    bool inv_fail = false;  // Signals failure of equation solving
    try {
      step = -H->SolveForx(g,SYM_POSDEF,p.EquationSolverTol(),p.EquationSolverMaxIter());
      ncf = cfo.cf(p.Par()+step);
    }
    catch(...) {
      inv_fail = true;
    }
    if (!inv_fail && (success = (ncf < p.CF()))) {              // If last step successful
      olambda = 0.0;                                          // Pristine Hessian, so no need to undo old lambda
      p.SetPar(p.Par()+step);                                 // Set attempt as new parameters
      p.SetLambda(p.Lambda()/10.0);                           // Decrease nudge factor
      // Check for convergence based on small decrease of cf
      if (zero_cf_diff_conv(p.CF(),ncf,p.FractionalCFTolerance())) {
        p.SetCF(ncf); p.SetStatus(NL_CFCONV); return(p.Status()); 
      }
      p.SetCF(ncf);                                           // Store value of cost-function
    }
    else {                                                    // If last step was unsuccesful
      olambda = p.Lambda();                                   // Returning to same H, so must undo old lambda
      p.SetLambda(10.0*p.Lambda());                           // Increase nudge factor
      p.SetCF(p.CF());                                        // Push another copy of best cost function value thus far
      // Check for convergence based on _really_ large lambda
      if (p.Lambda() > p.LambdaConvergenceCriterion()) {
        p.SetStatus(NL_LCONV); return(p.Status());
      }
    }
  }
  // Getting here means we did too many iterations
  p.SetStatus(NL_MAXITER);

  return(p.Status());
}

// Main routine for Nelder-Mead simplex optimisation

NonlinOut amoeba(const NonlinParam& p, 
		 const NonlinCF& cfo)
{
  // cout << "Initialsing simplex" << endl; cout.flush();
  Simplex smplx(p.Par(),cfo,p.GetAmoebaStart());
  p.SetCF(smplx.BestFuncVal());

  while (p.NextIter()) {
    // Check for convergence based on fractional difference 
    // between best and worst points in simplex.
    // cout << "New iteration: Checking for convergence" << endl; cout.flush();
    if (zero_cf_diff_conv(smplx.WorstFuncVal(),smplx.BestFuncVal(),p.FractionalCFTolerance())) {
      p.SetStatus(NL_CFCONV);
      return(p.Status());
    }
    
    // cout << "Attempting reflexion" << endl; cout.flush();
    double newf = smplx.Reflect();   // Attempt reflexion
    // Extend into an expansion if reflexion very successful
    if (newf <= smplx.BestFuncVal()) { 
      // cout << "Reflexion succesful: attempting expansion" << endl; cout.flush();
      smplx.Expand(); // Attempt expansion
    }
    else if (newf >= smplx.SecondWorstFuncVal()) {
      // cout << "New value worse than second worst: attempting contraction" << endl; cout.flush();
      double worst_fval = smplx.WorstFuncVal();
      newf = smplx.Contract();     // Do a contraction towards plane of "better" points
      if (newf >= worst_fval) {    // Didn't work. Contract towards best point
	// cout << "Contraction unsuccesful: contracting towards best point" << endl; cout.flush();
	smplx.MultiContract();
      }
    } 
    smplx.UpdateRankIndicies();
    p.SetCF(smplx.BestFuncVal());
    p.SetPar(smplx.BestPar());
  }
  // If we're here we've exceeded max number of iterations
  p.SetStatus(NL_MAXITER);
  return(p.Status());
}

// Main routine for gradient-descent optimisation. It is
// included mainly as a debugging tool for when the more
// advanced methods fail and one wants to pinpoint the
// reasons for that. 

NonlinOut grades(const NonlinParam& np, const NonlinCF& cfo)
{
  // Set up initial values
  np.SetCF(cfo.cf(np.Par()));
  ColumnVector g = -cfo.grad(np.Par());

  while (np.NextIter()) {
    // Check for convergence based on zero gradient
    if (zero_grad_conv(np.Par(),g,np.CF(),np.FractionalGradientTolerance())) {
      np.SetStatus(NL_GRADCONV); return(np.Status());
    }  
    // Bracket minimum along g
    pair<double,double> lp, mp;                                                                      // Leftmost and middle point of bracket
    pair<double,double> rp = bracket(np.Par(),g,cfo,np.FractionalParameterTolerance(),1.0,&lp,&mp);  // Rightmost point of bracket
    if (rp == lp) {                                                                                  // If no smaller point along g
      np.SetStatus(NL_PARCONV); return(np.Status());                                                 // Assume it is because we are at minimum
    }
    // Find minimum along g between lp and rp
    pair<double,double> minp;                               // Minimum along g
    LinOut lm_status = linmin(np.Par(),g,cfo,1.0,lp,mp,rp,
                              np.LineSearchFractionalParameterTolerance(),
                              np.LineSearchMaxIterations(),&minp);
    // Check for problems with line-search
    if (lm_status == LM_MAXITER) {np.SetStatus(NL_LM_MAXITER); return(np.Status());} // Ouch!
    // Set new cf value and parameters
    np.SetPar(np.Par() + minp.first*g);
    // Check for convergence based on small decrease of cost-function
    if (zero_cf_diff_conv(np.CF(),minp.second,np.FractionalCFTolerance())) {
      np.SetCF(minp.second); 
      np.SetStatus(NL_CFCONV); 
      return(np.Status());
    }
    // Check for convergence based on neglible move in parameter space
    else if (zero_par_step_conv(minp.first*g,np.Par(),np.FractionalParameterTolerance())) {
      np.SetCF(minp.second); 
      np.SetStatus(NL_PARCONV); 
      return(np.Status());
    }
    else {  // If no covergence
      np.SetCF(minp.second);
      g = -cfo.grad(np.Par());
    }
  }
  // If we get here we have used too many iterations
  np.SetStatus(NL_MAXITER);
    
  return(np.Status());
}

// Main routine for conjugate-gradient optimisation. The 
// implementation follows that of Numerical Recipies 
// reasonably closely.

NonlinOut congra(const NonlinParam& np, const NonlinCF& cfo)
{
  // Set up initial values
  np.SetCF(cfo.cf(np.Par()));
  ColumnVector r = -cfo.grad(np.Par());
  ColumnVector p = r;

  while (np.NextIter()) {
    // Check for convergence based on zero gradient
    if (zero_grad_conv(np.Par(),r,np.CF(),np.FractionalGradientTolerance())) {
      np.SetStatus(NL_GRADCONV); return(np.Status());
    }  
    // Bracket minimum along p
    pair<double,double> lp, mp;                                                                      // Leftmost and middle point of bracket
    pair<double,double> rp = bracket(np.Par(),p,cfo,np.FractionalParameterTolerance(),1.0,&lp,&mp);  // Rightmost point of bracket
    if (rp == lp) {                                                                                  // If no smaller point along p
      np.SetStatus(NL_PARCONV); return(np.Status());                                                 // Assume it is because we are at minimum
    }
    // Find minimum along p between lp and rp
    pair<double,double> minp;                               // Minimum along p
    LinOut lm_status = linmin(np.Par(),p,cfo,1.0,lp,mp,rp,
                              np.LineSearchFractionalParameterTolerance(),
                              np.LineSearchMaxIterations(),&minp);
    // Check for problems with line-search
    if (lm_status == LM_MAXITER) {np.SetStatus(NL_LM_MAXITER); return(np.Status());} // Ouch!
    // Set new cf value and parameters
    np.SetPar(np.Par() + minp.first*p);
    // Check for convergence based on small decrease of cost-function
    if (zero_cf_diff_conv(np.CF(),minp.second,np.FractionalCFTolerance())) {np.SetCF(minp.second); np.SetStatus(NL_CFCONV); return(np.Status());}
    // Check for convergence based on neglible move in parameter space
    else if (zero_par_step_conv(minp.first*p,np.Par(),np.FractionalParameterTolerance())) {np.SetCF(minp.second); np.SetStatus(NL_PARCONV); return(np.Status());}
    else {  // If no covergence
      np.SetCF(minp.second);
      if (((np.NIter())%np.NPar()) == 0) {                          // Explicitly reset directions after npar iterations
        r = -cfo.grad(np.Par());
        p = r;
      }
      else {
        ColumnVector oldr = r;
        r = -cfo.grad(np.Par());
        if (np.ConjugateGradientUpdate() == CG_FR) {              // Get conjugate direction Fletcher-Reeves flavour
          p = r + (DotProduct(r,r)/DotProduct(oldr,oldr)) * p;
	}
        else if (np.ConjugateGradientUpdate() == CG_PR) {         // Get conjugate direction Polak-Ribiere flavour
          p = r + (DotProduct(r-oldr,r)/DotProduct(oldr,oldr)) * p;
        }
      }
    }
  }
  // If we get here we have used too many iterations
  np.SetStatus(NL_MAXITER);
    
  return(np.Status());
}

// Main routine for scaled conjugate-gradient optimisation. The
// idea of the algorithm is similar to that of Levenberg-
// Marquardt. In the LM algorithm the search direction is a 
// "compromise" between the Newton direction and the gradient
// direction, where the compromise depends on a factor lambda.
// A large lambda means that it is close to the gradient and a 
// small lambda that it is close to the Newton direction. The 
// value of lambda is updated each iteration depending on the
// success of the last step. In this method the compromise is
// between the "conjugate gradient" direction and the gradient
// direction. The variable names follow the (excellent!) paper 
// by Martin Möller (1993) Neural Networks 6:525-533.
// I have tried to follow the notation he uses in his paper, thus
// enabling that to be the "documentation" for the routine below.

NonlinOut sccngr(const NonlinParam& np, const NonlinCF& cfo)
{
  // Set up initial values
  np.SetCF(cfo.cf(np.Par()));           // Current value for cost-function (E in Moller 92).
  double sigma = 1.0e-2;                // Step-length when estmating H*p from g(w+sigma*p)-g(w)
  double lambda_bar = 0.0;              // Update for lambda if approximate hessian not positive definite
  ColumnVector r = -cfo.grad(np.Par()); // Negative gradient
  ColumnVector p = r;                   // Search direction
  bool success = true;                  // True if previous step was successful
  double delta = 0.0;                   // Used to check pos def of H in loop below
  ColumnVector s(np.NPar());            // Used as approximation to H*p in loop below 

  while (np.NextIter()) {
    double p2 = DotProduct(p,p);                            // p'*p, Temporary variable to save some time
    if (success == true) {                                  // If last step led to reduction of cost-function
      double sigma_k = sigma/std::sqrt(p2);                      // Normalised step-length when estimating H*p
      // cout << "np.NIter() = " << np.NIter() << ", p2 = " << p2 << ", sigma_k = " << sigma_k << endl;
      s = (cfo.grad(np.Par()+sigma_k*p) + r) / sigma_k;     // Approximation to H*p
      delta = DotProduct(p,s);                              // Approximation to p'*H*p
    }
    s += (np.Lambda()-lambda_bar)*p;                      // Equivalent to adding (l-lb)*I to H
    delta += (np.Lambda()-lambda_bar)*p2;                 // If <0 then H+(l-lb)*I not positive definite
    if (delta <= 0) {                                     // If it H is not positive definite
      s += (np.Lambda() - 2.0*(delta/p2)) * p;            // Make H more diagonal dominant to ensure pos def
      lambda_bar = 2.0*(np.Lambda() - delta/p2);
      delta  = np.Lambda()*p2 - delta;
      np.SetLambda(lambda_bar);
    }
    double mu = DotProduct(p,r);
    double alpha = mu/delta;                              // Step size in direction p
    double tmp_cf = cfo.cf(np.Par()+alpha*p);             // Value of cost-function at attempted new point

    // cout << "np.NIter() " << np.NIter() << ", delta = " << delta << ", mu = " << mu << ", alpha = " << alpha << endl;

    /*
    char fname[100]; 
    sprintf(fname,"scg_debug_gradient_%02d.txt",np.NIter());
    print_newmat(r,fname);
    sprintf(fname,"scg_debug_step_%02d.txt",np.NIter());
    ColumnVector  step(p); step *= alpha;
    print_newmat(step,fname);
    */
    
    
    double Delta = 2.0*delta*(np.CF()-tmp_cf) / (mu*mu);  // > 0 means attempted step reduced cost-function
    if (Delta >= 0) {                                     // If step reduces cost-function
      np.SetCF(tmp_cf);                                   // Update lowest observed value of cost-function
      np.SetPar(np.Par() + alpha*p);                      // Update best set of parameters
      lambda_bar = 0.0;
      success = true;
      if ((np.NIter()%np.NPar()) == 0) {                  // If npar iterations since last resetting of directions
        r = -cfo.grad(np.Par());                          // Reset search direction to negative gradient
        p = r;
      }
      else {
        ColumnVector oldr = r;
        r = -cfo.grad(np.Par());
        double beta = (DotProduct(r,r)-DotProduct(oldr,r)) / mu;
        // cout << "np.NIter() = " << np.NIter() << ", beta = " << beta << endl;
        p = r + beta*p;                              // New search direction
      } 
      if (Delta > 0.75) {                            // If attempted step was \emph{REALLY} good
	np.SetLambda(np.Lambda()/2.0);
      }
    }
    else {                                           // If step doesn't reduce cost-function
      lambda_bar = np.Lambda();
      success = false;
    }
    if (Delta < 0.25) {                              // If step reduced cost-function only "a little" (or not at all)
      np.SetLambda(4.0*np.Lambda());
    }
    if (zero_grad_conv(np.Par(),r,np.CF(),np.FractionalGradientTolerance())) {  // If gradient is (practically) zero
      np.SetStatus(NL_GRADCONV); return(np.Status());
    }
  }
  // If we get here we have exceeded allowed # of iterations
  np.SetStatus(NL_MAXITER);
  return(np.Status());
}
 
// Main routine for variable-metric optimisation. This implements
// the variable-metric optimisation with the BFGS or DFP updating
// schemes. The implementation details are mostly quite close to
// those described in Numerical Recipies in C.

NonlinOut varmet(const NonlinParam& p, const NonlinCF& cfo)
{
  // Get scale factor to ensure a relative scale beteween
  // parameters and cost-function such that fast and robust
  // convergence is acheieved.

  double sf = cfo.sf();                      // Suggestion by "user"
  ColumnVector grad = sf*cfo.grad(p.Par());  // Gradient of const-function
  if (p.VariableMetricAutoScale()) {
    sf = scale_factor(p.Par(),-grad,cfo,p.LineSearchMaxIterations(),sf);         // Refinement by "me"
    if (sf == 0.0) {                                                             // No minimum in indicated direction
      p.SetStatus(NL_PARCONV);                                                   // Assume this means we are already at minimum
      return(p.Status());
    }
    grad = (sf/cfo.sf()) * grad;
  }

  VarmetMatrix  iH(p.NPar(),VM_OPT,p.VariableMetricUpdate());   // Inverse Hessian
  p.SetCF(sf*cfo.cf(p.Par()));                                // Current value of cost-function
  ColumnVector  pdir = -(iH*grad);                            // Direction to search in

  double        lambda = 0.0;     // Step-length returned by linsrch
  double        newcf = 0.0;      // New value for cost-function
  ColumnVector  newpar(p.NPar()); // New point in parameter space

  while (p.NextIter()) {
    // Do a line-search to find a new point in parameter space
    LinOut status = linsrch(pdir,p.Par(),grad,cfo,p.CF(),sf,p.LineSearchMaxIterations(),
                            p.LineSearchMaxStep(),p.VariableMetricAlpha(),
                            p.LineSearchFractionalParameterTolerance(),&lambda,&newcf,&newpar);
    // Check for convergence/problems based on outcome of linsrch
    if (status == LM_MAXITER) {p.SetStatus(NL_LM_MAXITER); return(p.Status());}
    else if (status == LM_LAMBDA_NILL) { // This means we might be heading uphill and should restart
      if (p.NextRestart()) { // If we have spare restarts
        p.SetCF(p.CF());      // Another copy of old value
        p.SetPar(p.Par());    // Another copy of old values
        iH.reset();           // Back to being unity matrix
        pdir = -grad;
        continue;
      }
      else {  
        p.SetStatus(NL_PARCONV); return(p.Status());
      }
    }
    // Test for convergence based on distance between points in parameter space
    ColumnVector dpar = newpar - p.Par();
    p.SetPar(newpar);
    p.SetCF(newcf);
    // cout << "p.FractionalParameterTolerance() = " << p.FractionalParameterTolerance() << endl;
    // cout << "P.Par() = " << p.Par() << endl;
    // cout << "dpar = " << dpar << endl;
    if (zero_par_step_conv(p.Par(),dpar,p.FractionalParameterTolerance())) {p.SetStatus(NL_PARCONV); return(p.Status());}
    // Get gradient at new point
    ColumnVector newgrad = sf*cfo.grad(p.Par());
    // Test for convergence based on "zero" gradient
    if (zero_grad_conv(p.Par(),newgrad,p.CF(),p.FractionalGradientTolerance())) {p.SetStatus(NL_GRADCONV); return(p.Status());}
    // Update estimate of inverse Hessian
    iH.update(dpar,newgrad-grad);
    // Update parameters and get new direction to go in
    grad = newgrad;
    pdir = -(iH*grad);      // N.B. no unary - op for iH, parenthesis necessary
  }

  // If we get here we have exceeded the allowed # of iterations
  p.SetStatus(NL_MAXITER);

  return(p.Status());  
}  

LinOut linsrch(// Input
               const ColumnVector&  dir,    // Search direction
               const ColumnVector&  p0,      // Current parameter values
               const ColumnVector&  grad,    // Gradient at p0
               const NonlinCF&      cfo,     // Cost-function object
               double               f0,      // Current cost-function value
               double               sf,      // Scale factor for cost-function
               double               maxiter, // Max # of iterations
               double               sm,      // Stepmax
               double               alpha,   // Alpha (sorry).
               double               ptol,    // Tolerance in parameter space
               // Output
               double               *lambda, // Resulting step length
               double               *of,     // Value of cost-function on output
               ColumnVector         *np)     // New parameters
{
  const double lmin = 0.1;             
  const double lmax = 0.5; 

  // First make sure that the step-length suggested
  // by pdir isn't completely unreasonable.

  double totstep=std::sqrt(DotProduct(dir,dir));
  ColumnVector pdir(dir);
  if (totstep > sm) {pdir *= sm/totstep;}

  // Calculate expected rate of change in the direction 
  // given by pdir.

  double fp0 = DotProduct(grad,pdir);

  // Calculate smallest meaningful lambda given what is
  // smallest meaningful change in parameter value.

  double almin=0.0;
  for (int i=0; i<p0.Nrows(); i++) {
    almin = std::max(almin,std::abs(pdir.element(i))/std::max(std::abs(p0.element(i)),1.0));
  }
  almin = ptol / almin;

  // First try a step of full lambda
  
  *lambda = 1.0;                        // Start with that
  (*np) =  p0 + (*lambda)*pdir;         // First new parameters to try
  double f2 = sf * cfo.cf(*np);         // Cost-function value for par

  // See if that does it (fat chance!)
  if (f2 < f0 + alpha*(*lambda)*DotProduct(grad,(*np)-p0)) {*of = f2; return(LM_CONV);}

  // Calculate Quadratic based on f(0), f'(0)
  // and f(1) and find minimum of that quadratic.

  *lambda = - fp0 / (2.0*(f2-f0-fp0));  // Minumum of f(lambda)
  // Make sure new lambda is 0.1*old_l < lambda < 0.5*old_l
  *lambda = std::max(lmin,*lambda);
  *lambda = std::min(lmax,*lambda);
  (*np) =  p0 + (*lambda)*pdir;         // Second set of new parameters to try
  double f1 = sf * cfo.cf(*np);              // Cost-function value for par

  // Now we will start fitting cubics to f(0), f'(0),
  // f(lambda_1) and f(lambda_2) where lambda_1 is the
  // lambda most recently tested and where lambda_2 is
  // the lambda tested second to last.

  double       l2 = 1.0;     // Second to last lambda
  double       l1 = *lambda; // Last lambda
  Matrix       X(2,2);
  ColumnVector y(2);       
  
  for (int iter=0; iter<maxiter; iter++) {
    // See if present lambda might be too small
    if (*lambda < almin) {*of = f1; return(LM_LAMBDA_NILL);}
    // See if present value is acceptable
    if (f1 < f0 + alpha*(*lambda)*DotProduct(grad,(*np)-p0)) {*of = f1; return(LM_CONV);}
    // Find parameter values for cubic and square on lambda
    X << std::pow(l1,3.0) << std::pow(l1,2.0) << std::pow(l2,3.0) << std::pow(l2,2.0);
    y << f1-fp0*l1-f0 << f2-fp0*l2-f0;
    ColumnVector b = X.i()*y;
    // Find value for lambda that yield minimum of cubic
    *lambda = (-b.element(1) + std::sqrt(std::pow(b.element(1),2.0) - 3.0*b.element(0)*fp0)) / (3.0*b.element(0));
    // Make sure new lambda is 0.1*old_l < lambda < 0.5*old_l
    *lambda = std::max(lmin*l1,*lambda);
    *lambda = std::min(lmax*l1,*lambda);
    // Get new function value and update parameters
    f2 = f1;
    (*np) = p0 + (*lambda)*pdir;
    f1 = sf * cfo.cf(*np);
    l2 = l1;
    l1 = *lambda;
  }

  // If we are here we have exceeded # of iterations
  
  *of = f1;
  return(LM_MAXITER);
}

// Will try and find a scale factor for the cost-function such that
// the step length (lambda) for the first iteration of the variable-
// metric method is ~0.25. Empricially I have found that such a scaling
// that yields a first step length in the range 0.1 -- 0.5 will yield 
// robust and fast convergence. N.B. though that I have only tested that
// for non-linear reg, and it is concievable that it is different for
// other applications with fewer parameters.
 
double scale_factor(const ColumnVector&  p,       // Current parameter values
                    const ColumnVector&  pdir,    // Search direction
                    const NonlinCF&      cfo,     // Cost-function object
                    int                  maxiter, // Max # of iterations
                    double               sf)      // Scale factor.
{
  const double        dl = 0.25;   // Desired Lambda
  const double        ftol = 0.01; // Fractional tolerance for minimum

  // Start out by finding an upper bound for lambda such that
  // a minimum is guaranteed to be between 0 and rp
  
  pair<double,double> lp;
  pair<double,double> mp;
  pair<double,double> rp = bracket(p,pdir,cfo,ftol,sf,&lp,&mp);
  if (rp == mp) { // If there is no minimum in the indicated direction
    return(0.0);
  }

  // Now find a minimum with a fractional accuracy of ~1%

  pair<double,double> minpoint;
  
  if (linmin(p,pdir,cfo,sf,lp,mp,rp,ftol,maxiter,&minpoint) == LM_MAXITER) {
    throw NonlinException("Failed to find minimum along search direction");
  }

  sf *= minpoint.first/dl;
  return(sf);
}

// Will find the minimum of the cost-function as a function of 
// lambda. This routine will find the minimum to a fractional
// tolerance of lambda. This is NOT practical to use for the
// Variable Metric minimisation (too slow), but is used for
// the conjugate-gradient method and for finding the initial 
// scaling between parameters and cost-function for the 
// variable metric method.

LinOut linmin(// Input
              const ColumnVector&   p,      // Current parameter values
              const ColumnVector&   pdir,   // Search direction
              const NonlinCF&       cfo,    // Cost-function object
              double                sf,     // Scale factor for cost-function
              pair<double,double>   lp,     // Left point
              pair<double,double>   mp,     // Point somewhere in interval
              pair<double,double>   rp,     // Right point
              double                ftol,   // Fractional tolerance
              int                   maxiter,// Max # of iterations
              // Output
              pair<double,double>   *x)     // Best point
{
  const double         gold = 0.382;// Golden section
  pair<double,double>  test;        // New point to test
  pair<double,double>  w = mp;      // Second best point
  pair<double,double>  v = mp;      // Last value of second best point
  double               step = 0.0;
  double               ostep = 0.0; // Length of 2nd to last step taken
  double               d = 0.0;     // Length of last step taken
  ColumnVector         y(3);        // Used for fitting parabolic
  Matrix               X(3,3);      // Used for fitting parabolic
  *x = mp;                          // Initialise "best" point

  for (int i=0; i<maxiter; i++) {
    double midp = (rp.first+lp.first)/2.0;                       // Midpoint of bracketing points
    double tol = 2.0*ftol*std::abs(x->first)+MISCMATHS::EPS;          // Std::Absolute tolerance
    if  (std::abs(x->first-midp) <= (tol-0.5*(rp.first-lp.first))) {  // Convergence check
      return(LM_CONV);
    }
    // Try parabolic fit, but not before third iteration
    double tmp = 10.0*std::sqrt(MISCMATHS::EPS);
    if (std::abs(ostep) > tol/2.0 &&           // If second to last step big enough
        std::abs(x->first-w.first) > tmp && 
        std::abs(x->first-v.first) > tmp && 
        std::abs(w.first-v.first) > tmp) {     // And points not degenerate
      step = ostep;
      ostep = d;
      y << x->second << w.second << v.second;
      X << std::pow(x->first,2.0) << x->first << 1.0 <<
	   std::pow(w.first,2.0) << w.first << 1.0 <<
	   std::pow(v.first,2.0) << v.first << 1.0;
      ColumnVector b = X.i() * y;
      if (b.element(0) < 4*MISCMATHS::EPS ||                   // If on line or going for maximum
          (test.first = -b.element(1)/(2.0*b.element(0))) <= lp.first 
          || test.first >= rp.first ||                         // If outside bracketed interval
          std::abs(test.first-x->first) > 0.5*step) {               // Or if step too big (indicates oscillation)
        // Take golden step into larger interval
        if (rp.first-x->first > x->first-lp.first) {           // If right interval larger
          test.first = x->first + gold * (rp.first - x->first);
        }
        else {
          test.first = x->first - gold * (x->first - lp.first); 
        }
      }
    }
    else { // Take golden step into larger interval   
      if (x->first < midp) {            // If right interval larger
        ostep = rp.first - x->first;
        test.first = x->first + gold * (rp.first - x->first);
      }
      else {
        ostep = x->first - lp.first;
        test.first = x->first - gold * (x->first - lp.first); 
      }
    }
    d = test.first - x->first;                              // Signed length of step
    test.second = sf*cfo.cf(p+test.first*pdir);                // Evaluate cf at new point
    // Now we have a new point, and we need to figure out what to do with it
    if (test.second <= x->second) {  // If it beats the best step
      if (test.first > x->first) {lp = *x;}
      else {rp = *x;}
      v = w; w = *x; *x = test;
    }
    else {
      if (test.first < x->first) {lp = test;}
      else {rp = test;}
      if (test.second <= w.second || w.first == x->first) {
        v = w; w = test;
      }
      else if (test.second <= v.second || v.first == x->first || v.first == w.first) {
        v = test;
      }
    }
  }
  // If we are here we have used too many iterations
  return(LM_MAXITER); // Error status    
}
              
// Will return a value lambda such that a function minimum is guaranteed to
// lie somewhere in the interval between p and p+lambda*pdir. The second value
// of the returned pair is the cost-function value at the point.

pair<double,double> bracket(// Input 
                            const ColumnVector& p,      // Current parameter values
                            const ColumnVector& pdir,   // Search direction
                            const NonlinCF&     cfo,    // Cost-function object
                            double              ptol,   // Relative tolerance for parameter values
                            double              sf,     // Scale factor of cost-function
                            // Output
                            pair<double,double> *p_0,   // Cost-function value at p
                            pair<double,double> *p_m)   // Point between p_0 and p_l
{
  pair<double,double>  p_l;
  const double gr = 0.618034;
  const double maxstep = 100.0;
  p_0->first = 0.0;
  double cf0 = sf*cfo.cf(p);
  double l1 = 1.0;
  double cf1 = sf*cfo.cf(p+l1*pdir);

  // Find maximum relative component of search direction

  double test = 0.0;
  for (int i=0; i<pdir.Nrows(); i++) {test = std::max(test,std::abs(pdir.element(i))/std::max(p.element(i),1.0));}

  // Do a crude initial search for order of magnitude

  while (!isfinite(cf1)) {
    l1 *= 0.1;
    cf1 = sf*cfo.cf(p+l1*pdir);
  }

  // Get third point to get us started
  double l2 = 0.0;
  double cf2 = 0.0;
  if (cf1 < cf0) {l2 = (2.0 + gr)*l1; cf2 = sf*cfo.cf(p+l2*pdir);}
  else {l2 = l1; cf2 = cf1; l1 = gr * l2; cf1 = sf*cfo.cf(p+l1*pdir);}
  // Check if we already have a bracket
  if (cf1 < cf0 && cf1 < cf2) {
    p_l.first = l2; p_l.second = cf2;
    p_m->first = l1; p_m->second = cf1;
    p_0->second = cf0; 
    return(p_l);
  }

  Matrix        X(2,2);
  ColumnVector  y(2);
  double        lt = 0.0;
  double        cft = 0.0;

  while (!(cf1 < cf0 && cf1 < cf2)) {  // If minimum still not bracketed
    if (l2*test < ptol) {              // If interval ridicously small
      p_l = *p_0;
      *p_m = *p_0;
      return(p_l);
    }
    // Let's see if a parabolic might help us
    if (std::abs(l2-l1) > 10.0*std::sqrt(MISCMATHS::EPS)) {
      X << std::pow(l1,2.0) << l1 << std::pow(l2,2.0) << l2;
      y << cf1 << cf2;
      ColumnVector b = X.i()*y;
      if (b.element(0) > 4.0*MISCMATHS::EPS) {           // Check they are not on a line and not for maximum
        lt = - (b.element(1) / (2.0 * b.element(0))); // Tentative point
        if (lt > 0 && lt < l2) {                      // If in range of previous points
          cft = sf*cfo.cf(p+lt*pdir);
          if (cft < cf0 && cft < cf2) {l1=lt; cf1=cft; continue;}
          else if (cft > cf0 && lt < l1 && cft < cf1) {l2=l1; cf2=cf1; l1=lt; cf1=cft; continue;}
          else if (cft > cf0 && lt > l1 && cft < cf2) {l2=lt; cf2=cft; continue;}
        }
        else if (lt > 0 && lt < maxstep*l2) {       // If expansion in allowed range
          cft = sf*cfo.cf(p+lt*pdir);
          l1 = l2; cf1 = cf2; 
          l2 = lt; cf2 = cft;
          continue;
        }
      }
    }
    // If we are here the parabolic was of no use
    if (cf2 < cf0) { // We need to expand
      lt = (2.0 + gr)*l2;
      cft = sf*cfo.cf(p+lt*pdir);
      l1 = l2; cf1 = cf2;
      l2 = lt; cf2 = cft;
    }
    else { // We need to contract
      lt = gr * l1;
      cft = sf*cfo.cf(p+lt*pdir);
      l2 = l1; cf2 = cf1;
      l1 = lt; cf1 = cft;
    }     
  }  

  // If we are here we know that there is a minimum
  // somewhere between 0 and l2;

  p_0->second = cf0;
  p_m->first = l1;
  p_m->second = cf1;
  p_l.first = l2;
  p_l.second = cf2;
  return(p_l);   
}  

// Utility routines that checks for convergence based on various criteria

// Based on zero (neglible) gradient

bool zero_grad_conv(const ColumnVector&   par,
                    const ColumnVector&   grad,
                    double                cf,
                    double                gtol)
{
  double test = 0.0;     // test will be largest relative component of gradient
  for (int i=0; i<par.Nrows(); i++) {
    test = std::max(test,std::abs(grad.element(i))*std::max(std::abs(par.element(i)),1.0));
  }
  test /= std::max(cf,1.0);   // Protect against near-zero values for cost-function

  return(test < gtol);
}

// Based on zero (neglible) decrease in cost-function

bool zero_cf_diff_conv(double cfo,
                       double cfn,
                       double cftol)
{
  return(2.0*std::abs(cfo-cfn) <= cftol*(std::abs(cfo)+std::abs(cfn)+MISCMATHS::EPS));
}

// Based on zero (neglible) step in parameter space

bool zero_par_step_conv(const ColumnVector&   par,
                        const ColumnVector&   step,
                        double                ptol)
{
  double test = 0.0;
  for (int i=0; i<par.Nrows(); i++) {
    test = std::max(test,std::abs(step.element(i))/std::max(std::abs(par.element(i)),1.0));
  }
  return(test < ptol);
}

// Utility routines that allow the user to check accuracy of their own grad and hess functions

pair<ColumnVector,ColumnVector> check_grad(const ColumnVector&  par,
                                           const NonlinCF&      cfo)
{
  pair<ColumnVector,ColumnVector> rv;
  
  rv.first = cfo.NonlinCF::grad(par);
  rv.second = cfo.grad(par);

  return(rv);
}

pair<boost::shared_ptr<BFMatrix>,boost::shared_ptr<BFMatrix> > check_hess(const ColumnVector& par,
                                                                          const NonlinCF&     cfo)
{
  pair<boost::shared_ptr<BFMatrix>,boost::shared_ptr<BFMatrix> > rv;
  
  rv.first = cfo.NonlinCF::hess(par);
  rv.second = cfo.hess(par);

  return(rv);
}     
           

void print_newmat(const NEWMAT::GeneralMatrix&  m,
                  std::string                   fname)
{
  if (!fname.length()) {
    cout << endl << m << endl;
  }
  else {
    try {
      std::ofstream  fout(fname.c_str());
      fout << setprecision(10) << m;
    }
    catch(...) {
      std::string  errmsg("print_newmat: Failed to write to file " + fname);
      throw NonlinException(errmsg);
    }
  }
}
  
} // End namespace MISCMATHS
