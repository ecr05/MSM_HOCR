/*! \file Simplex.cpp
    \brief Contains definitions of Simplex class that can be used for Nelder-Mead simplex minimisation.

    \author Jesper Andersson
    \version 1.0b, Oct., 2013.
*/
// Contains definitions of Simplex class that can 
// be used for Nelder-Mead simplex minimisation.
//
// Simplex.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2013 University of Oxford 
//

#include <iostream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include "newmat.h"
#include "miscmaths.h"
#include "nonlin.h"
#include "Simplex.h"

using namespace MISCMATHS;

/****************************************************************//**
*
* Constructs a Simplex object given a vector of starting guesses
* and a cost-function object derived from the NonlinCF base class.
* The simplex is created by placing n points at a distance 0.1p
* from the nx1 start guess. If p is zero along any dimension the
* corresponding point is placed at unity distance instead.
* \param p Starting guess for the parameters.
* \param cf Cost-function object of a class derived from the virtual
*        NonlinCF base class.
*
********************************************************************/ 
Simplex::Simplex(const NEWMAT::ColumnVector& p,
		 const MISCMATHS::NonlinCF&  cf)
: _cf(cf), _sp(p)
{
  NEWMAT::ColumnVector l(_sp.Nrows());
  for (int i=0; i<l.Nrows(); i++) {
    l(i+1) = (p(i+1)) ? 0.1*p(i+1) : 1.0;
  }
  setup_simplex(l);
  UpdateRankIndicies();
}

/****************************************************************//**
*
* Constructs a Simplex object given a vector of starting guesses,
* a cost-function object derived from the NonlinCF base class and
* a vector of perturbations used to build the simplex.
* \param p Starting guess for the parameters.
* \param cf Cost-function object of a class derived from the virtual
*        NonlinCF base class.
* \param l Perturbations used to create the nodes of the simplex.
*          On initialisation the ith point of the simplex will be
*          smplx[i] = p; smplx[i](i) += l(i);
*
********************************************************************/ 
Simplex::Simplex(const NEWMAT::ColumnVector& p,
		 const MISCMATHS::NonlinCF&  cf,
		 const NEWMAT::ColumnVector& l) 
: _cf(cf), _sp(p) 
{ 
  if (l.Nrows() != _sp.Nrows()) throw ;
  setup_simplex(l); 
  UpdateRankIndicies();
}

/****************************************************************//**
*
* Minimises the simplex, i.e. it will find the set of parameters
* that minimmises the NonlinCF cost-function that the Simplex
* object was constructed with.
* Suggested use:
* \code
* Simplex my_simplex(my_guess_par,my_cost_func);
* if (my_simplex(my_tol,1000)) {
*   ColumnVector optimal_par = my_simplex.BestPar();
* }
* else { cout << "bugger" << endl; }
* \endcode
* 
* \param ftol Fractional cost-function tolerance for convergence.
* \param miter Maximum allowed number of iterations.
*
********************************************************************/ 
bool Simplex::Minimise(double ftol, 
		       unsigned int miter)
{
  UpdateRankIndicies(); // Make sure it is ready for use

  for (unsigned int i=0; i<miter; i++) {
    if (HasConverged(ftol)) return(true); // Check for convergence 
    
    double newf = Reflect();   // Attempt reflexion
    // Extend into an expansion if reflexion very successful
    if (newf <= BestFuncVal()) { 
      Expand(); // Attempt expansion
    }
    else if (newf >= SecondWorstFuncVal()) {
      double worst_fval = WorstFuncVal();
      newf = Contract();     // Do a contraction towards plane of "better" points
      if (newf >= worst_fval) {    // Didn't work. Contract towards best point
	MultiContract();
      }
    } 
    UpdateRankIndicies();
  }
  return(false);
}
  
double Simplex::Reflect()
{
  calculate_reflexion_point(_wrsti); // Updates _rp
  NEWMAT::ColumnVector newp = 2.0*_rp - _smx[_wrsti];
  double newf = _cf.cf(newp);
  if (newf < _fv[_wrsti]) {
    _smx[_wrsti] = newp;
    _fv[_wrsti] = newf;
  }
  return(newf);
}

double Simplex::Expand()
{
  NEWMAT::ColumnVector newp = 2.0*_smx[_wrsti] - _rp;
  double newf = _cf.cf(newp);
  if (newf < _fv[_wrsti]) {
    _smx[_wrsti] = newp;
    _fv[_wrsti] = newf;
  }
  return(newf);
}

double Simplex::Contract()
{
  NEWMAT::ColumnVector newp = 0.5 * (_smx[_wrsti] + _rp);
  double newf = _cf.cf(newp);
  if (newf < _fv[_wrsti]) {
    _smx[_wrsti] = newp;
    _fv[_wrsti] = newf;
  }
  return(newf);
}

void Simplex::MultiContract()
{
  for (unsigned int i=0; i<_smx.size(); i++) {
    if (i != _bsti) {
      _smx[i] = 0.5 * (_smx[i] + _smx[_bsti]);
      _fv[i] = _cf.cf(_smx[i]);
    }
  }
  return;
}

void Simplex::UpdateRankIndicies()
{
  double minv = std::numeric_limits<double>::max(); 
  double maxv = -std::numeric_limits<double>::max(); 
  for (unsigned int i=0; i<_fv.size(); i++) {
    if (_fv[i] < minv) { minv = _fv[i]; _bsti = i; }
    if (_fv[i] > maxv) { maxv = _fv[i]; _wrsti = i; }
  }
  maxv = -std::numeric_limits<double>::max(); 
  for (unsigned int i=0; i<_fv.size(); i++) {
    if (i != _wrsti) {
      if (_fv[i] > maxv) { maxv = _fv[i]; _nwsti = i; }
    }
  }
  return;
}

void Simplex::setup_simplex(const NEWMAT::ColumnVector& l)
{
  _smx.resize(_sp.Nrows()+1);
  _fv.resize(_smx.size());
  _smx[0] = _sp;
  _fv[0] = _cf.cf(_smx[0]);
  for (int i=1; i<=_sp.Nrows(); i++) {
    _smx[i] = _sp;
    _smx[i](i) += l(i);
    _fv[i] = _cf.cf(_smx[i]);
  }
  return;
}
void Simplex::calculate_reflexion_point(unsigned int ii)
{
  if (_rp.Nrows() != _sp.Nrows()) _rp.ReSize(_sp.Nrows());
  _rp = 0.0;
  for (unsigned int i=0; i<_smx.size(); i++) {
    if (i != ii) _rp += _smx[i];
  }
  _rp /= static_cast<double>(_rp.Nrows());
}
