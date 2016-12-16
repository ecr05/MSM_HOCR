/*! \file Simplex.h
    \brief Contains declaration of Simplex class that can be used for Nelder-Mead simplex minimisation.

    \author Jesper Andersson
    \version 1.0b, Oct., 2013.
*/
// Contains declaration of Simplex class that can 
// be used for Nelder-Mead simplex minimisation.
//
// Simplex.h
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2013 University of Oxford 
//

#ifndef Simplex_h
#define Simplex_h

#include <iostream>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>
#include "newmat.h"
#include "miscmaths.h"
#include "nonlin.h"

namespace MISCMATHS {

/****************************************************************//**
*
* \brief Class used for implementing Nelder-Mead minimisation.
*
* Implements a class for implementing Nelder-Mead downhill slope
* simplex minimisation. It implements the full minimisation through
* its member function Minimise, but it also provides a set of utility
* functions for anyone wanting to provide a slightly different 
* implementation.
*
********************************************************************/ 
class Simplex
{
public:
  Simplex(const NEWMAT::ColumnVector& p,
	  const MISCMATHS::NonlinCF&   cf);
  Simplex(const NEWMAT::ColumnVector& p,
	  const MISCMATHS::NonlinCF&   cf,
	  const NEWMAT::ColumnVector& l); 
  ~Simplex() {}
  /// Minimises the costunction until the fractional difference between points in simplex is < ftol
  bool Minimise(double ftol, unsigned int miter);
  /// Checks if the fractional difference between points in simplex is < ftol
  bool HasConverged(double ftol) const { 
    return(2.0*std::abs(WorstFuncVal()-BestFuncVal()) < ftol*(std::abs(BestFuncVal())+std::abs(WorstFuncVal()))); 
  }
  /// Returns the number of parameters
  unsigned int NoPar() const { return(static_cast<unsigned int>(_sp.Nrows())); }
  /// Returns the "best" (lowest function value) parameters in the simplex
  const NEWMAT::ColumnVector& BestPar() const { return(_smx[_bsti]); }
  /// Returns the "best" (lowest) function value in the simplex
  double BestFuncVal() const { return(_fv[_bsti]); }
  /// Returns the 2nd to worst (highest) function value in the simplex
  double SecondWorstFuncVal() const { return(_fv[_nwsti]); }
  /// Returns the worst (highest) function value in the simplex
  double WorstFuncVal() const { return(_fv[_wrsti]); }
  /// Reflects the worst point through the average of the remaining points
  double Reflect();
  /// Expands upon a previous reflexion
  double Expand();
  /// Contracts the worst point half-way towards the average of the remaining points
  double Contract();
  /// Contracts all except the best point half-way towards the best point
  void MultiContract();
  /// Update the indicies for best, worst and 2nd to worst points.
  void UpdateRankIndicies();
private:
  const MISCMATHS::NonlinCF&                           _cf;
  const NEWMAT::ColumnVector                           _sp;
  std::vector<NEWMAT::ColumnVector>                    _smx;
  std::vector<double>                                  _fv;
  unsigned int                                         _bsti;     // Best (lowest function value) point
  unsigned int                                         _wrsti;    // Worst (highest function value) point
  unsigned int                                         _nwsti;    // Second worst (2nd highest function value) point
  NEWMAT::ColumnVector                                 _rp;       // Latest reflexion point

  void setup_simplex(const NEWMAT::ColumnVector& l);
  void calculate_reflexion_point(unsigned int ii);
};

} // End namespace MISCMATHS

#endif // End #ifndef Simplex_h
