/*  f2z.h

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(__f2z_h)
#define __f2z_h
  
#include <iostream>
#include <fstream>
#include "newmatap.h"
#include "newmatio.h"
#include "base2z.h"
//#include "miscmaths.h"

using namespace NEWMAT;

namespace MISCMATHS {

  class F2z : public Base2z
    {
    public:
      static F2z& getInstance();
      ~F2z() { delete f2z; }
      
      float convert(float f, int d1, int d2);

      static void ComputeFStats(const ColumnVector& p_fs, int p_dof1, int p_dof2, ColumnVector& p_zs);
      static void ComputeFStats(const ColumnVector& p_fs, int p_dof1, const ColumnVector& p_dof2, ColumnVector& p_zs);
      static void ComputeFStats(const ColumnVector& p_fs, const ColumnVector& p_dof1, const ColumnVector& p_dof2, ColumnVector& p_zs);

    private:
      F2z() : Base2z()
	{}
      
      const F2z& operator=(F2z&);
      F2z(F2z&);
      
      bool issmalllogp(float logp);
      bool islargef(float t, int d1, int d2, float &logp);
      float largef2logp(float t, int d1, int d2);

      static F2z* f2z;

    };
 
  inline F2z& F2z::getInstance(){
    if(f2z == NULL)
      f2z = new F2z();
  
    return *f2z;
  }

}

#endif
