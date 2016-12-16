/*  t2z.h

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(__t2z_h)
#define __t2z_h
  
#include <iostream>
#include <fstream>
#include "newmatap.h"
#include "newmatio.h"
#include "base2z.h"

using namespace NEWMAT;

namespace MISCMATHS {

  class T2z : public Base2z
    {
    public:
      static T2z& getInstance();
      ~T2z() { delete t2z; }
      
      float convert(float t, int dof);
      float converttologp(float t, int dof);
   
      static void ComputePs(const ColumnVector& p_vars, const ColumnVector& p_cbs, int p_dof, ColumnVector& p_ps);
      static void ComputeZStats(const ColumnVector& p_vars, const ColumnVector& p_cbs, int p_dof, ColumnVector& p_zs);
      static void ComputeZStats(const ColumnVector& p_vars, const ColumnVector& p_cbs, const ColumnVector& p_dof, ColumnVector& p_zs);

    private:
      T2z() : Base2z()
	{}
      
      const T2z& operator=(T2z&);
      T2z(T2z&);
      
      bool issmalllogp(float logp);
      bool islarget(float t, int dof, float &logp);
      float larget2logp(float t, int dof);

      static T2z* t2z;

    };
 
  inline T2z& T2z::getInstance(){
    if(t2z == NULL)
      t2z = new T2z();
  
    return *t2z;
  }



  class Z2t
    {
    public:
      static Z2t& getInstance();
      ~Z2t() { delete z2t; }
      
      float convert(float t, int dof);

    private:
      Z2t()
	{}
      
      const Z2t& operator=(Z2t&);
      Z2t(Z2t&);

      static Z2t* z2t;

    };

  inline Z2t& Z2t::getInstance(){
    if(z2t == NULL)
      z2t = new Z2t();
  
    return *z2t;
  }

}

#endif

