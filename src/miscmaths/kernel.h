/*  kernel.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

/*  CCOPYRIGHT  */

// General kernel interpolation class

#if !defined(kernel_h)
#define kernel_h

#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include "newmat.h"
      
using namespace NEWMAT;
using namespace std;

namespace MISCMATHS {

  /////////////////////////////////////////////////////////////////////////

  // Interpolation kernel storage class

  class kernelstorage
    {
    private:
      // NB: all widths are kernel half-widths (i.e. x \in [ -w, +w ] )
      int p_widthx;
      int p_widthy;
      int p_widthz;
      ColumnVector p_kernelx;
      ColumnVector p_kernely;
      ColumnVector p_kernelz;

      // stop all forms of creation except the constructors below
      kernelstorage();
      const kernelstorage& operator=(kernelstorage&);
      kernelstorage(kernelstorage&);


    public:
      float *storex;
      float *storey;
      float *storez;

      kernelstorage(const ColumnVector& kx, const ColumnVector& ky,
		    const ColumnVector& kz, int wx, int wy, int wz)	
        {
	  p_kernelx = kx; p_kernely = ky; p_kernelz = kz;
	  p_widthx = wx;  p_widthy = wy;  p_widthz = wz;
	  storez = new float[2*wz+1];
	  storey = new float[2*wy+1];
	  storex = new float[2*wx+1];
 	}

      ~kernelstorage()
      { 
	delete storex;
	delete storey;
	delete storez;
      }
      
      class comparer
	{
	public:
	  bool operator()(const kernelstorage* k1, 
			  const kernelstorage* k2) const
	    {
	      // comparison of sizes and values (toleranced)
	      if ( (k1->p_widthx!=k2->p_widthx) || 
		   (k1->p_widthy!=k2->p_widthy) || 
		   (k1->p_widthz!=k2->p_widthz) )
		return false;
	      if ( ( (k1->p_kernelx - k2->p_kernelx).MaximumAbsoluteValue()
		     > 1e-8 * k1->p_kernelx.MaximumAbsoluteValue() ) ||
		   ( (k1->p_kernely - k2->p_kernely).MaximumAbsoluteValue() 
		     > 1e-8 * k1->p_kernely.MaximumAbsoluteValue() ) ||
		   ( (k1->p_kernelz - k2->p_kernelz).MaximumAbsoluteValue() 
		     > 1e-8 * k1->p_kernelz.MaximumAbsoluteValue() ) )
		return false;
	      return true;
	    }
	};

      friend class comparer;

      int widthx() const { return p_widthx; }
      int widthy() const { return p_widthy; }
      int widthz() const { return p_widthz; }
      const ColumnVector& kernelx() const { return p_kernelx; }
      const ColumnVector& kernely() const { return p_kernely; }
      const ColumnVector& kernelz() const { return p_kernelz; }

    };


  /////////////////////////////////////////////////////////////////////////////

  class kernel
    {
    private:
      static set<kernelstorage*, kernelstorage::comparer> existingkernels;
      kernelstorage* storedkernel;

    public:
      kernel() { storedkernel = 0; }

      const kernel& operator=(const kernel& source)
      {
	// am allowed to copy pointers if other class either
	//  always exists or manages reference counts and self-deletes
	this->existingkernels = source.existingkernels;
	this->storedkernel = source.storedkernel;
	// signal storedkernel has an extra reference
	//   and that old storedkernel has one less reference
	return *this;
      }

      kernel(const kernel& source)
      {
	this->operator=(source);
      }

      virtual ~kernel() 
      { 
	// signal storedkernel it has one less reference
      }
      
      
      void setkernel(const ColumnVector& kx, const ColumnVector& ky,
		     const ColumnVector& kz, int wx, int wy, int wz)
      {		  
	// see if already in list:
	storedkernel = new kernelstorage(kx,ky,kz,wx,wy,wz);
	set<kernelstorage*, kernelstorage::comparer>::iterator 
	  it = existingkernels.find(storedkernel);
	if (it==existingkernels.end()) {		  
	  existingkernels.insert(storedkernel);
	  // signal that this is the first reference for storedkernel
	} else {
	  delete storedkernel;
	  storedkernel = *it;
	  // signal that *it has another reference now
	}
      }

      const kernelstorage* kernelvals() { return storedkernel; }
      
  };


  /////////////////////////////////////////////////////////////////////////

  //////// Support functions /////////
  
  float kernelval(float x, int w, const ColumnVector& kernel);
  float sincfn(float x);
  float hanning(float x, int w);
  float blackman(float x, int w);
  float rectangular(float x, int w);
  ColumnVector sinckernel1D(const string& sincwindowtype, int w, int n);
  kernel sinckernel(const string& sincwindowtype, int w, int nstore);
  kernel sinckernel(const string& sincwindowtype, 
		    int wx, int wy, int wz, int nstore);
  float extrapolate_1d(const ColumnVector& data, const int index);
  float interpolate_1d(const ColumnVector& data, const float index);
  float kernelinterpolation_1d(const ColumnVector& data, float index, const ColumnVector& userkernel, int width);
  float kernelinterpolation_1d(const ColumnVector& data, float index);
  float kernelinterpolation_1d(RowVector data, float index);
  float hermiteinterpolation_1d(const ColumnVector& data, int p1, int p4, float t);
}

#endif

