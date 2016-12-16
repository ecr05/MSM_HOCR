/*  base2z.h

    Mark Woolrich & Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(__base2z_h)
#define __base2z_h

#include <iostream>
#include <fstream>

namespace MISCMATHS {
  
  class Base2z
    {
    public:
      Base2z()
	{}
      
      virtual ~Base2z() { delete base2z; }      
      
      float convertlogp2z(float logp); 
      float logp2largez(float logp);
      float logbeta(float v, float w);

      virtual bool issmalllogp(float logp) = 0;

    private:
      
      const Base2z& operator=(Base2z&);
      Base2z(Base2z&);
      
      static Base2z* base2z;

    };

}

#endif

