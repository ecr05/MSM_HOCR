/*  histogram.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 4.0 (c) 2007, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/1112. */

#if !defined(_histogram_h)
#define _histogram_h

#include <iostream>
#include <fstream>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"

using namespace NEWMAT;

namespace MISCMATHS {
 
  class Histogram
    {
    public:
      Histogram(){};
      const Histogram& operator=(const Histogram& in){
	sourceData=in.sourceData; calcRange=in.calcRange; histMin=in.histMin; histMax=in.histMax; bins=in.bins; histogram=in.histogram; CDF=in.CDF; datapoints=in.datapoints; exclusion=in.exclusion;
	return *this;
      }

      Histogram(const Histogram& in){ *this=in;}

      Histogram(const ColumnVector& psourceData, int numBins)
	: sourceData(psourceData), calcRange(true), bins(numBins){}

      Histogram(const ColumnVector& psourceData, float phistMin, float phistMax, int numBins) 
	: sourceData(psourceData), calcRange(false), histMin(phistMin), histMax(phistMax), bins(numBins){}
      
      void set(const ColumnVector& psourceData, int numBins) {	
	sourceData=psourceData; calcRange=true; bins=numBins;
      }

      void set(const ColumnVector& psourceData, float phistMin, float phistMax, int numBins) {	
	sourceData=psourceData; calcRange=false; histMin=phistMin; histMax=phistMax; bins=numBins;
      }

      void generate();
      void generate(ColumnVector exclusion_values);
      void generateCDF();

      float getHistMin() const {return histMin;}
      float getHistMax() const {return histMax;}
      void setHistMax(float phistMax) {histMax = phistMax;}
      void setHistMin(float phistMin) {histMin = phistMin;}
      void setexclusion(ColumnVector exclusion_values) {exclusion =exclusion_values;}
      void smooth();

      int integrateAll() {return sourceData.Nrows();}

      const ColumnVector& getData() {return histogram;}
      void setData(const ColumnVector& phist) { histogram=phist;}

      int integrateToInf(float value) const { return integrate(value, histMax); }
      int integrateFromInf(float value) const { return integrate(histMin, value); }
      int integrate(float value1, float value2) const;
    
      void match(Histogram &);
      
      float mode() const;

      int getBin(float value) const;
      float getValue(int bin) const;
      float getPercentile(float perc);

      inline int getNumBins() const {return bins;}
      inline ColumnVector getCDF() const {return CDF;}
      inline ColumnVector getsourceData()const {return sourceData;}
    protected:

    private:

      ColumnVector sourceData;
      ColumnVector histogram;
      ColumnVector exclusion;
      ColumnVector CDF;

      bool calcRange;

      float histMin;
      float histMax;

      int bins; // number of bins in histogram
      int datapoints;
    };

  inline int Histogram::getBin(float value) const
    {
      float binwidth=(histMax-histMin)/bins;
      return Max(1, Min((int)((((float)bins)*((float)(value-(histMin-binwidth))))/((float)(histMax-histMin))),bins));
    }
  
  inline float Histogram::getValue(int bin) const
    {
      return (bin*(histMax-histMin))/(float)bins + histMin;
    }

  inline void Histogram::generateCDF() 
    { 
      
      CDF.ReSize(bins);
      
      
      CDF(1)=histogram(1)/datapoints;
  
      for (int i=2;i<=bins;i++)
	CDF(i)=CDF(i-1)+ histogram(i)/datapoints;

     
     
    }
}

#endif





