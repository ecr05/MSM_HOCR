/*  kertest.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

/*  CCOPYRIGHT  */

#include "kernel.h"

using namespace MISCMATHS;
using namespace NEWMAT;

int main (int argc,char** argv)
{
  cerr << "Test program for 1D kernel interpolation" << endl;

  ColumnVector data(10);
  data << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1;
  ColumnVector newVec = data;
  
  cerr << "Input: " << data << endl;

  for (int index = 0; index <= 10; index++)
    newVec[index] = kernelinterpolation_1d(data, index+0.5);

  cerr << "Result: " << newVec << endl;
  
  return 0;
}
