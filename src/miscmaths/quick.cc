/*  quick.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#define WANT_STREAM
#define WANT_MATH

#include "miscmaths.h"
#include "t2z.h"
//#include "libprob.h"

using namespace MISCMATHS;

int main(int argc, char *argv[])
{
  try{

    Matrix X = read_vest("/usr/people/woolrich/matlab/vbbabe/data/design2.mat").t();
 
    ColumnVector Y = read_vest("/usr/people/woolrich/matlab/vbbabe/data/sdf2.mat").t();
 
    ColumnVector m_B;
    SymmetricMatrix ilambda_B;

    glm_vb(X, Y, m_B, ilambda_B, 30);

    write_ascii_matrix(m_B,"/usr/people/woolrich/matlab/vbbabe/data/m_B");
    
  }
  catch(Exception p_excp) 
    {
      cerr << p_excp.what() << endl;
    }
  return 0;
}

