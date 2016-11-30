/*  estimate_metric_distortion.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/* this program is designed to downsample freesurfer label files to be used in combination with the SPH6.vtk or other downsampled meshes*/

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include "newmat.h"
#include "newmatio.h"
#include "utils/options.h"
#include "newimage/newimageall.h"


#include "MeshReg/meshreg.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "estimate_metric_distortion  <original mesh> <transformed mesh > < output> <-option>   " << endl;
  cout << " Compares triangle areas before and after registration - can be used to assess strength of regularisation" << endl;
  cout << " options: " << endl;
  cout << " -target project result onto a target image (requires argument)" << endl;
  cout << " -abs take absolute value " << endl;
  
}


int main(int argc, char **argv){

  
  newmesh ORIG, TRANSFORMED,TARG,DISTMAP;
  resampler R;
  string outdir;
  boost::shared_ptr<NEWMESH::newmesh> EXCL;
  int ok;
  bool _targ=false,_abs=false;
  if(argc < 4){

    Usage();
    exit(0);
  }

 
  ORIG.load(argv[1]);
  argc--; 
  argv++;
  TRANSFORMED.load(argv[1]);
  argc--; 
  argv++;
  outdir=argv[1];
  argc--; 
  argv++;

  TARG=TRANSFORMED;

  while (argc > 1) {
    ok = 0;
   
    if((ok == 0) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      _targ=true;
      TARG.load(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }if((ok == 0) && (strcmp(argv[1], "-abs") == 0)){
      argc--;
      argv++;
      _abs=true;
      ok = 1;
    }
    
    else{cout << " option doesn't exist " << endl; exit(1);}
  }


 
  DISTMAP=TRANSFORMED;
  ColumnVector distortion=getarealdistortion(ORIG,TRANSFORMED);
  double meandist=0;
  double val=0;
  for(int i=0;i<DISTMAP.nvertices();i++){
    if(_abs==true)val=abs(distortion(i+1));
    else  val=distortion(i+1);
    DISTMAP.set_pvalue(i,val);
    meandist+=val;
  }
  meandist/=DISTMAP.nvertices();
  cout << " mean=" << meandist << endl;
  DISTMAP.save(outdir + "distortion.func");

  if(_targ){ double MVD=Calculate_MVD(TRANSFORMED);  
  R.resample_scalar(DISTMAP,TARG,asin(MVD/(2*RAD)),EXCL);
  DISTMAP.save(outdir + "distortion_on_target.func");
  }
}
