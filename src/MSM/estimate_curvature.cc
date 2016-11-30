/*  estimate_curvature.cc

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
  cout << "estimate_curvature  <anatomical mesh>  <output base> <-option>   " << endl;
  cout << " estimates mean curvature by fiting a least-squares quadratic patch to the local neighborhood of a vertex f(x,y) = ax^2 + by^2 + cxy + dx + ey + f (inspired by patchcuvature.m the matlab function " << endl;
  cout << " options: " << endl;
  cout << " -fit-radius default 2 " << endl; 
  cout << " " << endl;
}


int main(int argc, char **argv){

  
  newmesh ORIG,SPHERE;
  string outdir;
  double fit_radius=2, smoothval=2.0;
  boost::shared_ptr<RELATIONS> REL;
  bool _calcrel=false,_smooth=false;
  bool _outputrho=false;
  int ok;

  if(argc < 3){

    Usage();
    exit(0);
  }

 
  ORIG.load(argv[1]);
  argc--; 
  argv++;
  outdir=argv[1];
  argc--; 
  argv++;


  while (argc > 1) {
    ok = 0;
   
    if((ok == 0) && (strcmp(argv[1], "-fit-radius") == 0)){
       argc--; 
       argv++;
       fit_radius=atof(argv[1]);
       argc--; 
       argv++;
       ok=1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-sphere") == 0)){
       argc--; 
       argv++;
       SPHERE.load(argv[1]);
       _calcrel=true;
       argc--; 
       argv++;
       ok=1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-smooth") == 0)){
       argc--; 
       argv++;
       _smooth=true;
       smoothval=atof(argv[1]);
       argc--; 
       argv++;
       ok=1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-rho") == 0)){
       argc--; 
       argv++;
       _outputrho=true;
       ok=1;
    } else{cout << " option doesn't exist " << endl; exit(1);}
  }

 if(_calcrel){
    REL=boost::shared_ptr<RELATIONS >(new RELATIONS (SPHERE,SPHERE,3*asin(fit_radius/RAD))); 
    REL->update_RELATIONS(SPHERE);
  }
  cout << " estimate curvature " << endl;
  mean_curvature(fit_radius,ORIG,REL);

  if(_smooth){
    if(!_calcrel) {cout << " need sphere for smoothing currently " << endl; exit(1);}
    resampler R;
    boost::shared_ptr<BFMatrix> scalardata;
    boost::shared_ptr<newmesh> EXCL;
    scalardata =boost::shared_ptr<BFMatrix> (new FullBFMatrix (1,SPHERE.nvertices()));
    R.set_method("GAUSSIAN");

    for (int i=0;i<ORIG.nvertices();i++)
      scalardata->Set(1,i+1,ORIG.get_pvalue(i));

    R.resampledata(SPHERE,SPHERE,EXCL,scalardata,smoothval,REL)  ;

    for (int i=0;i<ORIG.nvertices();i++)
      ORIG.set_pvalue(i,scalardata->Peek(1,i+1));

  }
  char filename[1000];
  sprintf(filename,"%s-meancurvature.func",outdir.c_str());
  ORIG.save(filename);
  if(_outputrho){

    for (int i=0;i<ORIG.nvertices();i++)
      ORIG.set_pvalue(i,1/ORIG.get_pvalue(i));

    sprintf(filename,"%s-rho.func",outdir.c_str());
    ORIG.save(filename);
  }

 
}
