/*  intensity_normalise.cc

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
  cout << "intensity_normalise  <mesh_source> <mesh_target> <source data>  <target data> <output> <-option>   " << endl;
  cout << "-option " << endl;
  cout << " -exclin exclude zeros for input mesh" << endl;
  cout << " -exclref exclude zeros for reference mesh" << endl;

  cout << " " << endl;
}


int main(int argc, char **argv){

  boost::shared_ptr<BFMatrix> DATAIN; // holds generic BFMATRIX data which can be sparse or full matrices
  boost::shared_ptr<BFMatrix> DATAREF;
  boost::shared_ptr<NEWMESH::newmesh> EXCL_IN, EXCL_REF;
int ok=0;
  newmesh ORIG, TARGET;
  bool _exclude_in=false, _exclude_ref=false, _scale=false;
  string output;
 
  if(argc < 5){

    Usage();
    exit(0);
  }

 
  ORIG.load(argv[1]);
  argc--; 
  argv++;
  TARGET.load(argv[1]);
  argc--; 
  argv++;
  set_data(argv[1],DATAIN,ORIG);
  argc--; 
  argv++;
  set_data(argv[1],DATAREF,TARGET);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;

  while (argc > 1) {
    ok = 0;
	if((ok == 0) && (strcmp(argv[1], "-exclin") == 0)){
      argc--;
      argv++;
      _exclude_in=true;
      ok = 1;
    }
	else if((ok == 0) && (strcmp(argv[1], "-exclref") == 0)){
      argc--;
      argv++;
      _exclude_ref=true;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-scale") == 0)){
      argc--;
      argv++;
      _scale=true;
      ok = 1;
    }
    else{cout << " option doesn't exist " << endl; exit(1);}
  }
  
   if(_exclude_in){
    EXCL_IN= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(ORIG,DATAIN->AsMatrix(),0,0))); 
   }
   if(_exclude_ref){
    EXCL_REF= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(TARGET,DATAREF->AsMatrix(),0,0))); 
   }
  
  multivariate_histogram_normalization(*DATAIN,*DATAREF,EXCL_IN,EXCL_REF,_scale);
  ORIG.set_pvalues(DATAIN->AsMatrix());
  ORIG.save(output);
  	
}
