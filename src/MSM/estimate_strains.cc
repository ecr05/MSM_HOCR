/*  estimate_strains.cc

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
  cout << "estimate_strains  <original anatomical mesh>  <transformed mesh > <outdir> <-option>   " << endl;
  cout << " estimates principal strains tangent to surface" << endl;
  cout << " options: " << endl;
  cout << " -fit-radius default 2 " << endl; 
  cout << " -sphere X.surf.gii supply input sphere for approximate neighbourhood searching (will speed up calculation) " << endl; 
   cout <<"-bulk bulk modulus (default 10) " << endl; 
  cout << "-shear shear modulus (default 10) " << endl; 

  cout << " " << endl;
}


int main(int argc, char **argv){

  clock_t start = clock(); 
   
  newmesh ORIG, TRANSFORMED, STRAINS,SPHERE;
  boost::shared_ptr<Matrix>  VECS; VECS= boost::shared_ptr<Matrix> (new Matrix);
  string outdir;
  double fit_radius=0;
  double MU=0.1,KAPPA=10;
  boost::shared_ptr<RELATIONS> REL;
  boost::shared_ptr<newmesh > SPHERETMP;
  int ok;
  bool _targ=false,_calcrel=false;
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
     else if((ok == 0) && (strcmp(argv[1], "-bulk") == 0)){
       argc--; 
       argv++;
       KAPPA=atof(argv[1]);
       argc--; 
       argv++;
       ok=1;
   }
     else if((ok == 0) && (strcmp(argv[1], "-shear") == 0)){
       argc--; 
       argv++;
       MU=atof(argv[1]);
       argc--; 
       argv++;
       ok=1;
    } else{cout << " option doesn't exist " << endl; exit(1);}
  }

  /* if(_calcrel){
    REL=boost::shared_ptr<RELATIONS >(new RELATIONS (SPHERE,SPHERE,3*asin(fit_radius/RAD))); 
    REL->update_RELATIONS(SPHERE);
  }
  */
  cout << " calculate strains " << MU << " " << KAPPA <<  endl;
  ORIG.estimate_normals(); TRANSFORMED.estimate_normals();
  if(SPHERE.nvertices()==ORIG.nvertices()) {
    SPHERETMP=boost::shared_ptr<newmesh >(new newmesh (SPHERE));  
  }
  
  if(fit_radius>0)
    STRAINS=calculate_strains(fit_radius,ORIG,TRANSFORMED,VECS,REL);
  else
    STRAINS=calculate_triangular_strains(ORIG,TRANSFORMED,MU,KAPPA);

  
  STRAINS.save(outdir+"strains.func");

  /*  ofstream strain1,strain2;

  strain1.open(outdir+"U1");
  strain2.open(outdir+"U2");

  cout << " output U1 and U2 " << endl; 
  if(VECS->Nrows()==ORIG.nvertices()){
    for(int i=1;i<=ORIG.nvertices();i++){
      for(int j=1;j<=3;j++)
	strain1 << (*VECS)(i,j) << " " ;
      strain1 << endl;

      for(int j=4;j<=6;j++)
	strain2 << (*VECS)(i,j) << " " ;
      strain2 << endl;
  }
  }

  strain1.close(); strain2.close();
  cout<<"Time elapsed: " << ((double)clock() - start) / CLOCKS_PER_SEC << endl;
  */
}
