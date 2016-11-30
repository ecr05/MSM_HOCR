/*  surfconvert.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include <time.h>
#include <stdio.h>
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"


using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "surfconvert  < surfaces> <output path> <data> " << endl;
  cout << "options:  " << endl;
  cout << "--data add data (requires argument)  " << endl;
  cout << "--rescale (requires argument)  " << endl;

}

int main( int argc, char **argv )
{
  newmesh SURFACE;
  int ok=0;
  double rescale=1.0;
  string outname;
  
  if(argc < 2){

    Usage();
    exit(0);
  }

  cout << "argv[1" << argv[1] << endl;
  SURFACE.load(argv[1]);
  argc--; 
  argv++;
  cout << "argv[1" << argv[1] << endl;
  outname=argv[1];
  argc--; 
  argv++;
 
  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "--data") == 0)){
      argc--;
      argv++;
      SURFACE.load(argv[1],false,true);
      argc--; 
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "--centre") == 0)){
      argc--;
      argv++;
      recentre(SURFACE);
     
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "--rescale") == 0)){
      argc--;
      argv++;
      rescale=atof(argv[1]);
      argc--; 
      argv++;
      ok = 1;
    } else{      cout << argv[1] << " option doesn't exist " << endl; exit(1);}
  }

  if(rescale>1.0){
    for(int i=0;i<SURFACE.nvertices();i++){
      Pt p=SURFACE.get_coord(i);
      p.normalize(); p=p*rescale;
      SURFACE.set_coord(i,p);
    }
  }
  SURFACE.save(outname);
}
