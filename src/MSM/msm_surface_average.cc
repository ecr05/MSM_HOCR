/*  msm_surface_average.cc

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
#include "utils/options.h"
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"


using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "msm_surface_average  <asciii list of data> <output> <-options>   " << endl;
  cout << "options:  " << endl;
  cout << "-weighted X supply a list of distances from the mean" << endl;
  cout << "-sigma X  sigma values " << endl;
}


int main(int argc, char **argv){

 
  NEWMESH::newmesh IN, AVERAGE;
  
  string output;
  vector<string> INlist,distances;
  bool _weighted=false;
  double _sigma;
  double weight;
  double sumweighted=0;
  int ok=0;
 

  if(argc < 2){

    Usage();
    exit(0);
  }

  INlist=read_ascii_list(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  
  

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-weighted") == 0)){
      argc--;
      argv++;
      _weighted=true;
      distances=read_ascii_list(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;
      argv++;
      _sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else{cout << argv[1] << " option doesn't exist " << endl; exit(1);}
  }

  AVERAGE.load(INlist[0]);  

  for(int i=0;i<AVERAGE.nvertices();i++){
    Pt p;
    AVERAGE.set_coord(i,p);
  }
  double mean=0;
  //for (unsigned int i=0;i<INlist.size();i++){
	//mean+=atof(distances[i].c_str());
//}
///mean/=INlist.size();
  /////////////////////////////////////////////////////
  for (unsigned int i=0;i<INlist.size();i++){
    cout << i << " " << INlist[i] << endl;

    IN.load(INlist[i]);
  //  cout << " recentre " << endl;
   // recentre(IN);
    if(_weighted){
      if(_sigma==0){ cout << "must set sigma for weighted" <<endl; exit(1);}
      weight=1/(_sigma*(sqrt(2*PI)))*exp(-0.5*(((atof(distances[i].c_str())-mean))/_sigma)*(((atof(distances[i].c_str())-mean))/_sigma));
      sumweighted+=weight;
      cout << "distance " <<  distances[i] << " weight " << weight << endl;
    }else{ weight=1;sumweighted+=1;}

    
    for( int n=0;n<IN.nvertices();n++){
      Pt p=AVERAGE.get_coord(n);
      Pt p2=IN.get_coord(n);
    
      Pt newpt;
      newpt.X=p.X+weight*p2.X;      
      newpt.Y=p.Y+weight*p2.Y;
      newpt.Z=p.Z+weight*p2.Z;
      AVERAGE.set_coord(n,newpt);
    }

  }

  for( int n=0;n<AVERAGE.nvertices();n++){
      Pt p=AVERAGE.get_coord(n);
      p.X=p.X/sumweighted;
      p.Y=p.Y/sumweighted;
      p.Z=p.Z/sumweighted;
      AVERAGE.set_coord(n,p);
  }

  //recentre(AVERAGE);
  AVERAGE.save(output+"AVERAGE.surf");
 
}
