/*  msm_makeicos.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "newmat.h"
#include "newmesh/meshfns.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



using namespace NEWMAT;
using namespace NEWMESH;

void Usage()
{ cout << " msm_makeicos <ico resolution> " << endl;
 
}


int main(int argc, char **argv){

  
  newmesh ico;
  int res;
  char filename[1000];

  if(argc < 1){

    Usage();
    exit(0);
  }

 
  
  res=atoi(argv[1]);
  argc--; 
  argv++;
  
 
  cout << " Make ico " << endl;
  ico.make_mesh_from_icosa(res); true_rescale(ico,RAD); 
  sprintf(filename,"ico-%d.surf.gii",res);
  ico.save(filename);
}
