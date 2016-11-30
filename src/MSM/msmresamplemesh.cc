/*  msmresample.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "newmat.h"
#include "newmesh/meshfns.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



using namespace NEWMAT;
using namespace NEWMESH;

void Usage()
{ cout << " msmresamplemesh <in_anat> <in_sphere> <ico resolution> <output base> -options " << endl;
  cout << " -data X supply data for resampling " << endl;
  cout << " -adap resample adaptively " << endl;

}


int main(int argc, char **argv){

  
  newmesh in_anat,in_sphere,ico;
  newmesh ANAT_res;
  string output;
  double res;
  char filename[1000];
  resampler R; R.set_method("ADAP_BARY");
  int ok;
  bool _resampledata=false;
  bool _adap=false;
  if(argc < 3){

    Usage();
    exit(0);
  }

 
  in_anat.load(argv[1]);
  argc--; 
  argv++;
  in_sphere.load(argv[1]);
  argc--; 
  argv++;
  res=atoi(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-data") == 0)){
      argc--;
      argv++;    
      _resampledata=true;
      in_sphere.load(argv[1],false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-adap") == 0)){
      argc--;
      argv++;    
      _adap=true;
      ok = 1;
    }else{cout << " option doesn't exist " << endl; exit(1);}
  }
  cout << " Make ico " << endl;
  ico.make_mesh_from_icosa(res); true_rescale(ico,RAD); 
  
 
  Matrix TRANSLATE=recentre(in_sphere);
  true_rescale(in_sphere,RAD);
  //recentre anat ///
  Pt mean;
  /* if(TRANSLATE.NormFrobenius() > EPSILON){
  
 // ColumnVector P_in(4);
  
  cout << " translate anat " << endl;
  for (int i=0;i< in_anat.nvertices();i++)
      mean+=in_anat.get_coord(i);

    mean/=in_anat.nvertices();
    cout << "inant  before mean" << mean.X << " mean.Y " << mean.Y << " mean.Z " << mean.Z <<  endl;
   
    mean*=0;
    for ( vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i= in_anat.vbegin(); i!=in_anat.vend(); i++){
       NEWMESH::Pt p = (*i)->get_coord();
       NEWMESH::Pt p2;
       if(p.norm()){
	 P_in(1) = p.X; P_in(2) = p.Y; P_in(3) = p.Z; P_in(4) = 1;  
	 P_in = TRANSLATE * P_in;
	 p2.X=P_in(1); p2.Y=P_in(2);  p2.Z=P_in(3);
	 (*i)->set_coord(p2);
       }
      
     }	
    
    for (int i=0;i< in_anat.nvertices();i++)
      mean+=in_anat.get_coord(i);

    mean/=in_anat.nvertices();
    cout << "inant  mean" << mean.X << " mean.Y " << mean.Y << " mean.Z " << mean.Z <<  endl;
  
  }
  */
  ANAT_res=mesh_resample(in_anat,in_sphere,ico,_adap);
  sprintf(filename,"%s-anat.surf",output.c_str());
  mean*=0;
  for (int i=0;i< ANAT_res.nvertices();i++)
      mean+=ANAT_res.get_coord(i);

    mean/=ANAT_res.nvertices();
    cout << "anat res  mean" << mean.X << " mean.Y " << mean.Y << " mean.Z " << mean.Z <<  endl;
  ANAT_res.save(filename);
  sprintf(filename,"%s-regular_sphere.surf",output.c_str());
  ico.save(filename);

  if(_resampledata){
    R.resample_scalar(in_sphere,ico,1);
    sprintf(filename,"%s-resampled_data.func",output.c_str());
    in_sphere.save(filename);

  }
}
