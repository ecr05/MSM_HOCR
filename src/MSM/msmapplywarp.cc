/*  msmapplywarp.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/* this program is designed to downsample freesurfer label files to be used in combination with the SPH6.vtk or other downsampled meshes*/

#include "newmat.h"
#include "newmesh/meshfns.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



using namespace NEWMAT;
using namespace NEWMESH;

void Usage()
{ cout << " msmapplywarp <to-be-transformedmesh> <outputname> -options  " << endl;
  cout << " Projects the to-be-transformed mesh though a transformation defined by an original mesh and its deformed counterpart" << endl;
  cout << " It is optional to supply the undeformed (original) mesh, but if no original mesh is supplied the algorithm will assume that the warp is prescribed by an icospheric template" << endl;
  cout << " -options " << endl;
  cout << " -original X.surf.gii " << endl;  
  cout << " -deformed DEFORMED.surf.gii (MUST be supplied in order to warp data)" << endl;
  cout << " -anat TARGET_SPHERE.surf.gii TARGET_ANAT.surf (2 inputs!). This will effectively project the INPUT anatomical mesh through the spherical warp." << endl;   
  cout << " -nospherical don't save spherical warp " << endl;  
  cout << " -affine (estimate affine transformation between input (-original) and deformed (-deformed) meshes and apply this to the to-be-transformed-mesh " << endl;  
  cout << " -readaffine X; where X  is an affine transformation matrix " << endl;  
  cout << " -writeaffine write out affine transformation matrix" << endl;  


}

void get_areas(newmesh &M){
  double val;
  ColumnVector meanarea(M.nvertices()); meanarea=0;

  for ( vector<NEWMESH::Triangle>::const_iterator i=M.tbegin() ; i!=M.tend(); i++){
    NEWMESH::Pt v1 = (*i).get_vertex_coord(0),  v2 = (*i).get_vertex_coord(1), v3 = (*i).get_vertex_coord(2);
	  
    val=computeArea(v1, v2, v3);
    meanarea((*i).get_vertex_no(0)+1)+=(1/3)*val; 
    meanarea((*i).get_vertex_no(1)+1)+=(1/3)*val; 
    meanarea((*i).get_vertex_no(2)+1)+=(1/3)*val; 
	
  }
  for (int i=0;i<M.nvertices();i++){
    M.set_pvalue(i,meanarea(i+1));
  }


}
int main(int argc, char **argv){

  
  newmesh to_be_deformed,to_be_deformed2, initial,final;
  newmesh TARGETSPHERE, TARGETANAT;
  newmesh ANAT_TRANS;
  Matrix _affinemat;
  string output;
  
  boost::shared_ptr<RELATIONS > _rel;
  bool _isregular=true, _projectanat=false;
  bool _deformed=false, _areas=false;
  bool _affine=false, _supplyaffine=false, _writeaffine=false;
  char filename[1000];
  int ok;

  if(argc < 3){

    Usage();
    exit(0);
  }

 
  to_be_deformed.load(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  recentre(to_be_deformed);
  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-original") == 0)){
      argc--;
      argv++;    
      initial.load(argv[1]);
      recentre(initial);
      _isregular=false;
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-deformed") == 0)){
      argc--;
      argv++;    
      final.load(argv[1]);
      recentre(final);
      _deformed=true;
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-anat") == 0)){
      argc--;
      argv++;   
      TARGETSPHERE.load(argv[1]); 
      recentre(TARGETSPHERE); 
      argc--;
      argv++; 
      TARGETANAT.load(argv[1]);
      argc--;
      argv++; 
      _projectanat=true;
     
      ok = 1;
    } 
    else if((ok == 0) && (strcmp(argv[1], "-affine") == 0)){
      argc--;
      argv++;    
      _affine=true;   
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-outputareas") == 0)){
      argc--;
      argv++;    
      _areas=true;   
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-readaffine") == 0)){
      argc--;
      argv++;
      _supplyaffine=true;
      _affinemat=read_ascii_matrix(argv[1]);    
      argc--;
      argv++;
      ok = 1;
    }
     else if((ok == 0) && (strcmp(argv[1], "-writeaffine") == 0)){
      argc--;
      argv++;
      _writeaffine=true;
      ok = 1;
    }
    else{ cout <<argv[1] << "option npot present " << endl; exit(1);} 
  }
  
  if(!_deformed && !_supplyaffine) {
    if(!_projectanat){ cout <<_supplyaffine <<  "Cannot resample without a target and anatomical mesh. Please supply a target and anatomical mesh (-anat option) for resampling or a deformed mesh (-deformed option) for warping" << endl; exit(1);}
    ANAT_TRANS=projectmesh(to_be_deformed,TARGETSPHERE,TARGETANAT);
    sprintf(filename,"%s_anatresampled.surf.gii",output.c_str());
    ANAT_TRANS.save(filename);
    if(_areas){
      get_areas(ANAT_TRANS);
	
	
      sprintf(filename,"%s_anatresampled_areas.func.gii",output.c_str());
      ANAT_TRANS.save(filename);

    }

  }
  else{
 
	if(_supplyaffine){
		cout << _affinemat << endl;
		affine_transform(to_be_deformed,_affinemat);
        
	}else{
	  newmesh ICO;

	  if(_isregular){
	    int ico=final.get_ico_resolution();
	    
	    
	    ICO.make_mesh_from_icosa(ico); true_rescale(ICO,RAD);
 
	  
	  }else ICO=initial;

    if(_affine){
      cout << " affine transformation" << endl;
      Matrix affinetrans=affine_transform(to_be_deformed,ICO,final);
      if(_writeaffine){
	sprintf(filename,"%s_affinewarp.txt",output.c_str());
	write_ascii_matrix(filename,affinetrans);
      }	
    }
    else
      barycentric_mesh_interpolation(to_be_deformed,ICO,final);
    
   
    if(_projectanat){
      ANAT_TRANS=projectmesh(to_be_deformed,TARGETSPHERE,TARGETANAT);
      sprintf(filename,"%s_projectedwarp.surf.gii",output.c_str());
      ANAT_TRANS.save(filename);
      
      if(_areas){
	get_areas(ANAT_TRANS);
	
	
	sprintf(filename,"%s_projectedwarp_areas.func.gii",output.c_str());
	ANAT_TRANS.save(filename);

      }
      
    }
	
  }
	sprintf(filename,"%s_warp.surf.gii",output.c_str());
	to_be_deformed.save(filename);
  }
}
