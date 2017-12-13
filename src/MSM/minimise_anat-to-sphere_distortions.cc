/*  minimise_anat-to-spahre_distortions.cc

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
#include <string>

#include "MeshReg/meshreg.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "minimise_anat-to-sphere_distortion  <anat> <sphere> < output> <-option>   " << endl;
  cout << " Relaxes sphere vertex locations to minimise distorions" << endl;
  cout << "options" << endl;
  cout << "-levels " << endl;
  cout << "-resolutions_per_level (argument is a comma separated list) " << endl;
  cout << "-iters  (argument is a comma separated list, needs to have as many items as resolutions) " << endl;
  cout << "-setexp " << endl;
  cout << "-bulk bulk modulus parameter (default 10)" << endl;
  cout << "-shear shear modulus parameter (default 0.1)" << endl;

}


vector<int> read_stringstream(string params){
  vector<int> list;
  int i;
  std::stringstream ss(params);
  while (ss >> i)
    {
      list.push_back(i);

      if (ss.peek() == ',')
	ss.ignore();
    }
  return list;
}

int main(int argc, char **argv){


	  newmesh ANAT, SPHERE;
	  newmesh ANATres, SPHEREres,ico;
	  newmesh controlgrid, transformed_controlgrid;
	  double energy=0,newenergy=0;
	  double mu=0.1,kappa=10,rexp=1.0;

	  boost::shared_ptr<SRegDiscreteModel> MODEL;

	  string outdir;
	  int ok,resolutionlevels=1,iter;
	  string iters;
	  string resolutions;
	  double area_anat,area_sphere,maxnorm=0;
	  std::vector<int> anat_res(1,5),iters_res(1,5);
	  bool _debug=true;

	  if(argc < 4){

		Usage();
		exit(0);
	  }


	  ANAT.load(argv[1]);
	  argc--;
	  argv++;
	  SPHERE.load(argv[1]);
	  argc--;
	  argv++;
	  outdir=argv[1];
	  argc--;
	  argv++;


	  while (argc > 1) {
		ok = 0;
		cout << argv[1] << endl;

		if((ok == 0) && (strcmp(argv[1], "-levels") == 0)){
		   argc--;
		   argv++;
		   resolutionlevels=atoi(argv[1]);
		   cout << " levels " << resolutionlevels << endl;

		   argc--;
		   argv++;
		   ok=1;
		}
		else if((ok == 0) && (strcmp(argv[1], "-iters") == 0)){
		   argc--;
		   argv++;
		   iters=argv[1];
		   iters_res=read_stringstream(iters);
		   cout << " iterations " << iters_res[0] << endl;
		   argc--;
		   argv++;
		   ok=1;
		}
		else if((ok == 0) && (strcmp(argv[1], "-setexp") == 0)){
		   argc--;
		   argv++;
		   rexp=atof(argv[1]);
		   cout << " exp " << rexp << endl;
		   argc--;
		   argv++;
		   ok=1;
		}
		else if((ok == 0) && (strcmp(argv[1], "-bulk") == 0)){
		   argc--;
		   argv++;
		   kappa=atof(argv[1]);
		   cout << " kappa " << kappa << endl;
		   argc--;
		   argv++;
		   ok=1;
		}
		else if((ok == 0) && (strcmp(argv[1], "-shear") == 0)){
		   argc--;
		   argv++;
		   mu=atof(argv[1]);
		   cout << " mu " << mu << endl;
		   argc--;
		   argv++;
		   ok=1;
		}
		else if((ok == 0) && (strcmp(argv[1], "-resolutions_per_level") == 0)){
		   argc--;
		   argv++;
		   resolutions=argv[1];
		   anat_res=read_stringstream(resolutions);
		   argc--;
		   argv++;
		   ok=1;
		}
		else{cout << " option doesn't exist " << endl; exit(1);}
	  }



	  char filename[1000];
	  newmesh STRAINS,SPHEREnew;
	  double totstrains=0,tristrain=0;
	  Reduction HOCR_mode=HOCR;
	  true_rescale(SPHERE,RAD);
	  SPHEREnew=SPHERE;
	  SPHEREnew=calculate_triangular_strains(ANAT,SPHEREnew,mu,kappa);
	  sprintf(filename,"%sSTRAINShiINIT.func",outdir.c_str());

	  SPHEREnew.save(filename);
	  for(int it=0;it<resolutionlevels;it++){
		cout << " Initialising level " << it+1 << " " << anat_res[it] << endl;

		ico.make_mesh_from_icosa(anat_res[it]); true_rescale(ico,RAD);
		SPHEREres=ico;
		ANATres=mesh_resample(ANAT,SPHERE,SPHEREres);
		STRAINS=ANATres;

		if(it==(resolutionlevels-1)){sprintf(filename,"%sSPHEREINIT.LR.surf",outdir.c_str());SPHEREres.save(filename);
			sprintf(filename,"%sANAT.LR.surf",outdir.c_str());ANATres.save(filename);}

		totstrains=0;
		for(int i=0;i<STRAINS.ntriangles();i++){
		  tristrain=calculate_triangular_strain(i,ANATres,SPHEREres,mu,kappa);
		  totstrains+=tristrain;
		  //cout << i << " STRAINS.get_pvalue(i) " << tristrain <<endl;
		}

		MODEL=boost::shared_ptr<RegularisationDiscreteModel>(new MetricDistortionDiscreteModel(anat_res[it]+2,mu,kappa,rexp));
		MODEL->reset_meshspace(ANATres);
		cout << it << " init strains " << totstrains<< endl;

		MODEL->Initialize(SPHEREres);
		if(it>0){  MODEL->warp_CPgrid(controlgrid,transformed_controlgrid) ; SPHEREres=MODEL->get_CPgrid(); }

		if(_debug) MODEL->set_debug();
		iter=0;
		while(iter < iters_res[it]){
		  cout << " iter " << iter << " iters_res[it] " << endl;
		  MODEL->setupCostFunction();
		  int *Labels=MODEL->getLabeling();
		  newenergy=Fusion::optimize(MODEL,HOCR_mode,_debug);

		  MODEL->applyLabeling();
		  if(iter>0 && iter%2==0 && newenergy>=energy) { cout << iter << "energy" << energy << " newenergy " << newenergy << " level has converged" << (iter-1) %2  << endl; break;}
		  for (int i = 0; i <SPHEREres.nvertices(); i++){
		  cout << i << " _iter " << iter << " label " << Labels[i] <<  endl;
		  }
		  //sprintf(filename,"SPHEREres-%d-iter-%d.surf",it,iter);
		  transformed_controlgrid=MODEL->get_CPgrid();
		  //   unfold(transformed_controlgrid); --- is unfold going to be necessary???
		  //MODEL->reset_CPgrid(transformed_controlgrid);
		  //transformed_controlgrid.save(filename);

		  totstrains=0;
		  for(int i=0;i<STRAINS.ntriangles();i++){
			tristrain=calculate_triangular_strain(i,ANATres,transformed_controlgrid,mu,kappa);
			totstrains+=tristrain;
		//	cout << i << " STRAINS.get_pvalue(i) " << tristrain <<endl;
		  }
		  cout << it <<" " << iter <<  " current strains " << totstrains<< endl;

		  energy=newenergy;
		  iter++;
		}

		transformed_controlgrid=MODEL->get_CPgrid();
		controlgrid=ico;

	  }
	  sprintf(filename,"%sSPHERE.LR.surf",outdir.c_str());
	  transformed_controlgrid.save(filename);
	  transformed_controlgrid=calculate_triangular_strains(ANATres,transformed_controlgrid,mu,kappa);
	  sprintf(filename,"%sSTRAINS.LR.func",outdir.c_str());
	  transformed_controlgrid.save(filename);
	  barycentric_mesh_interpolation(SPHEREnew,controlgrid,transformed_controlgrid);
	  sprintf(filename,"%sSPHERE.surf",outdir.c_str());
	  SPHEREnew.save(filename);
	  SPHEREnew=calculate_triangular_strains(ANAT,SPHEREnew,mu,kappa);
	  sprintf(filename,"%sSTRAINShi.func",outdir.c_str());

	  SPHEREnew.save(filename);
	  SPHEREnew=SPHERE;

}
