/*  groupmeshreg.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "groupmeshreg.h"

#include <time.h>

namespace MESHREG {

 void GroupMeshReg::Initialize_level(int level)
  {
    cout << " In Initialize level " <<  DATAlist.size() <<endl;
    MeshModify::Initialize();
    vector<double>  sigma(MESHES.size(),_sigma_in[level]);
    NEWMESH::newmesh CONTROL; // original low res icosphere mesh
    CONTROL.make_mesh_from_icosa(_gridres[level] ); true_rescale(CONTROL,RAD);
    FEAT=boost::shared_ptr<featurespace>(new featurespace(DATAlist));
    FEAT->set_smoothing_parameters(sigma);

    FEAT->set_cutthreshold(_threshold); // will also generate exclusion masks at the same mesh resolution as datagrid
    FEAT->logtransform(_logtransform);// if true logtransforms AND normalises
    FEAT->varnorm(_varnorm);// variance normalises
    FEAT->intensitynormalize(_IN, _scale); // matches the intensities of the source to to the target (will rescale all to the top feature of the target if scale is true)
    FEAT->resamplingmethod(_dataInterpolator);
    FEAT->is_sparse(_issparse);
    SPH_orig=FEAT->Initialize(_genesis[level],MESHES,_exclude);  /// downsamples and smooths data, creates and exclusion mask if exclude is true
    cout << " after FEAT INIT " << endl;

    if(cost[level]=="AFFINE"){
      throw  MeshReg_error("GroupMeshREG ERROR:: affine not currently available");
    }else{

      cout << " 1 " << endl;
      bool multivariate=true;
      if(FEAT->get_dim()==1){multivariate=false; if (_simval[level]==4) {throw  MeshReg_error("MeshREG ERROR:: simval option 4 (alphaMI) is not suitable for univariate costfunctions");}}

      if(multivariate==true)  throw  MeshReg_error("GroupMeshREG ERROR:: multivariate not currently available");
     
      PARAMETERS.insert(parameterPair("multivariate",multivariate));  // check whether data is multivariate (and NOT using patch based alignment) and add this to the parameter set

      cout << " create MODEL " << _debug << endl;
      MODEL=boost::shared_ptr<SRegDiscreteModel>(new GroupDiscreteModel(PARAMETERS));
      cout << " after create Model " << endl;
      if(_debug) MODEL->set_debug();
      // if(_L1path=="" && _quartet==true)  throw  MeshReg_error("GroupMeshREG ERROR::quartet version requires matlab path");
      //else MODEL->set_L1path(_L1path);
      cout << " set feat " << endl;
      MODEL->set_featurespace(FEAT);
      cout << " set meshspace " << endl;

      MODEL->set_meshspace(SPH_orig,SPH_orig,MESHES.size());
      MODEL->Initialize(CONTROL);
    
      
    }

  }

  void GroupMeshReg::Evaluate()
  {

    cout << " in Evaluate " <<_level<< endl;
    newmesh OLDREG;

    for(int n=0;n<MESHES.size();n++){
     
      if(_level==1) ALL_SPH_REG.push_back(project_CPgrid(SPH_orig,OLDREG)); // first project data grid through any predefined transformation or, from transformation from previous resolution level
      else {
	OLDREG=ALL_SPH_REG[n];
	ALL_SPH_REG[n]=project_CPgrid(SPH_orig,OLDREG,n);
      }
    }
   
    cout << "  run discrete opt " <<  ALL_SPH_REG.size() << endl;
      // sprintf(filename,"SPHreg-EVal%d.surf", _level); SPH_reg.save(filename);
    run_discrete_opt(ALL_SPH_REG);
    
    if(_verbose) cout << "exit main algorithm     " << endl;
  
  }

  void GroupMeshReg::Transform(const string &filename){
    char fullpath[1000];
    // char buffer[1000];

    if(_verbose) cout << " Transform Group" << endl;
    
    for(int n=0;n<MESHES.size();n++){
      barycentric_mesh_interpolation(MESHES[n],SPH_orig,ALL_SPH_REG[n]); 
   
      sprintf(fullpath,"%ssphere-%d.reg%s",filename.c_str(), n,_surfformat.c_str());
      cout <<n << " " << SPH_orig.nvertices() << " " << ALL_SPH_REG[n].nvertices() << " fullpath " << fullpath <<  endl;

      MESHES[n].save(fullpath);
    }

  }

  void GroupMeshReg::saveTransformedData(const double &sigma, const string &filename){
    if(_verbose) cout << " save transformed data " << endl;

    resampler R; R.set_method("ADAP_BARY"); 
    boost::shared_ptr<RELATIONS > REL;
    double ang; ang=R.guess_angular_spacing(TEMPLATE.nvertices());
    char fullpath[1000];

 
    for(int n=0;n<MESHES.size();n++){
      sprintf(fullpath,"%stransformed_and_reprojected-%d%s",filename.c_str(),n,_dataformat.c_str());
	
      boost::shared_ptr<BFMatrix> DATA;
   

      REL = boost::shared_ptr<RELATIONS > ( new RELATIONS(MESHES[n],TEMPLATE,ang));
      REL->update_RELATIONS(MESHES[n]);
  
      set_data(DATAlist[n],DATA,MESHES[n]);
    
      R.resampledata(MESHES[n],TEMPLATE,DATA,0.0,REL);

      boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(DATA);
      TEMPLATE.set_pvalues(DATA->AsMatrix());
      TEMPLATE.save(fullpath);
    }
  }

  ////// iterates over discrete optimisations/////
    void GroupMeshReg::run_discrete_opt(vector<NEWMESH::newmesh> &source){
    resampler R;  R.set_method("ADAP_BARY");
    int iter=1;
    
    NEWMESH::newmesh transformed_controlgrid,targetmesh; 
    vector<newmesh> controlgrid;
    boost::shared_ptr<RELATIONS> controlgrid_neighbourhood;
    myparam::iterator it;
    int res,_itersforlevel;
    int numNodes;
    double energy=0,newenergy=0;
    char filename[1000];
    ofstream Energyout; 

    it=PARAMETERS.find("CPres");res=boost::get<int>(it->second);
    it=PARAMETERS.find("iters");_itersforlevel=boost::get<int>(it->second);
   
    numNodes=MODEL->getNumNodes();
    targetmesh=MODEL->get_TARGET();
      
  
    while(iter <= _itersforlevel){
      controlgrid.clear();
      cout << " reset meshspace " << endl;
      for(int n=0;n<source.size();n++){
	cout << n << " " << source[n].nvertices() << endl;
	MODEL->reset_meshspace(source[n],n); // source mesh is updated and control point grids are reset
	//	sprintf(filename,"sourcebeforereg-level%d-it%d-m%d.surf.gii",_level, iter,n); source[n].save(filename);
	controlgrid.push_back(MODEL->get_CPgrid(n)); 
	//sprintf(filename,"controlbeforereg-it%d-m%d.surf.gii",iter,n); controlgrid[n].save(filename);

      }
   
      cout << " MODEL->setupCostFunction(); " <<endl;
      MODEL->setupCostFunction();
     
      int *Labels=MODEL->getLabeling();
     
#ifdef HAS_HOCR

	Reduction HOCR_mode;

	if(_discreteOPT.compare(0,4,"HOCR")==0)
	  HOCR_mode=HOCR;
	else if(_discreteOPT=="ELC"){
	  HOCR_mode=ELC_HOCR;
	}else if(_discreteOPT=="ELC_approx"){
	  HOCR_mode=ELC_APPROX;
	}else {throw  MeshReg_error("discrete optimisation mode is not available");}
	
	//	tick_count start = tick_count::now();

	cout << _discreteOPT << " optimise 2 " << endl;
	srand (time(NULL));
	newenergy=Fusion::optimize(MODEL,HOCR_mode,_verbose);

#endif	//	MODEL->saveSTRAINS(iter);
    

	if(iter>1 && ((iter-1) % 2==0) && newenergy>=energy) { cout << iter << " level has converged" << (iter-1) %2  << endl; break;}

 
	for (int i = 0; i <numNodes; i++){
	  if(_verbose)
	    cout << i << " _iter " << iter <<" _iters " << _itersforlevel <<  " label " << Labels[i] <<  endl;
	}
   
	MODEL->applyLabeling();
	controlgrid_neighbourhood=MODEL->get_cp_neighbourhood();

	// apply these choices in order to deform the CP grid 
	for(int n=0;n<source.size();n++){
	  transformed_controlgrid=MODEL->get_CPgrid(n);     
	  cout << n << " transformed_controlgrid.nvertices() " << transformed_controlgrid.nvertices() << " controlgrid[n].nvertices() " << controlgrid[n].nvertices() << " source vertices " << " " << source[n].nvertices() << " " <<  controlgrid_neighbourhood->Ncols() << endl;
	  barycentric_mesh_interpolation(source[n],controlgrid[n],transformed_controlgrid,controlgrid_neighbourhood);
      
	  unfold(transformed_controlgrid);
	  MODEL->reset_CPgrid(transformed_controlgrid,n); // source mesh is updated and control point grids are reset   
	//}

      //sprintf(filename,"controlgridtrans-%d-%d.surf",iter,res);
      //transformed_controlgrid.save(filename);
      //sprintf(filename,"sourcetrans-%d-%d.surf",iter,res);
      //source.save(filename);

	  unfold(source[n]);
	}
	energy=newenergy;
	iter++;
    }

      

  }
}
   

 
