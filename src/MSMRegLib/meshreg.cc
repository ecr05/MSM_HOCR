/*  meshreg.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "meshreg.h"
#include <time.h>

namespace MESHREG {


  void MeshReg::run_multiresolutions(const int &levels, const double &sigma, const string &parameters){
#ifdef  HAS_TBB
    tick_count start = tick_count::now();
#else 
    clock_t start; start=clock();
#endif
  
  

   srand (time(NULL));

   parse_reg_options(parameters);
   if(levels!=0) _resolutionlevels=levels;

  

   for(int i=0;i<_resolutionlevels;i++){
     _level=i+1;
     if(_verbose) cout << " Initialising level " << i+1 << endl;
     ///// set parameters /////////////////
     fix_parameters_for_level(i);
     Initialize_level(i);
     Evaluate();
     if(cost[i]=="AFFINE") affinecf.reset();  // clear memory following affine initilisalisation
   }
  
   if(_verbose){
#ifdef  HAS_TBB
     tick_count end = tick_count::now();
     cout<<"Time elapsed: " << (end - start).seconds()  << endl;
#else
     cout<<"Time elapsed: " << ( clock() - start ) / (double) CLOCKS_PER_SEC  << endl;
#endif
   }
   /////////////////// WRITE OUT DATA ACCORDING TO FORMAT /////////////////////////
  
   Transform(_outdir);   
   
   

   saveSPH_reg(_outdir);
   saveTransformedData(sigma,_outdir);
  


  }

  
  
  void MeshReg::Initialize_level(int level)
  {
    MeshModify::Initialize();
    Matrix CombinedWeight;
    vector<double>  sigma(2,0);
    sigma[0]=_sigma_in[level]; sigma[1]=_sigma_ref[level];
    //  MESHES[0].save("inputinsavetransformed_orig.surf.gii"); 
    //MESHES[1].save("refinsavetransformed_orig.surf.gii");  
     

    //////////////  INITIALIZE FEATURESPACE ////////////////////////////////////
    if(_usetraining){
      FEAT=boost::shared_ptr<featurespace>(new featurespace(CMfile_in,DATAlist));

    }
    else
      FEAT=boost::shared_ptr<featurespace>(new featurespace(CMfile_in,CMfile_ref));
   

    FEAT->set_smoothing_parameters(sigma);

    FEAT->set_cutthreshold(_threshold); // will also generate exclusion masks at the same mesh resolution as datagrid
    FEAT->logtransform(_logtransform);// if true logtransforms AND normalises
    FEAT->varnorm(_varnorm);// variance normalises
    FEAT->intensitynormalize(_IN, _scale,_cut); // matches the intensities of the source to to the target (will rescale all to the top feature of the target if scale is true)
    FEAT->resamplingmethod(_dataInterpolator);
    FEAT->is_sparse(_issparse);
    SPH_orig=FEAT->Initialize(_genesis[level],MESHES,_exclude);  /// downsamples and smooths data, creates and exclusion mask if exclude is true
    MVD=Calculate_MVD(SPH_orig);
    boost::shared_ptr<NEWMESH::newmesh> IN_EXCL=FEAT->get_input_excl();
    boost::shared_ptr<NEWMESH::newmesh> REF_EXCL=FEAT->get_reference_excl();
  
    SPHin_CFWEIGHTING=downsample_cfweighting(_sigma_in[level],SPH_orig,IN_CFWEIGHTING,IN_EXCL); /// COST FUNCTION WEIGHTINGS ARE IN ADDITION TO EXCLUSION MASKS - EXCLUSION MASKS ARE BINARY
    SPHref_CFWEIGHTING=downsample_cfweighting(_sigma_ref[level],SPH_orig,REF_CFWEIGHTING,REF_EXCL);
   
    if(cost[level]=="AFFINE"){ /// currently the affine alignment uses a different class to the discrete, but will hopefully implement an affine discrete optimisation soon
      affinecf = boost::shared_ptr<affineMeshCF>(new affineMeshCF(SPH_orig,SPH_orig,FEAT)); 
      affinecf->set_parameters(PARAMETERS);
      if(_simval[level]!=1) affinecf->set_simmeasure(1); 
      affinecf->Initialize();
      _isaffine=true;
    }
    if(cost[level]=="DISCRETE"){
      if(_simval[level]==3) cout << " warning NMI similarity measure does not take into account cost function weights or exclusion values" << endl; 
      bool multivariate;
      if(FEAT->get_dim()==1){multivariate=false; if (_simval[level]==4) {throw  MeshReg_error("MeshREG ERROR:: simval option 4 (alphaMI) is not suitable for univariate costfunctions");}}
      else multivariate=true;

      PARAMETERS.insert(parameterPair("multivariate",multivariate));  // check whether data is multivariate (and NOT using patch based alignment) and add this to the parameter set

      MODEL=boost::shared_ptr<SRegDiscreteModel>(new NonLinearSRegDiscreteModel(PARAMETERS));

      if(_usetraining){
	MODEL->set_L1path(_L1path);
      }
      if(_debug) MODEL->set_debug();

      MODEL->set_featurespace(FEAT,_concattraining);
      MODEL->set_meshspace(SPH_orig,SPH_orig);
      NEWMESH::newmesh CONTROL; // original low res icosphere mesh
      CONTROL.make_mesh_from_icosa(_gridres[level] ); true_rescale(CONTROL,RAD);

      if(_anat){
	vector<vector<int > > ANAT_face_neighbourhood;
	vector<map<int,double> > ANAT_to_CP_baryweights;
	newmesh aICO=resample_anatomy(CONTROL,ANAT_to_CP_baryweights, ANAT_face_neighbourhood, level);
	newmesh ANAT_target=mesh_resample(ref_anat,MESHES[1],aICO);

	MODEL->set_anatomical_meshspace(aICO,ANAT_target,aICO,ANAT_orig);
	MODEL->set_anatomical_neighbourhood(ANAT_to_CP_baryweights,ANAT_face_neighbourhood);
	
      }else{if(_regmode==5){throw  MeshReg_error("STRAINS based regularisation requires anatomical meshes");}
	else if(_regmode==4){throw  MeshReg_error("You have specified angular based penalisation of anatomical warp, for which you require anatomical meshes");}}
     
    
      MODEL->Initialize(CONTROL);
      
   
      
      
      _isaffine=false;
    }
   
         
    
  }

  newmesh  MeshReg::resample_anatomy(newmesh control_grid, vector<map<int,double> > &baryweights, vector<vector<int > > & ANAT_to_CPgrid_neighbours,int level){

    if(in_anat.nvertices()!=MESHES[0].nvertices() || ref_anat.nvertices()!=MESHES[1].nvertices())  {throw  MeshReg_error("MeshREG ERROR:: input/reference anatomical mesh resolution is inconconsistent with input/reference spherical mesh resolution ");}
    
    newmesh ANAT_ico=control_grid;
    resampler R;
   
    ///// save face neighbourhood relationships during retesselation ////
    vector<vector<int> >  FACE_neighbours, FACE_neighbours_tmp;
    vector<int> tmp;
    Pt v0,v1,v2;
    int id0,id1,id2;

    if( (_anatres[level]-_gridres[level]) > 0) {
      for (int i=0;i<(_anatres[level]-_gridres[level]) ;i++){
	ANAT_ico.retessellate(FACE_neighbours_tmp);
	if(i>0) {
	  ANAT_to_CPgrid_neighbours.clear();
	  for(unsigned int j=0;j<FACE_neighbours.size();j++){
	    ANAT_to_CPgrid_neighbours.push_back(tmp);
	    for(unsigned int k=0;k<FACE_neighbours[j].size();k++){
	      ANAT_to_CPgrid_neighbours[j].insert(ANAT_to_CPgrid_neighbours[j].begin(),FACE_neighbours_tmp[FACE_neighbours[j][k]].begin(),FACE_neighbours_tmp[FACE_neighbours[j][k]].end());
	      
	    }
	   
	  }
	  FACE_neighbours=ANAT_to_CPgrid_neighbours;
	}else 
	  {FACE_neighbours=FACE_neighbours_tmp;
	  }
	
	FACE_neighbours_tmp.clear();
      }
     
      if(ANAT_to_CPgrid_neighbours.size()==0) { ANAT_to_CPgrid_neighbours=FACE_neighbours;  } // i.e. for one increase in resolution use result directly from restesselation

     
    }else{
      for (int i=0;i<control_grid.ntriangles() ;i++){
	ANAT_to_CPgrid_neighbours.push_back(tmp);
	ANAT_to_CPgrid_neighbours[i].push_back(i);
      }
    }
    true_rescale(ANAT_ico,RAD);
    
    /// now get barycentric weights 
    baryweights.resize(ANAT_ico.nvertices(),map<int,double> ());

    for(unsigned int i=0;i<ANAT_to_CPgrid_neighbours.size();i++){
      id0=control_grid.get_triangle_vertexID(i,0);      id1=control_grid.get_triangle_vertexID(i,1);   id2=control_grid.get_triangle_vertexID(i,2);
      v0=control_grid.get_triangle_vertex(i,0);      v1=control_grid.get_triangle_vertex(i,1);   v2=control_grid.get_triangle_vertex(i,2);
      
      for(unsigned int j=0;j<ANAT_to_CPgrid_neighbours[i].size();j++){
	for(int k=0;k<3;k++){
	  Pt ci=ANAT_ico.get_triangle_vertex(ANAT_to_CPgrid_neighbours[i][j],k);
	    int id=ANAT_ico.get_triangle_vertexID(ANAT_to_CPgrid_neighbours[i][j],k);
	    baryweights[id]=get_barycentric_weights(v0,v1,v2,ci,id0,id1,id2);
	    //  ANAT_ico.set_pvalue(id,i);
	}
      }
    }
    
    ANAT_orig=mesh_resample(in_anat,MESHES[0],ANAT_ico);
    
   
    return ANAT_ico;
  }

  Matrix  MeshReg::downsample_cfweighting(const double &sigma , NEWMESH::newmesh SPH, boost::shared_ptr<NEWMESH::newmesh> CFWEIGHTING, boost::shared_ptr<NEWMESH::newmesh> EXCL){
    resampler R;
    Matrix DATA;
    R.set_method("NN");
    
    
    if(EXCL.get()){
      if(CFWEIGHTING.get()){
	DATA=CFWEIGHTING->get_pvalues();
      }else{  
	CFWEIGHTING=boost::shared_ptr<NEWMESH::newmesh>(new newmesh(*EXCL)); // if no CF weighting mask is supplied, then CF weighting is just equal to the EXCLUSION mask 
	DATA=CFWEIGHTING->get_pvalues();
      }
      
      /// resample onto reference mesh and print
      R.resampledata(*CFWEIGHTING,SPH,EXCL,DATA,sigma);  /// resample onto the data grid with smoothing
   
    
    }
    else if(CFWEIGHTING.get()){
      DATA=CFWEIGHTING->get_pvalues();
      R.resampledata(*CFWEIGHTING,SPH,EXCL,DATA,0);    
    }
    else{
      DATA.ReSize(1,SPH.nvertices()); DATA=1;
 
    }
    return DATA;
  }

  ///////////////////////////////////// MAIN FUNCTION //////////////////////////////////
  
  void MeshReg::Evaluate() {

    ///Initiallise deformation mesh
   
    SPH_reg=project_CPgrid(SPH_orig,SPH_reg); // first project data grid through any predefined transformation or, from transformation from previous resolution level
	

    if(_isaffine){
      affinecf->update_source(SPH_reg);
      SPH_reg=affinecf->run();
    }
    else{
      // sprintf(filename,"SPHreg-EVal%d.surf", _level); SPH_reg.save(filename);
      run_discrete_opt(SPH_reg);
     
    }
    if(_verbose) cout << "exit main algorithm     " << endl;
    
  }
  
  
  void MeshReg::Transform(const string &filename){

    newmesh INtransformed;
    char fullpath[1000];
    // char buffer[1000];
    
    if(_verbose) cout << " Transform " << endl;
    if(!_meshInterpolator.compare("BARY")){ barycentric_mesh_interpolation(MESHES[0],SPH_orig,SPH_reg); }
    else{
      if(SPH_orig.nvertices()<MESHES[0].nvertices()){// TPS interpolation, matrices can be singular if there aren't enough samples range is controlled by angle 2*asin(MVD/RAD)
	upsample_transform_RBF(MESHES[0],SPH_orig,SPH_reg,2*asin(MVD/RAD));
      }else if(SPH_orig.nvertices()>MESHES[0].nvertices()){
	project(MESHES[0],SPH_orig,SPH_reg,2*asin(MVD/RAD)); // for above reason local affine projections work better than TPS in this instance
	unfold(MESHES[0]);
      }else{  MESHES[0]=SPH_reg;}
    }
    sprintf(fullpath,"%ssphere.reg%s",filename.c_str(), _surfformat.c_str());

    MESHES[0].save(fullpath);
    
  }

  void MeshReg::saveTransformedData(const double &sigma, const string &filename){
    if(_verbose) cout << " save transformed data " << endl;
    char fullpath[1000];

    sprintf(fullpath,"%stransformed_and_reprojected%s",filename.c_str(),_dataformat.c_str());
    resampler R;
    boost::shared_ptr<RELATIONS > REL;
    double ang;

    boost::shared_ptr<BFMatrix> DATA;
    boost::shared_ptr<BFMatrix> DATAREF;
    boost::shared_ptr<NEWMESH::newmesh> IN_EXCL;
    boost::shared_ptr<NEWMESH::newmesh> REF_EXCL;

    if(sigma<1e-3){   R.set_method("ADAP_BARY"); 
      if(MESHES[1].nvertices() > MESHES[0].nvertices()){ang=R.guess_angular_spacing(MESHES[0].nvertices()); }
	else{ang=R.guess_angular_spacing(MESHES[0].nvertices()); }
      }
      else{ R.set_method("GAUSSIAN"); ang=2*asin(sigma/(RAD));}
    
    REL = boost::shared_ptr<RELATIONS > ( new RELATIONS(MESHES[1],MESHES[0],ang));
    REL->update_RELATIONS(MESHES[1]);
  
    if(_usetraining){
      set_data(CMfile_in,DATA,MESHES[0]);
      R.resampledata(MESHES[0],MESHES[1],IN_EXCL,DATA,sigma,REL);

    }
    else{
      /// binarize costfunction weights for use as exclusion masks during resampling ////
    
    
    
      set_data(CMfile_in,DATA,MESHES[0]);
      set_data(CMfile_ref,DATAREF,MESHES[1]);
  
      
      if(_exclude){ // no longer necessary as EXCL masks aren't downsampled anymore? Could save from initialisation
	IN_EXCL= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(MESHES[0],DATA->AsMatrix(),_threshold[0],_threshold[1]))); 
	REF_EXCL= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(MESHES[1],DATAREF->AsMatrix(),_threshold[0],_threshold[1]))); 
      
      }

      if(_IN){ /// INTENSITY NORMALIZE
	if(_verbose) cout << " intensity normalise " << endl;
	multivariate_histogram_normalization(*DATA,*DATAREF,IN_EXCL,REF_EXCL);    
	MESHES[0].set_pvalues(DATA->AsMatrix());
	//MESHES[0].save("intensity_normalised.func.gii");
      }
    
      
      /// resample onto reference mesh and print
    
      R.resampledata(MESHES[0],MESHES[1],IN_EXCL,DATA,sigma,REL);

    }

   
  
    newmesh TRANSFORMED=MESHES[1];
    boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(DATA);
    
    if(pin){
      
      TRANSFORMED.set_pvalues(DATA->AsMatrix());
      // cout << "write " << endl;
      TRANSFORMED.save(fullpath);
    }else DATA->Print(fullpath);
   
    if(_anat){
      char filename2[1000];
    
      newmesh ANAT_TRANS=projectmesh(MESHES[0],MESHES[1],ref_anat);
      sprintf(filename2,"%sanat.reg.surf",_outdir.c_str());
      ANAT_TRANS.save(filename2);
      boost::shared_ptr<Matrix>  VECS; VECS= boost::shared_ptr<Matrix> (new Matrix);
      in_anat.estimate_normals(); ANAT_TRANS.estimate_normals();
      newmesh STRAINSmesh=calculate_strains(2,in_anat,ANAT_TRANS,VECS);calculate_triangular_strains(in_anat,ANAT_TRANS);
      sprintf(filename2,"%sSTRAINS.func",_outdir.c_str());

    
      STRAINSmesh.save(filename2);
    }
  }
   
  ////////////////////// PROJECT TRANSFORMATION FROM PREVIOUS LEVEL TO UPSAMPLED SOURCE ////////////////////////////////////
  NEWMESH::newmesh MeshReg::project_CPgrid(NEWMESH::newmesh SPH_in,NEWMESH::newmesh  REG, int num){ /// num indices which warp for group registration
 
   
   if(_level==1){
      if(transformed_mesh.nvertices()>0){ 
      //// project into alignment with transformed mesh
	if(transformed_mesh==MESHES[num]){
	  cout << " WARNING:: transformed mesh has the same coordinates as the input mesh " << endl;
	}else{
	  // unfold(transformed_mesh);
	  barycentric_mesh_interpolation(SPH_in,MESHES[num],transformed_mesh);
	  if(MODEL.get()) MODEL->warp_CPgrid(MESHES[num],transformed_mesh,num); // for tri clique model control grid is continously deformed
	  // SPH_in.save("SPH_in.surf");
	}
      }
      
   }else {
      ///// following first round always start by projecting Control and data grids through warp at previous level
     
     /////////// PROJECT CPgrid into alignment with warp from previous level
     NEWMESH::newmesh icotmp;
    
     int ico=REG.get_ico_resolution();
     icotmp.make_mesh_from_icosa(ico); true_rescale(icotmp,RAD);    
     
     
     ///// projecct datagrid though warp defined for the high resolution meshes (the equivalent to if registration is run one level at a time )
     newmesh inorig=MESHES[0];
     newmesh incurrent=MESHES[0];
      
     barycentric_mesh_interpolation(incurrent,icotmp,REG); 
     barycentric_mesh_interpolation(SPH_in,inorig,incurrent); 
     if(MODEL.get())	 MODEL->warp_CPgrid(inorig,incurrent); 
	 
     if(_debug){
       char filename[1000];
       sprintf(filename,"%ssphere.regLR.Res%d.surf" ,_outdir.c_str(),_level);
       incurrent.save(filename);
     }
     

    // if(SPH_orig.nvertices()==REG.nvertices()){SPH_in=REG;} //// if resolution of data grid has not changed just initialise equal to previous warp
     //else{
     // double MVDtmp=Calculate_MVD(icotmp);
	
      // if(SPH_orig.nvertices()<REG.nvertices()){cout << "Warning  The next datamesh grid is lower resolution than the previous "  << endl; 
	 //barycentric_mesh_interpolation(SPH_in,icotmp,REG);}
     //  else{
	 //if(!_meshInterpolator.compare("BARY")){
	 //  barycentric_mesh_interpolation(SPH_in,icotmp,REG);
	// }
	 //else{
	 
	// upsample_transform_RBF(SPH_in,icotmp,REG,2*asin(MVDtmp/RAD));
	// }
      // }
	
       
   //  }
   }
  
   unfold(SPH_in);
   // SPH_in.save("SPH_in2.surf");
	 
    return SPH_in;
     
  }

  ////// iterates over discrete optimisations/////
  void MeshReg::run_discrete_opt(NEWMESH::newmesh &source){
    resampler R;  R.set_method("ADAP_BARY");
    int iter=1;
    
    NEWMESH::newmesh controlgrid,transformed_controlgrid,targetmesh; 
    boost::shared_ptr<RELATIONS> controlgrid_neighbourhood;
    myparam::iterator it;
    Matrix ResampledRefWeight;
    Matrix CombinedWeight;
    newmesh sourcetmp=source;
    int _itersforlevel;
    int numNodes;
    double energy=0,newenergy=0;
   
      
    it=PARAMETERS.find("iters");_itersforlevel=boost::get<int>(it->second);
   
    numNodes=MODEL->getNumNodes();
    targetmesh=MODEL->get_TARGET();
    controlgrid=MODEL->get_CPgrid();   
  
    while(iter <= _itersforlevel){

      //////// resample and combine the reference cost function weighting with the source
      ResampledRefWeight=SPHref_CFWEIGHTING;
        
     
      R.resampledata(targetmesh,source,ResampledRefWeight,0); 
	
      CombinedWeight=combine_costfunction_weighting(SPHin_CFWEIGHTING,ResampledRefWeight); // elementwise multiplication
    
    
      MODEL->reset_meshspace(source); // source mesh is updated and control point grids are reset
      MODEL->setupCostFunctionWeighting(CombinedWeight);
      MODEL->setupCostFunction();

      int *Labels=MODEL->getLabeling();
      if(_verbose) cout << "run optimisation" << endl;

      if(_discreteOPT=="FastPD"){
	MODEL->computeUnaryCosts();
	MODEL->computePairwiseCosts();
	

	FastPD opt(MODEL, 100 );
	newenergy = opt.run();

	boost::shared_ptr<RELATIONS> reltmp=MODEL->get_cp_neighbourhood();
	opt.getLabeling(Labels);

      }else {
#ifdef HAS_HOCR

	Reduction HOCR_mode;

	if(_discreteOPT.compare(0,4,"HOCR")==0)
	  HOCR_mode=HOCR;
	else if(_discreteOPT=="ELC"){
	  HOCR_mode=ELC_HOCR;
	}else if(_discreteOPT=="ELC_approx"){
	  HOCR_mode=ELC_APPROX;
	}else {throw  MeshReg_error("discrete optimisation mode is not available");}
	

	srand (time(NULL));
	newenergy=Fusion::optimize(MODEL,HOCR_mode,_verbose);


#endif
      }

     
      if(iter>1 && ((iter-1) % 2==0) && (energy -newenergy < 0.001) ) 
     	{ if(_verbose) cout << iter << " level has converged" << (iter-1) %2  << " newenergy " << newenergy <<  " energy " << energy <<  " energy- newenergy " <<  energy -newenergy <<   endl; break;}
      else {if(_verbose) cout <<  " newenergy " << newenergy <<  " energy " << energy <<  " energy- newenergy " <<  energy -newenergy <<  endl;}
	
  
      //
      //Run optimization
      //
      //Get labeling result

      if(_verbose) cout << "Output Control Point Label assignment:" << endl;
       for (int i = 0; i <numNodes; i++){
	if(_verbose)
	  cout << i << " _iter " << iter <<" _iters " << _itersforlevel <<  " label " << Labels[i] <<  endl;
      }
      //Get initial energy
      //   double m_initial_energy = opt.getInitialEnergy();
       //  double energybedorelabeling=energy->evaluateTotalCostSum();
      /// get labelling choices from the FastPD optimiser
      MODEL->applyLabeling();
      // apply these choices in order to deform the CP grid 
      transformed_controlgrid=MODEL->get_CPgrid();     
      controlgrid_neighbourhood=MODEL->get_cp_neighbourhood();

      /// use the control point updates to warp the source mesh
      if(!_meshInterpolator.compare("BARY")){
	
	barycentric_mesh_interpolation(source,controlgrid,transformed_controlgrid,controlgrid_neighbourhood);
      }
      else{
	upsample_transform_RBF(source,controlgrid,transformed_controlgrid,*controlgrid_neighbourhood);
      }
      //if(iter==2) exit(1);
      // if(_discreteOPT.compare(0,4,"HOCR")==0){ 
      controlgrid=transformed_controlgrid; /// higher order frameowrk continuous deforms the CP grid whereas the original FW resets the grid each time
      //controlgrid.save("transformed_controlgrid.surf");

      //cout << " unfold " << endl;
      unfold(transformed_controlgrid);
     // cout << " reset control grid " << endl;
      MODEL->reset_CPgrid(transformed_controlgrid); // source mesh is updated and control point grids are reset   
//transformed_controlgrid.save("transformed_controlgridunfold.surf");
//source.save("source.surf");


      unfold(source);
     // source.save("sourceunfold.surf");

      energy=newenergy;
      iter++;
    }

     

  }

}
   

 
