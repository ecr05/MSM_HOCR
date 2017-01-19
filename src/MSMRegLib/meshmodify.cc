/*  meshmodify.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "meshmodify.h"

//// NOTE SIM KERNEL IS TRANSPOSED RELATIVE TO ALEKS CODE TO CORRESPEND TO TRANSPOSED SPARSE MATRICES - input subjects represent columns
namespace MESHREG {

  void MeshModify::parse_reg_options(const string &parameters)
  {
    string title="msm configuration parameters ";
    string examples="";
    OptionParser options(title,examples);

    vector<string> costdefault;
    Option<vector<string> >  optimizer(string("--opt"),costdefault,
				       string("optimisation method. Choice of: AFFINE,DISCRETE (default)"),
				       false,requires_argument);

    vector<int> intdefault;
    Option<vector<int> > simval(string("--simval"), intdefault,
				string("code for determining which similarty measure is used to assess cost during registration: options are 1) SSD; 2) pearsons correlation (default); 3) NMI;)"),
				false, requires_argument);
 
    Option< vector<int> > iterations(string("--it"), intdefault,
				     string("number of iterations at each resolution (default -–it=3,3,3)"),
				     false, requires_argument);
    
    vector<float> floatdefault;
    
    Option< vector<float> > sigma_in(string("--sigma_in"),floatdefault,
				     string("smoothing parameter for input image (default --sigma_in=2,2,2)"),
				     false, requires_argument);
    
    Option< vector<float> > sigma_ref(string("--sigma_ref"),  floatdefault,
				      string("Sigma parameter - smoothing parameter for reference image  (set equal to sigma_in by default)"),
				      false, requires_argument);
    
    Option<vector<float> >  lambda(string("--lambda"),  floatdefault,
				   string("Lambda parameter - controls contribution of regulariser "),
				   false, requires_argument);
    
    Option<vector<int> >  datagrid(string("--datagrid"),intdefault,
				   string("DATA grid resolution (default --datagrid=5,5,5). If parameter = 0 then the native mesh is used."),
				   false, requires_argument);
    
  Option<vector<int> >  cpgrid(string("--CPgrid"),intdefault,
			       string("Control point grid resolution (default --CPgrid=2,3,4)"),
			       false, requires_argument);
  
  Option< vector<int> >  sampgrid(string("--SGgrid"),intdefault,
				  string("Sampling grid resolution (default = 2 levels higher than the control point grid)"),
				  false, requires_argument); 
  Option< vector<int> >  anatgrid(string("--anatgrid"),intdefault,
				  string("Anatomical grid resolution (default = 2 levels higher than the control point grid)"),
				  false, requires_argument); 
  
  
  Option< vector<int> > alpha_knn(string("--aKNN"),intdefault,
				  string("Number of neighbours for estimation of kNN graph and alpha entropy measure (default --aKNN=5,5,5)"),
				  false, requires_argument); 
  
  
  vector<float> cutthresholddefault(2,0); cutthresholddefault[1]=0.0001;
  Option< vector<float> > cutthreshold(string("--cutthr"),cutthresholddefault,
				       string("Upper and lower thresholds for defining cut vertices (default --cutthr=0,0)"),
				       false, requires_argument);
  
  
  Option<string>  meshinterpolationmethod(string("--mInt"), "BARY",
					  string("Method used for mesh interpolations, options: TPS or BARY (default)"),
					  false,requires_argument);
  
  Option<string>  datainterpolationmethod(string("--dInt"), "ADAP_BARY",
					  string("Method used for data interpolations, options: GAUSSIAN or ADAP_BARY (default)"),
					  false,requires_argument);
  
#ifdef HAS_HOCR
  Option<int>  regulariseroption(string("--regoption"), 1,
		  string("Choose option for regulariser form lambda*weight*pow(cost,rexp). Where cost can be PAIRWISE or TRI-CLIQUE based. Options are: 1) PAIRWISE - penalising diffences in rotations of neighbouring points (default); 2) TRI_CLIQUE Angle deviation penalty (for spheres); 3) TRI_CLIQUE: Strain-based (for spheres);  4) TRI_CLIQUE Angle deviation penalty (for anatomy); 5) TRI_CLIQUE: Strain-based (for anatomy)"),
				 false, requires_argument);
  Option<string>  doptimizer(string("--dopt"),"FastPD",
			     string("discrete optimisation implementation. Choice of: FastPD (default), HOCR (will reduce to QBPO for pairwise), ELC, ELC_approx"), 
			     false,requires_argument,false);
  
  Option<bool>  tricliquelikeihood(string("--triclique"), false,
				   string("estimate similarity for triangular patches (rather than circular)"),
				   false, no_argument);				   
  Option<float>  shear(string("--shearmod"), 0.4,
		       string("shear modulus (default 0.4); for use with --regoptions 3 "),
			 false,requires_argument);
  Option<float>  bulk(string("--bulkmod"), 1.6,
		      string("bulk mod (default 1.6); for use with --regoptions 3 "),
		      false,requires_argument);
  
  Option<float>  grouplambda(string("--glambda_pairs"), 1,
			     string("scaling for pairwise term in greoup alignment"),
			     false,requires_argument,false);
  Option<float>  kexponent(string("--k_exponent"), 2,
		      string("exponent inside strain equation (default 2)"),
		      false, requires_argument);
#endif
  
  Option<float>  expscaling(string("--scaleexp"), 1,
			    string("Scaling for weight exponent (default 1.0)"),
			    false,requires_argument,false);
  
  Option<float>  regulariserexp(string("--regexp"), 2.0,
				string("Regulariser exponent 'rexp' (default 2.0)"),
				false,requires_argument);
  
  Option<bool>  distweight(string("--weight"), false,
			   string("weight regulariser cost using areal distortion weighting"),
			   false, no_argument);
  Option<bool>  anorm(string("--anorm"), false,
			   string("norm regulariser cost using mean angle (for HCP compatibility)"),
			   false, no_argument);
  Option<bool>  rescale_labels(string("--rescaleL"), false,
				   string("rescale label grid rather than using barycentres"),
				   false, no_argument);
				   
  Option<float>  maxdist(string("--maxdist"), 4,
			 string("Set areal distortion threshold (default 4); for use with --regoptions 2 "),
			 false,requires_argument,false);
  
  Option<float>  pottsenergy(string("--potts"),0.0,
			     string("Use potts model for mrf regulariser (strain only thus far). Supply threshold value"),
			     false, requires_argument, false);
  
  Option<float>  controlptrange(string("--cprange"), 1.0,
				string("Range (as % control point spacing) of data samples (default 1) "),
				false,requires_argument,false);
  
  
  Option<bool>  logtransform(string("--log"), false,
			     string("log transform and normalise the data"),
			     false, no_argument);
  Option<bool>  intensitynormalize(string("--IN"), false,
				   string("Normalize intensity ranges using histogram matching "),
				   false, no_argument); 
  Option<bool>  intensitynormalizewcut(string("--INc"), false,
				   string("Normalize intensity ranges using histogram matching excluding cut"),
				   false, no_argument); 				   
  Option<bool>  variancenormalize(string("--VN"), false,
				  string("Variance normalize data "),
				  false, no_argument); 
  
  Option<bool>  scaleintensity(string("--scale"), false,
			       string("Scale intensity ranges of a features in multivariate data to be equal to that of the first (useful for multimodal contrasts)"),
			       false, no_argument); 
  Option<bool>   exclude(string("--excl"), false,
			 string("Ignore the cut when resampling the data"),
			 false, no_argument);
  
  Option<bool>   quartet(string("-Q"), false,
			 string("Estimate quartet low rank cost for group reg"),
			 false, no_argument,false);
  
  ///// AFFINE OPTIONS /////////////////////////
  Option<float>  affinestepsize(string("--stepsize"), 0.01,
				string("gradient stepping for affine optimisation (default 0.01)"),
				false, requires_argument);
  
  Option<float>  gradsampling(string("--gradsampling"), 0.5,
			      string("Determines the finite distance spacing for the affine gradient calculation (default 0.5)"),
			      false,requires_argument);
  Option<int>  threads(string("--numthreads"), 1,
		       string("number of threads for tbb (default 1)"),
		       false,requires_argument);
  int nonoptarg;
  
  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(optimizer);
    options.add(simval);
    options.add(iterations);
    options.add(sigma_in);
    options.add(sigma_ref);
    options.add(lambda);
    options.add(datagrid);
    options.add(cpgrid);
    options.add(sampgrid);
    options.add(anatgrid);
    options.add(alpha_knn);  
    options.add(cutthreshold);
    options.add(meshinterpolationmethod);
    options.add(datainterpolationmethod);
#ifdef HAS_HOCR
    options.add(regulariseroption);
    options.add(doptimizer);
    options.add(tricliquelikeihood);
    options.add(shear);
    options.add(bulk);
    options.add(grouplambda);
    options.add(kexponent);
#endif
    options.add(expscaling);
    options.add(regulariserexp); 
    options.add(distweight);
    options.add(anorm);
    options.add(rescale_labels); 
    options.add(maxdist);
    options.add(pottsenergy);
    options.add(controlptrange);
    options.add(logtransform);
    options.add(intensitynormalize);
    options.add(intensitynormalizewcut);

    options.add(variancenormalize);
    options.add(scaleintensity);
    options.add(exclude);
    options.add(quartet);
    options.add(affinestepsize);
    options.add(gradsampling);
    options.add(threads);

    if(parameters=="usage"){
      options.usage();
      exit(2);
    }
    else if(parameters!=""){
      nonoptarg = options.parse_config_file(parameters);
    }
   
  }
  catch(Utilities::X_OptionError& e) {

    options.usage();
    cerr << "\n" << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
  catch(std::exception &e) {
    cerr << "\n" << e.what() << "\n";
    exit(EXIT_FAILURE);
  }

  if( ! options.check_compulsory_arguments()){
      options.usage();
      
      exit(2);
    }     

  if(parameters==""){
    /// if no config is supplied use sulc config (as of Sept 2014) as default 
    cost.resize(3,"DISCRETE"); cost.insert(cost.begin(),"AFFINE");  _resolutionlevels=cost.size();
    _lambda.resize(4,0); _lambda[1]=0.1;_lambda[2]=0.2;_lambda[3]=0.3;
    _simval.resize(4,2); _simval[0]=1;
    _sigma_in.resize(4,2); _sigma_in[2]=3;_sigma_in[3]=2;
    _sigma_ref.resize(4,2); _sigma_ref[2]=1.5;_sigma_ref[3]=1;
    _iters.resize(4,3); _iters[0]=50;
    _alpha_kNN.resize(cost.size(),5);
    _gridres.resize(4,0); _gridres[1]=2;_gridres[2]=3;_gridres[3]=4; 
    _anatres.resize(4,0); _anatres[1]=4;_anatres[2]=5;_anatres[3]=6; 

    _genesis.resize(4,4); _genesis[2]=5;_genesis[3]=6;
    _sampres.resize(4,0); _sampres[1]=4;_sampres[2]=5;_sampres[3]=6;
   
  }
  else{
    cost=optimizer.value();
    _lambda=lambda.value(); 
    
    // now check for assignments and else set defaults
    if(simval.set()){_simval=simval.value();}
    else _simval.resize(cost.size(),2);
    
    if(iterations.set()){_iters=iterations.value();}
    else _iters.resize(cost.size(),3);
    
    if(sigma_in.set()){_sigma_in=sigma_in.value();}
    else _sigma_in.resize(cost.size(),2);
    
    if(sigma_ref.set()){_sigma_ref=sigma_ref.value();}
    else _sigma_ref=_sigma_in;

    if(datagrid.set()){_genesis=datagrid.value();}
    else _genesis.resize(cost.size(),5);

    if(cpgrid.set()){_gridres=cpgrid.value();}
    else {_gridres.resize(cost.size(),2);
      for (unsigned int i=1;i<cost.size();i++) _gridres[i]=_gridres[i-1]+1;}

   if(anatgrid.set()){_anatres=anatgrid.value();}
    else {_anatres.resize(cost.size(),2);
      for (unsigned int i=0;i<cost.size();i++) _anatres[i]=_gridres[i]+2; }

   if(sampgrid.set()){_sampres=sampgrid.value();}
   else{ _sampres.resize(cost.size());
      for (unsigned int i=0;i<cost.size();i++) _sampres[i]=_gridres[i]+2;}
   
   if(alpha_knn.set()){_alpha_kNN=alpha_knn.value();}
    else _alpha_kNN.resize(cost.size(),5);
  }
  
  _resolutionlevels=cost.size();
  _logtransform=logtransform.value(); _scale=scaleintensity.value();
#ifdef HAS_HOCR
  if(grouplambda.set()) {_set_group_lambda=true;
    _pairwiselambda=grouplambda.value();}
  _regmode=regulariseroption.value(); 
  _discreteOPT=doptimizer.value();
  _tricliquelikeihood=tricliquelikeihood.value();
  _shearmod=shear.value();
  _bulkmod=bulk.value();
  _k_exp=kexponent.value();
  
#else
  _regmode=1;
  _discreteOPT="FastPD";
#endif

   if(intensitynormalizewcut.set()){
	_IN=intensitynormalizewcut.value();
	_cut=true;}
	else{
	_IN=intensitynormalize.value();
	_cut=false;
   }	
  _varnorm=variancenormalize.value(); _exclude=exclude.value(); _quartet=quartet.value();
  _meshInterpolator=meshinterpolationmethod.value();  _dataInterpolator=datainterpolationmethod.value(); _weight=distweight.value(); _regoption2norm=anorm.value(); _potts=pottsenergy.value();
  _threshold=cutthreshold.value(); _regscaling=expscaling.value(); _regexp=regulariserexp.value(); _maxdist=maxdist.value(); _cprange=controlptrange.value();  _affinestepsize=affinestepsize.value(); _affinegradsampling=gradsampling.value(); _numthreads=threads.value();
  _rescale_labels=rescale_labels.value();
 
  if(_verbose){
    cout << "cost.size() " << cost.size() << " " << cost << endl;
    cout << " iters " << _iters << endl;
    cout << "lambda " << _lambda << endl;
    cout << " sigmain " << _sigma_in << endl;
    cout << " sigmaref " << _sigma_ref << endl;
    cout << " datagrid " << _genesis << endl;
    cout << " cp gridres " << _gridres << endl;
    cout << " sp grid " << _sampres << endl;
    cout << " alpha_knn " << _alpha_kNN << endl;
    
    cout << " in parse options " << " _simval " << _simval << " interations " << _iters << " input sigma " << _sigma_in <<" ref sigma " << _sigma_ref << " lambda " << _lambda  << "  _genesis "  << _genesis << " _gridres " <<_gridres  << " _sampres " <<   _sampres <<   " opt " << cost <<  " _scale " << _scale << " mesh interpolator " << _meshInterpolator << " discrete implementation " << _discreteOPT << " regoption " <<  _regmode << " potts " << _potts << endl; 
    
  }

  
  if(_regmode>1 && _discreteOPT=="FastPD"){throw  MeshReg_error("MeshREG ERROR:: you cannot run higher order clique regularisers with fastPD ");}
  if((int) _threshold.size()!=2){throw  MeshReg_error("MeshREG ERROR:: the cut threshold does not contain a limit for upper and lower threshold (too few inputs)");}
  if((int) _simval.size()!=_resolutionlevels) {cout << _simval.size() << " " << _resolutionlevels << endl;throw  MeshReg_error("MeshREG ERROR:: config file parameter list lengths are inconsistent: --simval");}
  if((int)_iters.size()!=_resolutionlevels) {cout << "2 " << endl; throw  MeshReg_error("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--it");}
  if((int)_sigma_in.size()!=_resolutionlevels) {cout << "3 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file  parameter list lengths are inconsistent: --sigma_in");}
  if((int)_sigma_ref.size()!=_resolutionlevels) {cout << "4 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--sigma_ref");}
  if((int)cost.size()!=_resolutionlevels) {cout << "5" << endl;throw  MeshReg_error("MeshREG ERROR:: config file parameter list lengths are inconsistent:--opt");}
  if((int)_lambda.size()!=_resolutionlevels) {cout << "6 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--lambda");}
  if((int)_genesis.size()!=_resolutionlevels) {cout << "7 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--datagrid");}
  if((int)_gridres.size()!=_resolutionlevels) {cout << "8 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file parameter list lengths are inconsistent:--CPgrid");}
  if((int)_sampres.size()!=_resolutionlevels) {cout << "9 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file parameter list lengths are inconsistent:--SGres");}
  if((int)_alpha_kNN.size()!=_resolutionlevels) {cout << "10 " << endl;throw  MeshReg_error("MeshREG ERROR:: config file parameter list lengths are inconsistent:--aKNN");}

 
}
  
  void MeshModify::fix_parameters_for_level(const int & i){
    PARAMETERS.clear();
    
    PARAMETERS.insert(parameterPair("lambda",_lambda[i]));
    PARAMETERS.insert(parameterPair("lambda_pairs",_pairwiselambda));
    PARAMETERS.insert(parameterPair("set_lambda_pairs",_set_group_lambda));
    PARAMETERS.insert(parameterPair("iters",_iters[i]));
    PARAMETERS.insert(parameterPair("simmeasure",_simval[i]));
    PARAMETERS.insert(parameterPair("sigma_in",_sigma_in[i]));
    PARAMETERS.insert(parameterPair("CPres",_gridres[i]));
    PARAMETERS.insert(parameterPair("SGres",_sampres[i]));
    PARAMETERS.insert(parameterPair("anatres",_anatres[i]));
    PARAMETERS.insert(parameterPair("quartet",_quartet));
    PARAMETERS.insert(parameterPair("regularisermode",_regmode));
    PARAMETERS.insert(parameterPair("TriLikelihood",_tricliquelikeihood));
    PARAMETERS.insert(parameterPair("rescalelabels",_rescale_labels));
    PARAMETERS.insert(parameterPair("maxdistortion",_maxdist));
    PARAMETERS.insert(parameterPair("shearmodulus",_shearmod));
    PARAMETERS.insert(parameterPair("bulkmodulus",_bulkmod));
    PARAMETERS.insert(parameterPair("pottsthreshold",_potts));
    PARAMETERS.insert(parameterPair("range",_cprange));
    PARAMETERS.insert(parameterPair("exponent",_regexp));
    PARAMETERS.insert(parameterPair("weight",_weight));
    PARAMETERS.insert(parameterPair("anorm",_regoption2norm));
    PARAMETERS.insert(parameterPair("scaling",_regscaling));
    PARAMETERS.insert(parameterPair("kNN",_alpha_kNN[i]));
    PARAMETERS.insert(parameterPair("verbosity",_verbose));
    PARAMETERS.insert(parameterPair("outdir",_outdir));
    PARAMETERS.insert(parameterPair("stepsize",_affinestepsize));
    PARAMETERS.insert(parameterPair("gradsampling",_affinegradsampling));
    PARAMETERS.insert(parameterPair("numthreads",_numthreads));
    PARAMETERS.insert(parameterPair("kexponent",_k_exp));

    
  }
  

  void MeshModify::check(){

  if(((MESHES[0].get_coord(0)).norm()- RAD) > 1e-5){ throw  MeshReg_error("MeshModify ERROR:: input mesh radius has not been normalised to RAD=100");}
  if(IN_CFWEIGHTING.get() && (IN_CFWEIGHTING->get_coord(0).norm()-RAD)>1e-5 ){throw  MeshReg_error("MeshModify ERROR:: input exclusion mesh radius has not been normalised to RAD=100");}
  if(REF_CFWEIGHTING.get() && (REF_CFWEIGHTING->get_coord(0).norm()-RAD)>1e-5 ){throw  MeshReg_error("MeshModify ERROR::reference exclusion mesh radius has not been normalised to RAD=100");}
 
}

  void MeshModify::Initialize()
  {
    
    check();
  }
  
 
  void MeshModify::set_output_format(string type){
    
    if(type=="GIFTI"){
      _surfformat=".surf.gii";
      _dataformat=".func.gii";
   }
   else if (type=="ASCII" || type=="ASCII_MAT" ){
     _surfformat=".asc";
     if(type=="ASCII")
       _dataformat=".dpv";
     else
       _dataformat=".txt"; // for multivariate
   }
   else{
     _surfformat=".vtk";
     _dataformat=".txt";
   }
    
  }
  
  //// each Matrix can have a different number of rows, for example if the user supplies a costfuncion weighting mask for multivariate features only on the reference, and also sets exclusion weighting then the reference weighting will be multivariate and the source weighting will be univariate - this combines into one mask on the source grid (therefore resampling much be reimplemented at every registration step)
  Matrix MeshModify::combine_costfunction_weighting(const Matrix & sourceweight, const Matrix & resampledtargetweight){

    Matrix NEW;
    int nrows;
    if(sourceweight.Nrows()>=resampledtargetweight.Nrows()){ NEW=sourceweight; nrows=resampledtargetweight.Nrows();}
    else{ NEW=resampledtargetweight; nrows=sourceweight.Nrows();}

    for (int i=1;i<=NEW.Ncols();i++){
      for (int j=1;j<=nrows;j++){
	NEW(j,i)=(sourceweight(j,i)+resampledtargetweight(j,i))/2;
      }

    }
    
    return NEW;
  }
  
}
