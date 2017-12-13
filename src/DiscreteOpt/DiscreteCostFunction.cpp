/*  DiscreteCostFunction.cpp

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "DiscreteCostFunction.h"
#include <algorithm>

namespace DISCRETEOPT{

  //================================BASECLASS===========================================================================//
  DiscreteCostFunction::DiscreteCostFunction()// m_is_memory(false), m_is_memory_triplets(false)
    : m_num_nodes(0), m_num_labels(0), m_num_pairs(0), m_num_triplets(0), m_num_quartets(0), unarycosts(0), paircosts(0),pairweights(0), _pairs(0), _triplets(0), _quartets(0)//, tripletcosts(0), quartetcosts(0), pairweights(0), tripletweights(0),
	  {
	  }

  //====================================================================================================================//
  DiscreteCostFunction::~DiscreteCostFunction()
  {
		if(unarycosts) { delete[] unarycosts; unarycosts = 0; } // replace pointers with special pointers?
		if(paircosts) { delete[] paircosts; paircosts = 0; }
		// if(quartetcosts) { delete[] quartetcosts; quartetcosts = 0; }
		if(pairweights) { delete[] pairweights; pairweights = 0; }
		//  if(tripletcosts) { delete[] tripletcosts; tripletcosts = 0; }
		// if(tripletweights) { delete[] tripletweights; tripletweights = 0; }
  }

  //====================================================================================================================//
  void DiscreteCostFunction::initialize(int numNodes,int numLabels, int numPairs, int numTriplets, int numQuartets)
  {

		if (m_num_nodes != numNodes || m_num_labels != numLabels)
		  {
		if (unarycosts){ delete [] unarycosts; }
		unarycosts = new double[numNodes*numLabels];
		  }

		if (m_num_pairs != numPairs || m_num_labels != numLabels)
		 {

		   if (paircosts) delete [] paircosts;
		   if (pairweights) delete [] pairweights;
		   paircosts = new double[numPairs*numLabels*numLabels];
		   pairweights = new double[numPairs*numLabels*numLabels]; // not used? delete?
		 }

		m_num_nodes = numNodes;
		m_num_labels = numLabels;
		m_num_pairs = numPairs;
		m_num_triplets = numTriplets;
		m_num_quartets = numQuartets;


		std::fill(unarycosts,unarycosts+m_num_labels*m_num_nodes,0.0f);
		std::fill(paircosts,paircosts+m_num_labels*m_num_labels*m_num_pairs,0.0f);
		std::fill(pairweights,pairweights+m_num_labels*m_num_labels*m_num_pairs,0.0f);

  }

  //====================================================================================================================//
  void DiscreteCostFunction::reset()
  {
		if(unarycosts) std::fill(unarycosts,unarycosts+m_num_labels*m_num_nodes,0.0f);
		if(paircosts) std::fill(paircosts,paircosts+m_num_labels*m_num_labels*m_num_pairs,0.0f);
		if(pairweights) std::fill(pairweights,pairweights+m_num_pairs,0.0f);

  }

  //====================================================================================================================//
  void DiscreteCostFunction::initPairWeights(const double *weights) // for msm distortion weights we modify pairwise costs rather than adding a weight
  {
		if(!pairweights) pairweights = new double[m_num_pairs];

		if(weights)
		  {
		for(int p = 0; p < m_num_pairs; ++p)
		  {
			pairweights[p] = weights[p];
		  }
		  }
		else
		  {
		std::fill(pairweights,pairweights+m_num_pairs,1.0f);
		  }
  }

  //====================================================================================================================//
  double DiscreteCostFunction::evaluateTotalCostSumZeroLabeling()
  {
		int label = 0;
		double cost_sum_unary = 0.0f;
		double cost_sum_pairwise = 0.0f;
		double cost_sum_triplet = 0.0f;
		double cost_sum_quartet = 0.0f;

		for(int i = 0; i < m_num_nodes; ++i){
		  cost_sum_unary += computeUnaryCost(i,label);
		  }
		for(int p = 0; p < m_num_pairs; ++p)
		  {
		cost_sum_pairwise += computePairwiseCost(p,label,label);
		  }
		for(int t = 0; t < m_num_triplets; ++t)
		  {
		cost_sum_triplet += computeTripletCost(t,label,label,label);
		  }

	   for(int q = 0; q < m_num_quartets; ++q)
		  {

		cost_sum_quartet += computeQuartetCost(q,label,label,label,label);
		  }

	   if(_debug) cout << m_num_quartets << " cost_sum_unary " << cost_sum_unary << " cost_sum_pairwise " << cost_sum_pairwise  << " cost_sum_triplet " << cost_sum_triplet << " cost_sum_quartet " << cost_sum_quartet << " total " <<  cost_sum_unary + cost_sum_pairwise + cost_sum_triplet + cost_sum_quartet << " m_num_triplets " << m_num_triplets <<   endl;

		return cost_sum_unary + cost_sum_pairwise + cost_sum_triplet + cost_sum_quartet;
  }

  //====================================================================================================================//
  double DiscreteCostFunction::evaluateTotalCostSum(const int *labeling, const int *pairs, const int *triplets, const int *quartets)
  {
		double cost_sum_unary = 0.0f;
		double cost_sum_pairwise = 0.0f;
		double cost_sum_triplet = 0.0f;
		double cost_sum_quartet = 0.0f;
		for(int i = 0; i < m_num_nodes; ++i){
		  cost_sum_unary += computeUnaryCost(i,labeling[i]);
		}

		for(int p = 0; p < m_num_pairs; ++p)
		  {
		const int index = p*2;
		cost_sum_pairwise += computePairwiseCost(p,labeling[pairs[index]],labeling[pairs[index+1]]);
		  }

		for(int t = 0; t < m_num_triplets; ++t)
		  {
		const int index = t*3;
		cost_sum_triplet += computeTripletCost(t,labeling[triplets[index]],labeling[triplets[index+1]],labeling[triplets[index+2]]);
		//cout << t << " " << cost_sum_triplet << " abeling[triplets[index]] " << labeling[triplets[index]] << " " << labeling[triplets[index+1]] << " " << labeling[triplets[index+2]] <<  endl;
		if (cost_sum_triplet!=cost_sum_triplet) cout << " cost_sum_triplet!=cost_sum_triplet " << endl;
		 }

		for(int q = 0; q < m_num_quartets; ++q)
		  {
		const int index = q*4;
		cost_sum_quartet += computeQuartetCost(q,labeling[quartets[index]],labeling[quartets[index+1]],labeling[quartets[index+2]],labeling[quartets[index+3]]);
		  }

		if(_verbosity)  cout << m_num_quartets << " cost_sum_unary " << cost_sum_unary << " cost_sum_pairwise " << cost_sum_pairwise  << " cost_sum_triplet " << cost_sum_triplet << " cost_sum_quartet " << cost_sum_quartet <<" total " <<  cost_sum_unary + cost_sum_pairwise + cost_sum_triplet + cost_sum_quartet << " m_num_triplets 2 " << m_num_triplets <<  endl;
		return cost_sum_unary + cost_sum_pairwise + cost_sum_triplet + cost_sum_quartet;
  }

  //====================================================================================================================//
  double DiscreteCostFunction::evaluateUnaryCostSum(const int *labeling)
  {
		double cost_sum_unary = 0.0f;
		for(int i = 0; i < m_num_nodes; ++i)
		  cost_sum_unary += computeUnaryCost(i,labeling[i]);// unarycosts[labeling[i]*m_num_nodes+i];
		return cost_sum_unary;
  }

  //====================================================================================================================//
  double DiscreteCostFunction::evaluatePairwiseCostSum(const int *labeling, const int *pairs)
  {
		double cost_sum_pairwise = 0.0f;
		for(int p = 0; p < m_num_pairs; ++p)
		  {
		const int index = p*2;
		cost_sum_pairwise += computePairwiseCost(p,labeling[pairs[index]],labeling[pairs[index+1]]);
		  }
		return cost_sum_pairwise;
  }

  //====================================================================================================================//
  double DiscreteCostFunction::evaluateTripletCostSum(const int *labeling, const int *triplets)
  {
		double cost_sum_triplet = 0.0f;
		for(int t = 0; t < m_num_triplets; ++t)
		  {
		const int index = t*3;
		cost_sum_triplet += computeTripletCost(t,labeling[triplets[index]],labeling[triplets[index+1]],labeling[triplets[index+2]]);
		  }
		return cost_sum_triplet;
  }

  //====================================================================================================================//
  double DiscreteCostFunction::evaluateQuartetCostSum(const int *labeling, const int *quartets)
  {
		double cost_sum_quartet = 0.0f;
		for(int q = 0; q < m_num_quartets; ++q)
		  {
		const int index = q*3;
		cost_sum_quartet += computeQuartetCost(q,labeling[quartets[index]],labeling[quartets[index+1]],labeling[quartets[index+2]],labeling[quartets[index+3]]);
		  }
		return cost_sum_quartet;
  }
 //================================GENERIC SURFACE CLASS===========================================================================//
  SRegDiscreteCostFunction::SRegDiscreteCostFunction() :_simmeasure(2)
  {
    _reglambda=0; _debug=false; _MEANANGLE=0; _verbosity=false; _outdir=""; _concat=false; _sigma=1;
  }

 //====================================================================================================================//
  void SRegDiscreteCostFunction::initialize(int numNodes,  int numLabels, int numPairs, int numTriplets, int numQuartets)
  {

		if(_TARGET.nvertices()==0 || _SOURCE.nvertices()==0) {    throw  DISCRETEOPTHOCRException("CostFunction::You must supply source and target meshes.");}
		// _AREALDIST.ReSize(_SOURCE.nvertices()); _AREALDIST=0.0;
		if(_HIGHREScfweight.Ncols()!=_SOURCE.nvertices()) _HIGHREScfweight.ReSize(1, _SOURCE.nvertices());
		if(_HIGHREScfweight.Nrows()!=1 && _HIGHREScfweight.Nrows()!=FEAT->get_dim()){throw  DISCRETEOPTHOCRException("DiscreteModel ERROR:: costfunction weighting has dimensions incompatible with data");}
		DiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);


  }

  void SRegDiscreteCostFunction::set_parameters(myparam & ALLPARAMS){

		myparam::iterator it;
		it=ALLPARAMS.find("lambda");_reglambda=boost::get<float>(it->second);
		it=ALLPARAMS.find("range");_controlptrange=boost::get<float>(it->second);
		it=ALLPARAMS.find("CPres");_RES=boost::get<int>(it->second);
		it=ALLPARAMS.find("anatres");_aRES=boost::get<int>(it->second);
		it=ALLPARAMS.find("simmeasure");_simmeasure=boost::get<int>(it->second); sim.set_simval(_simmeasure);
		it=ALLPARAMS.find("verbosity");_verbosity=boost::get<bool>(it->second);
		it=ALLPARAMS.find("outdir");_outdir=boost::get<string>(it->second);
		it=ALLPARAMS.find("regularisermode");_rmode=boost::get<int>(it->second);
		it=ALLPARAMS.find("sigma_in");_sigma=boost::get<float>(it->second);
		it=ALLPARAMS.find("numthreads");_threads=boost::get<int>(it->second);

  }

  void SRegDiscreteCostFunction::initialize_regulariser(){

		if(_aSOURCE.nvertices()>0 && _rmode>=3){


		  anatMVD=2*asin(Calculate_MaxVD(_aICO)/RAD);
		  _anatrel.Initialize(_aICO,_TARGEThi,anatMVD);

		}

  }

  newmesh SRegDiscreteCostFunction::project_anatomical(){
		newmesh TRANS;
		newmesh _aICOtrans=_aICO;
		barycentric_mesh_interpolation(_aICOtrans,_ORIG,_SOURCE);
		TRANS=projectmesh(_aICOtrans,(_TARGEThi),(_aTARGET));

		return TRANS;
  }


  /// this is for debug

  void SRegDiscreteCostFunction::reset_anatomical(const string &outdir, const int &iter){


		newmesh STRAINSmesh;
		ColumnVector strainstmp(_aSOURCE.ntriangles());
		double perc=0;
		_iter=iter;
		if(_aSOURCE.nvertices()>0 ) {

		  _aSOURCEtrans=project_anatomical();
		  MAXstrain=0;

		  for (int i=0;i<_aSOURCE.ntriangles();i++){
		strainstmp(i+1)=calculate_triangular_strain(i,_aSOURCE,_aSOURCEtrans,_mu,_kappa);
		if(strainstmp(i+1)>MAXstrain)
		  MAXstrain=strainstmp(i+1);
		  }

		  if(MAXstrain>1e-8){
		Histogram strainHist(strainstmp,256);

		perc=strainHist.getPercentile(0.95);
		  }

		  if(iter==1) strain95=perc;

		}

  }
  //====================================================================================================================//
  //================================AFFINE SURFACE CLASS===========================================================================//
  //AffineSRegDiscreteCostFunction::AffineSRegDiscreteCostFunction()
  // : _iters(50)
  //{
  //}

  //====================================================================================================================//
  //================================Non Linear SURFACE CLASS===========================================================================//
  NonLinearSRegDiscreteCostFunction::NonLinearSRegDiscreteCostFunction()
  {
    _expscaling = 1;
    _k_exp = 2.0;
    _maxdist=4;  _rexp=2;
     _kNN=5; _rmode=1;
     _mu=0.4; _kappa=1.6;
  }

  void NonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets,int numQuartets)
  {
		///// initialise CF weightings
		SRegDiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);

  }

  void NonLinearSRegDiscreteCostFunction::set_parameters(myparam & ALLPARAMS){
    myparam::iterator it;
		it=ALLPARAMS.find("maxdistortion");_maxdist=boost::get<float>(it->second);
		it=ALLPARAMS.find("exponent");_rexp=boost::get<float>(it->second);
		it=ALLPARAMS.find("weight");_dweight=boost::get<bool>(it->second);
		it=ALLPARAMS.find("anorm");_anorm=boost::get<bool>(it->second);
		it=ALLPARAMS.find("scaling");_expscaling=boost::get<float>(it->second);
		it=ALLPARAMS.find("shearmodulus");_mu=boost::get<float>(it->second);
		it=ALLPARAMS.find("bulkmodulus");_kappa=boost::get<float>(it->second);
		it=ALLPARAMS.find("kexponent");_k_exp=boost::get<float>(it->second);

		it=ALLPARAMS.find("kNN");_kNN=boost::get<int>(it->second);
		it=ALLPARAMS.find("pottsthreshold");_pottsthreshold=boost::get<float>(it->second);

		SRegDiscreteCostFunction::set_parameters(ALLPARAMS);


  }


  double NonLinearSRegDiscreteCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC){

		vector<int> id;
		Pt normal;
		map<int,Pt>  v;
		Pt v0,v1,v2;
		double cost=0;
		double weight=1.0;
		double distortion=1;
		double likelihood=0; ///////////// WILL ONLY BE ESTIMATED FOR HOCR method
		char filename[1000];

		map<int,Pt>  transformed_points;

		if(triplet==0 &&  _debug){ sumlikelihood=0; sumregcost=0; }
			   //  sprintf(filename,"CPgridintrip-%d.surf",_iter);
			//	_CPgrid.save(filename) ;
			  //   sprintf(filename,"SOURCEintrip-%d.surf",_iter);

			//	_SOURCE.save("filename");
				//_SOURCE.save("SOURCEintrip.func"); }
		if(!_triplets){throw  DISCRETEOPTHOCRException("DiscreteModel ERROR:: must run settripletss() prior to computTripleCost ");}


		v0=_CPgrid.get_coord(_triplets[3*triplet]);
		v1=_CPgrid.get_coord(_triplets[3*triplet+1]);
		v2=_CPgrid.get_coord(_triplets[3*triplet+2]);




		id.push_back(_triplets[3*triplet]);
		id.push_back(_triplets[3*triplet+1]);
		id.push_back(_triplets[3*triplet+2]);


		v[id[0]]=(*ROTATIONS)[id[0]]*_labels[labelA];
		v[id[1]]=(*ROTATIONS)[id[1]]*_labels[labelB];
		v[id[2]]=(*ROTATIONS)[id[2]]*_labels[labelC];

		likelihood=triplet_likelihood(triplet,id[0],id[1],id[2],v[id[0]],v[id[1]],v[id[2]]);


		Triangle TRI(v[id[0]],v[id[1]],v[id[2]],0);

		Triangle TRI_noDEF(v0,v1,v2,0);
		NEWMESH::Pt norm_new= TRI.normal();
		NEWMESH::Pt norm_old = TRI_noDEF.normal();
		/// only estimate cost if it doesn't cause folding
		if((norm_new|norm_old) < 0){ cost=1;  weight+=1e6 ;
			//cout << triplet << " triplet fold " <<  labelA  << " " << labelB << " " <<  labelC << endl;
			}
		else{


		  switch(_rmode) {
		  case 2:
		{
		  vector<double> deformed_angles, orig_angles;
		  Pt o_v0,o_v1,o_v2;
		  double val=0,tweight=1;
		  double diff=0;
		  Triangle TRI_ORIG;

		  o_v0=_ORIG.get_coord(_triplets[3*triplet]);
		  o_v1=_ORIG.get_coord(_triplets[3*triplet+1]);
		  o_v2=_ORIG.get_coord(_triplets[3*triplet+2]);
		  TRI_ORIG.set(o_v0,o_v1,o_v2,0);
		  deformed_angles=TRI.get_angles();
		  orig_angles=TRI_ORIG.get_angles();
		  distortion=log2(TRI.area()/TRI_ORIG.area());

		  if(_dweight && distortion!=1){
			tweight=exp(abs(distortion-1));

		  }else{tweight=1;}

		  for (int i=0;i<3;i++){

			val=deformed_angles[i]-orig_angles[i];
			diff+=val*val;

		  }

		  cost=tweight*(sqrt(diff));
		  if(_anorm){ weight*=(1/_MEANANGLE);}

		  break;
		}
		  case 3:{

		 ////
		  Pt o_v0,o_v1,o_v2;
		  Triangle TRI_ORIG;

		  o_v0=_ORIG.get_coord(_triplets[3*triplet]);
		  o_v1=_ORIG.get_coord(_triplets[3*triplet+1]);
		  o_v2=_ORIG.get_coord(_triplets[3*triplet+2]);
		  TRI_ORIG.set(o_v0,o_v1,o_v2,0);


		  cost=calculate_triangular_strain(TRI_ORIG, TRI, _mu, _kappa, boost::shared_ptr<ColumnVector>(), _k_exp);
		  // if(write) cout << triplet << " " << labelA <<  " " << labelB << " " <<  labelC << " " << cost << " (v0- v[id[0]]).norm() " << (v0- v[id[0]]).norm() << endl;


		break;
		  }
		  case 4:
		{
		  vector<double> deformed_angles, orig_angles;
		  double val=0;
		  double diff=0;
		  double tweight=1,distortion;
		  Triangle TRIorig,TRItrans;
		  map<int,bool> moved2;
		  map<int,Pt>  transformed_points;

		  for (unsigned int n=0;n<NEARESTFACES[triplet].size();n++){

			// for each face with each CP triangle
			TRIorig = _aSOURCE.get_triangle(NEARESTFACES[triplet][n]);
			TRItrans=deform_anatomy(triplet,n,v,moved2,transformed_points);

			deformed_angles=TRItrans.get_angles();
			orig_angles=TRIorig.get_angles();
			distortion=log2(TRItrans.area()/TRIorig.area());
			if(_dweight && distortion!=1){
			  tweight=exp(abs(distortion-1));

			}else{tweight=1;}



			for (int i=0;i<3;i++){

			  val=deformed_angles[i]-orig_angles[i];
			  diff+=val*val;

			}

			diff=tweight*sqrt(diff);
			cost+=diff;



		  }

		  cost=cost/NEARESTFACES[triplet].size();

		  break;
		}
		  case 5:
		{

		  Triangle TRIorig,TRItrans;
		  map<int,bool> moved2;
		  map<int,Pt>  transformed_points;


		  for (unsigned int n=0;n<NEARESTFACES[triplet].size();n++){
		  // for each face with each CP triangle
			TRIorig = _aSOURCE.get_triangle(NEARESTFACES[triplet][n]);
			TRItrans=deform_anatomy(triplet,n,v,moved2,transformed_points);

			cost+=calculate_triangular_strain(TRIorig, TRItrans, _mu, _kappa, boost::shared_ptr<ColumnVector>(), _k_exp);


		  }


		  cost=cost/NEARESTFACES[triplet].size();



		  break;
		}
		  default:
		throw  DISCRETEOPTHOCRException("DiscreteModel computeTripletCost regoption does not exist");



		  }
		}


		if(abs(cost)<1e-8){ cost=0;  }
		if(_debug){ sumlikelihood+=likelihood;
		sumregcost+=weight*_reglambda*MISCMATHS::pow(cost,_rexp);
		}
		cost=likelihood+ weight*_reglambda*MISCMATHS::pow(cost,_rexp);


		return cost; // normalise to try and ensure equivalent lambda for each resolution level
  }

  double  NonLinearSRegDiscreteCostFunction::computePairwiseCost(int pair, int labelA, int labelB){

		Matrix fullrotA,fullrotB;
		Matrix R1,R2,R_diff,logR,R_diff2,logR2;
		double theta_MVD=2*asin(MVDmax/(2*RAD));
		double cost=0;

		if(!_pairs){throw  DISCRETEOPTHOCRException("DiscreteModel ERROR:: must run setPairs() prior to computePairwiseCost ");}

		R1=estimate_rotation_matrix(_CPgrid.get_coord(_pairs[2*pair]),(*ROTATIONS)[_pairs[2*pair]]*_labels[labelA]); /// with additional rotation estimated for this label
		R2=estimate_rotation_matrix(_CPgrid.get_coord(_pairs[2*pair+1]),(*ROTATIONS)[_pairs[2*pair+1]]*_labels[labelB]);

		R_diff=R1.t()*R2;

		double theta=acos((R_diff.Trace()-1)/2);
		double weight=1;
		double totalcost=0.0;

		/// check for flipping ////
		Pt v0,v1;
		v0=_CPgrid.get_coord(_pairs[2*pair]);
		v1=_CPgrid.get_coord(_pairs[2*pair+1]);


		_CPgrid.set_coord(_pairs[2*pair],(*ROTATIONS)[_pairs[2*pair]]*_labels[labelA]);
		_CPgrid.set_coord(_pairs[2*pair+1],(*ROTATIONS)[_pairs[2*pair+1]]*_labels[labelB]);



		/////////////////////////////////////////


		if(fabs(1-(R_diff.Trace()-1)/2) > EPSILON){


		  for(vector<int>::const_iterator j=_CPgrid.tIDbegin(_pairs[2*pair]); j!=_CPgrid.tIDend(_pairs[2*pair]);j++){
		  Triangle TRI_trans=_CPgrid.get_triangle(*j);
		  Triangle TRI_orig=_oCPgrid.get_triangle(*j);
		  Pt norm_orig=TRI_orig.normal();
		  Pt norm_trans=TRI_trans.normal();
		  if((norm_orig|norm_trans) < 0){  cost=1;  weight+=1e6 ;  break;}
		  // cout <<pair <<  " (norm_orig|norm_trans) < 0 " << endl;
		  }


		  cost=(sqrt(2)*theta)/(theta_MVD);


		  if(_rexp==1)
		totalcost=weight*_reglambda*cost;
		  else
		totalcost=weight*_reglambda*MISCMATHS::pow(cost,_rexp) ;
		}

		_CPgrid.set_coord(_pairs[2*pair],v0);
		_CPgrid.set_coord(_pairs[2*pair+1],v1);

		return totalcost;
  }

  void  NonLinearSRegDiscreteCostFunction::computePairwiseCosts(const int *pairs)
  {
		int pair=0;
		Matrix fullrotA,fullrotB,fullrotA_2,fullrotB_2;
		Matrix R1_B,R2_B,R_diff,logR,R_diff2,logR2;

		while(pair<m_num_pairs){


		  for(unsigned int j=0;j<_labels.size();j++){

			for(unsigned int k=0;k<_labels.size();k++){

		  paircosts[pair*m_num_labels*m_num_labels+k*m_num_labels+j]=computePairwiseCost(pair,j,k);
		  pairweights[pair*m_num_labels*m_num_labels+k*m_num_labels+j]=1;


		}

		  }
		  // coudl use wcosts to downweight points that have no intensity gradient
		  pair++;

		}


  }

 void  NonLinearSRegDiscreteCostFunction::computeUnaryCosts(){

		/// for each control point resample data into blocks each influencing a single control point
		for( int j=0;j<m_num_labels;j++){

		  //// calculates similarity
		  for (int k=0;k<_CPgrid.nvertices();k++){
		unarycosts[j*m_num_nodes+k]=computeUnaryCost(k,j);
		  }

		}

  }


  Triangle NonLinearSRegDiscreteCostFunction::deform_anatomy(const int &trip,const int &n, map<int,Pt>  &vertex, map<int,bool> & moved, map<int,Pt>  &transformed){
		Triangle TRItrans;


		Pt newPt;
		int tindex;
		map<int,double> weight;

		for(int i=0;i<3;i++){ // for each point in face


		  tindex=_aSOURCE.get_triangle_vertexID(NEARESTFACES[trip][n],i);
		  if(moved.find(tindex) == moved.end()) {
		moved[tindex]=true;
		newPt.X=0; newPt.Y=0;newPt.Z=0;

		for (std::map<int,double>::iterator it= _ANATbaryweights[tindex].begin(); it!= _ANATbaryweights[tindex].end(); ++it){

		  newPt+=vertex[it->first]*it->second;
		}



		weight=R.get_barycentric_weight_for_ind(tindex,anatMVD,newPt,_TARGEThi, _anatrel);

		newPt.X=0;  newPt.Y=0;  newPt.Z=0;
		for (std::map<int,double>::iterator it= weight.begin(); it!= weight.end(); ++it){

		  newPt+=_aTARGET.get_coord(it->first)*it->second;

		}
		transformed[tindex]=newPt;

		TRItrans.set_vertice(i, newPt);

		  }else
		TRItrans.set_vertice(i, transformed[tindex]);
		}

		return TRItrans;
  }

  void NonLinearSRegDiscreteCostFunction::resample_weights(){ ///////// TAKE DATA TERM WEIGHTING AS MAXIMUM OF WEIGHTS ACROSS ALL DIMENSIONS //////////////

	    R.set_method("ADAP_BARY");
		double maxweight;
		AbsoluteWeights.ReSize(_SOURCE.nvertices()); AbsoluteWeights=0;

		for (int k=1;k<=_SOURCE.nvertices();k++){

		  maxweight=0;
		  for (int j=1;j<=_HIGHREScfweight.Nrows();j++){
		if(_HIGHREScfweight(j,k)>maxweight)
		  maxweight=_HIGHREScfweight(j,k);
		  }
		  AbsoluteWeights(k)=maxweight;
		}

		////////////////////////// resample to CPGRID /////////////////////////////////////
		R.resampledata(_SOURCE,_CPgrid,AbsoluteWeights,2*asin(MVDmax/RAD),_controlrel);



  }


  void NonLinearSRegDiscreteCostFunction::get_target_data(const int& node,const Matrix &PtROTATOR){//, const RELATIONS &cpneighbours
		Pt v0,v1,v2,CP,tmp;
		int n0,n1,n2;
		int ind,sourceind;

		reset_target_data(node);
		ind=-1;


		for(unsigned int i=0;i<_sourceinrange[node-1].size();i++){
		  sourceind=_sourceinrange[node-1][i];

		  tmp=PtROTATOR*_SOURCE.get_coord(sourceind-1);


		  get_data_index(node-1,sourceind-1,ind);

		  vector< pair<float,int> > NEIGHBOURS=_targetrel->update_RELATIONS_for_ind(tmp);

		  if(!_TARGET.return_closest_points(NEIGHBOURS[0].second,v0,v1,v2,tmp,n0,n1,n2)){

		for (unsigned int s=1;s<NEIGHBOURS.size();s++){

		  if(_TARGET.return_closest_points(NEIGHBOURS[s].second,v0,v1,v2,tmp,n0,n1,n2)){
			break;
		  }
		  if(s== NEIGHBOURS.size()-1){ cout << ind << " can't find face intersections " << " " << NEIGHBOURS[0].second <<  endl;
			n0=n1=n2=NEIGHBOURS[0].second;
		  }
		}

		  }

		  resample_target_data(v0,v1,v2,tmp,n0,n1,n2,ind);


		}

  }

  void NonLinearSRegDiscreteCostFunction::get_target_data(const int& triplet,const Pt& new_CP0,const Pt& new_CP1,const Pt& new_CP2){
		Pt v0,v1,v2,SP,tmp2,tmp;
		Pt CP0,CP1,CP2;
		int n0,n1,n2;
		int ind,sourceind;

		reset_target_data(triplet);
		 ind=-1;

		CP0=_CPgrid.get_coord(_triplets[3*(triplet-1)]);
		CP1=_CPgrid.get_coord(_triplets[3*(triplet-1)+1]);
		CP2=_CPgrid.get_coord(_triplets[3*(triplet-1)+2]);


		for(unsigned int i=0;i<_sourceinrange[triplet-1].size();i++){
		  sourceind=_sourceinrange[triplet-1][i];
		  SP=_SOURCE.get_coord(sourceind-1);
		  projectPoint(CP0,CP1,CP2,SP);
		  tmp= barycentric(CP0,CP1,CP2,SP,new_CP0,new_CP1,new_CP2);

		  tmp.normalize(); ////////// ASSUMING SPHERE"S ORIGIN IS 0,0,0 (which  we can project from the face to the surface of the sphere this way
		  tmp=tmp*RAD;
		  get_data_index(triplet-1,sourceind-1,ind);

		  vector< pair<float,int> > NEIGHBOURS=_targetrel->update_RELATIONS_for_ind(tmp);


		  if(!_TARGET.return_closest_points(NEIGHBOURS[0].second,v0,v1,v2,tmp,n0,n1,n2)){

		for (unsigned int s=1;s<NEIGHBOURS.size();s++){

		  if(_TARGET.return_closest_points(NEIGHBOURS[s].second,v0,v1,v2,tmp,n0,n1,n2)){
			break;
		  }
		  if(s==NEIGHBOURS.size()-1){ cout << ind << " can't find face intersections " << " " << NEIGHBOURS[0].second <<  endl;
			n0=n1=n2=NEIGHBOURS[0].second;
		  }
		}

		  }


		  resample_target_data(v0,v1,v2,tmp,n0,n1,n2,ind);



		}

  }




  //====================================================================================================================//
  //================================ UNIVARIATE Non Linear SURFACE CLASS===========================================================================//
  void UnivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets,int numQuartets)
  {
		NonLinearSRegDiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);
		_sourcedata.clear(); _sourcedata.resize(_CPgrid.nvertices(), vector<double> ());
		_sourceinrange.clear(); _sourceinrange.resize(_CPgrid.nvertices(), vector<int> ());
		_targetdata.clear(); _targetdata.resize(_CPgrid.nvertices(), vector<double> ());
		_weights.clear(); _weights.resize(_CPgrid.nvertices(), vector<double> ());

  }


  void UnivariateNonLinearSRegDiscreteCostFunction::get_source_data(){ //const RELATIONS &cpneighbours
		// get source data
		int sourceind,ind;//index;
		double distance,weight;
		ind=-1;

		for (unsigned int i=0;i<_sourcedata.size();i++) _sourcedata[i].clear();
		for (int k=1;k<=_sourcerel.Ncols();k++){


		for(int i=1;i<=_sourcerel.Nrows(k);i++){
		sourceind=_sourcerel(i,k);
		if(within_controlpt_range(k-1,sourceind-1,distance)){

		  get_data_index(k-1,sourceind-1,ind);
		  _sourceinrange[k-1].push_back(sourceind);
		  for(int d=1;d<=FEAT->get_dim();d++){
			  _sourcedata[ind].push_back(FEAT->get_input_val(d,sourceind));
			 // cout << i << " " << k << " " << ind << " " << sourceind << " FEAT->get_input_val(d,sourceind) " << FEAT->get_input_val(d,sourceind) << " _sourcerel.Nrows(k) " << _sourcerel.Nrows(k) <<  endl;
			  if(_HIGHREScfweight.Nrows()>=d) weight=_HIGHREScfweight(d,sourceind); else weight=1;

			  _weights[ind].push_back(weight);
			}
		  }
		}


		}
		resample_weights();


  }


  double UnivariateNonLinearSRegDiscreteCostFunction::computeUnaryCost(int node, int label){
		Matrix ROT;
		Pt newpt,CP;


		newpt=(*ROTATIONS)[node]*_labels[label];   // this rotates the label set over the control point
		CP=_CPgrid.get_coord(node);
		ROT=estimate_rotation_matrix(CP,newpt); // this Matrix will rotate control point k to the label
		get_target_data(node+1,ROT);

		double cost=AbsoluteWeights(node+1)*sim.get_sim_for_min(_sourcedata[node], _targetdata[node],_weights[node]);// weights similarity by cost function weighting


		if(_simmeasure==1 || _simmeasure==2)
		  return cost;
		else
		  return -cost;
  }

  void UnivariateNonLinearSRegDiscreteCostFunction::resample_target_data(const Pt &v0, const Pt &v1, const Pt &v2, const Pt &tmp, const int &n0, const int &n1, const int &n2, int &index){
		for(int d=1;d<=FEAT->get_dim();d++){
		  double targint=barycentric(v0,v1,v2,tmp,FEAT->get_ref_val(d,n0 + 1),FEAT->get_ref_val(d,n1 + 1),FEAT->get_ref_val(d,n2 + 1));
		  _targetdata[index].push_back(targint);
		}


  }


 //====================================================================================================================//
  //================================Multivariate Non Linear SURFACE CLASS===========================================================================//

  void MultivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets, int numQuartets)
  {
		NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels,numPairs,numTriplets,numQuartets);
		_sourcedata.clear(); _sourcedata.resize(_SOURCE.nvertices(), vector<double> ());
		_sourceinrange.clear(); _sourceinrange.resize(_CPgrid.nvertices(), vector<int> ());
		_targetdata.clear(); _targetdata.resize(_SOURCE.nvertices(), vector<double> ());
		_weights.clear(); _weights.resize(_SOURCE.nvertices(), vector<double> ());

  }

  void MultivariateNonLinearSRegDiscreteCostFunction::get_source_data(){
		// get source data
		double distance,weight;
		int ind=0;
		for (unsigned int i=0;i<_sourcedata.size();i++) _sourcedata[i].clear();

		for (int k=0;k<_SOURCE.nvertices();k++){

		  for (int j=1;j<=_controlrel->Nrows(k+1);j++){
		ind=(*_controlrel)(j,k+1);
		if(within_controlpt_range(ind-1,k,distance)){
		  _sourceinrange[ind-1].push_back(k+1);
		}
		  }

		  for(int d=1;d<=FEAT->get_dim();d++){
		_sourcedata[k].push_back(FEAT->get_input_val(d,k+1));

		if(_HIGHREScfweight.Nrows()>=d) weight=_HIGHREScfweight(d,k+1); else weight=1;
		_weights[k].push_back(weight);
		  }

		}

		resample_weights();

  }


  double MultivariateNonLinearSRegDiscreteCostFunction::computeUnaryCost(int node, int label){
		Matrix ROT;
		Pt newpt,CP;
		double cost	=0;

		newpt=(*ROTATIONS)[node]*_labels[label];   // this rotates the label set over the control point
		CP=_CPgrid.get_coord(node);
		ROT=estimate_rotation_matrix(CP,newpt); // this Matrix will rotate control point k to the label
		get_target_data(node+1,ROT);

		for (unsigned int i=0;i< _sourceinrange[node].size();i++){

		  cost+=sim.get_sim_for_min(_sourcedata[_sourceinrange[node][i]-1], _targetdata[_sourceinrange[node][i]-1],_weights[_sourceinrange[node][i]-1]);// weights similarity by cost function weighting

		}

		if(_sourceinrange[node].size()>0) cost/=_sourceinrange[node].size();

		if(_simmeasure==1 || _simmeasure==2)
		  return AbsoluteWeights(node+1)*cost;
		else
		  return -AbsoluteWeights(node+1)*cost;

  }


  void MultivariateNonLinearSRegDiscreteCostFunction::resample_target_data(const Pt &v0, const Pt &v1, const Pt &v2, const Pt &tmp, const int &n0, const int &n1, const int &n2, int &index){

	  	  barycentric2(v0,v1,v2,tmp,n0+1,n1+1,n2+1,*FEAT->get_reference_data(),_targetdata[index]);

  }


  void HOUnivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets,int numQuartets)
  {
		NonLinearSRegDiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);
		_sourcedata.clear(); _sourcedata.resize(_CPgrid.ntriangles(), vector<double> ());
		_sourceinrange.clear(); _sourceinrange.resize(_CPgrid.ntriangles(), vector<int> ());
		_targetdata.clear(); _targetdata.resize(_CPgrid.ntriangles(), vector<double> ());
		_weights.clear(); _weights.resize(_CPgrid.ntriangles(), vector<double> ());


  }


  void HOUnivariateNonLinearSRegDiscreteCostFunction::get_source_data(){ //const RELATIONS &cpneighbours
		// get source data
		int sourceind;//index;
		double weight;
		R.set_method("ADAP_BARY");

		for (unsigned int i=0;i<_sourcedata.size();i++) _sourcedata[i].clear();

		for (int k=1;k<=_sourcerel.Ncols();k++){


		  for(int i=1;i<=_sourcerel.Nrows(k);i++){
		sourceind=_sourcerel(i,k);


		_sourceinrange[k-1].push_back(sourceind);

		for(int d=1;d<=FEAT->get_dim();d++){
		  _sourcedata[k-1].push_back(FEAT->get_input_val(d,sourceind));

		  if(_HIGHREScfweight.Nrows()>=d) weight=_HIGHREScfweight(d,sourceind); else weight=1;

		  _weights[k-1].push_back(weight);
		}
		  }

		}



		resample_weights();
  }

  double HOUnivariateNonLinearSRegDiscreteCostFunction::triplet_likelihood(const int &triplet,const int & CP_id0,const int & CP_id1,const int &CP_id2,const Pt & CP_def0,const Pt &CP_def1,const Pt &CP_def2){
		get_target_data(triplet+1,CP_def0,CP_def1,CP_def2);

		double cost=(1.0/3)*(AbsoluteWeights(CP_id0+1)+AbsoluteWeights(CP_id1+1)+AbsoluteWeights(CP_id2+1))*sim.get_sim_for_min(_sourcedata[triplet], _targetdata[triplet],_weights[triplet]);// weights similarity by cost function weighting


		if(_simmeasure==1 || _simmeasure==2)
		  return cost;
		else
		  return -cost;

  }


  void HOMultivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes,int numLabels, int numPairs, int numTriplets, int numQuartets)
    {
		  NonLinearSRegDiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);
		  _sourcedata.clear(); _sourcedata.resize(_SOURCE.nvertices(), vector<double> ());
		  _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.ntriangles(), vector<int> ());
		  _targetdata.clear(); _targetdata.resize(_SOURCE.nvertices(), vector<double> ());
		  _weights.clear(); _weights.resize(_SOURCE.nvertices(), vector<double> ());

    }

    void HOMultivariateNonLinearSRegDiscreteCostFunction::get_source_data(){
			// get source data
			double weight;

			for (unsigned int i=0;i<_sourcedata.size();i++) _sourcedata[i].clear();

			for (int k=1;k<=_sourcerel.Ncols();k++){
				for(int i=1;i<=_sourcerel.Nrows(k);i++){
					int sourceind=_sourcerel(i,k);
					_sourceinrange[k-1].push_back(sourceind);
				}
			}


			for (int k=0;k<_SOURCE.nvertices();k++){
				for(int d=1;d<=FEAT->get_dim();d++){
					_sourcedata[k].push_back(FEAT->get_input_val(d,k+1));

					if(_HIGHREScfweight.Nrows()>=d) weight=_HIGHREScfweight(d,k+1); else weight=1;
					_weights[k].push_back(weight);
				}
			}

			resample_weights();

    }

  double HOMultivariateNonLinearSRegDiscreteCostFunction::triplet_likelihood(const int &triplet,const int & CP_id0,const int & CP_id1,const int &CP_id2,const Pt & CP_def0,const Pt &CP_def1,const Pt &CP_def2){
		get_target_data(triplet+1,CP_def0,CP_def1,CP_def2);
		double cost=0;

		for (unsigned int i=0;i< _sourceinrange[triplet].size();i++){

		  cost+=sim.get_sim_for_min(_sourcedata[_sourceinrange[triplet][i]-1], _targetdata[_sourceinrange[triplet][i]-1],_weights[_sourceinrange[triplet][i]-1]);
		}
		if(_sourceinrange[triplet].size()>0) cost/=_sourceinrange[triplet].size();
		if(_simmeasure==1 || _simmeasure==2)
		  return (1.0/3)*(AbsoluteWeights(CP_id0+1)+AbsoluteWeights(CP_id1+1)+AbsoluteWeights(CP_id2+1))*(cost);
		else
		  return -(1.0/3)*(AbsoluteWeights(CP_id0+1)+AbsoluteWeights(CP_id1+1)+AbsoluteWeights(CP_id2+1))*(cost);

		return cost;
  }


  RegularisationDiscreteCostFunction::RegularisationDiscreteCostFunction(double shear, double bulk, double exp)
    {  _rexp=exp;	   _mu=shear; _kappa=bulk;}

  void RegularisationDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets, int numQuartets)
  {
	  	DiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);
  }

  double RegularisationDiscreteCostFunction::computeTripletRegularisation(int triplet, int labelA, int labelB, int labelC){
		vector<int> id;
		Pt normal;
		vector<Pt> v(_CPgrid.nvertices(), Pt());
		Pt v0,v1,v2;
		double W;




		boost::shared_ptr<ColumnVector> strains; /// for debugging
		//if(triplet==0){ sumlikelihood=0; sumregcost=0;}

		if(_debug) strains=boost::shared_ptr<ColumnVector>(new ColumnVector (2));
		v0=_CPgrid.get_triangle_vertex(triplet,0);
		v1=_CPgrid.get_triangle_vertex(triplet,1);
		v2=_CPgrid.get_triangle_vertex(triplet,2);

		id.push_back(_CPgrid.get_triangle_vertexID(triplet,0));
		id.push_back(_CPgrid.get_triangle_vertexID(triplet,1));
		id.push_back(_CPgrid.get_triangle_vertexID(triplet,2));


		v[id[0]]=(*ROTATIONS)[id[0]]*_labels[labelA];
		v[id[1]]=(*ROTATIONS)[id[1]]*_labels[labelB];
		v[id[2]]=(*ROTATIONS)[id[2]]*_labels[labelC];

		_SOURCE.set_coord(id[0],v[id[0]]);
		_SOURCE.set_coord(id[1],v[id[1]]);
		_SOURCE.set_coord(id[2],v[id[2]]);

		//cost=calculate_triangular_strain(TRI_ORIG, TRI, _mu, _kappa, boost::shared_ptr<ColumnVector>(), _k_exp);

		W=calculate_triangular_strain(triplet,_SOURCE,_CPgrid,_mu,_kappa,strains);
//W=calculate_triangular_strain(triplet,_ANAT,_SPHERE,_mu,_kappa,strains);
		W=MISCMATHS::pow(W,_rexp);


		_CPgrid.set_coord(id[0],v0);
		_CPgrid.set_coord(id[1],v1);
		_CPgrid.set_coord(id[2],v2);



		//sumregcost+=W;

		return W; // normalise to try and ensure equivalent lambda for each resolution level
  }

}
