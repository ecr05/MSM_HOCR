/*  DiscreteGroupCostFunction.cpp

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "DiscreteGroupCostFunction.h"
#include <algorithm>

namespace DISCRETEOPTHOCR{
  
  //====================================================================================================================//
  //================================GroupSURFACE CLASS===========================================================================//
  GroupDiscreteCostFunction::GroupDiscreteCostFunction() //: L1data(0)
  {
    _rexp=2;
    _mu=0.1; _kappa=10; 
    _reglambda=0; _debug=false; _MEANANGLE=0; _verbosity=false; _outdir=""; _concat=false;
    _lambdapairs=1;
    VERTICES_PER_SUBJ=0; TRIPLETS_PER_SUBJ=0;

    /* if (!(eng = engOpen("matlab -nodisplay -nodesktop -nojvm"))){
      cout << " eng not open " << endl; throw  DISCRETEOPTHOCRException("DiscreteModel ERROR:: Can't start MATLAB engine"); 
      engSetVisible(eng, 0);
    }
    */
  }   
  /*
  GroupDiscreteCostFunction::~GroupDiscreteCostFunction()
  {

    // free arrays
    
    if(L1data!=0){
      for (int i =0 ; i< num_subjects; i++ ) 
	delete[] L1data[i];

      delete[]L1data;
    }
    
    engClose(eng);
   
    
  }
  */
 
  void GroupDiscreteCostFunction::set_parameters(myparam & ALLPARAMS){
    myparam::iterator it;
    

    it=ALLPARAMS.find("lambda");_reglambda=boost::get<float>(it->second);
    it=ALLPARAMS.find("lambda_pairs");_lambdapairs=boost::get<float>(it->second);
    it=ALLPARAMS.find("set_lambda_pairs");_setpairs=boost::get<bool>(it->second);
    it=ALLPARAMS.find("CPres");_RES=boost::get<int>(it->second);
    it=ALLPARAMS.find("simmeasure");_simmeasure=boost::get<int>(it->second); sim.set_simval(_simmeasure);
    it=ALLPARAMS.find("verbosity");_verbosity=boost::get<bool>(it->second);
    it=ALLPARAMS.find("outdir");_outdir=boost::get<string>(it->second);
    it=ALLPARAMS.find("shearmodulus");_mu=boost::get<float>(it->second);
    it=ALLPARAMS.find("bulkmodulus");_kappa=boost::get<float>(it->second);
    it=ALLPARAMS.find("sigma_in");_sigma=boost::get<float>(it->second);
    it=ALLPARAMS.find("quartet");_quadcost=boost::get<bool>(it->second);

    cout << " _reglambd " << _reglambda << " _lambdapairs " << _lambdapairs << " setpair " << _setpairs << endl;
  }

  void GroupDiscreteCostFunction::set_relations(const boost::shared_ptr<RELATIONS> &CONTROL,const boost::shared_ptr<RELATIONS> &TARG){ 
    _controlrel=CONTROL; _targetrel=TARG; 
    _sourcerel=_controlrel->invert_relationsTR( _CONTROLMESHES[0],_DATAMESHES[0]);  
  }

  void GroupDiscreteCostFunction::get_spacings(){
    Pt CP;
    double dist;
    ColumnVector vMAXmvd(_CONTROLMESHES[0].nvertices());
    SPACINGS.clear();

    for(int n=0; n<num_subjects;n++){
      vMAXmvd=0;
      for (int k=0;k<_CONTROLMESHES[n].nvertices();k++){
   
	CP=_CONTROLMESHES[n].get_coord(k); 
	for (vector<int>::const_iterator it=_CONTROLMESHES[n].nbegin(k);it !=_CONTROLMESHES[n].nend(k);it++){

	  dist= 2*RAD*asin((CP-_CONTROLMESHES[n].get_coord(*it)).norm()/(2*RAD));
	     

	if(dist > vMAXmvd(k+1))
	  vMAXmvd(k+1)=dist;
	
	}
      }
      SPACINGS.push_back(vMAXmvd);
    }

  }

  void GroupDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets,int numQuartets)
  {
    cout << " GroupDiscreteCostFunction::initialize" << numQuartets << endl;

    DiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);

    ///////////////////// identify TEMPLATE grid cell that must be sampled from for each point //////////
    if(!_targetrel.get()) {throw  DISCRETEOPTHOCRException(" GroupDiscreteCostFunction::initialize  must initialise RELATIONS in MODEL before running this function");}
    if(_CONTROLMESHES.size()==0){throw  DISCRETEOPTHOCRException(" GroupDiscreteCostFunction::initialize  must initialise  meshes before running this function");}

    define_template_patches();
    resample_to_template();
    get_spacings();
    resample_patches();

    //   L1data=0;
    if(_quadcost==true && _setpairs==false){ _lambdapairs=(float)m_num_quartets/(float)m_num_pairs; cout << m_num_quartets << " true " << m_num_pairs << " " << (float)m_num_quartets/(float)m_num_pairs <<  endl;}
    cout << " INIT " << _reglambda << " _lambdapairs " << _lambdapairs << " setpair " <<  _setpairs << endl;
   
  }

  void GroupDiscreteCostFunction::define_template_patches(){ /// make patches cover overlap with all between subject neighbours
    int ind,mesh_ID,num;
    vector<vector<int> > TOTALCELLS(m_num_nodes,vector<int> ());
    vector<vector<int> > TEMPLATE_CELLS; 
    vector<int> PTS;
    char filename[1000];
    TEMPLATEPTS.clear();
    cout <<" define_template_patches() " << endl;
    if(!_pairs)  {throw  DISCRETEOPTHOCRException(" GroupDiscreteCostFunction::initialize  must set pairs before function");}

    for (int n=0;n<m_num_nodes;n++){ //(int n=53;n<=53;n++){//
      mesh_ID=floor(n/VERTICES_PER_SUBJ);
      ind=n-mesh_ID*VERTICES_PER_SUBJ;
      Pt p=_CONTROLMESHES[mesh_ID].get_coord(ind);
      // cout << n << " mesh_ID " << mesh_ID << " *MVD_LR " << MVD_LR << " p.X " << p.X <<" p.Y " << p.Y <<" p.Z " << p.Z << endl;
      TEMPLATE_CELLS.push_back(_targetrel->return_cell_group(p,1.75*asin(MVD_LR/RAD)));
      TEMPLATEPTS.push_back(vector<int>());
      //cout << TEMPLATEPTS.size() << " " << TEMPLATEPTS[n].size() << endl;

    }
    // exit(1);
    /////////////// make patch constant across all node pairs
    for(int pair=0;pair<m_num_pairs;pair++){
      for (int i=0;i<2;i++){
	for(int j=0;j<TEMPLATE_CELLS[_pairs[2*pair+i]].size();j++){
	  TOTALCELLS[_pairs[2*pair+i]].push_back(TEMPLATE_CELLS[_pairs[2*pair+i]][j]);
	}
      }
    }

    for (int n=0;n<m_num_nodes;n++){
      sort( TOTALCELLS[n].begin(), TOTALCELLS[n].end() );
      TOTALCELLS[n].erase( unique( TOTALCELLS[n].begin(), TOTALCELLS[n].end() ), TOTALCELLS[n].end() );
      //   cout << n << " TOTALCELLS[n].size() " << TOTALCELLS[n].size() << endl;
      
      for(int i=0;i< TOTALCELLS[n].size();i++){
	PTS=_targetrel->get_cell_members(TOTALCELLS[n][i]);
	TEMPLATEPTS[n].insert(TEMPLATEPTS[n].end(),PTS.begin(),PTS.end());
	//	cout << i << " " << n << "PTS.size " << PTS.size() << "TEMPLATEPTS.size() " << TEMPLATEPTS.size() <<endl;
     
      }

    }
    /*
    for (int n=0;n<10;n++){
      for(int i=0;i<  _TEMPLATE.nvertices();i++)
	_TEMPLATE.set_pvalue(i,0);
      
      for(int i=0;i< TEMPLATEPTS[n].size();i++)
	_TEMPLATE.set_pvalue(TEMPLATEPTS[n][i],1);
     
      sprintf(filename,"TEMPLATE_PATCH2-%d-CPRES%d.func.gii",n,_RES);
      _TEMPLATE.save(filename);
    
    }
    */
  }
      
    

    



  void GroupDiscreteCostFunction::resample_to_template(){

    resampler R; R.set_method("ADAP_BARY"); 
    char filename[1000];
    RESAMPLEDDATA.clear();
    cout << " in resample data " << _RES << endl;
    for(int n=0; n<num_subjects;n++){
      _targetrel->update_RELATIONS(_DATAMESHES[n]);
      RESAMPLEDDATA.push_back(FEAT->get_data_matrix(n));
      R.resampledata(_DATAMESHES[n],_TEMPLATE,RESAMPLEDDATA[n],0.0,_targetrel);
      //sprintf(filename,"dataresampledtotemplate-%d-res%d-iter%d.func.gii",n,_RES,_iter);
      //_TEMPLATE.set_pvalues(RESAMPLEDDATA[n]);_TEMPLATE.save(filename);
      // sprintf(filename,"DATAMESHES-n%d-res%d-iter%d.surf.gii",n,_RES,_iter);_DATAMESHES[n].save(filename);
    }
    
  }

  //rpca version 
  /* double GroupDiscreteCostFunction::computeQuartetCost(int quartet, int labelA, int labelB, int labelC,int labelD){
    double cost=0;
    vector<int> indices(4,0);
    vector<int> rows;
    int mesh_ID, index=0;
    map<int,int> testsubjects;
    //  ofstream out; out.open("rows.txt");
    bool found;

    indices[0]=labelA+_quartets[4*quartet]*_labels.size();
    indices[1]=labelB+_quartets[4*quartet+1]*_labels.size();
    indices[2]=labelC+_quartets[4*quartet+2]*_labels.size();
    indices[3]=labelD+_quartets[4*quartet+3]*_labels.size();

    // if(quartet%10==0) cout << quartet << " computeQuartetCost " << labelA << " " << labelB << " " << labelC << " " << labelD << " num_quartets/num_pairs " << ((float)m_num_quartets/(float)m_num_pairs) <<  endl;
    for(int i=0;i<4;i++){
      mesh_ID=floor(_quartets[4*quartet+i]/VERTICES_PER_SUBJ);
      testsubjects[mesh_ID]=indices[i];
      // cout << i << " meshID "<<  mesh_ID << " index " << indices[i] << " patch size " << PATCHDATA[indices[i]].size() <<  endl;
    }

    for (map<int, float>::const_iterator iter = PATCHDATA[indices[0]].begin(); iter != PATCHDATA[indices[0]].end(); ++iter){
      found=true;
      for(int i=1;i<4;i++){
	if(PATCHDATA[indices[i]].find(iter->first) == PATCHDATA[indices[i]].end()){
	  found=false;
	}
      }
      if(found){ rows.push_back(iter->first);
	//	out << iter->first << endl;
      }
    }

    if(L1data!=0){
      //  cout << " L1data!=0" << L1data << endl;
      for (int i =0 ; i< num_subjects; i++ ) 
	delete[] L1data[i];
      delete[]L1data;
    

    }else {cout << "L1  is NULL " <<  endl;}
   
    L1data=new double*[PATCHDATA[indices[0]].size()];
    // Matrix outdata(rows.size(),num_subjects);outdata=0;
    // L1data = new double[num_subjects];
    //cout << " estimate L1 data " << rows.size() << endl;
    for(int n=0;n<num_subjects;n++){

	L1data[n] = new double[rows.size()];

	for(int r=0;r<rows.size();r++){
	  //// first add all the test subjects then add remaining subjects?
	  if(testsubjects.find(n) != testsubjects.end()){
	    //    if(r==0) cout <<n <<  " subject found " << testsubjects.find(n)->second << endl;
	    L1data[n][r]=PATCHDATA[testsubjects[n]].find(rows[r])->second;
	  }else{
	    L1data[n][r]=RESAMPLEDDATA[n](1,rows[r]);
	  }
	  //  outdata(r+1,n+1)=L1data[n][r] ;
	}
	
    }
    //write_ascii_matrix(outdata,"L1data");
    cost=sim.rpca(L1data,eng,num_subjects,rows.size());
    
    return cost;
  }
  */

  // pca version
  double GroupDiscreteCostFunction::computeQuartetCost(int quartet, int labelA, int labelB, int labelC,int labelD){
    double cost=0;
    vector<int> indices(4,0);
    vector<int> rows;
    int mesh_ID, index=0;
    map<int,int> testsubjects;
    //  ofstream out; out.open("rows.txt");
    bool found;

    indices[0]=labelA+_quartets[4*quartet]*_labels.size();
    indices[1]=labelB+_quartets[4*quartet+1]*_labels.size();
    indices[2]=labelC+_quartets[4*quartet+2]*_labels.size();
    indices[3]=labelD+_quartets[4*quartet+3]*_labels.size();

    // if(quartet%10==0) cout << quartet << " computeQuartetCost " << labelA << " " << labelB << " " << labelC << " " << labelD << " num_quartets/num_pairs " << ((float)m_num_quartets/(float)m_num_pairs) <<  endl;
    for(int i=0;i<4;i++){
      mesh_ID=floor(_quartets[4*quartet+i]/VERTICES_PER_SUBJ);
      testsubjects[mesh_ID]=indices[i];
      // cout << i << " meshID "<<  mesh_ID << " index " << indices[i] << " patch size " << PATCHDATA[indices[i]].size() <<  endl;
    }

    for (map<int, float>::const_iterator iter = PATCHDATA[indices[0]].begin(); iter != PATCHDATA[indices[0]].end(); ++iter){
      found=true;
      for(int i=1;i<4;i++){
	if(PATCHDATA[indices[i]].find(iter->first) == PATCHDATA[indices[i]].end()){
	  found=false;
	}
      }
      if(found){ rows.push_back(iter->first);
	//	out << iter->first << endl;
      }
    }
   
    Matrix Groupdata(rows.size(),num_subjects); Groupdata=0;
    // L1data = new double[num_subjects];
    // cout << " estimate PCA " << rows.size() << endl;
    if(rows.size()>0){
      for(int n=0;n<num_subjects;n++){

	for(int r=0;r<rows.size();r++){
	  //// first add all the test subjects then add remaining subjects?
	  if(testsubjects.find(n) != testsubjects.end()){
	    //    if(r==0) cout <<n <<  " subject found " << testsubjects.find(n)->second << endl;
	    Groupdata(r+1,n+1)=PATCHDATA[testsubjects[n]].find(rows[r])->second;
	  }else{
	    Groupdata(r+1,n+1)=RESAMPLEDDATA[n](1,rows[r]+1);
	  }
	  //  outdata(r+1,n+1)=L1data[n][r] ;
	}
	
      }

      DiagonalMatrix eigenvals;
      SVD(Groupdata,eigenvals);

      for(int i=1;i<=eigenvals.Nrows(); i++){
	cost+=eigenvals(i);
	//if(quartet%10==0 && labelA==0 && labelB==0 && labelC==0 && labelD==0) cout << quartet << " i " << i << " eigenvals(i) " << eigenvals(i) << " " << cost << " rows.size() " << rows.size() <<  endl;
      }
      //cost/=rows.size();
    }else cost=MAXSIM;
    return _lambdapairs*cost;
  }
  double GroupDiscreteCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC){
    double cost,weight=1;
    int meshID;
    Pt v0,v1,v2;
    Pt n_v0,n_v1,n_v2;
    vector<int> id;
    Pt o_v0,o_v1,o_v2;
    int m0,m1,m2;
    //  cout << triplet << " in compute triplet " << TRIPLETS_PER_SUBJ << endl;
    /////////// GET VERTICES AND INDICES FOR THIS MESH
    meshID=floor(triplet/TRIPLETS_PER_SUBJ);
    m0=_triplets[3*triplet]-meshID*VERTICES_PER_SUBJ;
    m1=_triplets[3*triplet+1]-meshID*VERTICES_PER_SUBJ;
    m2=_triplets[3*triplet+2]-meshID*VERTICES_PER_SUBJ;

    // cout << m0 << " " << m1 << " " << m2 << " num_triplets " << m_num_triplets << " meshID " << meshID <<  " triplet 1 " << _triplets[3*triplet] << " " << _triplets[3*triplet+1] << " " << _triplets[3*triplet+2] << " labelA " << labelA << " " << labelB << " " << labelC << "_labels.szie() " << _labels.size() <<  endl;

    v0=_CONTROLMESHES[meshID].get_coord(m0);
    v1=_CONTROLMESHES[meshID].get_coord(m1);
    v2=_CONTROLMESHES[meshID].get_coord(m2);
  
    o_v0=_ORIG.get_coord(m0);
    o_v1=_ORIG.get_coord(m1);
    o_v2=_ORIG.get_coord(m2);
    
    //id.push_back(_triplets[3*triplet]);
    //id.push_back(_triplets[3*triplet+1]);
    //id.push_back(_triplets[3*triplet+2]);

    //cout << "v0.X " << v0.X << " " << v0.Y << " " << v0.Z << " " << "v1.X " << v1.X << " " << v1.Y << " " << v1.Z << " " << "v0.X " << v2.X << " " << v2.Y << " " << v2.Z << " " << endl;
    //cout << "o_v0.X " << o_v0.X << " " << o_v0.Y << " " << o_v0.Z << " " <<"o_v1.X " << o_v1.X << " " << o_v1.Y << " " << o_v1.Z << " " <<"o_v2.X " << o_v2.X << " " << o_v2.Y << " " << o_v2.Z << " " << endl;

    n_v0=(*ROTATIONS)[_triplets[3*triplet]]*_labels[labelA];
    n_v1=(*ROTATIONS)[_triplets[3*triplet+1]]*_labels[labelB];
    n_v2=(*ROTATIONS)[_triplets[3*triplet+2]]*_labels[labelC];

    // cout << " n_v0.X  " << n_v0.X << " " << n_v0.Y << " " << n_v0.Z << " " << n_v1.X << " " <<  n_v1.Y << " " <<  n_v1.Z << " " <<  n_v2.X << " " <<n_v2.Y << " " <<n_v2.Z << " size " << ROTATIONS->size() << endl;
    
    //cout << (*ROTATIONS)[_triplets[3*triplet]] << endl;
    //cout << (*ROTATIONS)[_triplets[3*triplet+1]] << endl;
    //cout << (*ROTATIONS)[_triplets[3*triplet+2]] << endl;

    //////////ESTIMATE STRAIN REGULARISATION /////////////

    Triangle TRI(n_v0,n_v1,n_v2,0);
    Triangle TRI_ORIG(o_v0,o_v1,o_v2,0);  
    Triangle TRI_noDEF(v0,v1,v2,0);

    NEWMESH::Pt norm_new= TRI.normal();
    NEWMESH::Pt norm_old = TRI_noDEF.normal();

   
    if((norm_new|norm_old) < 0){ weight+=1e6 ; }

    cost=weight*_reglambda*MISCMATHS::pow(calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa),_rexp);
    // cout << weight << " cost " << cost <<  " mu " << _mu << " kappa " << _kappa <<" labelA " << labelA << " " << labelB << " " << labelC << "_labels.szie() " << endl;
   
    return cost; // normalise to try and ensure equivalent lambda for each resolution level
  }


  double  GroupDiscreteCostFunction::computePairwiseCost(int pair, int labelA, int labelB){

    ///// first set template sample grid to be within fixed range  or p and q
    double cost;
    int index1,index2;
    Pt newpt,CP;

    cost=MAXSIM;

    ////////////// Get template patch
   
    int ind1=0, ind2=162;
   
    index1=labelA+_pairs[2*pair]*_labels.size();
    index2=labelB +_pairs[2*pair+1]*_labels.size();
    
   
    if(PATCHDATA[index1].size()>0 && PATCHDATA[index1].size()>0){
      cost=sim.corr(PATCHDATA[index1],PATCHDATA[index2]);
    }
    
    if(cost==MAXSIM || cost!=cost){
      cout << " pairwise cost " << cost << " MAXSIM " << MAXSIM << " PATCHDATA[index1].size( " << PATCHDATA[index1].size() << " " << PATCHDATA[index2].size() << " pair " << pair << " labelA " << labelA << " " << labelB << " _pairs[2*pair] " << _pairs[2*pair] << " " << _pairs[2*pair+1] << " index1 " << index1 << " index2 " << index2 <<    endl;

    }
   
    /* if(pair<10){
      //  sprintf(filename,"sourceresampled-%d.surf.gii",ID);
      //newsource.save(filename);
      //sprintf(filename,"sourceresampled-%d.func.gii",ID);
      // if(labelA==0 && labelB==0) cout << pair << " patch data size " << PATCHDATA[index1].size() << " " << PATCHDATA[index2].size() << " cost " << cost <<  endl;

      Matrix DATAMAT(2,_TEMPLATE.nvertices()); DATAMAT=0;
      char filename[1000];
      //newsource.save(filename);
      for (map<int, float>::const_iterator iter = PATCHDATA[index1].begin(); iter != PATCHDATA[index1].end(); ++iter)
	DATAMAT(1,iter->first+1)=iter->second;
       
      for (map<int, float>::const_iterator iter = PATCHDATA[index2].begin(); iter != PATCHDATA[index2].end(); ++iter)
	DATAMAT(2,iter->first+1)=iter->second;
     
      sprintf(filename,"pairresampled-to-template-%d-labels%d-%d-cpres-%d-iter=%d.func.gii",pair,labelA,labelB,_RES,_iter);
      _TEMPLATE.set_pvalues(DATAMAT); _TEMPLATE.save(filename);
    }
    */
    return cost;
  }

  void GroupDiscreteCostFunction::resample_patches(){

    ///// first set template sample grid to be within fixed range  or p and q
    double cost;
    int mesh_ID;
    int node;
    double dist;
    map<int,float> nodedata;
    int labnum;
    int num;
    vector<bool> withinCPrange(m_num_nodes*_labels.size(),false);
       
    cost=MAXSIM;
    PATCHDATA.clear();
    
    ////////////// Get template patch

    
  
    for (int n=0;n<m_num_nodes;n++){
      mesh_ID=floor(n/VERTICES_PER_SUBJ);
      node=n-mesh_ID*VERTICES_PER_SUBJ;
      for(int lab=0;lab<_labels.size();lab++){
	PATCHDATA.push_back(map<int,float>());
    
	dist=(_labels[lab]-_labels[0]).norm();
	if(dist <= 0.5*SPACINGS[mesh_ID](node+1))
	  withinCPrange[lab+n*_labels.size()]=true;
      }
    }

    // boost::thread_group group;
    //for (int i = 0; i <_threads; ++i)
    //  group.create_thread(boost::bind(&GroupDiscreteCostFunction::resampler_worker_function, this,0,m_num_nodes*_labels.size(),withinCPrange));

    resampler_worker_function(0,m_num_nodes*_labels.size(),withinCPrange);
    
  }

  void GroupDiscreteCostFunction::resampler_worker_function(const int &begin, const int &end,const vector<bool> &INrange){
    int mesh_ID;
    int node,num;
    int labnum;
    Pt CP, newpt;

    for (int n=begin;n<end;n++){
      num=floor(n/_labels.size());
      mesh_ID=floor(num/VERTICES_PER_SUBJ);
      node=num-mesh_ID*VERTICES_PER_SUBJ;
      
      labnum=n-num*_labels.size();

      CP=_CONTROLMESHES[mesh_ID].get_coord(node);
      newpt=(*ROTATIONS)[num]*_labels[labnum];   // this rotates the label set over the control point

      if(INrange[labnum+num*_labels.size()]==true){
		
	//// move source mesh according to label rotation and resample on template
	  
	cout << " lab+num*_labels.size() " << labnum+num*_labels.size() << " PATCHDATA size " << PATCHDATA.size() <<  endl;
	PATCHDATA[labnum+num*_labels.size()]=resample_onto_template(node,mesh_ID,newpt,TEMPLATEPTS[num]);
	  
      }else {cout << " (newpt-CP).norm() > 0.5*SPACINGS[meshID](node+1) " << endl; }
      if(n%1000==0) cout << " resample patch: " << n << " " <<  labnum << " " << PATCHDATA[labnum+num*_labels.size()].size() << endl;

    }
  }

  map<int,float>  GroupDiscreteCostFunction::resample_onto_template(const int &node,const int &ID, const Pt &newCP,const vector<int> &PTS){//, const RELATIONS 
    map<int,float>   sampledata;
    map<int,float>    weights;
    double g_weight,dist;
    int sourceind, id0,id1,id2;
    Pt tmp,SP;
    Pt CP0,CP1,CP2;
    Pt newCP0,newCP1,newCP2;
    
    char filename[1000];
    int trind=0;
    bool _print=false;

    
    for(int i=0;i<  PTS.size();i++){
      sampledata[PTS[i]]=0;
      weights[PTS[i]]=0;
   }

    
    //cout << node << "_sourcerel " << _sourcerel.Ncols()<< " CP norm " <<  newCP.norm() << endl;

    for (vector<int>::const_iterator j =_CONTROLMESHES[ID].tIDbegin(node); j!=_CONTROLMESHES[ID].tIDend(node); j++){
      trind++;
      NEWMESH::Triangle tri=_CONTROLMESHES[ID].get_triangle(*j);
      CP0 = tri.get_vertex_coord(0);  CP1 = tri.get_vertex_coord(1); CP2 = tri.get_vertex_coord(2);
      id0=tri.get_vertex_no(0);  id1=tri.get_vertex_no(1);  id2=tri.get_vertex_no(2);
   
      if(id0==node){
	newCP0=newCP; newCP1=CP1;newCP2=CP2;
      }
      else if(id1==node){
	newCP0=CP0; newCP1=newCP;newCP2=CP2;

      }else if(id2==node){
	newCP0=CP0; newCP1=CP1;newCP2=newCP;

      }else{cout << node << " triangle node not found " << " id0" << id0 << " " << id1 << " " << id2 << endl; exit(1);}
   
      for(int i=1;i<=_sourcerel.Nrows(*j+1);i++){
	sourceind=_sourcerel(i,*j+1);
	///// move source using barycentric interpolation
	SP=_DATAMESHES[ID].get_coord(sourceind-1);
	projectPoint(CP0,CP1,CP2,SP);
	tmp= barycentric(CP0,CP1,CP2,SP,newCP0,newCP1,newCP2);
     
	tmp.normalize(); ////////// ASSUMING SPHERE"S ORIGIN IS 0,0,0 (which  we can project from the face to the surface of the sphere this way
	tmp=tmp*RAD;
	

	/*if(_debug && _print){
	  if(ID==0 &&  trind==1 && i==1 ){ //  || i%10 ==0
	    cout << node << " " << i <<  " ID " << ID << " _sourcerel.Nrows(k) " << _sourcerel.Nrows(*j+1) << "sourcept.X " <<SP.X << " " << SP.Y << " " << SP.Z << " sourcept.norm() " << SP.norm() <<  "CP0.norm() " <<  CP0.norm() <<  " CP1.norm() " << CP1.norm() << " " << CP2.norm() <<   endl;	  
	  }
	}
	*/
	////////////// resample onto template using gaussian 
	vector< pair<float,int> > NEIGHBOURS=_targetrel->update_RELATIONS_for_ind(tmp);
	//	if(_RES==2 && node==0) cout << " neighbours.size() " << NEIGHBOURS.size() << " " << _sigma << endl;
	for (int s=0;s<NEIGHBOURS.size();s++){
	  dist=(tmp-_TEMPLATE.get_coord(NEIGHBOURS[s].second)).norm();
	  g_weight=exp(-(dist*dist)/(2*_sigma*_sigma));
	  //	  if(_RES==3 && s<5 && node==0) cout << node << " ID " << ID << " g_weight " << g_weight << " dist " << dist << endl;

	  sampledata[NEIGHBOURS[s].second]+=g_weight*FEAT->get_data_val(1,sourceind,ID);
	  weights[NEIGHBOURS[s].second]+=g_weight;
	}
      }
    }

    for(int i=0;i<PTS.size();i++){
      if(weights[PTS[i]]==0)
	sampledata[PTS[i]]=RESAMPLEDDATA[ID](1,PTS[i]+1);
      else
	sampledata[PTS[i]]/=weights[PTS[i]];

    }
    
    
    
    

    return sampledata;

  }
};



 
