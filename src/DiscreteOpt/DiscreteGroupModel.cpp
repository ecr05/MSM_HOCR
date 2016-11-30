/*  DiscreteGroupModel.cpp

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "DiscreteGroupModel.h"
#include <algorithm>
#include <iostream>


namespace DISCRETEOPTHOCR{


  ///////////// BETWEEN MESH SIMILARITY IS EITHER PAIRIWISE (simple model) OR QUADS

  void GroupDiscreteModel::initialize_pairs(){
  
    m_num_pairs=0;
    for(int n=0;n<m_num_subjects;n++){
      for(int i=0;i<m_CONTROLMESHES[n].nvertices();i++){ 
	for(int n2=0;n2<m_num_subjects;n2++){
	  if(n2>n){
	    m_num_pairs++;
	  }
	}
      }
    }

   

  }

  void GroupDiscreteModel::initialize_quartets(){

    m_num_quartets=0; /// number of combinations *nvertices
    subjects.clear();
    int k=4;
    
    for (int i = 0; i < m_num_subjects; i++) ////m_num_subjects
      subjects.push_back(i);

    cout << " initialise combinations " << subjects.size() <<  endl;
    
    
    if(k>m_num_subjects) throw  DISCRETEOPTHOCRException("GroupDiscreteModel::initialize_quartets() not enough examples!");
    for(int i=0;i<m_CONTROLMESHES[0].nvertices();i++){ 

      do
	{
	  m_num_quartets++;
	}
      while(next_combination(subjects.begin(),subjects.begin() + k,subjects.end()));
    }
  }

  void GroupDiscreteModel::estimate_pairs(){

    /////// NEED to find closest vertices BETWEEN IMAGE PAIRS USING TEMPLATE SPACE
    int pair=0;
    pairs  = new int [m_num_pairs*2];
    // ofstream out1; out1.open("pairsorig");

    for(int n=0;n<m_num_subjects;n++){
      for(int i=0;i<m_CONTROLMESHES[n].nvertices();i++){ 
	for(int n2=0;n2<m_num_subjects;n2++){
	  if(n2>n){
	    int node_ids[2] = { i+n*control_grid_size, BETWEEN_SUBJECT_PAIRS[n][n2][i]};
	     //  cout << " node_ids[0] " << node_ids[0] << " " << node_ids[1] << endl;
	     if(2*pair+1> m_num_pairs*2) { cout <<2*pair+1 << " m_num_pairs  " <<  m_num_pairs << " " << m_num_pairs*2 << endl; throw  DISCRETEOPTHOCRException("GroupDiscreteModel::estimate_pairs() index exceeds array size!");}
	     if((node_ids[1] >= (n-1)*m_CONTROLMESHES[n].nvertices()) && (node_ids[1] < n*m_CONTROLMESHES[n].nvertices())){throw  DISCRETEOPTHOCRException("GroupDiscreteModel::estimate_pairs() between mesh index pairs are from the same mesh!");}
	     if((node_ids[0] > m_num_subjects*control_grid_size) || (node_ids[1] > m_num_subjects*control_grid_size)){throw  DISCRETEOPTHOCRException("GroupDiscreteModel::estimate_pairs() indices exceed the total number of control points!");}
	     sort_nodeids(node_ids,2);
	     pairs[2*pair]=node_ids[0];
	     pairs[2*pair+1]=node_ids[1];
	     //   cout << " pair " << pair << "  " << node_ids[0] << " " << node_ids[1] << endl;
	     pair++;
	   }
	 }
   
      }
    }
   

  }

  void GroupDiscreteModel::estimate_quartets(){
    
    quartets  = new int [m_num_quartets*4];
    cout << " estimate quartets " << endl;
    estimate_combinations(4,quartets);
   
  }

  void GroupDiscreteModel::estimate_combinations(int k, int *combinations){
   
    int ind=0;
    
    for(int i=0;i<m_CONTROLMESHES[0].nvertices();i++){ 

      do
       {
	 vector<int> node_ids; 
	 node_ids.push_back(i+subjects[0]*control_grid_size);
	 for (int n = 1; n < k; ++n)
	   {
	     node_ids.push_back(BETWEEN_SUBJECT_PAIRS[subjects[0]][subjects[n]][i]);
	    
	   }
	 std::sort(node_ids.begin(), node_ids.end());
	 std::cout << i << " " ;
	 for (int n = 0; n < k; ++n){
	   combinations[k*ind+n]=node_ids[n];
	   std::cout << node_ids[n] << " " ;
	 }
	 std::cout << "\n";
	 ind++;
       }
      while(next_combination(subjects.begin(),subjects.begin() + k,subjects.end()));
    }
    cout << " ind " << ind << " " <<  m_num_quartets << endl;
  }
 
  
  void GroupDiscreteModel::estimate_triplets(){ ////// for regularisation

    m_num_triplets=m_num_subjects*m_CONTROLMESHES[0].ntriangles();
   
    triplets  = new int [m_num_triplets*3];
    cout << " in estimate triplets " << m_num_triplets << endl;
    for(int n=0;n<m_num_subjects;n++){
      for(int i=0;i<m_CONTROLMESHES[n].ntriangles();i++){
	int node_ids[3] = {m_CONTROLMESHES[n].get_triangle_vertexID(i,0)+n*m_CONTROLMESHES[n].nvertices(),m_CONTROLMESHES[n].get_triangle_vertexID(i,1)+n*m_CONTROLMESHES[n].nvertices(), m_CONTROLMESHES[n].get_triangle_vertexID(i,2)+n*m_CONTROLMESHES[n].nvertices()};
	sort_nodeids(node_ids,3);
	if(3*i+n*m_CONTROLMESHES[n].ntriangles()+2 > m_num_triplets*3 ) {throw  DISCRETEOPTHOCRException("GroupDiscreteModel::estimate_triplets index exceeds array size");}
	if((node_ids[0] >= m_num_subjects*m_CONTROLMESHES[n].nvertices() ) || (node_ids[1] >= m_num_subjects*m_CONTROLMESHES[n].nvertices() )|| (node_ids[2] >= m_num_subjects*m_CONTROLMESHES[n].nvertices() )) {throw  DISCRETEOPTHOCRException("GroupDiscreteModel::estimate_triplets node index exceeds mesh size");}
       	//cout <<i<< " " << n <<  " node_ids[0] " << node_ids[0] <<  " 3*triplet+2 " <<  3*i+n*3*m_CONTROLMESHES[n].ntriangles()+2 << " node_ids[0] +n*m_CONTROLMESHES[n].nvertices() " << node_ids[0]  << " " << node_ids[1] << " " << node_ids[2]  << " m_CONTROLMESHES[n].ntriangles() " << m_CONTROLMESHES[n].ntriangles() << " " <<  n*m_CONTROLMESHES[n].ntriangles() << endl;
	triplets[3*i+n*3*m_CONTROLMESHES[n].ntriangles()]=node_ids[0];
	triplets[3*i+n*3*m_CONTROLMESHES[n].ntriangles()+1]=node_ids[1];
	triplets[3*i+n*3*m_CONTROLMESHES[n].ntriangles() +2]=node_ids[2];
	
      }
    }
    // exit(1);
  }

  void GroupDiscreteModel::get_between_subject_pairs(){
    cout << " GroupDiscreteModel::get_between_subject_pairs " << endl;
    vector<boost::shared_ptr<Mpoint> > ALLPOINTS=m_CONTROLMESHES[0].get_points();
    int ind,mesh_ID;
    int found=0;
    vector<int> indices(control_grid_size,-1);   //// Between subjects pairs is a per subject matrix of subjectsxvertex indicators of closest vertices for each node
    double ang=2*asin(MVD/RAD),angtmp;
    BETWEEN_SUBJECT_PAIRS.clear();

    cout << BETWEEN_SUBJECT_PAIRS.size() << endl;
    for(int n=0; n<m_num_subjects;n++){
      vector<vector< int > > mat(m_num_subjects,indices);
      BETWEEN_SUBJECT_PAIRS.push_back(mat);
    }

    //////// create one giant mesh for all subject control grids	
    for (int n=1;n<m_num_subjects;n++){
      vector<boost::shared_ptr<Mpoint> > tmppoints=m_CONTROLMESHES[n].get_points();
      ALLPOINTS.insert(ALLPOINTS.end(),tmppoints.begin(),tmppoints.end());
    }
    cout << " ALLPOINTSsize " << ALLPOINTS.size() << " " << m_TEMPLATE_LR.nvertices() << endl;
      
    m_TEMPLATE_LR_ALL_RELATIONS=boost::shared_ptr<RELATIONS >(new RELATIONS (m_TEMPLATE_LR,ALLPOINTS,ang)); 

    //// for all subjects and all vertices find the closest between mesh neighbours 
    for(int n=0;n<m_num_subjects;n++){      
      angtmp=ang;
      int i=1;
      while (i<=control_grid_size) {
	found=0;
	m_TEMPLATE_LR_ALL_RELATIONS->update_RELATIONS_for_ind(i,m_CONTROLMESHES[n], angtmp);
	//	cout << i << " m_TEMPLATE_LR_ALL_RELATIONS->Nrows(i) " << m_TEMPLATE_LR_ALL_RELATIONS->Nrows(i) << " " <<  m_TEMPLATE_LR_ALL_RELATIONS->Ncols() << endl;
	for (int r=1;r<=m_TEMPLATE_LR_ALL_RELATIONS->Nrows(i);r++){
	  ind=(*m_TEMPLATE_LR_ALL_RELATIONS)(r,i)-1;
	  mesh_ID=floor(ind/control_grid_size);
	  if(mesh_ID!=n && BETWEEN_SUBJECT_PAIRS[n][mesh_ID][i-1]==-1) { BETWEEN_SUBJECT_PAIRS[n][mesh_ID][i-1]=ind;
	    found++;
	    //  cout << r << " " << ind << " mesh ID " <<  mesh_ID << " " << n << " " <<  ALLPOINTS.size() << " found " << found <<" m_num_subjects-1 " << m_num_subjects-1 <<  endl;

	    if(found==m_num_subjects-1) break;
	  }
	  

	} 
	if(found!=m_num_subjects-1){ angtmp=2*ang; //// check that neighbours have been found between all other meshes
	  cout << " GroupDiscreteModel::get_between_subject_pairs() not between subjects neighbour pairs found for subject " << n << endl; break;}
	else {i++;}
      }
      
    }

   
    // exit(1);
  }
    
 void GroupDiscreteModel::get_rotations(vector<Matrix> &ROT){ /// rotates sampling grid to each control point
    Matrix R,R2,R_n,R_n2,R_diff;
    Pt ci,CP,CP_n;
    ROT.clear();
    ci=m_samplinggrid.get_coord(m_centroid);


    for (int n=0;n<m_num_subjects;n++){

      for (int k=0;k<m_CONTROLMESHES[n].nvertices();k++){
	CP=m_CONTROLMESHES[n].get_coord(k);
	//// rotate sampling points to overlap with control point
	R=estimate_rotation_matrix(ci,CP);    
	ROT.push_back(R);
    
      }
    }
    
  }

  void GroupDiscreteModel::Initialize(const newmesh &CONTROLGRID){
    
    char filename[1000];

    m_num_nodes=CONTROLGRID.nvertices()*m_num_subjects;
    control_grid_size=CONTROLGRID.nvertices();
    m_CONTROLMESHES.clear();
    m_TEMPLATE_LR=CONTROLGRID;
   
    if(m_debug){

      sprintf(filename,"TEMPLATE_res%d.surf.gii",m_CPres); m_TEMPLATE.save(filename);
      sprintf(filename,"TEMPLATELR_res%d.surf.gii",m_CPres); m_TEMPLATE_LR.save(filename);

    }

    for(int n=0; n<m_num_subjects;n++)
      m_CONTROLMESHES.push_back(CONTROLGRID);
    
    cout << " GroupDiscreteModel::Initialize " << m_CONTROLMESHES.size() << " " << m_DATAMESHES.size() << endl;
    /////////////////////////// INITIALIZE REGULARISATION TRIPLETS /////////////////////
    initialize_pairs();
    if(_estquartet && m_num_subjects>=4) initialize_quartets();
    estimate_triplets();
    
   
    /////////////////// CALCULATE FIXED GRID SPACINGS////////
    MVD=Calculate_MVD(CONTROLGRID);
    m_maxs_dist=0.4*Calculate_MaxVD(CONTROLGRID);
    /// INITIALIAZE LABEL GRID/////////
    initLabeling();
    Initialize_sampling_grid();  

    //////////// INITIALIZE NEIGHBOURHOODS ////////////////
    m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS (m_DATAMESHES[0],CONTROLGRID,2*asin(MVD/RAD))); /// as source mesh moves with control grid the relationships are constant for all meshes
    m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS (m_DATAMESHES[0],m_TEMPLATE,2*asin(MVD/RAD)));
    // sprintf(filename,"%sSOURCE1-%d.surf",m_outdir.c_str(),m_iter); m_SOURCE.save(filename);
   
    cout << " update control RELATIONS " << endl;
    m_cp_neighbourhood->update_RELATIONS(m_DATAMESHES[0]);
    cout << " set meshes cf " << endl;
    if(costfct.get()){
      costfct->set_meshes(m_TEMPLATE,m_DATAMESHES[0],CONTROLGRID,m_DATAMESHES.size());

      costfct->set_relations(m_cp_neighbourhood,m_inputrel);
    }
    m_iter=1;
    
  }


 
  void GroupDiscreteModel::setupCostFunction()
  {

    char filename[1000];
    cout << " GroupDiscreteModel::setupCostFunction " << endl;
    resetLabeling(); // initialise label array to zero

    ////// use geodesic distances ////////
    for(int n=0;n<m_num_subjects;n++){
      costfct->reset_CPgrid(m_CONTROLMESHES[n],n);
      // sprintf(filename,"Controlgrid_set_up_it%d_res%d_num%d.surf.gii",m_iter, m_CPres,n);
      //m_CONTROLMESHES[n].save(filename);
    }
    costfct->set_iter(m_iter);
    //////////////////////// GET BETWEEN MESH GROUPINGS //////////////
    get_between_subject_pairs();
    
    estimate_pairs(); 
    if(_estquartet) estimate_quartets();    

    ////////////////// GET LABEL SPACE //////////
    get_rotations(m_ROT);  
    if(m_iter%2==0) m_labels=m_samples;
    else m_labels=m_barycentres;

    m_num_labels=m_labels.size();

    costfct->set_labels(m_labels,m_ROT);
    ///////////////// INIT ////////
    if(m_verbosity)cout << " initialize cost function 2" << m_iter << " m_num_triplets " << m_num_triplets << endl;
    costfct->setPairs(pairs);
    costfct->setTriplets(triplets);
    if(_estquartet) costfct->setQuartets(quartets);
    
    costfct->initialize(m_num_nodes,m_num_labels,m_num_pairs,m_num_triplets,m_num_quartets);
    
    
   
   
    
      
   if(m_verbosity) cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs <<  endl;
    m_iter++;
  }


  void GroupDiscreteModel::applyLabeling(int *dlabels){
    cout << " GroupDiscreteModel applylabelling " << endl;
    if(dlabels){
      for (int n=0;n<m_num_subjects;n++){
	cout << " " << m_CONTROLMESHES[n].nvertices() << " m_ROT.size() " <<  m_ROT.size() << endl;
	for (int i=0;i<m_CONTROLMESHES[n].nvertices();i++){
	  //// rotate sampling points to overlap with control point
	  //// transform grid point to new position given by label
	  Pt newpt=m_ROT[i +n*m_CONTROLMESHES[n].nvertices()]*m_labels[dlabels[i +n*m_CONTROLMESHES[n].nvertices()]];

	  m_CONTROLMESHES[n].set_coord(i,newpt);
	}

      }
    }
   
  } 


}
 
