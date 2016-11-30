/*  DiscreteGroupModel.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2014 University of Oxford  */

/*  CCOPYRIGHT  */

#pragma once

#include "DiscreteGroupCostFunction.h"
#include "DiscreteModel.h"
#include <vector>

 

namespace DISCRETEOPTHOCR{
 
 

  class GroupDiscreteModel: public SRegDiscreteModel
  {     
  protected:
    vector<NEWMESH::newmesh> m_DATAMESHES; // TARGET MESH
    vector<NEWMESH::newmesh> m_CONTROLMESHES; // TARGET MESH
   
    NEWMESH::newmesh m_TEMPLATE;   /// template data grid
    NEWMESH::newmesh m_TEMPLATE_LR; /// template control grid

    /////////////////////// NEIGHBOURHOOD INFO ////////////////////
    boost::shared_ptr<RELATIONS>  m_TEMPLATE_LR_ALL_RELATIONS; /// hold neighbourhood between control grids and to LR target
    vector<vector<vector<int> > >  BETWEEN_SUBJECT_PAIRS;
    vector<int> subjects;   
  
    int m_num_subjects;
    int control_grid_size;
  public:
    /**
     * Constructor.
     */
    GroupDiscreteModel(){};  
    ~GroupDiscreteModel(){};

    GroupDiscreteModel(myparam & P){  
      m_CPres=2; m_SGres=4; m_simmeasure=2; m_multivariate=false; m_verbosity=false; m_outdir="";m_debug=false;  
      set_parameters(P);
      cout << " create cf group 1 " << endl;
      costfct=boost::shared_ptr<GroupDiscreteCostFunction>(new GroupDiscreteCostFunction());
      m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS());
      m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS ()); 
      costfct->set_parameters(P);
     
    };

    void set_meshspace(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source, const int num=1){  m_TEMPLATE=target;
      m_DATAMESHES.clear();
      m_num_subjects=num;
      for (int i=0;i<num;i++)
	m_DATAMESHES.push_back(source);
    }

    void reset_meshspace(const NEWMESH::newmesh & source, int num=0){  
      char filename[1000];
      cout << " reset meshspace" << endl;
      m_DATAMESHES[num]=source; 
      //sprintf(filename,"SOURCE_set_up_it%d_res%d_num%d.surf.gii",m_iter, m_CPres,num);  m_DATAMESHES[num].save(filename);
      if(costfct.get()){costfct->reset_source(source,num);}
      else{throw  DISCRETEOPTHOCRException("Discrete Model:: You cannot reset the source mesh withou initialising the discrete costfunction");} 
    }
    
    void reset_CPgrid(const NEWMESH::newmesh & grid, int num=0){  
       m_CONTROLMESHES[num]=grid;
    }
  
    void warp_CPgrid(NEWMESH::newmesh & START, NEWMESH::newmesh & END, int num=0){ barycentric_mesh_interpolation(m_CONTROLMESHES[num],START,END); unfold(m_CONTROLMESHES[num]);	}


    void applyLabeling(){applyLabeling(labeling);} ;
    void applyLabeling(int *discreteLabeling);

    void initialize_quartets();
    void initialize_pairs();

    void estimate_pairs();
    void estimate_triplets();
    void estimate_quartets();
    void estimate_combinations(int, int*);
 

    void Initialize(const newmesh &);
    void Initialize(){};
    void get_rotations(vector<Matrix>  &);
    void setupCostFunction();

    //// functions for keeping track of between mesh relationships
    void get_between_subject_pairs();
    void resample_to_template();

    NEWMESH::newmesh get_CPgrid(int num=0){return m_CONTROLMESHES[num];}

  };


}
