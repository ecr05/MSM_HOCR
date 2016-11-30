/*  DiscreteGroupCostFunction.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2014 University of Oxford  */

/*  CCOPYRIGHT  */
#include <vector>
#include "DiscreteCostFunction.h"
//#include <boost/thread.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
//#include "tbb/tick_count.h"
//using namespace tbb;

typedef map<string, boost::variant<int, string, double, float, bool> > myparam;
typedef pair<string, boost::variant<int, string, double, float, bool> > parameterPair;


namespace DISCRETEOPTHOCR{

  class GroupDiscreteCostFunction: public SRegDiscreteCostFunction
  {
    
  public:
    //
    // Constructor.
    //
    GroupDiscreteCostFunction();
    // ~GroupDiscreteCostFunction();
    //// SET UP //////////////////

    void set_parameters(myparam & ); 
    /// neighbouhood info
    void set_relations(const boost::shared_ptr<RELATIONS> &CONTROL,const boost::shared_ptr<RELATIONS> &TARG);
    void get_spacings();
    void set_meshes(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source, const NEWMESH::newmesh & GRID, int num=1){_TEMPLATE=target; 

      _ORIG=source;
      num_subjects=num;
      VERTICES_PER_SUBJ=GRID.nvertices();
      TRIPLETS_PER_SUBJ=GRID.ntriangles();
      MVD_LR=Calculate_MVD(GRID);
      cout << " costfunction set meshes " << num << endl;
      for (int i=0;i<num;i++){
	_DATAMESHES.push_back(source);
	_CONTROLMESHES.push_back(GRID);
      }
    } 

    //// INITIALISATION //////////////////

    virtual void initialize(int numNodes, int numLabels, int numPairs, int numTriplets = 0,int numQuartets = 0); // quartets not used yet so no code for them below
    void define_template_patches();
    void resample_to_template(); /// resamples all data ontotemplate

    /* inline void set_matlab_path(string s){_matlabpath=s;   
      cout << " in group set matlab path " << _matlabpath.c_str() << endl;
      char filename[1000];
      sprintf(filename,"addpath('%s')",_matlabpath.c_str());
      engEvalString(eng,filename);
     
   }
    */
    //////////////// Updates ///////////////////////
    void reset_source(const NEWMESH::newmesh & source, int num=0){_DATAMESHES[num]=source;}
    virtual void reset_CPgrid(const NEWMESH::newmesh & grid,int num=0){_CONTROLMESHES[num]=grid;}


    double computeQuartetCost(int quartet, int labelA, int labelB, int labelC,int labelD);

    double computeTripletCost(int triplet, int labelA, int labelB, int labelC);
    double computePairwiseCost(int pair, int labelA, int labelB);

    void resample_patches();
    void resampler_worker_function(const int &, const int &,const vector<bool> &);
    map<int,float>  resample_onto_template(const int &,const int &,const Pt &, const vector<int> &);

  protected:

    //////////////////////// MESHES //////////////////////
    vector<NEWMESH::newmesh> _DATAMESHES; // TARGET MESH
    vector<NEWMESH::newmesh> _CONTROLMESHES; // TARGET MESH
    newmesh _TEMPLATE;
    /////////////////////// NEIGHBOURHOOD INFO ////////////////////
    vector<boost::shared_ptr<RELATIONS> > _CONTROLRELATIONS; // hold control grid neighbours of each source vertex
    vector<boost::shared_ptr<RELATIONS> > _TEMPLATERELATIONS; // hold target grid neighbours of each source vertex
    vector<vector<int> > TEMPLATEPTS; 
    vector<ColumnVector> SPACINGS;

    //////////////////////DATA//////////////////
    vector<Matrix>  RESAMPLEDDATA;
    vector<map<int,float> > PATCHDATA;
    //double **L1data;

    //Engine *eng; /// matlab engine

    double MVD_LR;

    float _lambdapairs;
    int num_subjects;
    int TRIPLETS_PER_SUBJ;
    int VERTICES_PER_SUBJ;

    bool _quadcost;
    bool _setpairs;

  };
  

}
