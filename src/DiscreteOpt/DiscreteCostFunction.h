/*  DiscreteCostFunction.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include <vector>
#include "similarities.h"

#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>

typedef map<string, boost::variant<int, string, double, float, bool> > myparam;
typedef pair<string, boost::variant<int, string, double, float, bool> > parameterPair;


namespace DISCRETEOPT{

  class DiscreteCostFunction
  {
  public:
    //typedef std::vector< std::vector<float> >		PairwiseMemory;
    //typedef std::vector< std::vector<float> >		TripletMemory;
    /**
     * Constructor.
     */
    DiscreteCostFunction();


    /**
     * Destructor.
     */
    virtual ~DiscreteCostFunction();

    /**
     * Resets all costs to zero.
     */
    void reset();

    /**
     * Sets the pairs index buffer.
     */
    void setPairs(int* p) { _pairs = p; }

    /**
     * Sets the pairs index buffer.
     */
    void setTriplets(int* p) { _triplets = p; }

    /**
     * Sets the quartets index buffer.
     */
    void setQuartets(int* p) { _quartets = p; }

    /**
     * Returns the unary costs look-up table.
     */
    double* getUnaryCosts() { return unarycosts; }

    /**
     * Returns the pairwise costs look-up table.
     */
    double* getPairwiseCosts() { return paircosts; }

    /**
     * Returns the triplet costs look-up table.
     */
    //double* getTripletCosts() { return tripletcosts; }

    /**
     * Returns the triplet costs look-up table.
     */
    // double* getQuartetCosts() { return quartetcosts; }

    /**
     * Computes the unary costs look-up table.
     */
    virtual void computeUnaryCosts(){};

    /**
     * Computes the unary potential for a the given node
     */
    virtual double computeUnaryCost(int node, int label){ return 0; };

    /**
     * Computes the pairwise costs look-up table.
     */
    virtual void computePairwiseCosts(const int *pairs){};

    /**
     * Computes the pairwise potential for a the given pair and labels.
     */
    virtual double computePairwiseCost(int pair, int labelA, int labelB){ return 0; };

    /**
     * Initializes the pair weights with given weights. If no weights are provided, all are set to one.
     */
    virtual void initPairWeights(const double *weights = 0);

    /**
     * Computes the triplet costs look-up table.
     */
    // virtual void computeTripletCosts() {};

    /**
     * Computes the triplet potential for a the given triplet and labels.
     */
    virtual double computeTripletCost(int triplet, int labelA, int labelB, int labelC) { return 0; }

    /**
     * Initializes the triplet weights with given weights. If no weights are provided, all are set to one.
     */
    ///virtual void initTripletWeights(const double *weights = 0);


    /**
     * Computes the quartet costs look-up table.
     */
    // virtual void computeQuartetCosts() {};

    /**
     * Computes the quartet potential for a the given triplet and labels.
     */
    virtual double computeQuartetCost(int quartet, int labelA, int labelB, int labelC, int labelD) { return 0; }


    /**
     * Evaluates the total cost for the zero labeling.
     */
    double evaluateTotalCostSumZeroLabeling();

    /**
     * Evaluates the total cost for the given labeling.
     */
    virtual double evaluateTotalCostSum(const int *labeling, const int *pairs, const int *triplets = 0,const int *quartets = 0);

    /**
     * Evaluates the sum of unary costs for a given labeling.
     */
    double evaluateUnaryCostSum(const int *labeling);

    /**
     * Evaluates the sum of pairwise costs for a given labeling.
     */
    double evaluatePairwiseCostSum(const int *labeling, const int *pairs);

    /**
     * Evaluates the sum of triplet costs for a given labeling.
     */
    double evaluateTripletCostSum(const int *labeling, const int *triplets);

    /**
     * Evaluates the sum of quartet costs for a given labeling.
     */
    double evaluateQuartetCostSum(const int *labeling, const int *quartets);

    /**
     * Enables the memory for the pairwise potentials computation.
     */
    // void enablePairwiseMemory();

    /**
     * Enables the memory for the triplet potentials computatio
     */
    // void enableTripletMemory();

    virtual void set_parameters(myparam & )=0;

    virtual void report(){};

  protected:
    /**
     * Initializes the cost function.
     * \param numNodes Number of nodes.
     * \param numLabels Number of labels.
     * \param numPairs Number of pairs.
     * \param numTriplets Number of triplets.
     */
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets = 0,int numQuartets = 0);

    int				m_num_nodes;			///< Number of nodes.
    int				m_num_labels;			///< Number of labels.
    int				m_num_pairs;			///< Number of pairs.
    int				m_num_triplets;			///< Number of triplets.
    int				m_num_quartets;			///< Number of quartets.


    double			*unarycosts;			///< Unary potentials look-up table.
    double			*paircosts;				///< Pairwise potentials look-up table.
    //  double			*tripletcosts;			///< Triplet potentials look-up table.
    //  double			*quartetcosts;			///< Quartet potentials look-up table.

    ///(NOT USED!!)
    double 			*pairweights;			///< Weightings for the pairwise potentials.
    //double 			*tripletweights;		///< Weightings for the triplet potentials.
    //

    int *_pairs;
    int *_triplets;
    int *_quartets;
    int *labels;							///< Labeling array.

    float _reglambda;  // scaling parameter for regulariser

    bool _debug;
    bool  _verbosity;
    bool _concat;
    string _outdir;
    string _matlabpath;
//bool write;
    //PairwiseMemory	m_paircosts_memory;		///< Pairwise potentials memory to avoid re-computations.
    //TripletMemory	m_tripletcosts_memory;	///< Triplet potentials memory to avoid re-computations.
    //bool			m_is_memory;			///< Indicates whether the pairwise potentials memory is enabled.
    //bool			m_is_memory_triplets;	///< Indicates whether the triplet potentials memory is enabled.
  };

  /// should be able to implement affine spherical registration using discrete labels represent fixed rotations
  class DummyCostFunction: public DiscreteCostFunction
  {
  public:
    DummyCostFunction(){m_num_labels=2;};

    void setUnaryCost(int node, double cost0, double cost1){

      unaryenergies.insert(pair<int, vector<double> >(node, vector<double>()));
      unaryenergies[node].push_back(cost0);
      unaryenergies[node].push_back(cost1);
    }

    void setPairwiseCost(int ind, double E00, double E01, double E10, double E11){
      pairenergies.insert(pair<int, vector<double> >(ind, vector<double>()));
      pairenergies[ind].push_back(E00);
      pairenergies[ind].push_back(E01);
      pairenergies[ind].push_back(E01);
      pairenergies[ind].push_back(E11);
    }

    double computePairwiseCost(int pair, int labelA, int labelB){
      double cost;
      if(labelA==0 && labelB==0)
	cost=pairenergies[pair][0];
      else if (labelA==0 && labelB==1)
	cost=pairenergies[pair][1];
      else if  (labelA==1 && labelB==0)
	cost=pairenergies[pair][2];
      else
	cost=pairenergies[pair][3];
      // cout << pair << " " << labelA << " " << labelB << " computePairwiseCost " << cost << endl;
      return cost;
    }

    void convertenergies(int numNodes,int numPairs, int numLabels){
      m_num_nodes=numNodes; m_num_labels=numLabels; m_num_pairs=numPairs;
      if (unarycosts){ delete [] unarycosts; }
      unarycosts = new double[numNodes*numLabels];
      for(int i=0;i<numLabels;i++){
	for(int j=0;j<numNodes;j++){
	  unarycosts[i*numNodes+j]=unaryenergies[j][i];
	}
      }
    }

    void reset(){unaryenergies.clear();  pairenergies.clear();
      // if (unarycosts){ delete [] unarycosts; }
    }
    void set_parameters(myparam & ){};
  protected:
    map<int, vector<double> > unaryenergies; /// maps of nodes  xlabels x vals
    map<int,vector<double> > pairenergies;
  };

  class SRegDiscreteCostFunction: public DiscreteCostFunction
  {

  public:
    //
    // Constructor.
    //
    SRegDiscreteCostFunction();

    ///////////////////// VARIABLE ASSIGNMENT //////////////////////////
    // data input and reference mesh, plus low resolution control point grid
    virtual void set_meshes(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source, const NEWMESH::newmesh & GRID, int num=1){_TARGET=target; _SOURCE=source; _ORIG=source;  _CPgrid=GRID; _oCPgrid=GRID;}
    void set_meshes(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source){_TARGET=target; _SOURCE=source;};

    void set_anatomical(const NEWMESH::newmesh &targetS,const NEWMESH::newmesh &targetA,const NEWMESH::newmesh & sourceS,const NEWMESH::newmesh & sourceA){_TARGEThi=targetS; _aTARGET=targetA;_aICO=sourceS, _aSOURCE=sourceA;}
    void set_anatomical_neighbourhood(const vector<map<int,double> > & weights, const vector<vector<int> > neighbourhood){_ANATbaryweights=weights; NEARESTFACES=neighbourhood; }

    ///// holds data for source and target mesh, used for look up during similarity estimation

    newmesh project_anatomical();
    void set_featurespace(const boost::shared_ptr<featurespace> &features, bool _concatenate=false){FEAT=features; _concat=_concatenate;}

    void initialize_regulariser();

    void set_initial_angles(const vector<vector<double> >  angles){

      double meanang=0;
      for (size_t i = 0; i < angles.size(); i++)
      {
    	for (size_t j = 0; j < angles[i].size(); j++)
    	{
	      meanang+=angles[i][j];
    	}
      }
      _MEANANGLE=meanang/(3.0*angles.size());
    }

    virtual void set_matlab_path(string s){_matlabpath=s;}
    ///
    void set_dataaffintyweighting(const Matrix &HRWeight){_HIGHREScfweight=HRWeight;}

    //// holds neighbourhood information
    virtual void set_relations(const boost::shared_ptr<RELATIONS> &CONTROL,const boost::shared_ptr<RELATIONS> &TARG){ _controlrel=CONTROL; _targetrel=TARG;
 _sourcerel=_controlrel->invert_relations(_CPgrid,_SOURCE);  }

    void set_parameters(myparam & );
    //MAXSEP is control grid vertex spacings
    virtual void set_spacings(const ColumnVector &spacings, const double MAX, int num=0){MAXSEP=spacings;MVDmax=MAX;}

    void debug(){_debug=true;} // for debuging
    /// label list changes between iterations for spherical optimisation (barycentres and vertices of regular sampling grid)
    void set_labels(const vector<Pt> labellist, const vector<Matrix> ROT=vector<Matrix>()){_labels=labellist;  ROTATIONS= boost::shared_ptr<vector<Matrix> >(new vector<Matrix> (ROT));  } // ROTATION MATRIX MAY ONLY BE REQUIRED FOR NON LINEAR FRAMEWORK

    inline void set_iter(int iter) {_iter=iter;}
    virtual void reset_source(const NEWMESH::newmesh & source, int num=0){_SOURCE=source; }

    void reset_anatomical(const string &, const int &);

    //void set_pairtrIDs(const vector<int> &pair_trIDs){ _trIDs=pair_trIDs;}
    virtual void reset_CPgrid(const NEWMESH::newmesh & grid, int num=0){_CPgrid=grid;}

    //////////////// GENERIC HELPER FUNCTIONS

    bool within_controlpt_range(const int &CPindex, const int &sourceindex){
      double dist=0.0;
      if(within_controlpt_range(CPindex,sourceindex,dist))return true;
      else return false;

    }

    bool within_controlpt_range(const int &CPindex, const int &sourceindex, double & dist){
      Pt CP=_CPgrid.get_coord(CPindex);
      dist= 2*RAD*asin((CP-_SOURCE.get_coord(sourceindex)).norm()/(2*RAD));
      if(dist<_controlptrange*MAXSEP(CPindex+1))return true;
      else return false;

    }
    virtual void initialize(int numNodes,int numLabels, int numPairs, int numTriplets = 0, int numQuartets = 0);

    virtual void resample_target_data(const Pt &, const Pt &, const Pt &, const Pt &, const int &, const int &, const int &, int &){};

    virtual void resample_weights(){};
    virtual void get_source_data(){};
    virtual void reset_target_data(int){};
    virtual double triplet_likelihood(const int &,const int &,const int &,const int &,const Pt &,const Pt &,const Pt &){ return 0; };

    NEWMESH::newmesh get_SOURCE(){return _SOURCE;}
    void report(){ if(_debug) cout << " sumlikelihood " << sumlikelihood << " sumregcost " << sumregcost <<endl; }
  protected:

    //////////////////////// MESHES //////////////////////
    NEWMESH::newmesh _TARGET; // TARGET MESH
    NEWMESH::newmesh _TARGEThi; // TARGET MESH
    NEWMESH::newmesh _SOURCE; // SOURCE MESH
    NEWMESH::newmesh _aTARGET; // ANATOMICAL TARGET MESH

    NEWMESH::newmesh _aSOURCE; // ANATOMICAL SOURCE MESH
    NEWMESH::newmesh _aICO; // ANATOMICAL SOURCE MESH

    NEWMESH::newmesh _aSOURCEtrans; // ANATOMICAL SOURCE MESH


    NEWMESH::newmesh _ORIG; // NON DEFORMED SOURCE
    NEWMESH::newmesh _oCPgrid; ///// NON DEFORMED CP GRID
    NEWMESH::newmesh _CPgrid; ///// CONTRO

    //////////////// CF WEIGHTINGS ////////////////////////////
    Matrix _HIGHREScfweight; /// source mesh cfweight

    ColumnVector MAXSEP;

    /////////////////////// NEIGHBOURHOOD INFO ////////////////////
    boost::shared_ptr<RELATIONS> _controlrel; // hold control grid neighbours of each source vertex
    boost::shared_ptr<RELATIONS> _targetrel; // hold target grid neighbours of each source vertex
    RELATIONS _sourcerel; // hold source grid neighbours of each CP(!!) vertex/triangle

    RELATIONS _anatrel;

    vector<vector<int> > _sourceinrange;
    vector<vector<int> > NEARESTFACES;
    //////////////////////// FEATURE SET //////////////////
    boost::shared_ptr<featurespace> FEAT; /// holds data
    sparsesimkernel<double> sim; // similarity object

    /////////////////// LABELLING PARAMETERS////////
    vector<Pt> _labels;
    boost::shared_ptr<vector<Matrix> > ROTATIONS; // rotates label set onto each control point

    double MVDmax; // max distance between CPs
    double anatMVD;
    double _MEANANGLE;
    float _controlptrange;
    double sumlikelihood;
    double sumregcost;

    float _mu; // shear modulus
    float _kappa; // bulk modulus
    float _pottsthreshold;
    float _rexp;
    float _sigma;

    ///////// USER DEFINED PARAMETERS /////////
    int _simmeasure;
    int _RES;
    int _aRES;
    int  _rmode;
    int _iter;
    int _threads;

    float _k_exp;
    vector<vector<double> > _sourcedata;
    vector<vector<double> > _targetdata;
    vector<vector<double> > _weights;

    vector<map<int,double> >  _ANATbaryweights;
    resampler R;
    double MAXstrain;
    double strain95;
  };


  /*
      class AffineSRegDiscreteCostFunction: public SRegDiscreteCostFunction
      {
      public:
      //
      // Constructor.
      //
      AffineDiscreteCostFunction();

      //
      // Computes the unary costs look-up table.

      void computeUnaryCosts() {};

      //
      // Computes the pairwise costs look-up table.

      void computePairwiseCosts() {};

      protected:

      };
    */
  class NonLinearSRegDiscreteCostFunction: public SRegDiscreteCostFunction
  {

  protected:

    ////////////// REGULARISER OPTIONS /////////////////////////
    float _maxdist;
    float _expscaling;
    bool _dweight;
    bool _anorm;

    int _kNN;
    int _currentlabelA;
    int _currentlabelB;
    int _currentlabelC;

    ColumnVector AbsoluteWeights;

    vector<Matrix> PreviousDeformations;
    vector<Matrix> ROTATE2LABEL;
    vector<Pt> _ORIGpositions;

  public:
    //
    // Constructor.
    //
    NonLinearSRegDiscreteCostFunction();

    virtual void initialize(int numNodes, int numLabels, int numPairs, int numTriplets = 0,int numQuartets = 0); // quartets not used yet so no code for them below

    void set_parameters(myparam & );

    //// compute unary costs
    void computeUnaryCosts();

    //// compute pairwise costs
    double computePairwiseCost(int pair, int labelA, int labelB);
    void computePairwiseCosts(const int *pairs) ;

    //compute triplet costs
    // void computeTripletCosts();
    double computeTripletCost(int triplet, int labelA, int labelB, int labelC);
    virtual double triplet_likelihood(const int &,const int &,const int &,const int &,const Pt &,const Pt &,const Pt &){ return 0;};
   /**
     * Computes the quartet costs look-up table.
     */
   // void computeQuartetCosts() {};

    /**
     * Computes the quartet potential for a the given triplet and labels.
     */
    double computeQuartetCost(int quartet, int labelA, int labelB, int labelC, int labelD) { return 0; }

    // for transforming anatomy
    Triangle deform_anatomy(const int&,const int&,map<int,Pt>  & ,map<int,bool> &, map<int,Pt>  &);


    void resample_weights();

    void get_target_data(const int &,const Matrix &);

    void get_target_data(const int&,const Pt&,const Pt&,const Pt&);



    virtual void get_data_index(const int & controlpoint, const int & sourcepoint, int & ind){}; // data is held in 2D vectors in each case, but for univariate case data is controlpointsxfeatures and for multivariate its sourcepointsxfeatures. This function determines the correct index for rows
  };

  class UnivariateNonLinearSRegDiscreteCostFunction: public NonLinearSRegDiscreteCostFunction
  {
  public:
    //
    // Constructor.
    //
    UnivariateNonLinearSRegDiscreteCostFunction(){};
    virtual void initialize(int numNodes, int numLabels, int numPairs,int numTriplets = 0, int numQuartets = 0); // quartets not used yet
    virtual void get_source_data();

    double computeUnaryCost(int node, int label);

    // void resample_source_data(const int &);
    void reset_target_data(int node){ _targetdata[node-1].clear();}
    virtual void resample_target_data(const Pt &, const Pt &, const Pt &, const Pt &, const int &, const int &, const int &, int &);
    void get_data_index(const int & controlpoint, const int & sourcepoint, int &ind){ind=controlpoint;}

  };

  class MultivariateNonLinearSRegDiscreteCostFunction: public NonLinearSRegDiscreteCostFunction
  {
  public:
    //
    // Constructor.
    //
    MultivariateNonLinearSRegDiscreteCostFunction(){};
    virtual void initialize(int numNodes,int numLabels, int numPairs, int numTriplets = 0,int numQuartets = 0); // quartets not used yet
    void get_source_data();

    double computeUnaryCost(int node, int label);

    //   void resample_source_data(const int & sourcept);
    void reset_target_data(int node){  for(unsigned int i=0;i<_sourceinrange[node-1].size();i++) { _targetdata[_sourceinrange[node-1][i]-1].clear();}} /// *****NOT THREAD SAFE ********

    virtual void resample_target_data(const Pt &, const Pt &, const Pt &, const Pt &, const int &, const int &, const int &, int &);
    void get_data_index(const int & controlpoint, const int & sourcepoint, int & ind){if(sourcepoint > (int) _sourcedata.size()){cout << "vector index too high " << endl;exit(1);}
      ind=sourcepoint;}

  };

  class alphaMINonLinearSRegDiscreteCostFunction: public MultivariateNonLinearSRegDiscreteCostFunction
  {
  public:
    //
    // Constructor.
    //
    alphaMINonLinearSRegDiscreteCostFunction(){};
    //   void resample_source_data(const int & sourcept);

  };

 class HOUnivariateNonLinearSRegDiscreteCostFunction: public UnivariateNonLinearSRegDiscreteCostFunction
  {
  public:
    //
    // Constructor.
    //
    HOUnivariateNonLinearSRegDiscreteCostFunction(){};
    void initialize(int numNodes,int numLabels, int numPairs,int numTriplets = 0, int numQuartets = 0); // quartets not used yet
    void get_source_data();

    double computeUnaryCost(int node, int label){return 0;};

    void set_relations(const boost::shared_ptr<RELATIONS> &CONTROL,const boost::shared_ptr<RELATIONS> &TARG){_controlrel=CONTROL; _targetrel=TARG;
      _sourcerel=_controlrel->invert_relationsTR(_CPgrid,_SOURCE);  }
    double triplet_likelihood(const int &,const int &,const int &,const int &,const Pt &,const Pt &,const Pt &);

  };

 class HOMultivariateNonLinearSRegDiscreteCostFunction: public MultivariateNonLinearSRegDiscreteCostFunction
  {
  public:
    //
    // Constructor.
    //
    HOMultivariateNonLinearSRegDiscreteCostFunction(){};
    void initialize(int numNodes,int numLabels, int numPairs,int numTriplets = 0, int numQuartets = 0); // quartets not used yet
    void get_source_data();

    double computeUnaryCost(int node, int label){return 0;};

    void set_relations(const boost::shared_ptr<RELATIONS> &CONTROL,const boost::shared_ptr<RELATIONS> &TARG){_controlrel=CONTROL; _targetrel=TARG;
      _sourcerel=_controlrel->invert_relationsTR(_CPgrid,_SOURCE); }

    double triplet_likelihood(const int &,const int &,const int &,const int &,const Pt &,const Pt &,const Pt &);
  };



 class RegularisationDiscreteCostFunction: public SRegDiscreteCostFunction
 {
 	 protected:

     //////////////////////// MESHES //////////////////////
    // NEWMESH::newmesh _SPHERE; // SOURCE MESH
     float _rexp;

   public:
     //
     // Constructor.
     //
     RegularisationDiscreteCostFunction(double shear, double bulk, double exp);


	 void initialize(int numNodes,int numLabels, int numPairs,int numTriplets = 0, int numQuartets = 0); // quartets not used yet
	 //void reset_CPgrid(const NEWMESH::newmesh & grid, int num=0){_SPHERE=grid;}
	 double computeTripletRegularisation(int triplet, int labelA, int labelB, int labelC);

 };

 class MetricDistortionDiscreteCostFunction: public RegularisationDiscreteCostFunction
  {

      protected:
      bool isinitialised;

  	  public:
		//
		// Constructor.
		//
		MetricDistortionDiscreteCostFunction(double shear, double bulk, double exp):RegularisationDiscreteCostFunction(shear, bulk,  exp){isinitialised=false;}


		double computeTripletCost(int triplet, int labelA, int labelB, int labelC){
			if (isinitialised==false)
				throw  DISCRETEOPTHOCRException("MetricDistortionDiscreteCostFunction::You must supply meshes.");
			return computeTripletRegularisation(triplet,labelA, labelB,labelC);

		}
};





}
