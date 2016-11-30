/*  ContinuosOpt.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  

    Some sections of code inspired by Alek Petrovic.*/
/*  CCOPYRIGHT  */
/// Class for gradient based optimisation - hopefully soon to be obselete just need to convert affine method to discrete ////

#ifndef ContinuosOpt_h
#define ContinuosOpt_h

#include "newmesh/featurespace.h"

//#include <FastPD/FastPD.h>

#ifdef HAS_HOCR
#include "Fusion.h"
#else
#include <FastPD/FastPD.h>
#endif

using namespace DISCRETEOPT;


#define CORLIM  1E-10  // when calculating corrzero consider any correlations below this value to be ZERO- NOT USED ANY MORE?

using namespace NEWMESH;

namespace MESHREG{
  
  class MESHREGException: public std::exception
    {
    private:
      std::string m_msg;
    public:
      MESHREGException(const std::string& msg) throw(): m_msg(msg) {}
      
      virtual const char * what() const throw() {
	return string("Fnirt: msg=" + m_msg).c_str();
      }
      
      ~MESHREGException() throw() {}
    };


  class MeshCF {    /// basic mesh cost function - uses the gradient of the similarity function, estimated using weighted least squares
   

  protected:
    NEWMESH::newmesh _TARGET; // TARGET MESH
    NEWMESH::newmesh _SOURCE; // SOURCE MESH
   
    boost::shared_ptr<RELATIONS> _rel; // saves mesh neighbourhood information

    Matrix _inweight; // exclusion/weighting mask for cost function masking
    Matrix _refweight;

    const boost::shared_ptr<featurespace> FEAT; /// holds data
    sparsesimkernel<double> sim; // similarity matrix
  
  
    double MVD; // mean vertex distance

    ///////// user defined parameters /////////
    int _simmeasure; 
    int _iters; // total iterations
    float _stepsize;
    float _spacing;
       //////////////////////////////////////////

    
    string m_outdir;
  
    vector<Tangs> BASIS; /////////////  tangent vector bases for all vertices

    bool  _verbosity;
    ////// GRAD DESCENT PARAMETERS /////////////////
    mutable ColumnVector current_sim;

  
    double min_sigma;
    double CF; 
    double totJP;
    double totJD;

   
  public:
   
     MeshCF(const NEWMESH::newmesh &        target, 
	   const NEWMESH::newmesh &        source,
	   const Matrix T_cfweight, 
	   const Matrix S_cfweight,
	   const boost::shared_ptr<featurespace> features)      
      : _TARGET(target),_SOURCE(source),_inweight(S_cfweight),_refweight(T_cfweight), FEAT(features)
    {
    
      ////  DEFAULT PARAMETRISATION /////
      _simmeasure=3;
      _iters=10000;
      _verbosity=false;
    }

    MeshCF (const NEWMESH::newmesh target, const NEWMESH::newmesh &source,
	    const boost::shared_ptr<featurespace> features):_TARGET(target), _SOURCE(source),FEAT(features){  
      _simmeasure=3;
      _iters=20;
      _verbosity=false;
      _inweight.ReSize(1, _SOURCE.nvertices()); _refweight.ReSize(1, _TARGET.nvertices());  // if no cf weighting then set weighting to 1 by default
      _inweight=1; _refweight=1;
      
    } 
    
    void set_parameters(myparam PAR){
      myparam::iterator it;
      it=PAR.find("iters");_iters=boost::get<int>(it->second);
      it=PAR.find("simmeasure");_simmeasure=boost::get<int>(it->second); 
      it=PAR.find("verbosity");_verbosity=boost::get<bool>(it->second);
      it=PAR.find("stepsize");_stepsize=boost::get<float>(it->second);
      it=PAR.find("gradsampling");_spacing=boost::get<float>(it->second);
      it=PAR.find("outdir");m_outdir=boost::get<string>(it->second);

    }

    ///////////////////// INITIALIZE ///////////////////
    virtual void Initialize();  
    void set_simmeasure(const int & simval){_simmeasure=simval;}
     
    //////////////// MAKE UPDATES ////////////////////
    
    void update_similarity();  
    void update_similarity(const int &, vector<int> &); 
    void update_source(const NEWMESH::newmesh& M) {_SOURCE=M; };
    
    //////////////// SIMILARITY GRADIENT ESTIMATION /////////////////////  
    ///// similarity gradient via weighted regression -
    ColumnVector WLS_simgradient(const Tangs &T, int ,const vector<int> &);  
    //// prepares data for sim gradient calculation
   
    ColumnVector Evaluate_SIMGradient(int i,const Tangs &T);
    
    virtual  NEWMESH::newmesh run()=0;
    ////////////////// RETURNING FUNCTIONS ///////////////////////

    NEWMESH::newmesh get_REG()const{return _SOURCE;}  
  };


  class affineMeshCF: public MeshCF 
  {


  public:
    
  affineMeshCF(const NEWMESH::newmesh & target, const NEWMESH::newmesh &source,const boost::shared_ptr<featurespace> features ): MeshCF(target,source,features){}
    virtual ~ affineMeshCF() {};
    
    void Initialize(){MeshCF::Initialize();}
    /////////////// AFFINE FUNCTIONS - uses image gradients and euler rotations - taken from Alek Petrovic's implementation.
    void Rotate_IN_mesh(const double &, const double &,const double &);
    double Affine_cost_mesh(const double & , const double & ,const double &);
    NEWMESH::newmesh run();
    ColumnVector return_transformation(){return REC;} 
   private:

    ColumnVector REC;
    affineMeshCF();
    virtual  affineMeshCF& operator=(const  affineMeshCF& inf) {throw MESHREGException("LM_MeshCF:: Assignment explicitly disallowed"); return(*this);}

  };
  

  
}

#endif
