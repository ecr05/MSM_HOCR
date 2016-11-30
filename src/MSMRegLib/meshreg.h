/*  meshreg.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/* Class for overseeing the running of each resolution level of the registration */
#ifndef meshreg_h
#define meshreg_h


#include "meshmodify.h"
#include "miscmaths/SpMat.h"


#define CORLIM  1E-10  // when calculating corrzero consider any correlations below this value to be ZERO
#define D_DIST 0.25   // sampling distance along tangent plane for computing derivatives - 0.25 is what Bruce uses for FreeSurfer

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

namespace MESHREG {

  struct spcoord
  {
    double theta;
    double phi;
    
    int ind;
  };
  
  struct SortPhi
  {
    bool operator() (const spcoord& m1, const spcoord& m2)
    { return m1.phi < m2.phi; }
  };
  
  struct SortAscTheta
  {
    bool operator() (const spcoord& m1, const spcoord& m2)
    { return m1.theta < m2.theta; }
  };
  
  struct SortDecTheta
  {
    bool operator() (const spcoord& m1, const spcoord& m2)
    {return m1.theta > m2.theta; }
  };
  
  
  class MeshReg: public MeshModify {
    
  private:
 
    newmesh icoref;
    newmesh ANAT_RES;
    
    bool _isaffine;

    boost::shared_ptr<affineMeshCF> affinecf;

  protected:
    int _level;
    boost::shared_ptr<SRegDiscreteModel> MODEL;
  public:
    
    // Constructors
    inline MeshReg(){_isaffine=false; _level=0;};
    
    //MeshReg(int);
    
    // Destructor
    inline ~MeshReg(){};
    
    //////////////// INITIALIZATION ////////////////////
    //// loops over steps/resolution levels
    void run_multiresolutions(const int &, const double &, const string &);
    ///// initialises featurespace and cost function model for each step (either a single resolution level of the discrete opt or an Affine step)
    void Initialize_level(int);
    /// downsample anatomy
    newmesh  resample_anatomy(newmesh, vector<map<int,double> > &, vector<vector<int > > &,int );
    ///// additionally resamples the costfunction weighting to match the downsampled data grids for that resolution leve;
    Matrix downsample_cfweighting(const double &, NEWMESH::newmesh, boost::shared_ptr<NEWMESH::newmesh>, boost::shared_ptr<NEWMESH::newmesh>);
    
    //////////////// RUN //////////////////////
    void Evaluate();
    void Transform(const string &);
    void saveTransformedData(const double &,const string &);

    ///// projects warp from a lower resolution level to the next one
    NEWMESH::newmesh project_CPgrid(NEWMESH::newmesh,NEWMESH::newmesh, int num=0);
    /////// runs all iterations of the discrete optimisation for each resolution level
    void run_discrete_opt(NEWMESH::newmesh &);
   
   
  };

  
  
  
}

  

#endif

