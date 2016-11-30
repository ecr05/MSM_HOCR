/*  groupmeshreg.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/* Class for overseeing the running of each resolution level of the registration */
#ifndef groupmeshreg_h
#define groupmeshreg_h


#include "meshreg.h"
#include "miscmaths/SpMat.h"


#define CORLIM  1E-10  // when calculating corrzero consider any correlations below this value to be ZERO
#define D_DIST 0.25   // sampling distance along tangent plane for computing derivatives - 0.25 is what Bruce uses for FreeSurfer

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

namespace MESHREG {
  
  
  class GroupMeshReg: public MeshReg {
    
  private:
 
    vector<newmesh> ALL_SPH_REG;

  public:
    
    // Constructors
    inline GroupMeshReg(){ };
    
    //MeshReg(int);
    
    // Destructor
    inline ~GroupMeshReg(){};
    
    
    void Initialize_level(int);
    void Evaluate();
    void Transform(const string &);
    void saveTransformedData(const double &,const string &);
    void run_discrete_opt(vector<NEWMESH::newmesh> &);

    inline void saveSPH_reg(const string &filename) const {
      char fullpath[1000];

      for(int i=0;i<ALL_SPH_REG.size();i++){
	sprintf(fullpath,"%ssphere-%d.LR.reg%s",filename.c_str(),i,_surfformat.c_str());
	cout <<i <<  " new save SPH " << fullpath << endl;
	    
	ALL_SPH_REG[i].save(fullpath);
	
      }
    
    }
  };

  
  
  
}

  

#endif

