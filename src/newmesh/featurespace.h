/*  featurespace.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/*  CLASS FOR PROCESSING DATA PRIOR TO MSM REGISTRATION */
/* ideally this class should probably be split up and the major components should go into newmesh for manipulation the mesh data for each mesh individually */

#if !defined(featurespace_h)
#define featurespace_h

#include <fstream>
#include <stdio.h>

#include "meshfns.h"
#include "utils/options.h"

using namespace Utilities;

namespace NEWMESH{

   enum StringValue { evNotDefined, 
		      evStringValue1, 
		      evStringValue2, 
		      evStringValue3, 
		      evEnd };
  // Map to associate the strings with the enum values
   static std::map<std::string, StringValue> s_mapStringValues;
  
   /* features - should this also inherit from newmat?  */
  class featurespace{
    
  private:
   
    NEWMESH::newmesh source;

    vector<boost::shared_ptr<BFMatrix> > DATA; // holds generic BFMATRIX data which can be sparse or full matrices
    vector<boost::shared_ptr<NEWMESH::newmesh> >  EXCL;  // exclusion masks for binary weighting of the data during resampling
   
    vector<string >  CMfile_in;  // path to data
    vector<double> _sigma_in;  // smoothing parameters for input and reference
    vector<float> _fthreshold;

    string inorig;   
    string reforig;
    string _resamplingmethod;
 
    bool _logtransform;  // will log transform and normalise
    bool _intensitynorm; // will histogram match
    bool _scale;  /// will rescale each feature to crudley match the distribution of the first in a multivariate distribution
    bool _issparse;  /// notes that data is sparse
    bool _varnorm;  // performs online variance normalisation - maybe replace with non online version called during logtransformandnormalise()
    bool _cut;
  public:

    featurespace(){};
    ~featurespace(){};
    featurespace(const string &datain, const string &dataref){
      _sigma_in.resize(2,5.0);_logtransform=false; _issparse=false;  _scale=false; _intensitynorm=false; _fthreshold.resize(2,0.0); 
      CMfile_in.push_back(datain); CMfile_in.push_back(dataref);
    };

    featurespace(const string &datain, const vector<string> &datareflist){
      _sigma_in.push_back(5);_logtransform=false; _issparse=false;  _scale=false; _intensitynorm=false; _fthreshold.resize(2,0.0); 
       CMfile_in.push_back(datain); 
       for(int i=0;i<(int) datareflist.size();i++){
	 CMfile_in.push_back(datareflist[i]);
	 _sigma_in.push_back(5);
       }

    };

    featurespace(const vector<string> &datalist){
      _logtransform=false; _issparse=false;  _scale=false; _intensitynorm=false; _fthreshold.resize(2,0.0); 
      for(int i=0;i<(int) datalist.size();i++){
	 CMfile_in.push_back(datalist[i]);
	 _sigma_in.push_back(5);
       }

    };
    /////////////////// INITIALIZE //////////////////////////////
    void set_smoothing_parameters(const vector<double> s){ 
      _sigma_in.clear();
      if(s.size()!=CMfile_in.size()){ 
	if(s.size()==1){ 
	  for (int i=0;i<(int) CMfile_in.size();i++)
	    _sigma_in.push_back(s[0]);	  
	}else throw  NEWMESHException("Mewmesh::featurespace:: smoothing sigma size incompatible with data dimensions");
      }else _sigma_in=s;};

    void set_cutthreshold(vector<float> & thr){_fthreshold=thr;};
    void logtransform(const bool &log){_logtransform=log;}
    void varnorm(const bool &norm){_varnorm=norm;}

    void is_sparse(const bool &sp){_issparse=sp;}
    void intensitynormalize(const bool & norm, const bool &scale, const bool& _exclcut){_intensitynorm=norm; _cut=_exclcut;};
    void resamplingmethod(string method){_resamplingmethod=method;}

    NEWMESH::newmesh Initialize(const int &, vector<NEWMESH::newmesh> &,const bool &);
    /////////////////// RESAMPLE DATA /////////////////////////
    void resample(const double &, boost::shared_ptr<BFMatrix> &,NEWMESH::newmesh &,const NEWMESH::newmesh &, boost::shared_ptr<NEWMESH::newmesh> &);
    void smooth(NEWMESH::newmesh &, boost::shared_ptr<BFMatrix> &, const double &);
    void varnorm(boost::shared_ptr<BFMatrix> &, boost::shared_ptr<NEWMESH::newmesh> &); // combine this and next function
    vector<double> online_variance_normalize(vector<vector<double> >  &); // mean centre and v ariance normalize

    //////////////////// ACCESS //////////////////////////////////////////////
    string get_path(const int i)const{return CMfile_in[i];}
    int get_dim()const{return DATA[0]->Nrows();}

    /////// PAIRWISE ACCESSS /////////////////////////////////////////////////////////////////////
    string get_input_path()const{return CMfile_in[0];}
    string get_reference_path()const{return CMfile_in[1];}
    double get_input_val(const int &i,const int &j)const{ return DATA[0]->Peek(i,j);}
    double get_ref_val(const int &i,const int &j)const{ return DATA[1]->Peek(i,j);}
    double get_data_val(const int &i,const int &j,const int n)const{ return DATA[n]->Peek(i,j);}
    void set_input_val(const int &i,const int &j, const double &val)const{ return DATA[0]->Set(i,j,val);}
    void set_ref_val(const int &i,const int &j, const double &val)const{ return DATA[1]->Set(i,j,val);}
    boost::shared_ptr<NEWMESH::newmesh>  get_input_excl()const{return EXCL[0];}
    boost::shared_ptr<NEWMESH::newmesh>  get_reference_excl()const{return EXCL[1];}

    boost::shared_ptr<BFMatrix>  get_input_data()const{
      return DATA[0];
    }

    boost::shared_ptr<BFMatrix>  get_data(const int i)const{
      return DATA[i];
    }
    Matrix  get_data_matrix(const int i)const{
       return DATA[i]->AsMatrix();
    }
    boost::shared_ptr<BFMatrix> get_reference_data()const{
      return DATA[1];
    }
  };
  
}
#endif
  
  
