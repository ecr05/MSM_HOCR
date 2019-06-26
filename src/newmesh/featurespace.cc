/*  featurespace.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "featurespace.h"

namespace NEWMESH {


  ////////////////////////// RESAMPLING OF DATA DIMENSION//////////////////////////////

  NEWMESH::newmesh featurespace::Initialize(const int &ico, vector<NEWMESH::newmesh> &IN,const bool & exclude){
   
    NEWMESH::newmesh icotmp;
    
    if(IN.size()!=CMfile_in.size()){ throw  NEWMESHException(" NEWMESH::featurespace::Initialize do not have the same number of datasets and surface meshes");	}
    else {DATA.resize(IN.size(),boost::shared_ptr<BFMatrix > ()); }
    
    if(ico>0){icotmp.make_mesh_from_icosa(ico); true_rescale(icotmp,RAD); }
     

    
    for (unsigned int i=0;i<IN.size();i++){
      boost::shared_ptr<BFMatrix> tmp=boost::shared_ptr<BFMatrix > ();
      bool isfunc=set_data(CMfile_in[i],DATA[i],IN[i],_issparse); 
      IN[i].set_pvalues(DATA[i]->AsMatrix());
    
      if(ico==0){icotmp=IN[i];} 
      // create exclusion mask 

 
      if(exclude || _cut){
		NEWMESH::newmesh excl_tmp=create_exclusion(IN[i],DATA[i]->AsMatrix(),_fthreshold[0],_fthreshold[1]);
		EXCL.push_back(boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(excl_tmp))); ///mask data according to some min (_fthreshold[0]) and max (_fthreshold[1]) threshold
	     
      }else{EXCL.push_back(boost::shared_ptr<NEWMESH::newmesh>()); }

      ///// downsample data to regular grid of resolution defined by "ico"
     
      resample(_sigma_in[i],DATA[i],icotmp,IN[i],EXCL[i]);
     
    }
    icotmp.set_pvalues(DATA[0]->AsMatrix());
  
    /// intensity normalise using histogram matching 
    if(_intensitynorm){ 
	  
      
      for (unsigned int i=1;i<IN.size();i++){
		multivariate_histogram_normalization(*DATA[i],*DATA[0],EXCL[i],EXCL[0],_scale);  // match input data feature distributions to equivalent in ref, rescale all to first feature in reference if _scale is 
   
      }
      
      if(!exclude){
		  for (unsigned int i=0;i<IN.size();i++){
		   EXCL[i]=boost::shared_ptr<NEWMESH::newmesh>();
		  }
      }
    }
    
    
    if(_logtransform){
      for (unsigned int i=0;i<IN.size();i++)
	log_transform_and_normalise(*DATA[i]);
    }

    if(_varnorm){
      for (unsigned int i=0;i<IN.size();i++){
	varnorm(DATA[i],EXCL[i]); 
       
      
      }
    }
   
    return icotmp;
 
  }

  ////////////////////////// RESAMPLING OF DATA DIMENSION//////////////////////////////


  void featurespace::resample(const double &sigma, boost::shared_ptr<BFMatrix> &DATAtmp, NEWMESH::newmesh &icotemp,const NEWMESH::newmesh &M, boost::shared_ptr<NEWMESH::newmesh> &EXCLtmp){
    
    resampler R;
    newmesh tmp=M;

    if(_resamplingmethod=="ADAP_BARY"){
      R.set_method("ADAP_BARY");  
      
      R.resampledata(tmp,icotemp,EXCLtmp,DATAtmp,sigma);

     
      if(sigma>0){
	R.set_method("GAUSSIAN");
	R.smooth_data(sigma,DATAtmp,icotemp);
	
      }   
    }
    else{
      R.set_method("GAUSSIAN");     
      R.resampledata(tmp,icotemp,EXCLtmp,DATAtmp,sigma);
    }
  }


  void featurespace::smooth(NEWMESH::newmesh &IN, boost::shared_ptr<BFMatrix> &tmp, const double &_sigma)
  {
    
    resampler R;
    R.set_method("GAUSSIAN");
    R.smooth_data(_sigma,tmp,IN);
   
  }

  void featurespace::varnorm(boost::shared_ptr<BFMatrix> & DATA, boost::shared_ptr<NEWMESH::newmesh> &EXCL){ /// need to condense code of next two functions into one
   

    vector<vector<double> > _data;
    vector<double> var1;
    vector<double> tmp;
    int ind;
    ind=0;
    
    
    for (unsigned int i=1;i<=DATA->Ncols();i++){
      if(!EXCL.get() || EXCL->get_pvalue(i-1)>0){
	_data.push_back(tmp);
	for(unsigned int k=1;k<=DATA->Nrows();k++){
	  _data[ind].push_back(DATA->Peek(k,i));
	}
	ind++;
      }
    } 
    ind=0;
  
    var1=online_variance_normalize(_data); /// can replace with normalise(_sourcedata)
   
    ind=0;

    for (unsigned int i=1;i<=DATA->Ncols();i++){
      if(!EXCL.get() ||EXCL->get_pvalue(i-1)>0){
	for(unsigned int k=1;k<=DATA->Nrows();k++){
	  DATA->Set(k,i,_data[ind][k-1])  ;  
	 
	}
	ind++;
      }
    }
       
  }


  /// calculate mean and variance simultaneously to speed up estimation
  vector<double> featurespace::online_variance_normalize( vector<vector<double> > &M){   // variance normalisation is susceptible to precision error especially when data points are close to the mean
  
    RowVector mean(M[0].size());
    vector<double> var(M[0].size(),0);
    double delta,maxvar;
    mean=0; maxvar=0;

    for (unsigned int j=0;j<M[0].size();j++){
      var[j]=0; mean(j+1)=0;
      for (unsigned int i=0;i<M.size();i++){
	delta=M[i][j]-mean(j+1);
	mean(j+1)=mean(j+1)+delta/(i+1);
	var[j]= var[j] +delta*(M[i][j]-mean(j+1));
      }
      var[j]= var[j]/(M.size()-1);

      for (unsigned int i=0;i<M.size();i++){
	M[i][j]-=mean(j+1);
	if(var[j]>0) M[i][j]=M[i][j]/sqrt(var[j]);
	if(var[j]>maxvar) maxvar=var[j];
      }

     }

    for (unsigned int i=0;i<M[0].size();i++){
	var[i]/=maxvar;
      }

    return var;
   
  }


 
    
}

  
  
