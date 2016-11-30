/* this program is designed to downsample freesurfer label files to be used in combination with the SPH6.vtk or other downsampled meshes*/

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include "utils/options.h"
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"


using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "msm_metric_average  <asciii list of data> <output> <-options>   " << endl;
  cout << "options:  " << endl;
  cout << "-target " << endl;
  cout << "-abs " << endl;
  cout << "-weighted distances" << endl;
  cout << "-sigma " << endl;
  cout << " -normalize  - normlalize intensity range to target (requires argument)" << endl;

}


int main(int argc, char **argv){


  boost::shared_ptr<NEWMESH::newmesh> EXCL_ref;
  SpMat<double> *DATA =new SpMat<double>();
  int ok;
  int _avdim=0;
  double newval;
  float _sigma;
  
  NEWMESH::newmesh totalcounts;
 
  NEWMESH::newmesh IN, TARG, HISTTARG,AVERAGE,VARIANCE,Zscore;
  
  string output,excl_in,excl_ref;
  vector<string> distances;

  vector<string> INlist,INsurflist;
  vector<newmesh> INPUTS;
  bool _isdatamatrix=false;
  string subtype;
  bool _exclude=false,_excludenorm=false,_normalize=false;
  bool _abs=false,_weighted=false;
  ColumnVector exclusion;

  if(argc < 2){

    Usage();
    exit(0);
  }

  INlist=read_ascii_list(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  newval=0.0;
  excl_in="";excl_ref="";
  EXCL_ref=boost::shared_ptr<NEWMESH::newmesh >(new NEWMESH::newmesh ()) ;
 

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      _isdatamatrix=1;
      cout << " normalize " << argv[1] << endl;
      TARG.load(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-inputsurfaces") == 0)){
      argc--;
      argv++;
      INsurflist=read_ascii_list(argv[1]);
      if(INsurflist.size()!=INlist.size()){ cout << " surf list is not the same length as data list" << endl; exit(1); }
      argc--;
      argv++;
      ok = 1;
    }
	else if((ok == 0) && (strcmp(argv[1], "-excl_in") == 0)){
      argc--;
      argv++;
      _exclude=1;
      excl_in=argv[1];
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-excl_ref") == 0)){
      argc--;
      argv++;
      _exclude=1;
      excl_ref=argv[1];
      EXCL_ref->load(excl_ref,false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-abs") == 0)){
      argc--;
      argv++;
      _abs=1;
      ok = 1;
    }
     else if((ok == 0) && (strcmp(argv[1], "-averagedim") == 0)){
      argc--;
      argv++;
      _avdim=atoi(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
   else if((ok == 0) && (strcmp(argv[1], "-weighted") == 0)){
      argc--;
      argv++;
      _weighted=true;
      distances=read_ascii_list(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;
      argv++;
      _sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-normalize") == 0)){
      argc--;
      argv++;
      _normalize=true;
      HISTTARG.load_gifti(argv[1],false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-excl") == 0)){
      argc--;
      argv++;
      _excludenorm=true;
      ok = 1;
    }
    else{cout << " option doesn't exist " << endl; exit(1);}
  }
  if(!_isdatamatrix)TARG.load(INlist[0]);
  if(excl_ref==""){
    *EXCL_ref=TARG;

    for (int i=0;i<EXCL_ref->nvertices();i++){
      EXCL_ref->set_pvalue(i,1);
    }
  }

  exclusion.ReSize(TARG.nvertices()); exclusion=1;

  /////////////// INITIALIZE //////////////////////
  int dim,vert;
  subtype = INlist[0].substr(INlist[0].size()-8, 4);
  if(_isdatamatrix){  /// set average either to dimensions of target mesh or, assuming all data have been resampled already to the size of first input mesh
    AVERAGE=TARG;
    vert=TARG.nvertices();
    if(INsurflist.size()==0) IN=AVERAGE;
  }
  else{
    AVERAGE.load(INlist[0]);
    vert=AVERAGE.nvertices();
  }
  
  if(subtype=="func" || subtype=="hape" ){ // define dimension of data
	IN.load(INlist[0],false,false);
	dim=IN.get_dimension();
  }else {
     DATA=new SpMat<double>(INlist[0]);
     if(DATA->Ncols()>DATA->Nrows()) dim=DATA->Nrows();
     else dim=DATA->Ncols();
  }
  

  /// initialise average and variances to zero 
  Matrix tmpData(dim,vert);  tmpData=0;
  AVERAGE.set_pvalues(tmpData);
  VARIANCE=AVERAGE;
  Zscore=AVERAGE;   
  totalcounts=AVERAGE;
  
  int atlasdim=AVERAGE.get_dimension();
  double sumweighted=0;
  vector<float> weight(INlist.size(),1);
  /////////////////////////////////////////////////////
  for (unsigned int i=0;i<INlist.size();i++){
	  cout << i << " " << INlist[i] << endl;
    if(subtype=="func"|| subtype=="hape"){
	  if(INsurflist.size()==INlist.size()) IN.load(INsurflist[i]);
      IN.load(INlist[i],false,false);
    }else if(_isdatamatrix){
      DATA=new SpMat<double>(INlist[i]);
      IN.set_pvalues(DATA->AsNEWMAT());
    }else
      IN.load(INlist[i]);

    if(IN.nvertices()!=TARG.nvertices() && TARG.nvertices()>0){
      //// resample data onto target
      resampler R; R.set_method("ADAP_BARY");
      Matrix tmpDATA=IN.get_pvalues();
      R.resampledata(IN,TARG,tmpDATA,1);
      IN=TARG;
      IN.set_pvalues(tmpDATA);
    }

   //// options to intensity normalise 
   if(_normalize)  {
    Matrix DATAREF=HISTTARG.get_pvalues();
    Matrix DATAIN=IN.get_pvalues();
    boost::shared_ptr<BFMatrix > BFIN,BFREF;
    boost::shared_ptr<NEWMESH::newmesh> EXCL_IN, EXCL_REF;
	BFIN = boost::shared_ptr<BFMatrix >(new FullBFMatrix (DATAIN));
    BFREF=boost::shared_ptr<BFMatrix > (new FullBFMatrix (DATAREF)); //

	if(_excludenorm){
		cout << " TARG.nvertices() " << TARG.nvertices() << " DATAREF.Ncols() " << DATAREF.Ncols() << endl;
		EXCL_IN= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(HISTTARG,DATAIN,0,0))); 
		EXCL_REF= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(HISTTARG,DATAREF,0,0))); 

	}
   
    multivariate_histogram_normalization(*BFIN,*BFREF,EXCL_IN,EXCL_REF,false);
    //IN.set_pvalues(BFIN->AsMatrix());
    //char filename[1000];
    //sprintf(filename,"IN_norm-%d.func",i);
    //IN.save(filename);

    
    }	
    
	INPUTS.push_back(IN);
    boost::shared_ptr<NEWMESH::newmesh> EXCL;
    int dim=IN.get_dimension();
    if(dim!=atlasdim){cout << i << " dim " << dim << " atlasdim " << endl; throw MeshReg_error("Feature dimensions do not agree");}
    
    if(_weighted){
      if(_sigma==0){ cout << "must set sigma for weighted" <<endl; exit(1);}
      weight[i]=exp(-0.5*((atof(distances[i].c_str()))/_sigma)*((atof(distances[i].c_str()))/_sigma));
     // sumweighted+=weight[i];
      cout << "distance " <<  distances[i] << " weight " << weight[i] << endl;
    }//else{ sumweighted=1;}
      cout << i << " weight " << weight[i] << endl;
    // cout << " add to total ... " << IN.nvertices() << endl;
    for(int d=0;d<atlasdim;d++){
      cout << d << " atlasdim " <<atlasdim << " dim " << dim << "EXCL_ref->nvertices() " << EXCL_ref->nvertices() << endl; 
      for (int j=0;j<AVERAGE.nvertices();j++){
	if(EXCL_ref->get_pvalue(j) && abs(IN.get_pvalue(j,d)) > EPSILON){
	  if(_abs) newval=AVERAGE.get_pvalue(j,d)+weight[i]*abs(IN.get_pvalue(j,d));
	  else newval=AVERAGE.get_pvalue(j,d)+weight[i]*IN.get_pvalue(j,d);
	  totalcounts.set_pvalue(j,totalcounts.get_pvalue(j,d)+weight[i]*1,d);
	  AVERAGE.set_pvalue(j,newval,d);
	}else{
	  if(_exclude)
	    exclusion(j+1)=0;
	}

      }
     
    }
    // cout <<" here "<< endl;
  }

  double sumav=0;
 
  for(int d=0;d<atlasdim;d++){  
    for (int i=0;i<AVERAGE.nvertices();i++){
      if(totalcounts.get_pvalue(i,d)==0) newval=0;
      else newval=AVERAGE.get_pvalue(i,d)/(totalcounts.get_pvalue(i,d));
      AVERAGE.set_pvalue(i,newval,d);
      if(_avdim==-1 || _avdim==d) sumav+=newval;
      //   cout << i << " newval " << newval << " AVERAGE.get_pvalue(d,i) " <<  AVERAGE.get_pvalue(d,i) << " totalcounts.get_pvalue(d,i) " << totalcounts.get_pvalue(d,i) << " sumav " << sumav << endl;
    }
  }
  

  for (unsigned int i=0;i<INlist.size();i++){
   
    for(int d=0;d<atlasdim;d++){
      for (int j=0;j<AVERAGE.nvertices();j++){
	if(EXCL_ref->get_pvalue(j) && exclusion(j+1)){
	  if(_abs) newval=weight[i]*(abs(INPUTS[i].get_pvalue(j,d))-AVERAGE.get_pvalue(j,d));
	  else  newval=weight[i]*(INPUTS[i].get_pvalue(j,d)-AVERAGE.get_pvalue(j,d));
	  newval=VARIANCE.get_pvalue(j,d)+newval*newval;
	  VARIANCE.set_pvalue(j,newval,d);
	}
      }

    }
  }
  double meanvar=0;
  for(int d=0;d<atlasdim;d++){
    for (int i=0;i<VARIANCE.nvertices();i++){
      if(EXCL_ref->get_pvalue(i) && totalcounts.get_pvalue(i,d) > 0 && exclusion(i+1)){
	newval=sqrt(VARIANCE.get_pvalue(i,d)/(totalcounts.get_pvalue(i,d)));
	VARIANCE.set_pvalue(i,newval,d);
	if(_avdim==-1 || _avdim==d) meanvar+=newval;
	//cout << i << " newval " << newval << " VARIANCE.get_pvalue(d,i) " <<  VARIANCE.get_pvalue(d,i) << " AVERAGE.get_pvalue(d,i) " <<  AVERAGE.get_pvalue(d,i) << " totalcounts.get_pvalue(d,i) " << totalcounts.get_pvalue(d,i) << " meanvar " << meanvar << endl;

      }
      else VARIANCE.set_pvalue(i,0,d);
    }
  }

  for(int d=0;d<atlasdim;d++){
    for (int i=0;i<VARIANCE.nvertices();i++){
      if(VARIANCE.get_pvalue(i,d)>0 && exclusion(i+1))
	newval=(abs(AVERAGE.get_pvalue(i,d)))/VARIANCE.get_pvalue(i,d);
      else
	newval=0;

      Zscore.set_pvalue(i,newval,d);
    }
  }

int totdim;
  if(_avdim==-1) totdim=atlasdim;
  else totdim=1;
  cout << " meanval " << sumav/(VARIANCE.nvertices()*totdim) << endl;

  cout << " meanvar " << meanvar/(VARIANCE.nvertices()*totdim) << endl;

  AVERAGE.save(output+"AVERAGE.func");
  VARIANCE.save(output+"STD.func");
  Zscore.save(output+"ZSCORE.func");
  delete DATA;
}
