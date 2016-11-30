/*  msmresamplemetric.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
/* this program is designed to downsample freesurfer label files to be used in combination with the SPH6.vtk or other downsampled meshes*/

#include "newmesh/newmesh.h"
#include "MeshReg/meshreg.h"
#include "miscmaths/SpMat.h"

using namespace MESHREG;

void Usage()
{
  cout << "msmresample  <inputmesh> <output> <-option>   " << endl;
  cout << "-option " << endl;
  cout << " -project - project final result back down onto surface (requires argument) " << endl;
  cout << " -labels - load labels file for input (inc. .func. and .shape, requires argument)" << endl;
  cout << " -datamat - multivariate data matrix (requires argument)" << endl; 
  cout << " -barycentric use barycentric interpolation" << endl;
  cout << " -adap_bary use adaptive barycentric interpolation" << endl;
  cout << " -linear - linear interpolation kernel with kernel size X (requires argument)" << endl;
  cout << " -gaussian - gaussian interpolation kernel with std deviation X (requires argument) " << endl;
  cout << " -normalize  - normlalize intensity range to target (requires argument)" << endl;
  cout << " -excl exclude the area of the cut form contributing to the resampling" << endl;
  cout << " -arealdist estimate areal distortion" << endl;
}


int main(int argc, char **argv){

  Matrix M;
  int ok,getwm;
  bool barycentric,adap_barycentric,_normalize;
  resampler resample;
  //MeshReg MR; // only used for read list -> move to resample?
  NEWMESH::newmesh in,SPH,REG,ICO,wm,HISTTARG;
  bool _outputlabel,_transform,_project,_exclude;
  bool _save_rel,_datamatrix, _datamatsp;
  string output;
  double thr;
  double sigma,rad;
  boost::shared_ptr<RELATIONS>  _rel;
  string rel_out;
  boost::shared_ptr<BFMatrix > datamat;
  bool labels=false;
  // BFMatrix *datamat;
  

  if(argc < 3){

    Usage();
    exit(0);
  }
 
  
  in.load(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;

  recentre(in);
  
  getwm=0; thr=0.0; sigma=1.0;rad=100.0;
  _transform=false;_project=false;_outputlabel=false;
  _datamatrix=false; _save_rel=false;  _datamatsp=false;
  _exclude=false; resample.set_method("NN");
  barycentric=false; adap_barycentric=false;  _normalize=false;
  
 
  while (argc > 1) {
    ok = 0;
    cout << argv[1] << endl;
   
    if((ok == 0) && (strcmp(argv[1], "-project") == 0)){
      argc--;
      argv++;
      _project=true;
      SPH.load(argv[1]);
      recentre(SPH);
      argc--;
      argv++;
      ok = 1;
    }else if((ok == 0) && (strcmp(argv[1], "-labels") == 0)){
      argc--;
      argv++;
      labels=true;
   
      in.load(argv[1],false,false);
      Matrix M=in.get_pvalues();
      datamat = boost::shared_ptr<BFMatrix >(new FullBFMatrix (M));

      argc--;
      argv++;
      ok = 1;
    }
  else if((ok == 0) && (strcmp(argv[1], "-datamat") == 0)){
      argc--;
      argv++;
      _datamatrix=true;
      Matrix M=read_ascii_matrix(argv[1]);
      datamat = boost::shared_ptr<BFMatrix >(new FullBFMatrix (M));
      in.set_pvalues(M);
      argc--;
      argv++;
      ok = 1;
    }
  else if((ok == 0) && (strcmp(argv[1], "-datamatsp") == 0)){
      argc--;
      argv++;
      _datamatsp=true;
      SpMat<double> M(argv[1]);
      cout << " after load sparsemat " << endl;
      datamat = boost::shared_ptr<BFMatrix >(new SparseBFMatrix<double> (M));
      in.set_pvalues(datamat->AsMatrix());
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-normalize") == 0)){
      argc--;
      argv++;
      _normalize=true;
      HISTTARG.load(argv[1],false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-wmmesh") == 0)){
      argc--;
      argv++;
      getwm=1;
      wm.load(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }   
    else if((ok == 0) && (strcmp(argv[1], "-excl") == 0)){
      argc--;
      argv++;
      _exclude=true;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-barycentric") == 0)){
      argc--;
      argv++;
      barycentric=true;
      resample.set_method("BARYCENTRIC");  
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-adap_bary") == 0)){
      argc--;
      argv++;
      adap_barycentric=true;
      resample.set_method("ADAP_BARY");  
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-linear") == 0)){
      argc--;
      argv++;
      resample.set_method("LINEAR");
      sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-gaussian") == 0)){
      argc--;
      argv++;
      resample.set_method("GAUSSIAN");
      sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    } 
    else{cout << " option doesn't exist " << endl; exit(1);}

   
  }


 
  boost::shared_ptr<NEWMESH::newmesh> EXCL_IN, EXCL_REF;
  Matrix DATAIN;
  if(_exclude){
    DATAIN=in.get_pvalues();
    EXCL_IN= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(in,DATAIN,0,0))); 
  }
 
  if(_normalize)  {
    Matrix DATAREF;
    boost::shared_ptr<BFMatrix > BFIN,BFREF;
    DATAREF=HISTTARG.get_pvalues();
    if(!datamat.get()){
      Matrix M=in.get_pvalues();
      datamat = boost::shared_ptr<BFMatrix >(new FullBFMatrix (M));

    }

    if(_exclude){
      cout << SPH.nvertices() << " " << DATAREF.Ncols() << endl;
      EXCL_REF= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(HISTTARG,DATAREF,0,0))); 
      
    }
    
    BFREF=boost::shared_ptr<BFMatrix > (new FullBFMatrix (DATAREF)); //
    multivariate_histogram_normalization(*datamat,*BFREF,EXCL_IN,EXCL_REF,false);
    in.set_pvalues(datamat->AsMatrix());

  }

 

  if(_project){
    for (int i=0;i<SPH.nvertices();i++)
      SPH.set_pvalue(i,0); 
    
    if(datamat.get()){
      resample.resampledata(in,SPH,EXCL_IN,datamat,sigma); 
     
           
    }
    else{
      resample.resample_scalar(in,SPH,sigma,EXCL_IN);
      
    }
  

    if(labels){
      SPH.set_pvalues(datamat->AsMatrix());
      SPH.save(output+".func");
    }else{
      datamat->Print(output);
   
    }
    
   

    if(getwm){
      
      for (int i = 0; i < SPH.nvertices(); i++)
	wm.set_pvalue(i,SPH.get_pvalue(i));
	
    }

 

    if(getwm)
      wm.save(output +"targetmesh");
  
  }
 
}
