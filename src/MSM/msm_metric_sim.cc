#include <time.h>
#include <stdio.h>
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"

using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "msm_metric_sim  <surface1> <surface2> <outbase>  " << endl;
  cout << "Input are func files assumes surfaces have been resampled to same mesh   " << endl;
  cout << " options: " << endl;
  cout << "-compare X compare correlation to this file " << endl;
}

int main( int argc, char **argv )
{
  newmesh SURFACE1, SURFACE2, OUTSURF,COMP;
  vector<double> data1,data2;
  string outputname;
  sparsesimkernel<double> sim;
  double correlation,sum=0.0;
  bool _compare=false;
  bool _exclude=false;
  int ok=0;

  if(argc < 3){

    Usage();
    exit(0);
  }

  SURFACE1.load(argv[1],false,false);
  argc--; 
  argv++;
  SURFACE2.load(argv[1],false,false);
  argc--; 
  argv++;
  outputname=argv[1];
  argc--; 
  argv++;

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-compare") == 0)){
      argc--;
      argv++;
      _compare=true;
      COMP.load(argv[1],false,true);
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
  }
    if(SURFACE1.npvalues()!=SURFACE2.npvalues() || SURFACE1.get_dimension()!=SURFACE2.get_dimension()) { cout << " surfaces have different number of vertices, or features abort! " << endl; exit(1);}

  boost::shared_ptr<NEWMESH::newmesh> EXCL_IN, EXCL_REF;
  if(_exclude){
	Matrix DATAIN,DATAREF;
    DATAIN=SURFACE1.get_pvalues();
    EXCL_IN= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(SURFACE1,DATAIN,0,0))); 
    DATAREF=SURFACE2.get_pvalues();
    EXCL_REF= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(SURFACE1,DATAREF,0,0))); 
  }else{
	  EXCL_IN= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(SURFACE1));
	  EXCL_REF= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(SURFACE2));
	  for(int i=0;i<SURFACE1.npvalues();i++){
		EXCL_IN->set_pvalue(i,1);
		EXCL_REF->set_pvalue(i,1);
	}
  }
  
  OUTSURF=SURFACE1;

  OUTSURF.clear_data();
  ColumnVector allcorr(SURFACE1.npvalues()); 


  if(SURFACE1.get_dimension()==1){
	for (int i=0;i<SURFACE1.npvalues();i++){	
		if(EXCL_IN->get_pvalue(i)>0 && EXCL_REF->get_pvalue(i) > 0){
			data1.push_back(SURFACE1.get_pvalue(i));
			data2.push_back(SURFACE2.get_pvalue(i));
      }
	}
	correlation=sim.corr(data1,data2);
	cout << " mean=" << correlation << endl;
  }
  else{
	for (int i=0;i<SURFACE1.npvalues();i++){
		data1.clear();data2.clear();
    
		for (int j=0;j<SURFACE1.get_dimension();j++){
		data1.push_back(SURFACE1.get_pvalue(i,j));
		data2.push_back(SURFACE2.get_pvalue(i,j));
      
		}
		correlation=sim.corr(data1,data2);
		allcorr(i+1)= correlation;
    
		sum+=correlation;
	}
  
  OUTSURF.set_pvalues(allcorr);
  double mean,var=0.0;
  mean=sum/SURFACE1.npvalues();
  for (int i=0;i<SURFACE1.npvalues();i++)
    var+=(OUTSURF.get_pvalue(i) -mean)*(OUTSURF.get_pvalue(i) -mean);

  var/=SURFACE1.npvalues();

  cout << " mean=" << mean << endl;
  cout << " var=" << var << endl;

  OUTSURF.save(outputname+".func.gii");
 }
  if(_compare){
    for (int i=0;i<SURFACE1.npvalues();i++){
      double diff=OUTSURF.get_pvalue(i)-COMP.get_pvalue(i);
      OUTSURF.set_pvalue(i,diff);
    }
  OUTSURF.save(outputname+"DIFF.func.gii");

  }

}
