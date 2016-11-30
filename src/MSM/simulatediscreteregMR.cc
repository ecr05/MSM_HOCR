/*  simulatediscreteregMR.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "MeshReg/meshreg.h"

#include <time.h>
#include <stdio.h>

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace MESHREG;
using namespace NEWMESH;

#define sigma 1

Pt get_square(NEWMESH::newmesh& IN){
  cout << " get sqaure " << endl;
  double theta,phi,mindist=RAD;
  NEWMESH::Pt ci,origin,closest;
  SpMat<double> data(1,IN.nvertices());
  double disttoorigin;
  resampler R;
  NEWMESH::newmesh tmpmesh=IN;
  R.set_method("GAUSSIAN");
  origin.X=RAD*sin(0.5*3.14)*cos(0);
  origin.Y=RAD*sin(0.5*3.14)*sin(0);
  origin.Z=RAD*cos(0.5*3.14);

  for (int i=0;i<IN.nvertices();i++){
    ci=IN.get_coord(i);
    phi=atan2(ci.Y,ci.X);
    theta=acos(ci.Z/RAD);
    disttoorigin=(ci-origin).norm();
    if(disttoorigin<mindist){
      mindist=disttoorigin;
      closest=ci;
    }
    if(theta>=0.3*3.14 && theta<=0.7*3.14){
      if(phi>=-0.2*3.14 && phi< 0.2*3.14)  {
	IN.set_pvalue(i,10);
      }
    }
  }

  R.resample_scalar(IN,tmpmesh,sigma);
  cout << " here " << endl;
  for (int i=0;i<IN.nvertices();i++){
    if(IN.get_pvalue(i)>1)
      data.Set(1,i+1,IN.get_pvalue(i));
    else 
      IN.set_pvalue(i,0);
  }
  data.Print("CMfile_reference");
  return closest;
}

void get_circle(NEWMESH::newmesh& IN, Pt p){
  cout << " get circle " << endl;
  double val,dist;
  NEWMESH::Pt ci;
  SpMat<double> data(1,IN.nvertices());
  NEWMESH::newmesh tmpmesh=IN;
  resampler R;
  R.set_method("GAUSSIAN");
  for (int i=0;i<IN.nvertices();i++){
    ci=IN.get_coord(i);
    dist=(ci-p).norm();
    val=std::ceil(100*exp(-(dist*dist)/(2*(sigma)*(sigma))));
    if(val>10){
     
      IN.set_pvalue(i,10);
    }
  }

  R.resample_scalar(IN,tmpmesh,0.1*sigma);
  for (int i=0;i<IN.nvertices();i++){
    if(IN.get_pvalue(i)>1)
      data.Set(1,i+1,IN.get_pvalue(i));
    else 
      IN.set_pvalue(i,0);
  }
 
  data.Print("CMfile_input");
}

void check_sim(NEWMESH::newmesh &IN){
  cout << " in check sim " << IN.nvertices() << endl;
  ///// creat a check patten by dividing theta range into 3 and phi range into 6 
  Matrix DATA(1,IN.nvertices()); DATA=0;
  Matrix checksize(3,2);
  checksize(1,1)=5;  checksize(1,2)=15;
  checksize(2,1)=7;  checksize(2,2)=21;
  checksize(3,1)=3;  checksize(3,2)=9;

  for (int row=1;row<=1;row++){
    Grid checks;

    checks.Initialize(IN,checksize(row,1),checksize(row,2));

    vector<vector<int> > cells=checks.return_cell2mesh();
    //  vector<int> val;
    int val=1;

    for(unsigned int i=0;i<cells.size();i++){
      //cout << i << " cells[i].size() " << cells[i].size() << endl;
      for (unsigned int j=0;j<cells[i].size();j++){
	//cout <<row << " " << j << " " <<  cells[i][j]+1 << " DATA.NCols() " << DATA.Ncols() << " i " << i << " cells.size() " << cells.size() << endl;
	DATA(row,cells[i][j]+1)=val;
      }
      if(val==4) val=1;
      else val++;
    }
  }
  // newmesh ICO;
  //ICO.load("/Users/emmar/MRIDATA/Phase2/Scalar/NATIVE/SURF/100307.L.myelin.native.surf.gii");
  
  //Matrix DATA2=ICO.get_pvalues(); 
  //resampler R;
  //R.set_method("GAUSSIAN");
  //R.resampledata(ICO,IN,DATA2,sigma); 
  //cout << IN.nvertices() << " DATA2.Ncols() " <<  DATA2.Ncols() << " DATA.Ncols() " <<  DATA.Ncols() << " DATA2.Nrows() " <<  DATA2.Nrows() << " DATA.Nrows() " <<  DATA.Nrows() << endl;
  //for (int i=1;i<IN.nvertices();i++){
  //cout << i << endl;
  // cout << "DATA2(1,i) " << DATA2(1,i) << endl;
  //cout << "DATA(2,i) " << DATA(2,i) << endl;

  //DATA(2,i)=DATA2(1,i);
  //}
  IN.set_pvalues(DATA);

  IN.save("INPUTchecks2.surf");IN.save("INPUTchecks2.func");
  
}


void initialize_data(NEWMESH::newmesh &IN, NEWMESH::newmesh &REF, const int res, const double ang, const double _sigma){
  /////////// initialize data /////////////////////////
  IN.make_mesh_from_icosa(res); true_rescale(IN,RAD); 
  REF=IN;
  NEWMESH:: newmesh ICO,REG,tmpM;
  SpMat<double> tmp(1,IN.nvertices()),tmp2(1,REF.nvertices());
  //  
  cout << "REF.mvertices() " << REF.nvertices() << endl;

  tmpM=IN;
  //Pt centre=get_square(REF);
  //cout << " centre.X " << centre.X <<  " centre.Y " << centre.Y <<  " centre.Z " << centre.Z << endl;
  //get_circle(IN,centre);
  
  check_sim(IN);
  
  REF=IN;

  
  cout << " init " << IN.nvertices() << " ref.nvertices() " <<REF.nvertices() <<  endl;
  ICO.load("/Users/emmar/MRIDATA/BadSubjects/205725/MNINonLinear/Native/205725.L.sphere.native.surf.gii");
  REG.load("/Users/emmar/MRIDATA/BadSubjects/205725/MNINonLinear/Native/205725.L.sphere.reg.native.surf.gii");
  // REG.load("/Users/emmar/MRIDATA/Myelindatanew/CP10051_v3/CP10051_v3.L.MyelinMap_on_sphere.reg.reg_LR.native.asc");
  double MVDinput=Calculate_MVD(IN);
  //cout <<  MVDinput << " sigma " <<  2*asin(MVDinput/(2*RAD))<< endl;

  tmpM=REF;
  cout << "before upsample " << endl;
  barycentric_mesh_interpolation(tmpM,ICO,REG);
  //upsample_transform_RBF(tmpM,ICO,REG,4*asin(ang*MVDinput/(2*RAD)));
  cout << "after upsample " << endl;
  resampler R;
  R.set_method("NN");
  Matrix DATA=IN.get_pvalues();
  R.resampledata(tmpM,IN,DATA,sigma); 
  //  R.resample_scalar(REF,tmpM,_sigma);
  
  REF.set_pvalues(DATA);
  REF.save("REFchecks2.func");
  cout << "REF.mvertices() " << REF.nvertices() << endl;
  
}

NEWMESH::newmesh create_conf_weighting(newmesh M){

  newmesh conf=M;
  Matrix DATA= M.get_pvalues();
  
  for (int i=1;i<=M.nvertices();i++)
    DATA(1,i)=1;
 
  // for (int i=1;i<=M.nvertices();i++){
  // if(DATA(2,i)>1.27) DATA(1,i)=1;
  // else DATA(1,i)=0;
  //}

  conf.set_pvalues(DATA);
  conf.save("conf_weighting.func.gii");

  return conf;
}

void smooth(NEWMESH::newmesh &M, const double & _sigma){

  resampler R;
  NEWMESH::newmesh target=M;

  R.set_method("GAUSSIAN");
  R.resample_scalar(M,target,_sigma);

}

int main( int argc, char **argv )
{
  NEWMESH::newmesh conf,INPUT,INPUTorig,REF;
  resampler resample;
  char filename[1000];
  int maxiters=10;
  vector<int> _iters;
  vector<double> _sigma;
  vector<double> ang;
  vector<int> _datares;
  MeshReg MR;
  ofstream out;
  int levels=3;

  out.open("simconf");
  out << "--simval=1,2,2,2" << endl;
  out << "--it=50,3,3,3" << endl;
  out << "--CPgrid=0,2,3,4" << endl;
  out << "--datagrid=5,5,5,6" << endl;
  out << " --lambda=0,0.005,0.003,0.002 " << endl;
  out << " --sigma_in=5,3,2,1 " << endl;
  out << " --aKNN=20,20,20,20 " << endl;
  out << " --opt=AFFINE,DISCRETE,DISCRETE,DISCRETE " << endl;
  //out << " --HO " << endl;

  out.close();

  _datares.push_back(6);
  _sigma.push_back(0);
  ang.push_back(1.2);

  initialize_data(INPUT,REF, _datares[0], ang[0],_sigma[0]);
  conf=create_conf_weighting(REF);

  INPUT.save("CMfile_input.func");
  INPUT.save("CMfile_input.surf");
  REF.save("CMfile_reference.func");
   
   //// INITIALISE VARIABLES //////////////////// 
  
  MR.set_verbosity(true);

  MR.set_input("CMfile_input.surf.gii");
  MR.set_reference("CMfile_input.surf.gii");
  MR.set_outdir("./");
  
  MR.set_CMpathin("CMfile_input.func.gii");
  MR.set_CMpathref("CMfile_reference.func.gii");

  // MR.set_reference_cfweighting("conf_weighting.func.gii");
  // MR.resample(true);
  MR.run_multiresolutions(levels,sigma,"simconf");  
    
  
   
}
  
  

  /////// nxt need potentials 
  
  
  
