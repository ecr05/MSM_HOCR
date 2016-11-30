/*  msm.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "MeshReg/meshreg.h"
#include "msmOptions.h"



#include "newran.h"
#include <time.h>

using namespace std;
using namespace NEWMESH;
using namespace NEWRAN;
using namespace MESHREG;

msmOptions* msmOptions::gopt = NULL;


 
int main(int argc, char *argv[]) {

  msmOptions& opts = msmOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  string CMpathin;

  opts.parse_command_line(argc,argv,logger);
  //opts.parse_command_line(argc,argv);
  MeshReg MR;

  //// INITIALISE VARIABLES //////////////////// 
  if(opts.printoptions.value()){ MR.print_config_options(); return 0; }
  // if(opts.version.value()){ cout << " MSM_XS version 0.0" << endl; return 0; }
  if(opts.verbose.value())    MR.set_verbosity(opts.verbose.value());
  if(opts.debug.value())    MR.set_debug(opts.debug.value());

  MR.set_input(opts.inputmesh.value());
  if(opts.referencemesh.value()=="")  MR.set_reference(opts.inputmesh.value());
  else MR.set_reference(opts.referencemesh.value());

  ////// add anatomical meshes for stress and strain ////////////
  if(opts.inputanatmesh.value()!=""){
    if(opts.referenceanatmesh.value()=="") { cout << " Error: must supply both anatomical meshes or none "<< endl; exit(1);}
    MR.set_anatomical(opts.inputanatmesh.value(),opts.referenceanatmesh.value());
  }

  MR.set_outdir(opts.outbase.value());
  MR.set_output_format(opts.outformat.value());
 

  if(opts.transformed_sphere.set()) MR.set_transformed(opts.transformed_sphere.value());  
  if(opts.cfweight_in.set())  MR.set_input_cfweighting(opts.cfweight_in.value());
  if(opts.cfweight_ref.set()) MR.set_reference_cfweighting(opts.cfweight_ref.value());

  if(opts.in_register.set()){
    //// if data is supplied at a different resolution to the input or transformed mesh it is necessary to resample i.e. HCP 32K to native
    //// the in_register sphere HAS to be the sphere on which the data is supplied and the input or transformed mesh MUST be in alignment with it.
    newmesh in_register, target;
    resampler R; R.set_method("ADAP_BARY");
    boost::shared_ptr<BFMatrix > DATA;
    char filename[1000];
    if(opts.transformed_sphere.set()) target.load(opts.transformed_sphere.value());
    else target.load(opts.inputmesh.value());

    in_register.load(opts.in_register.value());
    set_data(opts.CMmatrixin.value(),DATA,in_register);
    R.resampledata(in_register,target,DATA,0.1);
    sprintf(filename,"%sinput_in_register.func.gii",opts.outbase.value().c_str());
    target.set_pvalues(DATA->AsMatrix());
    target.save(filename);
    CMpathin=filename;
   
  }else CMpathin=opts.CMmatrixin.value();
 
  MR.set_CMpathin(CMpathin);
  MR.set_CMpathref(opts.CMmatrixref.value());
 

  MR.run_multiresolutions(opts.multiresolutionlevels.value(),opts.smoothoutput.value(),opts.parameters.value());
 
  if(opts.in_register.set()){
    /// if data is supplied at a lower mesh resolution, then resample final warp accordingly. This is typically a HCP formating issue
    newmesh ORIG, FINAL, in_register;
    char filename[1000];

    if(opts.transformed_sphere.set()) ORIG.load(opts.transformed_sphere.value());
    else ORIG.load(opts.inputmesh.value());

    in_register.load(opts.in_register.value());
    sprintf(filename,"%ssphere.reg%s",opts.outbase.value().c_str(),MR.get_surf_format().c_str());
   
    FINAL.load(filename);
    barycentric_mesh_interpolation(in_register,ORIG,FINAL); 
    sprintf(filename,"%ssphere.in_register.reg.%s",opts.outbase.value().c_str(),MR.get_surf_format().c_str());
    in_register.save(filename);

  }
}
