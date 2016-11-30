/*  msm_group.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "MeshReg/groupmeshreg.h"
#include "msmGOptions.h"
//#include "QPBO/QPBO.h"
//#include "HOCR/HOCR.h"

//#include "HOCR/QPBO.h"


#include "newran.h"
#include <time.h>

using namespace std;
using namespace NEWMESH;
using namespace NEWRAN;
using namespace MESHREG;

msmGOptions* msmGOptions::gopt = NULL;


 
int main(int argc, char *argv[]) {

  msmGOptions& opts = msmGOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  string CMpathin;

  opts.parse_command_line(argc,argv,logger);
  //opts.parse_command_line(argc,argv);
  GroupMeshReg MR;

  //// INITIALISE VARIABLES //////////////////// 
  if(opts.printoptions.value()){ MR.print_config_options(); return 0; }
  if(opts.version.value()){ cout << " MSM_XS version 0.0" << endl; return 0; }
  if(opts.verbose.value())    MR.set_verbosity(opts.verbose.value());
  if(opts.debug.value())    MR.set_debug(opts.debug.value());

  MR.set_inputs(opts.meshes.value()); 
  MR.set_template(opts.templatemesh.value()); 
  MR.set_data_list(opts.data.value());

  MR.set_outdir(opts.outbase.value());

  //  MR.set_matlab(opts.L1matlabpath.value()); 

  MR.run_multiresolutions(opts.multiresolutionlevels.value(),opts.smoothoutput.value(),opts.parameters.value());
 
}
