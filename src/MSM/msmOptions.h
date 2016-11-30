/*  msmOptions.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  

    Some sections of code inspired by A. Petrovic.
*/
/*  CCOPYRIGHT  */

#if !defined(msmOptions_h)
#define msmOptions_h

#include <string> 
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;
using namespace NEWMESH;

class msmOptions {

public:
  static msmOptions& getInstance();
  ~msmOptions() { delete gopt; }
  
 
  Option<bool>   help;
  Option<bool>   verbose;
  Option<bool>   printoptions;
  //  Option<bool> version;
  Option<bool> debug;
 

  Option<string> inputmesh; 
  Option<string> referencemesh;
  Option<string> inputanatmesh; 
  Option<string> referenceanatmesh; 
 

  Option<string> CMmatrixin;
  Option<string> CMmatrixref;
 

  Option<string> transformed_sphere; 
  Option<string> in_register; 
 
  Option<string> cfweight_in; 
  Option<string> cfweight_ref; 
  

  Option<string> outbase; 
  Option<string> outformat; 
  Option<string> parameters;  /// slowly replace most of these with the parameter file??
 
  Option<int> multiresolutionlevels;
  Option<float>  smoothoutput;

  bool parse_command_line(int argc, char** argv,Log& logger);
  
  vector<int> return_datarange(NEWMESH::newmesh);
 private:
  msmOptions();  
  const msmOptions& operator=(msmOptions&);
  msmOptions(msmOptions&);

  OptionParser options;
      
  static msmOptions* gopt;
  
};

 inline msmOptions& msmOptions::getInstance(){
 
   if(gopt == NULL)
     gopt = new msmOptions();
  
  return *gopt;
 }

inline msmOptions::msmOptions() :
		  help(string("-h,--help"), false,
		       string("display this message"),
		       false, no_argument),
		  verbose(string("-v,--verbose"), false,
		       string("switch on diagnostic messages"),
		       false, no_argument),
		  printoptions(string("-p,--printoptions"), false,
		       string("print configuration file options"),
		       false, no_argument),
		  // version(string("--version"), false,
		  //     string("report version informations"),
		  //     false, no_argument),
		  debug(string("--debug"), false,
		       string("run debugging or optimising options"),
			false, no_argument,false),
		  inputmesh(string("--inmesh"), string(""),
			    string("input mesh (available formats: VTK, ASCII, GIFTI). Needs to be a sphere"),
			    true , requires_argument),
		  referencemesh(string("--refmesh"), string(""),
				string("reference mesh (available formats: VTK, ASCII, GIFTI). Needs to be a sphere. If not included algorithm assumes reference mesh is equivalent input"),
				false , requires_argument),
		  inputanatmesh(string("--inanat"), string(""),
			    string("input anatomical mesh (available formats: VTK, ASCII, GIFTI). For Example, white, pial, midthickness  (must either supply both input and reference anatomical surfaces or none)"),
				false , requires_argument,false),
		  referenceanatmesh(string("--refanat"), string(""),
				string("reference mesh (available formats: VTK, ASCII, GIFTI). For Example, white, pial, midthickness (must either supply both input and reference anatomical surfaces or none)"),
				    false , requires_argument,false),
	
		  CMmatrixin(string("--indata"), string(""),
			     string("scalar or multivariate data for input - can be ASCII (.asc,.dpv,.txt) or GIFTI (.func.gii or .shape.gii) "),
			     false, requires_argument),
		  CMmatrixref(string("--refdata"), string(""),
			      string("scalar or multivariate data for reference - can be ASCII (.asc,.dpv,.txt) or GIFTI (.func.gii or .shape.gii) "),
			      false, requires_argument),
		 
		  transformed_sphere(string("--trans"), string(""),
			      string("\t Transformed source mesh (output of a previous registration). Use this to initiliase the current registration."),
			      false , requires_argument),
		  in_register(string("--in_register"), string(""),
			      string("\t Input mesh at data resolution. Used to resample data onto input mesh if data is supplied at a different resolution. Note this mesh HAS to be in alignment with either the input_mesh of (if supplied) the transformed source mesh. Use with supreme caution."),
			      false , requires_argument),
		  cfweight_in(string("--inweight"), string(""),
			      string("cost function weighting for input - weights data in these vertices when calculating similarity (ASCII or GIFTI). Can be multivariate provided dimension equals that of data "),
			      false , requires_argument),
		  cfweight_ref(string("--refweight"), string(""),
			       string("cost function weighting for reference - weights data in these vertices when calculating similarity (ASCII or GIFTI). Can be multivariate provided dimension equals that of data"),
			       false , requires_argument),
		  outbase(string("-o,--out"), string(""),
			  string("output basename"),
			  true, requires_argument),
		  outformat(string("-f,--format"), string("GIFTI"),
			  string("format of output files, can be: GIFTI, VTK, ASCII or ASCII_MAT (for full details of output file formats see MSM wiki)"),
			  false, requires_argument),
		  parameters(string("--conf"), string(""),
			     string("\tconfiguration file "),
			     false, requires_argument),
		  multiresolutionlevels(string("--levels"),0, 
					string("number of resolution levels (default = number of resolution levels specified by --opt in config file)"), 
					false, requires_argument),
		  smoothoutput(string("--smoothout"), 0, 
			       string("smooth tranformed output with this sigma (default=0)"), 
			       false, requires_argument),
		  
		  options("msm", "msm [options]\n")
{
  try {
   
        options.add(help); 
	options.add(verbose);
	options.add(printoptions);
	options.add(debug);
        options.add(inputmesh);
        options.add(referencemesh);
	options.add(inputanatmesh);
        options.add(referenceanatmesh);
	options.add(CMmatrixin);
	options.add(CMmatrixref);
	options.add(transformed_sphere);
	options.add(in_register);
	options.add(cfweight_in);
	options.add(cfweight_ref);
	options.add(outbase);
	options.add(outformat);
	options.add(parameters);
	options.add(multiresolutionlevels);
	options.add(smoothoutput);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     

     
}

inline bool msmOptions::parse_command_line(int argc, char** argv,Log& logger){
	 
	 
  for(int a = options.parse_command_line(argc, argv); a < argc; a++) ;
  if(!(printoptions.value()) ){
  if(help.value() || ! options.check_compulsory_arguments())
    {
      options.usage();
      //throw NEWMAT::Exception("Not all of the compulsory arguments have been provided");
      exit(2);
    }      
  
  logger.makeDir(outbase.value()+"logdir","MSM.log");
	  
 
  
  // do again so that options are logged
  for(int a = 0; a < argc; a++)
    logger.str() << argv[a] << " ";

  logger.str() << endl << "---------------------------------------------" << endl << endl;
  }  
  return true;
}



/* inline vector<int> msmOptions::return_datarange(NEWMESH::newmesh in){ */

/*   vector<int> V; */
/*   int range2; */
/*   int i; */
/*   int indexval=index.value(); */
 
/*   if(range.value()<1E-7){ */
/*     range2=in.nvertices(); */
/*     cout << " range2 " << range2 << endl;; */
/*   } */
/*   else{ */
/*     range2=range.value(); */
/*     cout << " here2 " << endl; */
/*   } */
/*   if(patch.value()==""){ */
/*       //      D.ReSize(range); */
/*       cout << " in return datarange" << range2 << " range.value() " << range.value() <<  endl; */

/*       // cout << " range " << range <<  " patch " <<opts.patch.value() << " V.Nrows() " << V.Nrows() << endl; */
/*       //shared(V,indexval,range2) private(i) */


/* #pragma omp parallel */

/*       { */
/*   // cout << " here " << endl; */
 
/* //cout << "omp_get_thread_num" << omp_get_thread_num() << endl; */
/* #pragma omp for nowait */
/*       for( i=0;i<range2;i++) */
/* 	V.push_back(indexval+i); // labels index from 1  */
/*       }     */
/*   } */
/*   else{ */
      
/*       Matrix patchm=read_ascii_matrix(patch.value());  */
      
/*       if(in.get_pvalue(0)){ */
/* 	for (int i=0;i<in.nvertices();i++) */
/* 	  in.set_pvalue(i,0); */
/*       } */

/*       cout << " label load not available for newmesh yet " << endl; */
/*       exit(0); */
/*       // in.load_fs_label(patch.value()); */
/*       // cout << "here " << endl; */
/*       int ind=0; */
/*       for(vector<boost::shared_ptr<Mpoint> >::const_iterator i=in.vbegin();i!=in.vend();i++){ */
/* 	//  cout << " here 1 " << (*i)->get_value() << endl; */
/* 	if(in.get_pvalue((*i)->get_no()) > 0){ */
	  
/* 	  V.push_back((*i)->get_no()+1); /// labels indexing runs from 1 !!! */
/* 	  ind++; */
/* 	  /// note this was in error before inasmuch as if I was using a patch I would count indexes from 0 but the datarange assumed indexing was running from 1 */
/* 	// 	cout << ind << " V(ind) " << V(ind) << endl; */
/*        } */
/*     } */
/*   } */

/*   return V; */

/* } */

#endif

