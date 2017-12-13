/*  Fusion.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */



#include <iostream>


#ifdef HAS_HOCR
#include "ELC/ELC.h"
#endif



#ifdef HAS_TBB

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"
using namespace tbb;
#endif

#if !defined(__DiscreteModel_h)
#define __DiscreteModel_h
#include "DiscreteOpt/DiscreteModel.h"
#endif

#include <FastPD/FastPD.h>

using namespace FPD;
using namespace detail;
using namespace DISCRETEOPT;


typedef	double REAL;
struct UnaryData	{ REAL buffer[2]; };
struct PairData		{ REAL buffer[4]; };
struct TripletData	{ REAL buffer[8]; };
struct QuartetData	{ REAL buffer[16]; };


namespace ELCReduce {
  template<typename T> class PBF;
}

template<typename T> class QPBO;


namespace MESHREG {

#ifdef HAS_TBB
  class computeTripletCosts{
    const boost::shared_ptr<DiscreteModel> my_energy;
    const int* my_trip;
    const int* my_labels;
    const int my_current_label;

    std::vector<TripletData> & my_DATA;

  public:
    // constructor copies the arguments into local storage
  computeTripletCosts(const boost::shared_ptr<DiscreteModel> & energy, const int *trip, const  int *labels, const int & label, std::vector<TripletData> & DATA) :
    my_energy(energy), my_trip(trip), my_labels(labels), my_current_label(label),my_DATA(DATA)
    {}
    // overload () so it does a vector multiply
    void operator() (const blocked_range<int> &r) const {
      for(size_t triplet=r.begin(); triplet!=r.end(); triplet++){
	const int nodeA = my_trip[triplet*3];
	const int nodeB = my_trip[triplet*3+1];
	const int nodeC = my_trip[triplet*3+2];

	my_DATA[triplet].buffer[0] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_labels[nodeB],my_labels[nodeC]);	//000
	my_DATA[triplet].buffer[1] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_labels[nodeB],my_current_label);			//001
	my_DATA[triplet].buffer[2] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_current_label,my_labels[nodeC]);			//010
	my_DATA[triplet].buffer[3] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_current_label,my_current_label);						//011
	my_DATA[triplet].buffer[4] = my_energy->computeTripletCost(triplet,my_current_label,my_labels[nodeB],my_labels[nodeC]);			//100
	my_DATA[triplet].buffer[5] = my_energy->computeTripletCost(triplet,my_current_label,my_labels[nodeB],my_current_label);						//101
	my_DATA[triplet].buffer[6] = my_energy->computeTripletCost(triplet,my_current_label,my_current_label,my_labels[nodeC]);						//110
	my_DATA[triplet].buffer[7] = my_energy->computeTripletCost(triplet,my_current_label,my_current_label,my_current_label);
      }
    }

  };

 class computeQuartetCosts{
    const boost::shared_ptr<DiscreteModel> my_energy;
    const int* my_quartet;
    const int* my_labels;
    const int my_current_label;

    std::vector<QuartetData> & my_DATA;

  public:
    // constructor copies the arguments into local storage
 computeQuartetCosts(const boost::shared_ptr<DiscreteModel> & energy, const int *quart, const  int *labels, const int & label, std::vector<QuartetData> & DATA) :
    my_energy(energy), my_quartet(quart), my_labels(labels), my_current_label(label),my_DATA(DATA)
    {}
    // overload () so it does a vector multiply
    void operator() (const blocked_range<int> &r) const {
      for(size_t quartet=r.begin(); quartet!=r.end(); quartet++){
	const int nodeA = my_quartet[quartet*4];
	const int nodeB = my_quartet[quartet*4+1];
	const int nodeC = my_quartet[quartet*4+2];
	const int nodeD = my_quartet[quartet*4+3];

	my_DATA[quartet].buffer[0]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[1]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[2]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[3]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_current_label,my_current_label);

	my_DATA[quartet].buffer[4]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[5]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[6]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[7]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_current_label,my_current_label);

	my_DATA[quartet].buffer[8]  = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[9]  = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[10] = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[11] = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_current_label,my_current_label);

	my_DATA[quartet].buffer[12] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[13] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[14] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[15] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_current_label,my_current_label);


      }
    }

  };

#endif

 enum Reduction
  {
    ELC_HOCR,
    ELC_APPROX,
    HOCR
  };
  class Fusion
  {
  public:
    /**
     * Constructor.
     */
    Fusion() {}
    /**
     * Destructor.
     */
    ~Fusion() {}

    template<typename OPTIMIZER> static void reduce_and_convert(ELCReduce::PBF<REAL>&,OPTIMIZER&, Reduction);
    /**
     * Runs the optimization.
     */
    static double optimize(boost::shared_ptr<DiscreteModel> energy, Reduction reductionMode, bool debug=false);
  };

}
