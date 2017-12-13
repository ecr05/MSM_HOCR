/*  Fusion.cpp

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  CCOPYRIGHT  */
#include "Fusion.h"
//#include "ELC/ELC.h"
//#ifdef HAS_QPBO
//#include "QPBO/QPBO.hpp"
//#else
//#include <FastPD/FastPD.h>
//#endif

#include <iostream>
#include <time.h>


#define NUM_SWEEPS 2
#define MAX_IMPROVEMENTS 0

using namespace ELCReduce;

namespace MESHREG {


  template<typename OPTIMIZER> void Fusion::reduce_and_convert(PBF<REAL>& pbf, OPTIMIZER& MODEL, Reduction mode)
  {
      switch(mode)
      {
      case ELC_HOCR:
      {
        PBF<REAL> qpbf;

        pbf.reduceHigher();
	pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce the remaining higher-order terms using HOCR adding auxiliary variables
	qpbf.convert(MODEL, qpbf.maxID()+1);
        pbf.clear();
        qpbf.clear();
        break;
      }
      case ELC_APPROX:
      {
        //PBF<REAL> qpbf;
	//qpbf = pbf;
	pbf.reduceHigherApprox();

        pbf.convert(MODEL, pbf.maxID()+1);
        pbf.clear();

        break;
      }
      case HOCR:
      {
        PBF<REAL> qpbf;
        pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce to Quadratic pseudo-Boolean function using HOCR.
        qpbf.convert(MODEL, qpbf.maxID()+1);
        pbf.clear();
        qpbf.clear();
        break;

      }

      }
    }

//====================================================================================================================//
  double Fusion::optimize( boost::shared_ptr<DiscreteModel> energy, Reduction reductionMode, bool verbose )
{
  const int numNodes    = energy->getNumNodes();
  const int numPairs    = energy->getNumPairs();
  const int numTriplets = energy->getNumTriplets();
  const int numQuartets = energy->getNumQuartets();
  //const int numthreads  = energy->getNumthreads();
  const int *pairs    = energy->getPairs();
  const int *triplets = energy->getTriplets();
  const int *quartets = energy->getQuartets();

  //const int numGraphNodes = numNodes + numTriplets + 5 * numQuartets;
  //const int numGraphEdges = numPairs + 3 * numTriplets + 17 * numQuartets;
  //  cout << "  numNodes  " <<  numNodes  << " numPairs " << numPairs << " " << numTriplets << endl;
  #ifdef HAS_QPBO
  if(verbose) cout << " has QPBO " << endl;
  QPBO<REAL> qpbo(numGraphNodes,numGraphEdges) ;
  #endif
  boost::shared_ptr<DiscreteModelDummy> FPDMODEL=boost::shared_ptr<DiscreteModelDummy>(new DiscreteModelDummy());

  //double avg_unlabeled_ratio = 0;

  const int numSweeps = NUM_SWEEPS;
  const int numLabels = energy->getNumLabels();
  int *labeling = energy->getLabeling();

  //energy->applytestLabeling(labeling,-1);
  double initEnergy = energy->evaluateTotalCostSum();

  double lastEnergy = initEnergy;
  double unlabeledAvg = 0;//unlabeledMax = 0, 
  double sumlabeldiff=0;


  for(int sweep = 0; sweep < numSweeps; ++sweep)
  {
    for(int label = 0; label < numLabels; ++label)
    {
      //cout << "1" << label <<  endl;
      sumlabeldiff=0;

      PBF<REAL> pbf;
      int improveCounter = 0;

      double ratioUnlabeled = 0;
      int nodesChanged = 0;



      std::vector<UnaryData> unary_data(numNodes);

      for(int node = 0; node < numNodes; ++node)
	{
	  unary_data[node].buffer[0] = energy->computeUnaryCost(node,labeling[node]);	//0
	  unary_data[node].buffer[1] = energy->computeUnaryCost(node,label);				//1

	  sumlabeldiff+=abs(label-labeling[node]);
	}

      //cout << " 2 " << endl;
      if(sumlabeldiff>0){
	for(int node = 0; node < numNodes; ++node)
	  {

	    pbf.AddUnaryTerm(node, unary_data[node].buffer[0], unary_data[node].buffer[1]);
	  }


	std::vector<PairData> pair_data(numPairs);

	for(int pair = 0; pair < numPairs; ++pair)
	  {
	    const int nodeA = pairs[pair*2];
	    const int nodeB = pairs[pair*2+1];

	      pair_data[pair].buffer[0] = energy->computePairwiseCost(pair,labeling[nodeA],labeling[nodeB]);	//00
	      pair_data[pair].buffer[1] = energy->computePairwiseCost(pair,labeling[nodeA],label);			//01
	      pair_data[pair].buffer[2] = energy->computePairwiseCost(pair,label,labeling[nodeB]);			//10
	      pair_data[pair].buffer[3] = energy->computePairwiseCost(pair,label,label);
					//11
	  }

	for(int pair = 0; pair < numPairs; ++pair)
	  {
	    const int nodeA = pairs[pair*2];
	    const int nodeB = pairs[pair*2+1];
	    int node_ids[2] = { nodeA, nodeB};

	    pbf.AddPairwiseTerm(node_ids[0], node_ids[1], pair_data[pair].buffer[0], pair_data[pair].buffer[1], pair_data[pair].buffer[2], pair_data[pair].buffer[3]);
	  }
	//cout << " 3 " << endl;
	std::vector<TripletData> triplet_data(numTriplets);
	//std::vector<TripletData> triplet_data_par(numTriplets);
#ifdef HAS_TBB

	tbb::task_scheduler_init init(numthreads);
	parallel_for( blocked_range<int>(0,numTriplets), computeTripletCosts(energy,triplets,labeling,label,triplet_data) );

#else

	for (int triplet = 0; triplet < numTriplets; ++triplet)
	  {
	    const int nodeA = triplets[triplet*3];
	    const int nodeB = triplets[triplet*3+1];
	    const int nodeC = triplets[triplet*3+2];

	    triplet_data[triplet].buffer[0] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],labeling[nodeC]);	//000
	    triplet_data[triplet].buffer[1] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],label);			//001
	    triplet_data[triplet].buffer[2] = energy->computeTripletCost(triplet,labeling[nodeA],label,labeling[nodeC]);			//010
	    triplet_data[triplet].buffer[3] = energy->computeTripletCost(triplet,labeling[nodeA],label,label);						//011
	    triplet_data[triplet].buffer[4] = energy->computeTripletCost(triplet,label,labeling[nodeB],labeling[nodeC]);			//100
	    triplet_data[triplet].buffer[5] = energy->computeTripletCost(triplet,label,labeling[nodeB],label);						//101
	    triplet_data[triplet].buffer[6] = energy->computeTripletCost(triplet,label,label,labeling[nodeC]);						//110
	    triplet_data[triplet].buffer[7] = energy->computeTripletCost(triplet,label,label,label);								//111

	  }

#endif
	//cout << " 4 " << endl;
	for (int triplet = 0; triplet < numTriplets; ++triplet)
	  {
	    const int nodeA = triplets[triplet*3];
	    const int nodeB = triplets[triplet*3+1];
	    const int nodeC = triplets[triplet*3+2];
	      int node_ids[3] = { nodeA, nodeB, nodeC };
	      pbf.AddHigherTerm(3, node_ids, triplet_data[triplet].buffer);

	  }


	std::vector<QuartetData> quartet_data(numQuartets);
#ifdef HAS_TBB

	parallel_for( blocked_range<int>(0,numQuartets), computeQuartetCosts(energy,quartets,labeling,label,quartet_data) );

#else

	for (int quartet = 0; quartet < numQuartets; ++quartet)
	  {
	    const int nodeA = quartets[quartet*4];
	    const int nodeB = quartets[quartet*4+1];
	    const int nodeC = quartets[quartet*4+2];
	    const int nodeD = quartets[quartet*4+3];
	    if(! (label==labeling[nodeA] && label==labeling[nodeB] && label==labeling[nodeC] && label==labeling[nodeD])){

	    quartet_data[quartet].buffer[0]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[1]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],labeling[nodeC],label);
	    quartet_data[quartet].buffer[2]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],label,labeling[nodeD]);
	    quartet_data[quartet].buffer[3]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],label,label);

	    quartet_data[quartet].buffer[4]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[5]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,labeling[nodeC],label);
	    quartet_data[quartet].buffer[6]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,label,labeling[nodeD]);
	    quartet_data[quartet].buffer[7]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,label,label);
	    quartet_data[quartet].buffer[8]  = energy->computeQuartetCost(quartet,label,labeling[nodeB],labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[9]  = energy->computeQuartetCost(quartet,label,labeling[nodeB],labeling[nodeC],label);
	    quartet_data[quartet].buffer[10] = energy->computeQuartetCost(quartet,label,labeling[nodeB],label,labeling[nodeD]);
	    quartet_data[quartet].buffer[11] = energy->computeQuartetCost(quartet,label,labeling[nodeB],label,label);
	    quartet_data[quartet].buffer[12] = energy->computeQuartetCost(quartet,label,label,labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[13] = energy->computeQuartetCost(quartet,label,label,labeling[nodeC],label);
	    quartet_data[quartet].buffer[14] = energy->computeQuartetCost(quartet,label,label,label,labeling[nodeD]);
	    quartet_data[quartet].buffer[15] = energy->computeQuartetCost(quartet,label,label,label,label);
	    }
	  }
#endif

	for (int quartet = 0; quartet < numQuartets; ++quartet)
	  {
	    const int nodeA = quartets[quartet*4];
	    const int nodeB = quartets[quartet*4+1];
	    const int nodeC = quartets[quartet*4+2];
	    const int nodeD = quartets[quartet*4+3];
	    if(! (label==labeling[nodeA] && label==labeling[nodeB] && label==labeling[nodeC] && label==labeling[nodeD])){

	      int node_ids[4] = { nodeA, nodeB, nodeC, nodeD };
	      pbf.AddHigherTerm(4, node_ids, quartet_data[quartet].buffer);
	    }
	  }
	double newEnergy;
#ifdef HAS_QPBO
 	qpbo.Reset();
	reduce_and_convert(pbf, qpbo, reductionMode);
	qpbo.MergeParallelEdges();
	qpbo.Solve();

	for(int node = 0; node < numNodes; ++node) {
	  if(labeling[node] != label)
	    {
	      if(qpbo.GetLabel(node) < 0) unlabeledNodes++;
	    }
	}

	//TRY QPBO-I to improve the solution
	if(unlabeledNodes > 0)
	  {
	    srand ( static_cast<unsigned int>(time(NULL)) );

	    int numTrials = MAX_IMPROVEMENTS;
	    const double ratioThresh = 0.3;
	    ratioUnlabeled = static_cast<double>(unlabeledNodes) / static_cast<double>(numNodes);
	    if (ratioUnlabeled < ratioThresh) numTrials = static_cast<int>(0.5+ratioUnlabeled*numTrials/ratioThresh);

	    if(MAX_IMPROVEMENTS > 0)
	      {
		for(int i = 0; i < numTrials; ++i)
		  {
		    if(qpbo.Improve()) improveCounter++;
		  }
	      }

	    if(ratioUnlabeled > unlabeledMax) unlabeledMax = ratioUnlabeled;
	    unlabeledAvg += ratioUnlabeled;
	  }





	for(int node = 0; node < numNodes; ++node)
	  {
	    if(labeling[node] != label)
	      {
		if(qpbo.GetLabel(node) == 1) { labeling[node] = label; nodesChanged++; }
	      }
	  }

	if(verbose){//energy->applytestLabeling(labeling,label);
			newEnergy = energy->evaluateTotalCostSum();
}

#else
	FPDMODEL->reset();
	//	cout << "HOCR fastPD " << endl;
	reduce_and_convert(pbf,*FPDMODEL, reductionMode);
	//	exit(1);
	FPDMODEL->initialise();
	int *Labels=FPDMODEL->getLabeling();
	FPD::FastPD opt(FPDMODEL, 100 );
	newEnergy = opt.run();
        opt.getLabeling(Labels);

	for(int node = 0; node < numNodes; ++node)
	  {
	    //  cout << node << " labeling[node] " << labeling[node] << endl;
	    if(labeling[node] != label)
	      {
		//	cout << "Labels[node]  " << Labels[node]  << endl;
		if(Labels[node] == 1) { labeling[node] = label; nodesChanged++; }
	      }
	  }
#endif


	//

	if(verbose){

	  //
	  energy->report();
	  std::cout << "  LAB " << label << ":\t" << lastEnergy << " -> " << newEnergy << " / " << ratioUnlabeled * 100 << "% UNL / " << nodesChanged / static_cast<double>(numNodes) * 100 << "% CHN / IMP: " << improveCounter << std::endl;
	  lastEnergy = newEnergy;
      }


	//	exit(1);
      }

      }
    }

    unlabeledAvg /= static_cast<double>(numLabels*numSweeps);
#ifndef PRINT_ENERGY
//cout << " evaluateTotalCost " << endl;
    //energy->applytestLabeling(labeling,numLabels+1);
    lastEnergy = energy->evaluateTotalCostSum();
#endif

    return lastEnergy;

	}
}
