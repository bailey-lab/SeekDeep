//
// Created by Nicholas Hathaway on 3/5/23.
//


#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/objects/seqContainers.h>
#include <njhseq/seqToolsUtils/seqToolsUtils.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSplitter.hpp>
#include <njhseq/objects/collapseObjects/collapser.hpp>
#include <njhseq/helpers.h>
#include <njhseq/objects/seqObjects/Clusters/clusterUtils.hpp>
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"

namespace njhseq {








struct KmerClusteringRatePars {

  double cutOff = 0.05;
  uint32_t repCutOff = 1;
  double freqCutOff = 0.005;
  std::string sizeCutOffStr = "0.05%,1";
  double idCutOff = .90;
  bool byScore = false;
  bool map = false;
  uint32_t numThreads = 1;
  bool doTies = false;
  bool tiesDuringMapping = true;
//	comparison passableErrors;
  uint32_t kmerStart = 2;
  uint32_t kmerStop = 5;
  uint32_t largeIndel = 10;
  uint32_t chiAllowableError = 0;
  bool writeInitialClusters = false;
  uint32_t largeGapSizeCutOff = 50;
  bool breakLargeIndelCons = false;
  std::string aSetRepName;
  uint32_t nodeSize = 50;
  bool visualize = false;
  bool noTies = false;
  std::string additionalOutLocationFile;
  bool recalcConsensus = false;

  bool verbose{false};

  uint32_t readLengthMinDiff{0};


  readDistGraph<double>::dbscanPars dbPars_{};

  bool useHDBScan = false;
  readDistGraph<double>::HDBScanInputPars hdbScanPars_;

  cluster::snpBreakoutPars breakoutPars;

  uint32_t subSamplingAmount_ = 2000;

  uint32_t subSamplingMinAmount_ = 500;


  KmerClusteringRatePars() {


    dbPars_.minEpNeighbors_ = 4;
    dbPars_.eps_ = std::numeric_limits<double>::min();
//		passableErrors.lqMismatches_ = 0;
//		passableErrors.hqMismatches_ = 0;
//		passableErrors.oneBaseIndel_ = 5;
//		passableErrors.twoBaseIndel_ = 5;
  }
};

std::vector<std::shared_ptr<seqWithKmerInfo>> createKmerReadVec(
    const SeqIOOptions & opts) {
  std::vector<std::shared_ptr<seqWithKmerInfo>> ret;
  seqInfo seq;
  SeqInput reader(opts);
  reader.openIn();
  while (reader.readNextRead(seq)) {
    ret.emplace_back(std::make_shared<seqWithKmerInfo>(seq));
    readVec::handelLowerCaseBases(ret.back(), opts.lowerCaseBases_);
  }
  return ret;
}








readDistGraph<double> genKmerAccerDistGraphWithDbSmartBuildDev(
    std::vector<std::shared_ptr<seqWithKmerInfo>> & reads,
    const KmerClusteringRatePars & pars,
    aligner &alignerObj) {
  if (pars.verbose) {
    std::cout << njh::bashCT::bold << "Computing kmer distances and Building Graph"
              << njh::bashCT::reset << std::endl;
  }
  if(pars.kmerStop < 3 ){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ": kmerStop is less than 3, should at least be 3 or greater " <<
       ", KmerStop: " << pars.kmerStop;
    throw std::runtime_error{ss.str()};
  }

  std::function<
      double(const std::shared_ptr<seqWithKmerInfo> &,
             const std::shared_ptr<seqWithKmerInfo> &)> disFun =
      [](const std::shared_ptr<seqWithKmerInfo> & read1,
         const std::shared_ptr<seqWithKmerInfo> & read2) {
        auto dist = read1->compareKmers(*read2);
        return dist.second;
      };
  std::unordered_map<uint32_t, std::vector<std::vector<double>>>distanceMaps;
  for (uint32_t k = 2; k < 4; ++k) {
    if(pars.verbose){
      std::cout << "K: " << k << std::endl;
      std::cout << "\tIndexing Kmers" << std::endl;
    }
    njh::stopWatch watch;
    allSetKmers(reads, k, false);

    if(pars.verbose){
      std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
      std::cout << "\tCalculating Distances" << std::endl;
    }
    watch.reset();
    distanceMaps[k] = getDistance(reads, pars.numThreads, disFun);
    if(pars.verbose){
      std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
    }
  }
  std::vector<std::vector<double>> distances = distanceMaps[2];
  for (const auto rowPos : iter::range(distances.size())) {
    for (const auto colPos : iter::range(distances[rowPos].size())) {
      /*double res = std::abs(
          distanceMaps[2][rowPos][colPos]
              - distanceMaps[3][rowPos][colPos] );*/
      double res = distanceMaps[2][rowPos][colPos] - distanceMaps[3][rowPos][colPos];
      if(res <= pars.dbPars_.eps_){
        distances[rowPos][colPos] = res;
      }else{
        distances[rowPos][colPos] = std::numeric_limits<double>::max();
      }
    }
  }
  for(const auto k : iter::range<uint32_t>(4,pars.kmerStop + 1)){
    if(pars.verbose){
      std::cout << "K: " << k << std::endl;
      std::cout << "\tIndexing Kmers" << std::endl;
    }
    njh::stopWatch watch;
    allSetKmers(reads, k, false);

    if(pars.verbose){
      std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
      std::cout << "\tCalculating Distances" << std::endl;
    }
    watch.reset();
    if(k > 5){
      for (const auto rowPos : iter::range(distances.size())) {
        for (const auto colPos : iter::range(distances[rowPos].size())) {
          std::vector<double> differences;
          for (uint32_t currentKLen = 2; currentKLen < pars.kmerStop; ++currentKLen) {
            /*differences.emplace_back(
             std::abs(
             distanceMaps[currentKLen][rowPos][colPos]
             - distanceMaps[currentKLen + 1][rowPos][colPos]));*/
            differences.emplace_back(
                distanceMaps[currentKLen][rowPos][colPos]
                - distanceMaps[currentKLen + 1][rowPos][colPos]);
          }
          if (std::numeric_limits<double>::max() != distances[rowPos][colPos]) {
            double res = vectorMean(differences);
            if(res <= pars.dbPars_.eps_){
              distances[rowPos][colPos] = res;
            }else{
              distances[rowPos][colPos] = std::numeric_limits<double>::max();
            }
          }
        }
      }
    }
    std::vector<std::vector<double>> currentKDists;
    std::vector<std::pair<uint32_t, uint32_t>> indices;
    for(const auto pos : iter::range(reads.size())){
      currentKDists.emplace_back(pos);
      for(const auto secondPos : iter::range(pos)){
        if(std::numeric_limits<double>::max() != distances[pos][secondPos]){
          indices.emplace_back(pos, secondPos);
        }
      }
    }
    if(pars.numThreads < 2 || pars.numThreads >= reads.size()){
      paritialDis(reads, indices, currentKDists, disFun);
    }else{
      std::vector<std::thread> threads;
      auto step = static_cast<uint32_t>(std::round(static_cast<double>(indices.size())/static_cast<double>(pars.numThreads)));
      std::vector<std::vector<std::pair<uint32_t, uint32_t>>> indsSplit;
      for(const auto tNum : iter::range(pars.numThreads - 1)){
        std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + tNum * step,
                                                         indices.begin() + (tNum + 1)*step};
        indsSplit.emplace_back(temp);
      }
      std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + (pars.numThreads - 1) * step,
                                                       indices.end()};
      indsSplit.emplace_back(temp);
      for(const auto tNum : iter::range(pars.numThreads)){
        threads.emplace_back(paritialDis<std::shared_ptr<seqWithKmerInfo>,double>, std::cref(reads),
                                      indsSplit[tNum], std::ref(currentKDists), disFun);
      }
      for(auto & t : threads){
        t.join();
      }
    }
    distanceMaps[k] = currentKDists;
    if(pars.verbose){
      std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
    }
  }



  readDistGraph<double> distanceGraph(reads);
  for (const auto rowPos : iter::range(distances.size())) {
    for (const auto colPos : iter::range(distances[rowPos].size())) {
      std::vector<double> differences;
      for (uint32_t k = 2; k < pars.kmerStop; ++k) {
        /*differences.emplace_back(
         std::abs(
         distanceMaps[k][rowPos][colPos]
         - distanceMaps[k + 1][rowPos][colPos]));*/
        differences.emplace_back(
            distanceMaps[k][rowPos][colPos]
            - distanceMaps[k + 1][rowPos][colPos]);
      }
      if (std::numeric_limits<double>::max() != distances[rowPos][colPos]) {
        double meanDecrease = vectorMean(differences);
        if (meanDecrease <= pars.dbPars_.eps_) {
          distanceGraph.addEdge(reads[rowPos]->seqBase_.name_,
                                reads[colPos]->seqBase_.name_, meanDecrease);
        }
      }
    }
  }
  //distanceGraph.turnOffEdgesAbove(pars.eps_);
  //turn off connections if they have large indels
  if(pars.readLengthMinDiff > 0){
    if (pars.verbose) {
      std::cout << njh::bashCT::bold << "Breaking connections with differences in length more than " << pars.readLengthMinDiff
                << njh::bashCT::reset << std::endl;
    }
    for (auto & e : distanceGraph.edges_) {
      if (e->on_) {
        auto seq1 = e->nodeToNode_.begin()->second.lock()->value_;
        auto seq2 = distanceGraph.nodes_[distanceGraph.nameToNodePos_[e->nodeToNode_.begin()->first]]->value_;
        if(uAbsdiff(len(*seq1), len(*seq2)) > pars.readLengthMinDiff){
          e->on_ = false;
        }
      }
    }
  }
  if (pars.breakLargeIndelCons) {
    if (pars.verbose) {
      std::cout << njh::bashCT::bold << "Breaking connections with large indels"
                << njh::bashCT::reset << std::endl;
    }
    for (auto & e : distanceGraph.edges_) {
      if (e->on_) {
        auto seq1 = e->nodeToNode_.begin()->second.lock()->value_;
        auto seq2 =
            distanceGraph.nodes_[distanceGraph.nameToNodePos_[e->nodeToNode_.begin()->first]]->value_;
        alignerObj.alignCacheGlobal(seq1, seq2);
//        alignerObj.alignCacheGlobalDiag(seq1, seq2);
        alignerObj.profilePrimerAlignment(seq1, seq2);
        bool foundLargeIndel = false;
        for (const auto & g : alignerObj.comp_.distances_.alignmentGaps_) {
          if (g.second.size_ >= pars.largeIndel) {
            foundLargeIndel = true;
            break;
          }
        }
        if (foundLargeIndel) {
          e->on_ = false;
        }
      }
    }
  }

  if (pars.useHDBScan) {
    distanceGraph.runHDBScan(pars.hdbScanPars_, distances);
  } else {
    distanceGraph.resetBestAndVistEdges();
    distanceGraph.resetVisitedNodes();
    distanceGraph.allDetermineLowestBest(pars.doTies);
    distanceGraph.removeOffEdges();
    distanceGraph.dbscan(pars.dbPars_);
    distanceGraph.assignNoiseNodesAGroup();
  }
  return distanceGraph;
}




template<typename T>
cluster getConsensusWithAReCalc(const std::vector<T> & reads, aligner & alignerObj,
                                const std::string & name){
  if(reads.empty()){
    return {};
  }
  cluster mainCluster(getSeqBase(*reads.begin()));
  if(reads.size() > 1){
    std::vector<cluster> inClusters;
    for(const auto readPos : iter::range<uint64_t>(1,reads.size())){
      inClusters.emplace_back(cluster(getSeqBase(reads[readPos])));
    }
    for(const auto & clus : inClusters){
      mainCluster.addRead(clus);
    }
  }

//	mainCluster.calculateConsensus(alignerObj, true);
//	//once a consensus has been built, build again against the built consensus
//	mainCluster.needToCalculateConsensus_ = true;
  baseCluster::calculateConsensusPars conPars(true);
  conPars.convergeAttempts = 20;
  conPars.convergeConsensus = true;
  mainCluster.calculateConsensus(alignerObj, conPars);
  mainCluster.setName(name);
  return mainCluster;
}



struct SimpleCollapsePars {

  comparison passableErrors_;
  uint32_t klenComp_ = 9;
  double initialKdistCutOff_ = .95;
  double kdistCutOff_ = 0.80;
  double kdistCutOffDrop_ = 0.10;
  bool forceOneComp_{false};
  bool debug_{false};

  SimpleCollapsePars(){
    passableErrors_.lqMismatches_ = 5;
    passableErrors_.hqMismatches_ = 0;
    passableErrors_.oneBaseIndel_ = 5;
    passableErrors_.twoBaseIndel_ = 5;
    passableErrors_.largeBaseIndel_ = 0.99;

  }
};

bool simpleCollapseWithParsWork(std::vector<cluster> & consensusReads,
                                aligner & alignerObj,
                                const SimpleCollapsePars & collapsePars) {

  if(consensusReads.size() <=1){
    return false;
  }
  readVecSorter::sort(consensusReads);
  //compare final consensus reads to see if any created very similar consensus and should be collapsed into one
  std::set<uint64_t> needToRemove;
  comparison noErrors;
//	uint32_t klenComp = 9;
//	double kdistCutOff = .95;
  std::vector<kmerInfo> kinfos;
  for(const auto & seq : consensusReads){
    //std::cout << seq.seqBase_.name_ << std::endl;
    kinfos.emplace_back(seq.seqBase_.seq_, collapsePars.klenComp_ , false);
  }
  //first do a comparison of no errors to collapse almost identical clusters
  for (const auto firstPos : iter::range(consensusReads.size())) {

    if (consensusReads[firstPos].remove) {
      continue;
    }
    for (const auto secondPos : iter::range(firstPos + 1,
                                            consensusReads.size())) {
      bool print = collapsePars.debug_;
//			if("AS2-S0-Sub0-mip0MID7G8-5ng-rep1.fastq.07_t1308.3" == consensusReads[firstPos].seqBase_.name_ && "AS2-S0-Sub0-mip0MID7G8-5ng-rep1.fastq.08_t41.3333" == consensusReads[secondPos].seqBase_.name_){
//				print = true;
//			}
      if (consensusReads[secondPos].remove) {
        continue;
      }
      //first check for an exact sequence match
      if (consensusReads[secondPos].seqBase_.seq_
          == consensusReads[firstPos].seqBase_.seq_) {
        consensusReads[firstPos].addRead(consensusReads[secondPos]);
        needToRemove.emplace(secondPos);
        consensusReads[secondPos].remove = true;
      } else {
        if(print){
          std::cout << __FILE__ << " " << __LINE__ << std::endl;
          std::cout << kinfos[firstPos].compareKmers(kinfos[secondPos]).second << std::endl;
        }
        if(kinfos[firstPos].compareKmers(kinfos[secondPos]).second < collapsePars.initialKdistCutOff_){
          continue;
        }
        //check to see if there are no errors
        //(basically if there is a little bit of overlap, could be dangerous for very different lengths)
//        alignerObj.alignCacheGlobalDiag(consensusReads[secondPos], consensusReads[firstPos]);
        alignerObj.alignCacheGlobal(consensusReads[secondPos], consensusReads[firstPos]);

        alignerObj.profileAlignment(consensusReads[secondPos],
                                    consensusReads[firstPos], false, true, false);
        if (noErrors.passErrorProfile(alignerObj.comp_)) {
          consensusReads[firstPos].addRead(consensusReads[secondPos]);
          needToRemove.emplace(secondPos);
          consensusReads[secondPos].remove = true;
        }
      }
    }
  }

  uint32_t nonRemoved = 0;
  for (const auto & consensusRead : consensusReads) {
    if (!consensusRead.remove) {
      ++nonRemoved;
    }
  }
  if(nonRemoved > 1){
    //now do a second compare with some errors allowed
    //decreacse kDist cut off
    double kdistCutOff = collapsePars.kdistCutOff_;
    std::vector<uint32_t> positions(consensusReads.size(), 0);
    njh::iota<uint32_t>(positions, 0);
    for (const auto & firstPos : iter::reversed(positions)) {
      if (consensusReads[firstPos].remove) {
        continue;
      }
      uint32_t comps = 0;
      bool run = true;
      while(run){
        for (const auto secondPos : iter::range(firstPos)) {
          if (consensusReads[secondPos].remove) {
            continue;
          }
//					bool print = false;
          bool print = collapsePars.debug_;
//			if("AS2-S0-Sub0-mip0MID7G8-5ng-rep1.fastq.07_t1308.3" == consensusReads[secondPos].seqBase_.name_ && "AS2-S0-Sub0-mip0MID7G8-5ng-rep1.fastq.08_t41.3333" == consensusReads[firstPos].seqBase_.name_){
//				print = true;
//			}
          if(print){
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            std::cout << "" << consensusReads[firstPos].seqBase_.name_ << " vs " << consensusReads[secondPos].seqBase_.name_ << std::endl;
            std::cout << "kdistCutOff: " << kdistCutOff << std::endl;
            std::cout << kinfos[firstPos].compareKmers(kinfos[secondPos]).second << std::endl;
          }
          if(kinfos[firstPos].compareKmers(kinfos[secondPos]).second < kdistCutOff){
            continue;
          }
          if(print){
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
            std::cout << "" << consensusReads[firstPos].seqBase_.name_ << " vs " << consensusReads[secondPos].seqBase_.name_ << std::endl;
            std::cout << "passed cut off of: " << kdistCutOff  << std::endl;
          }
          alignerObj.alignCacheGlobal(consensusReads[secondPos], consensusReads[firstPos]);
//          alignerObj.alignCacheGlobalDiag(consensusReads[secondPos], consensusReads[firstPos]);
          alignerObj.profileAlignment(consensusReads[secondPos],
                                      consensusReads[firstPos], false, true, false);
          if(print){
            std::ofstream outTest("outtest.fastq");
            alignerObj.alignObjectA_.seqBase_.outPutFastq(outTest);
            alignerObj.alignObjectB_.seqBase_.outPutFastq(outTest);
            std::cout << "passableErrors.passErrorProfile(alignerObj.comp_): " << njh::colorBool(collapsePars.passableErrors_.passErrorProfile(alignerObj.comp_)) << std::endl;
            if(collapsePars.passableErrors_.passErrorProfile(alignerObj.comp_)){
              std::cout << alignerObj.comp_.toJson() << std::endl;
            }
          }
          if (collapsePars.passableErrors_.passErrorProfile(alignerObj.comp_)) {
            consensusReads[firstPos].addRead(consensusReads[secondPos]);
            needToRemove.emplace(secondPos);
            consensusReads[secondPos].remove = true;
          }
        }
        run = false;
        if(collapsePars.forceOneComp_ && comps ==0 && kdistCutOff > 0){
          run = true;
          kdistCutOff -= collapsePars.kdistCutOffDrop_;
        }
      }
    }
  }



  //remove the reads that need to be removed
  for (const auto & pos : iter::reversed(needToRemove)) {
    consensusReads.erase(consensusReads.begin() + pos);
  }
  //calculate consensus again
  bool changed = true;
  uint32_t calcCount = 0;
  while(changed){
    ++calcCount;
    VecStr oldSeqs = readVec::getSeqs(consensusReads);
    clusterVec::allCalculateConsensus(consensusReads, alignerObj, true);
    readVec::allUpdateName(consensusReads);
    changed = false;
    for(const auto seqPos : iter::range(consensusReads.size())){
      if(oldSeqs[seqPos] != consensusReads[seqPos].seqBase_.seq_){
        changed = true;
        break;
      }
    }
  }
  return calcCount > 1;
}

void simpleCollapseWithPars(std::vector<cluster> & consensusReads,
                            aligner & alignerObj,
                            const SimpleCollapsePars & collapsePars){
  uint32_t rounds = 0;
  while(simpleCollapseWithParsWork(consensusReads, alignerObj, collapsePars)){
    ++rounds;
  }
}




int SeekDeepRunner::kmerClusteringRate(const njh::progutils::CmdArgs & inputCommands) {
  KmerClusteringRatePars pars;
  SimpleCollapsePars postCollapsePars;
  bool sortToTopBeforeSubGrouping = false;
  bool sortByErrorRateBeforeSubGrouping = false;
  uint32_t klenForPreSort = 4;

  bool checkIndelsWhenMapping = false;
  bool checkIndelsAgainstSNPsWhenMapping = false;
  bool doNotCheckIndelsAgainstSNPsWhenMapping = false;
  seqSetUp setUp(inputCommands);
  setUp.processDebug();
  setUp.processVerbose();
  pars.verbose = setUp.pars_.verbose_;
  setUp.setOption(checkIndelsWhenMapping, "--checkIndelsWhenMapping", "check Indels When Mapping");
  //setUp.setOption(checkIndelsAgainstSNPsWhenMapping, "--checkIndelsAgainstSNPsWhenMapping", "check Indels Against SNPs When Mapping");
  setUp.setOption(doNotCheckIndelsAgainstSNPsWhenMapping, "--doNotCheckIndelsAgainstSNPsWhenMapping", "don't check Indels Against SNPs When Mapping");
  checkIndelsAgainstSNPsWhenMapping = !doNotCheckIndelsAgainstSNPsWhenMapping;
  setUp.pars_.chiOpts_.parentFreqs_ = 2;
  setUp.setOption(sortToTopBeforeSubGrouping, "--sortToTopBeforeSubGrouping", "sort To Top Before Sub Grouping");
  setUp.setOption(sortByErrorRateBeforeSubGrouping, "--sortByErrorRateBeforeSubGrouping", "sort By Error Rate Before Sub Grouping");
  setUp.setOption(klenForPreSort, "--klenForPreSort", "k-len For Pre Sort");

  pars.readLengthMinDiff = 400;
  setUp.setOption(pars.readLengthMinDiff, "--readLengthMinDiff", "The Read Length Mininum Difference to form a connection");

  setUp.setOption(pars.subSamplingAmount_, "--subSamplingAmount", "Amount to sub-sample for each group to run initial clustering on to save on memory and speed, the input be divied in groups of this size and clustering ran on each group");
  setUp.setOption(pars.subSamplingMinAmount_, "--subSamplingMinAmount", "Minimum amount in a sub grouping if multiple groups");
  if(pars.subSamplingMinAmount_ > pars.subSamplingAmount_){
    setUp.failed_ = true;
    setUp.addWarning(njh::pasteAsStr("--subSamplingMinAmount:", pars.subSamplingMinAmount_, " can't be greater than --subSamplingAmount:", pars.subSamplingAmount_));
  }

  setUp.setOption(pars.visualize, "--visualize", "Visualize the initial clustering");
  setUp.setOption(pars.aSetRepName, "--repName", "Set the repname to this for all reads regardless of naming");
  setUp.setOption(pars.largeIndel, "--largeIndelSize", "Size for connections breaking to be considered large");
  setUp.setOption(pars.breakLargeIndelCons, "--breakLargeIndelCons", "Break Connections between nodes that have large indels(>" + estd::to_string(pars.largeIndel) + "bp)");
  setUp.setOption(pars.kmerStart, "--kmerLenStart", "Length for kmers to start at");
  setUp.setOption(pars.kmerStop, "--kmerLenStop", "Length for kmers to stop at (inclusive)");
  setUp.processComparison(postCollapsePars.passableErrors_);
  setUp.setOption(postCollapsePars.forceOneComp_, "--postCollapseForceOneComp", "Post Collapse Force One Comp");


//  setUp.setOption(pars.noTies, "--noTiesDuringMapping", "Don't do ties on mapping estimate");
//  pars.tiesDuringMapping = !pars.noTies;
  pars.tiesDuringMapping = false;
  setUp.setOption(pars.tiesDuringMapping, "--tiesDuringMapping", "distribute ties on mapping estimate");
  setUp.setOption(pars.doTies, "--doTies", "Do ties");
  setUp.setOption(pars.largeGapSizeCutOff, "--largeGapSizeCutOff", "During Mapping Exclude any reads have a gap larger than this value");
  setUp.setOption(pars.writeInitialClusters, "--writeInitialClusters", "Write the initial clusters before further collapsing");
  setUp.setOption(pars.chiAllowableError, "--chiAllowableError", "Allowable Error For Chi Overlap");
  setUp.setOption(setUp.pars_.chiOpts_.checkChimeras_, "--checkChimeras", "Check For Possible chimeric sequence");
  setUp.setOption(setUp.pars_.chiOpts_.parentFreqs_, "--parFreqs", "Chimeric parent frequency multiplier");
  setUp.setOption(pars.additionalOutLocationFile, "--additionalOut",
                  "A location file created by makeSampleDirectories");
  setUp.setOption(pars.numThreads, "--numThreads,-t", "Number of Threads to Use");
  setUp.setOption(pars.map, "--map", "Map the Input Reads Against Consensus Seqs");
  setUp.setOption(pars.byScore, "--byScore", "By Alignment Score");
  setUp.setOption(pars.repCutOff, "--repCutOff", "PCR Rep Count Cut Off");
  setUp.setOption(pars.sizeCutOffStr, "--sizeCutOff",
                  "Fraction of Reads for Cluster Size Cut Off, can give XX% for relative number or XX for absolute number");
  setUp.setOption(pars.freqCutOff, "--freqCutOff", "Mapped Frequency Cut off");

  setUp.setOption(pars.cutOff, "--cutOff", "Average decrease in kmer similarity score");
  pars.dbPars_.eps_ = pars.cutOff;
  setUp.setOption(pars.dbPars_.minEpNeighbors_, "--minNeighbors", "Number of minimum neighbors for dbscan");


//  bool doNotUseHDBScan = false;
//  setUp.setOption(doNotUseHDBScan, "--doNotUseHDBScan", "Don't use Hierarchical density based scan for initial cluster formation, just use classical DBScan");
//  pars.useHDBScan = !doNotUseHDBScan;
  setUp.setOption(pars.useHDBScan, "--useHDBScan", "use Hierarchical density based scan for initial cluster formation");

  if(pars.useHDBScan){
    pars.dbPars_.eps_ = std::numeric_limits<double>::max();
  }

  pars.hdbScanPars_.verbose = setUp.pars_.verbose_;
  pars.hdbScanPars_.debug = setUp.pars_.debug_;
  setUp.setOption(pars.hdbScanPars_.HDBSredetermineMaxEps, "--HDBSredetermineMaxEps", "HDBS redetermine Max Eps allowed in initial step based by setting it equal to mean of the non-same-group dist minus 2sd");
  bool doNotCountZeroNeighbors = false;
  setUp.setOption(doNotCountZeroNeighbors, "--HDBSdoNotCountZeroNeighbors", "HDBS when doing j-th nearest neighbor do Not Count Zero Neighbors");
  pars.hdbScanPars_.HDBScountZeroNeighbors = !doNotCountZeroNeighbors;
  setUp.setOption(pars.hdbScanPars_.HDBSmaxInitialEps, "--HDBSmaxInitialEps", "HDBS a hard cut off for max Initial Eps for initial DBSCAN step in H-DBSCAN");
  setUp.setOption(pars.hdbScanPars_.proposedClusters, "--HDBSproposedClusters", "HDBS proposed number of clusters");
  setUp.setOption(pars.hdbScanPars_.groupDiffToReCalc, "--HDBSgroupDiffToReCalc", "HDBS the difference in group counts to re-calculate centroid distances");

  setUp.setOption(pars.hdbScanPars_.HDBSCountSingletGroups, "--HDBSCountSingletGroups", "For HD DBscan count Singlet Groups, by default these are not included in towards the proposed group counts");
  pars.hdbScanPars_.numThreads = pars.numThreads;

  setUp.setOption(pars.idCutOff, "--idCutOff",
                  "Identity Cut Off to be Mapped to reads");
  setUp.pars_.ioOptions_.lowerCaseBases_ = "remove";
  setUp.pars_.ioOptions_.out_.outFilename_ = "output";
  setUp.processDefaultReader(true);
  setUp.setOption(setUp.pars_.colOpts_.kmerOpts_.kLength_, "--kmerLength,-k", "kLength");
  setUp.processDirectoryOutputName(true);
//  setUp.pars_.gapLeft_ = "5,1";
//  setUp.pars_.gapRight_ = "5,1";
//  setUp.pars_.gap_ = "5,1";
//  setUp.pars_.gapInfo_ = gapScoringParameters(5,1);
//
  setUp.pars_.gapLeft_ = "0,0";
  setUp.pars_.gapRight_ = "0,0";
  setUp.pars_.gap_ = "5,1";
  setUp.pars_.gapInfo_ = gapScoringParameters(5,1, 0,0, 0,0);

  setUp.processGap();

//  setUp.pars_.gapLeftRef_ = "5,1";
//  setUp.pars_.gapRightRef_ = "5,1";
//  setUp.pars_.gapRef_ = "5,1";
//  setUp.pars_.gapInfoRef_ = gapScoringParameters(5,1);

  setUp.pars_.gapLeftRef_ = "0,0";
  setUp.pars_.gapRightRef_ = "0,0";
  setUp.pars_.gapRef_ = "5,1";
  setUp.pars_.gapInfoRef_ = gapScoringParameters(5,1, 0,0, 0,0);
  setUp.processRefFilename(false);


  setUp.pars_.qScorePars_.qualThresWindow_ = 2;
  setUp.processQualThres();
  setUp.processScoringPars();
  setUp.processAlnInfoInput();
  setUp.setOption(pars.recalcConsensus, "--recalcConsensus", "Recalculate consensus after mapping reads to the consensus sequences");
  pars.breakoutPars.minSnps = postCollapsePars.passableErrors_.hqMismatches_ + 1;
  pars.breakoutPars.hardCutOff = pars.dbPars_.minEpNeighbors_;
  pars.breakoutPars.qScorePars = setUp.pars_.qScorePars_;
  //pars.breakoutPars.qScorePars.qualThresWindow_ = 0;
  setUp.setOption(pars.breakoutPars.qScorePars.qualThresWindow_, "--snpBreakoutParsQualThresWindow", "snpBreakout Pars Qual Thres Window");
  setUp.setOption(pars.breakoutPars.hardCutOff, "--snpBreakoutMinGroupSize", "A hard cut off for when breaking out clusters, clusters that larger(non-inclusive) than this are broken out");
  pars.breakoutPars.snpFreqCutOff = 0.025;
  pars.breakoutPars.hardSnpFreqCutOff = 0.025;
  setUp.setOption(pars.breakoutPars.snpFreqCutOff, "--snpFreqCutOff", "Cut off for when breaking out snp frequencies");
  setUp.setOption(pars.breakoutPars.hardSnpFreqCutOff, "--hardSnpFreqCutOff", "Hard SNP Freq Cut Off");
  setUp.setOption(pars.breakoutPars.minSnps, "--snpBreakoutMinSnps", "SNP Breakout Min Snps");

  double kDistCutOff = .90;
  uint32_t kLenComp = 6;
  uint32_t  threadChunkNum = 10;
  bfs::path colorFile  = "";
  setUp.setOption(colorFile, "--colorFile", "Name of a color file");
  setUp.setOption(threadChunkNum, "--threadChunkNum", "Thread Chunk Num");
  setUp.setOption(kDistCutOff, "--kDistCutOff", "For speeding up mapping use this k cut off to skip comparisons");
  setUp.setOption(kLenComp, "--kLenComp", "For speeding up mapping use this k len to calc distances to skip comparisons");
  setUp.finishSetUp(std::cout);

  setUp.startARunLog(setUp.pars_.directoryName_);
  setUp.rLog_.setCurrentLapName("Start");
  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.tab.txt",
                            false, false);

  if(setUp.pars_.verbose_){
    std::cout << "Gap Scoring Ref: " << std::endl;
    std::cout << setUp.pars_.gapInfoRef_.toJson() << std::endl;

    std::cout << "Gap Scoring: " << std::endl;
    std::cout << setUp.pars_.gapInfo_.toJson() << std::endl;

  }

  //read in reads and get the cluster size cut off
  //create the sequence with kmer info vector
  std::map<uint32_t, std::vector<std::shared_ptr<seqWithKmerInfo>>> inputReadsGrouped;
  std::vector<std::shared_ptr<seqWithKmerInfo>> allInputReads;
  uint32_t totalReadCnt = 0;
  uint64_t maxReadLength = 0;
  if(setUp.pars_.verbose_){
    std::cout << njh::bashCT::bold << "Reading" << std::endl;;

  }
  {
    seqInfo seq;
    SeqInput reader(setUp.pars_.ioOptions_);
    reader.openIn();
    while (reader.readNextRead(seq)) {
      readVec::getMaxLength(seq, maxReadLength);
      auto kSeq = std::make_shared<seqWithKmerInfo>(seq);
      readVec::handelLowerCaseBases(kSeq, setUp.pars_.ioOptions_.lowerCaseBases_);
      allInputReads.emplace_back(kSeq);
    }
    if(sortByErrorRateBeforeSubGrouping){
      readVecSorter::sortByTotalCountAE(allInputReads, true);
    }
    if(sortToTopBeforeSubGrouping){
      allSetKmers(allInputReads, klenForPreSort, false);
      seqWithKmerInfo top = *allInputReads.front();
      readVecSorter::sortReadVectorFunc(allInputReads,[&top](const std::shared_ptr<seqWithKmerInfo>& p1, const std::shared_ptr<seqWithKmerInfo> & p2){
        return top.compareKmers(*p1).second > top.compareKmers(*p2).second;
      });
    }

    for(const auto & inputSeq : allInputReads){
      inputReadsGrouped[totalReadCnt/pars.subSamplingAmount_].emplace_back(inputSeq);
      ++totalReadCnt;
    }
  }
  if(setUp.pars_.verbose_){
    std::cout << njh::bashCT::bold << "Done Reading" << std::endl;;
  }
  uint32_t clustersizeCutOff = processRunCutoff(pars.sizeCutOffStr, totalReadCnt);
  if(setUp.pars_.verbose_){
    std::cout << njh::bashCT::bold << "Read in " << totalReadCnt << " sequences and grouped into " << inputReadsGrouped.size() << " sub groups of size: " << pars.subSamplingAmount_ << std::endl;;
    std::cout << "Cluster Cut off " << clustersizeCutOff << njh::bashCT::reset
              << std::endl;
  }
  //add some meta data about file and analysis paths so latter trace back from finalClustering population can happen
  Json::Value metaData;
  auto analysisDirPath = njh::files::bfs::canonical(setUp.pars_.directoryName_);
  metaData["analysisDirPath"] = njh::json::toJson(analysisDirPath.string());
  auto fullPathToInput = njh::files::bfs::canonical(setUp.pars_.ioOptions_.firstName_);
  metaData["inputFile"] = njh::json::toJson(fullPathToInput);
  metaData["extractionDir"] = njh::json::toJson(fullPathToInput.parent_path());

  //get max length of the reads for the aligner
  if(maxReadLength < 500){
    maxReadLength = 700;
  }
  //get ref seqs
  std::vector<readObject> refSeqs;
  if (!setUp.pars_.refIoOptions_.firstName_.empty()) {
    refSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxReadLength);
  }

  if(setUp.pars_.verbose_){
    std::cout << "Creating Aligner" << std::endl;
  }
  //create aligner
  aligner alignerObj(maxReadLength,
                     setUp.pars_.gapInfo_,
                     setUp.pars_.scoring_,
                     KmerMaps(), setUp.pars_.qScorePars_,
                     setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                     setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
  if(setUp.pars_.verbose_){
    std::cout << "Reading In Alignments" << std::endl;
  }
  alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);



  setUp.rLog_.logCurrentTime("All by All Comparison and Bulding Graph");


  std::vector<std::vector<readObject>> outReads;
  std::vector<std::shared_ptr<readDistGraph<double>>> disGraphs;

  for(auto & seqs : inputReadsGrouped){
    if(inputReadsGrouped.size() > 1 && (seqs.second.size() <= 1 || seqs.second.size() < pars.subSamplingMinAmount_) ){
      continue;
    }
    if(setUp.pars_.verbose_){
      std::cout << "On group: " << seqs.first + 1 << "/" << inputReadsGrouped.size() << " with " << seqs.second.size() << " sequences" << std::endl;
    }
    auto distanceGraph = std::make_shared<readDistGraph<double>>(genKmerAccerDistGraphWithDbSmartBuildDev(seqs.second,pars,
                                                                                                          alignerObj));

    if(setUp.pars_.debug_){
      distanceGraph->printAdjByGroup(std::cout);
    }
    setUp.rLog_.logCurrentTime(njh::pasteAsStr("Traveling Through Graph - ", seqs.first));
    if(setUp.pars_.verbose_){
      std::cout << "Number of groups: " << distanceGraph->numberOfGroups_ << std::endl;
    }
    //create a vector of vectors for each group in the graph
    std::vector<std::vector<readObject>> currentOutReads(distanceGraph->numberOfGroups_);
    if(setUp.pars_.debug_){
      std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
      auto groupCounts = distanceGraph->getGroupCounts();
      std::cout << "group\t" << "groupCount" << std::endl;
      for(const auto & group : groupCounts){
        std::cout << group.first << "\t" << group.second << std::endl;
      }
    }
    for (const auto & n : distanceGraph->nodes_) {
      if(setUp.pars_.debug_){
        std::cout << "n->group_: " << n->group_ << std::endl;
      }
      currentOutReads[n->group_].emplace_back(*(n->value_));
    }
    if(setUp.pars_.debug_){
      std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    }
    addOtherVec(outReads, currentOutReads);
    disGraphs.emplace_back(distanceGraph);
    if(setUp.pars_.debug_){
      std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    }
  }

  if(setUp.pars_.debug_){
    std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
  }
  //sort the groups by the number of groups in each
  njh::sort(outReads,
            [](const std::vector<readObject> & vec1,
               const std::vector<readObject> & vec2) {return vec1.size() > vec2.size();});
  //create various output directories
  auto clusDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("clusters"));
  auto singlesDir = njh::files::makeDir(clusDir, njh::files::MkdirPar("belowSizeCutOff"));
  std::string repDir;
  if(pars.repCutOff > 1){
    repDir = njh::files::makeDir(clusDir, njh::files::MkdirPar("belowRepCutOff")).string();
  }
  if(setUp.pars_.debug_){
    std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
  }
  //create consensus sequence and filter clusters on replicate count and cluster size cut off
  std::vector<cluster> consensusReads;
  std::ofstream initialClusterInfoFile;
  openTextFile(initialClusterInfoFile,
               setUp.pars_.directoryName_ + "initialClusterInfo.tab.txt", ".tab.txt",
               false, false);
  initialClusterInfoFile << "kmerClusterId\treadCount\n";
  setUp.rLog_.logCurrentTime("Building Consensus Sequences");

  uint32_t totalClusterCount = 0;

  for (const auto vecPos : iter::range(outReads.size())) {
    if (outReads[vecPos].size() > clustersizeCutOff) {
      ++totalClusterCount;
    }else{
      SeqOutput::write(outReads[vecPos],
                       SeqIOOptions(singlesDir.string() + leftPadNumStr<uint32_t>(vecPos, outReads.size()),
                                    setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
    }
  }
  if(setUp.pars_.verbose_){
    std::cout << njh::bashCT::bold << "Building Consensus Sequences for " << totalClusterCount << " groups" << njh::bashCT::reset << std::endl;
  }
  for (const auto vecPos : iter::range(outReads.size())) {
    initialClusterInfoFile << vecPos << '\t'
                           << getPercentageString(readVec::getTotalReadCount(outReads[vecPos]),
                                                  totalReadCnt) << '\n';
  }
  {
    std::vector<uint32_t> readPositions(totalClusterCount);
    njh::iota<uint32_t>(readPositions, 0);
    concurrent::AlignerPool alignPool(alignerObj, pars.numThreads);
    alignPool.initAligners();
    njh::concurrent::LockableQueue<uint32_t> readPositionsQueue(readPositions);
    std::mutex consensusMut;
    njh::ProgressBar pbar(totalClusterCount);
    pbar.progColors_ = pbar.RdYlGn_;

    std::function<void()> buildConsensus = [&readPositionsQueue,&consensusMut,
        &consensusReads,&outReads,&alignerObj,
        &setUp,&pars,&clustersizeCutOff,
        &alignPool,&pbar,&repDir](){
      uint32_t vecPos = std::numeric_limits<uint32_t>::max();

      auto curThreadAligner = alignPool.popAligner();
      while(readPositionsQueue.getVal(vecPos)){
        //if below cluster cut off size throw out
        if (outReads[vecPos].size() > clustersizeCutOff) {
          //std::cout << "setUp.pars_.verbose_: " << njh::colorBool(setUp.pars_.verbose_) << std::endl;
          //std::cout << consensusReads.size() << consensusReads.size() << std::endl;

          cluster consensus;
          //if input is combined replicates do a replicate count,
          //only works if input reads named with rep name
          if (pars.repCutOff > 1) {
            std::set<std::string> repNames;
            for (const auto & read : outReads[vecPos]) {
              auto toks = tokenizeString(read.seqBase_.name_, ".");
              repNames.emplace(njh::replaceString(toks[0], "_Comp", ""));
            }
            if (repNames.size() < pars.repCutOff) {
              SeqOutput::write(outReads[vecPos],
                               SeqIOOptions(repDir + leftPadNumStr<uint32_t>(vecPos, outReads.size()),
                                            setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
            } else {
              //add consensus to consensus sequences
              readVecSorter::sortByTotalCountAE<readObject>(outReads[vecPos], true);
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
              consensus = getConsensusWithAReCalc(outReads[vecPos], *curThreadAligner,
                                                  leftPadNumStr<uint32_t>(vecPos, outReads.size()));
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
            }
          } else {
            //add consensus to consensus sequences
            readVecSorter::sortByTotalCountAE<readObject>(outReads[vecPos],true);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
            consensus = getConsensusWithAReCalc(outReads[vecPos], *curThreadAligner,
                                                leftPadNumStr<uint32_t>(vecPos, outReads.size()));
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
          }
          if(setUp.pars_.verbose_){
            pbar.outputProgAdd(std::cout, 1, true);
          }
          {
            std::lock_guard<std::mutex> conLock(consensusMut);
            consensusReads.emplace_back(consensus);
          }
        }
      }
      {
        std::lock_guard<std::mutex> conLock(consensusMut);
        alignerObj.alnHolder_.mergeOtherHolder(curThreadAligner->alnHolder_);
      }

    };
    //njh::concurrent::runVoidFunctionThreaded(buildConsensus, pars.numThreads);
    njh::concurrent::runVoidFunctionThreaded(buildConsensus, 1);


    readVecSorter::sort(consensusReads);
  }


  if(setUp.pars_.verbose_){
    std::cout << std::endl;
  }

  if(pars.writeInitialClusters){
    //write out initial clusters if needed, this is mostly a debugging option
    std::string initialDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("initialClusters")).string();
    std::string initialClustersDir = njh::files::makeDir(initialDir, njh::files::MkdirPar("clusters")).string();
    SeqOutput::write(consensusReads, SeqIOOptions(initialDir + "allConsensus",setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
    for(const auto & clus : consensusReads){
      clus.writeClustersInDir(initialClustersDir, setUp.pars_.ioOptions_);
    }
  }
  setUp.rLog_.logCurrentTime("Consensus Comparison");
  if(setUp.pars_.verbose_){
    std::cout << njh::bashCT::bold << "Comparing consensus sequences" << njh::bashCT::reset << std::endl;
  }
  //collapse clusters that ended up being the same exact sequence, or within a passable error profile
  simpleCollapseWithPars(consensusReads, alignerObj, postCollapsePars);
  if(setUp.pars_.verbose_){
    std::cout << njh::bashCT::bold << "Finished with " << consensusReads.size()
              << " Consensus Sequences" << njh::bashCT::reset << std::endl;
  }
  setUp.rLog_ << "Finished with " << consensusReads.size()
              << " Consensus Sequences"  << "\n";


  alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);


  if(pars.map){
    setUp.rLog_.logCurrentTime("Calculating Input Kmer info for comparison");
    allSetKmers(allInputReads, kLenComp, false);
  }


  std::vector<refMapContainer<seqInfo>> readCounts;
  std::string mappingDir;
  double unmappableAmount = 0.0;
  double indeterminateAmount = 0.0;
  double tiesAmount = 0.0;


  if (pars.map) {
    setUp.rLog_.logCurrentTime("Determining Variation");
    if (setUp.pars_.verbose_) {
      std::cout << njh::bashCT::bold << "Mapping Sequences" << njh::bashCT::reset << std::endl;
    }
    if (setUp.pars_.verbose_) {
      std::cout << njh::bashCT::bold << "Comparing Consensus Sequences" << njh::bashCT::reset << std::endl;
    }
    //create aligner with gap info from ref so a different gap scoring can be used for mapping determination
    aligner alignerObjMapper(maxReadLength, setUp.pars_.gapInfoRef_,
                             setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                             setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                             setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
    alignerObjMapper.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
    //set all reads' fraction
    njh::for_each(allInputReads, [&totalReadCnt](auto& seq) { getSeqBase(seq).setFractionByCount(totalReadCnt); });
    //create ref map containers to keep counts
    readCounts = refMapContainer<
        seqInfo>::createContainers<seqInfo, cluster>(consensusReads);
    //determine variation in final consensus sequences
    std::vector<refVariants> refVariationInfo;
    for (const auto refPos : iter::range(readCounts.size())) {
      refVariationInfo.emplace_back(readCounts[refPos].seqBase_);
    }
    for (const auto refPos : iter::range(readCounts.size())) {
      if (setUp.pars_.verbose_) {
        std::cout << njh::bashCT::bold << refPos + 1 << "/" <<  readCounts.size() << njh::bashCT::reset << std::endl;
      }
      for (const auto refSubPos : iter::range(readCounts.size())) {
        if (refPos == refSubPos) {
          continue;
        }
        if(uAbsdiff(len(readCounts[refSubPos].seqBase_), len(readCounts[refPos].seqBase_)) > pars.readLengthMinDiff){
          continue;
        }
        refVariationInfo[refPos].addVariant(readCounts[refSubPos].seqBase_,
                                            alignerObjMapper, false);
      }
    }
    std::vector<kmerInfo> conKInfos(consensusReads.size());

    setUp.rLog_.logCurrentTime("Calculating Consensus Kmer info for comparison");
    {
      std::vector<uint32_t> conReadPositions(consensusReads.size());
      njh::iota<uint32_t>(conReadPositions, 0);
      njh::concurrent::LockableQueue<uint32_t> conReadPosQueue(
          conReadPositions);
      auto fillInfos =
          [&consensusReads,&conKInfos,&conReadPosQueue,&kLenComp]() {
            uint32_t pos = std::numeric_limits<uint32_t>::max();
            while(conReadPosQueue.getVal(pos)) {
              conKInfos[pos] = kmerInfo(consensusReads[pos].seqBase_.seq_, kLenComp, false);
            }
          };
      std::vector<std::thread> threads;
      for (uint32_t tNum = 0; tNum < pars.numThreads; ++tNum) {
        threads.emplace_back(std::thread(fillInfos));
      }
      for (auto & t : threads) {
        t.join();
      }
    }


    setUp.rLog_.logCurrentTime("Mapping Sequences");
    if (setUp.pars_.verbose_) {
      std::cout << njh::bashCT::bold << "Mapping all sequences to Consensus Sequences" << njh::bashCT::reset << std::endl;
    }
    std::vector<std::pair<uint64_t, std::vector<uint64_t>>> ties;
    std::vector<seqInfo> tiesReads;
    std::mutex tiesMut;
    std::vector<seqInfo> unmappable;
    std::mutex unmappableMut;
    std::vector<seqInfo> indeterminate;
    std::mutex indeterminateMut;

    njh::ProgressBar pbar(allInputReads.size());
    pbar.progColors_ = pbar.RdYlGn_;
    //set the alignment best difference and function
    //the best alignment is considered within a certain range due such high error rates the best alignment is not
    //Necessary the one with the best score, it is the one that contains the segregating differences
    double difference = 0;
    std::function<double(const aligner &)> getScoreFunc;
    if (pars.byScore) {
      getScoreFunc = [](const aligner & alignerObject){
        return alignerObject.parts_.score_;
      };
      difference = 10;
    } else {
      getScoreFunc = [](const aligner & alignerObject){
        //return alignerObjMapper.comp_.distances_.eventBasedIdentity_;
        return alignerObject.comp_.distances_.eventBasedIdentityHq_;
      };
      difference = 0.001;
    }

    concurrent::AlignerPool alignPoolMap(alignerObjMapper, pars.numThreads);
    alignPoolMap.initAligners();
    alignPoolMap.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

    std::vector<uint32_t> readPositions(allInputReads.size());
    njh::iota<uint32_t>(readPositions, 0);
    njh::concurrent::LockableQueue<uint32_t> readPosQueue(readPositions);

    auto mapRead =
        [&pbar, &readPosQueue, &alignPoolMap, &pars, &getScoreFunc, &readCounts, &allInputReads,
            &tiesReads, &ties, &tiesMut, &unmappable, &unmappableMut, &indeterminate, &indeterminateMut,
            &difference, &refVariationInfo,&threadChunkNum, &setUp,
            &checkIndelsWhenMapping,&checkIndelsAgainstSNPsWhenMapping]() {
          std::vector<uint32_t> readPositionsChunk;
          auto curThreadAligner = alignPoolMap.popAligner();
          while(readPosQueue.getVals(readPositionsChunk, threadChunkNum)) {
            for(const auto readPos : readPositionsChunk){
              //best score to be used to get all the references that map within a certain distance to this latter
              double bestScore = 0;
              bool mappable = false;
              for (const auto refPos : iter::range(readCounts.size())) {

                /*if(conKInfos[refPos].compareKmers(reads[readPos]->kInfo_).second < kDistCutOff){
                  continue;
                }*/

                const auto & ref = readCounts[refPos];
                if(uAbsdiff(len(getSeqBase(ref)), len(allInputReads[readPos]->seqBase_)) > pars.readLengthMinDiff){
                  continue;
                }
//                curThreadAligner->alignCacheGlobalDiag(ref.seqBase_, allInputReads[readPos]->seqBase_);
                curThreadAligner->alignCacheGlobal(ref.seqBase_, allInputReads[readPos]->seqBase_);
                curThreadAligner->profilePrimerAlignment(ref.seqBase_,
                                                         allInputReads[readPos]->seqBase_);

                //check for very large gaps that wouldn't be expected no matter the technology
                bool passLargeGaps = true;
                for(const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                  if(g.second.size_ > pars.largeGapSizeCutOff) {
                    passLargeGaps = false;
                    break;
                  }
                }

                //check to see if the read is mappable by comparing to the percent identity cut off and has no large gaps
//                std::cout << std::this_thread::get_id() << std::endl;
//						std::cout << "curThreadAligner->comp_.distances_.eventBasedIdentity_: " << curThreadAligner->comp_.distances_.eventBasedIdentity_ << std::endl;
//						std::cout << "passLargeGaps: " << njh::colorBool(passLargeGaps) << std::endl;
                if (curThreadAligner->comp_.distances_.eventBasedIdentity_ >= pars.idCutOff && passLargeGaps) {
                  mappable = true;
                  double currentScore = getScoreFunc(*curThreadAligner);
                  if(currentScore > bestScore) {
                    bestScore = currentScore;
                  }
                }
              }
              // if mappable get the best references
              if (mappable) {
                std::vector<uint64_t> bestRefs;
                for (const auto refPos : iter::range(readCounts.size())) {
                  const auto & ref = readCounts[refPos];
                  curThreadAligner->alignCacheGlobal(ref.seqBase_,allInputReads[readPos]->seqBase_);
//                  curThreadAligner->alignCacheGlobalDiag(ref.seqBase_,allInputReads[readPos]->seqBase_);
                  curThreadAligner->profilePrimerAlignment(ref.seqBase_,
                                                           allInputReads[readPos]->seqBase_);
                  if (curThreadAligner->comp_.distances_.eventBasedIdentity_ >= pars.idCutOff) {
                    mappable = true;
                    double currentScore = getScoreFunc(*curThreadAligner);
                    bool passLargeGaps = true;
                    for(const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                      if(g.second.size_ > pars.largeGapSizeCutOff) {
                        passLargeGaps = false;
                        break;
                      }
                    }
                    //if the current score is a certain amount within the best score and doesn't have large gaps then put it
                    //into the vector for possible match
                    if (std::abs(currentScore - bestScore) < difference && passLargeGaps) {
                      bestRefs.emplace_back(refPos);
                    }
                  }
                }
//						std::cout << "bestRefs.size(): " << bestRefs.size() << std::endl;
                if (bestRefs.size() == 1) {
                  //if only matching one, then place with that one
                  readCounts[bestRefs.front()].addRead(allInputReads[readPos]->seqBase_);
                } else if (bestRefs.size() == 0) {
                  //if no mappable candidates found, put into un-mappable
                  std::lock_guard<std::mutex> unmLock(unmappableMut);
                  unmappable.emplace_back(allInputReads[readPos]->seqBase_);
                } else {
                  VecStr matchingRefNames;
                  for(const auto & refPos : bestRefs) {
                    matchingRefNames.emplace_back(readCounts[refPos].seqBase_.name_);
                  }
                  std::vector<uint64_t> matchingRefs;
                  for (const auto & refPos : bestRefs) {
                    const auto & ref = readCounts[refPos];
                    //get the current snp info for the current ref to the other matching refs
                    //auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 1);
                    auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 0);
                    curThreadAligner->alignCacheGlobal(ref.seqBase_, allInputReads[readPos]->seqBase_);
//                    curThreadAligner->alignCacheGlobalDiag(ref.seqBase_, allInputReads[readPos]->seqBase_);
                    curThreadAligner->profilePrimerAlignment(ref.seqBase_, allInputReads[readPos]->seqBase_);
                    //determine the snps for this current read to this ref
                    std::unordered_map<uint32_t, char> currentSnps;
                    for (const auto & m : curThreadAligner->comp_.distances_.mismatches_) {
                      currentSnps[m.second.refBasePos] = m.second.seqBase;
                    }
                    bool pass = true;
                    //iterate over segregating snps locations
                    for(const auto & loc : seqSnpPosBases) {
                      //determine if current read has a snp at location
                      auto search = currentSnps.find(loc.first);
                      if(search != currentSnps.end()) {
                        //if it does have a snp here, check to see if it is a known variant, if it is the read doesn't pass
                        if(njh::in(search->second, loc.second)) {
                          pass = false;
                          break;
                        }
                      }
                    }

                    if(pass && checkIndelsWhenMapping){
                      //check indels
                      //get the current insertions info for the current ref to the other matching refs
                      auto deletionsInOtherMatching = refVariationInfo[refPos].getVariantDeletionLociMap(matchingRefNames, 2, true);
                      auto insertionsInOtherMatching = refVariationInfo[refPos].getVariantInsertionLociMap(matchingRefNames, 2, true);

                      //determine the insertions for this current read to this ref
                      std::unordered_map<uint32_t, gap> currentInsertions;
                      std::unordered_map<uint32_t, gap> currentDeletions;
                      for (const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                        if(g.second.ref_){
                          currentInsertions.emplace(g.second.refPos_, g.second);
                        }else{
                          currentDeletions.emplace(g.second.refPos_, g.second);
                        }
                      }

                      //check deletions
                      //iterate over segregating deletions locations
                      for(const auto & loc : deletionsInOtherMatching) {
                        //determine if current read has a snp at location
                        auto search = currentDeletions.find(loc.first);
                        if(search != currentDeletions.end()) {
                          //if it does have a deletion here, check to see if it is a known variant, if it is the read doesn't pass
                          if(njh::contains(loc.second, search->second,  [](const gap & gap1, const gap& gap2){
                            return gap1.gapedSequence_ == gap2.gapedSequence_;
                          })) {
                            pass = false;
                            break;
                          }
                        }
                      }
                      //check insertions
                      //iterate over segregating insertions locations
                      for(const auto & loc : insertionsInOtherMatching) {
                        //determine if current read has a snp at location
                        auto search = currentInsertions.find(loc.first);
                        if(search != currentInsertions.end()) {
                          //if it does have a insertion here, check to see if it is a known variant, if it is the read doesn't pass
                          if(njh::contains(loc.second, search->second,  [](const gap & gap1, const gap& gap2){
                            return gap1.gapedSequence_ == gap2.gapedSequence_;
                          })) {
                            pass = false;
                            break;
                          }
                        }
                      }
                    }
                    if(pass && checkIndelsAgainstSNPsWhenMapping) {
                      auto seqSnpPos = refVariationInfo[refPos].getVariantSnpLoci(matchingRefNames,0);
                      for(const auto & loc : seqSnpPos) {
                        //determine if current read has a snp at location
                        if(currentSnps.find(loc) == currentSnps.end()) {
                          //if there is a gap here, don't believe it either
                          if(curThreadAligner->alignObjectB_.seqBase_.seq_[curThreadAligner->getAlignPosForSeqAPos(loc)] == '-') {
                            pass = false;
                            break;
                          }
                        }
                      }
                    }
                    //if all the segregating check out, put this ref in the matchingRefs
                    if(pass) {
                      matchingRefs.emplace_back(refPos);
                    }
                  }
                  if(matchingRefs.size() == 1) {
                    //if only one matching ref, add to it
                    readCounts[matchingRefs.front()].addRead(allInputReads[readPos]->seqBase_);
                  } else if(matchingRefs.size() == 0) {
                    //if none of the ref work out, put it in indeterminate
                    std::lock_guard<std::mutex> indeterminateLock(indeterminateMut);
                    indeterminate.emplace_back(allInputReads[readPos]->seqBase_);
                  } else {
                    //if there is more than one matching ref, put it in ties
                    std::lock_guard<std::mutex> tiesLock(tiesMut);
                    tiesReads.emplace_back(allInputReads[readPos]->seqBase_);
                    ties.emplace_back(readPos, matchingRefs);
                    const seqInfo & info = allInputReads[readPos]->seqBase_;
                    seqInfo tempObj(info.name_, info.seq_, info.qual_,
                                    info.cnt_ / matchingRefs.size());
                    tempObj.updateName();
                    if(pars.tiesDuringMapping){
                      for (const auto & best : matchingRefs) {
                        readCounts[best].addRead(tempObj);
                      } 
                    }
                  }
                }
              } else {
                std::lock_guard<std::mutex> unmLock(unmappableMut);
                unmappable.emplace_back(allInputReads[readPos]->seqBase_);
              }
            }
            if (setUp.pars_.verbose_) {
              pbar.outputProgAdd(std::cout, readPositionsChunk.size(), true);
            }
          }
        };


    std::vector<std::thread> threads;
    for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
      threads.emplace_back(std::thread(mapRead));
    }
    for(auto & t : threads){
      t.join();
    }


    mappingDir = njh::files::makeDir(setUp.pars_.directoryName_,
                                     njh::files::MkdirPar("mappingInfo")).string();
    SeqOutput::write(unmappable, SeqIOOptions(mappingDir + "unmappable", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
    SeqOutput::write(indeterminate, SeqIOOptions(mappingDir + "indeterminate", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
    std::vector<readObject> tiesObj;
    std::ofstream tiesInfo;
    openTextFile(tiesInfo, mappingDir + "tiesInfo", ".tab.txt",
                 setUp.pars_.ioOptions_.out_);
    tiesInfo << "readName\trefs\n";
    for (const auto & tie : ties) {
      tiesObj.emplace_back(allInputReads[tie.first]->seqBase_);
      tiesInfo << allInputReads[tie.first]->seqBase_.name_ << "\t";
      VecStr refNames;
      for (const auto & r : tie.second) {
        refNames.emplace_back(readCounts[r].seqBase_.name_);
      }
      tiesInfo << vectorToString(refNames, ",") << "\n";
    }
    SeqOutput::write(tiesObj, SeqIOOptions(mappingDir + "ties", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
    unmappableAmount = std::accumulate(unmappable.begin(),
                                       unmappable.end(), 0.0,
                                       [](double init, const seqInfo & seq) {return init + seq.cnt_;});
    indeterminateAmount = std::accumulate(indeterminate.begin(),
                                          indeterminate.end(), 0.0,
                                          [](double init, const seqInfo & seq) {return init + seq.cnt_;});
    tiesAmount  = std::accumulate(tiesReads.begin(),
                                  tiesReads.end(), 0.0,
                                  [](double init, const seqInfo & seq) {return init + seq.cnt_;});
    {
      if(pars.map){
        std::ofstream mapInfo;
        openTextFile(mapInfo, mappingDir + "initialMapInfo", ".tab.txt", setUp.pars_.ioOptions_.out_);
        std::set<std::string> repNames;
        for (const auto & read : allInputReads) {
          if("" != pars.aSetRepName){
            repNames.emplace(pars.aSetRepName);
          }else{
            auto toks = tokenizeString(read->seqBase_.name_, ".");
            repNames.emplace(njh::replaceString(toks[0], "_Comp", ""));
          }

        }
        mapInfo << "seqName\trefId\tftotalReadsMapped\ttotalMappedFraction";
        for (const auto & n : repNames) {
          mapInfo << "\t" << n << "_readsMapped\t" << n << "_mappedFraction";
        }
        if (setUp.pars_.refIoOptions_.firstName_ != "") {
          mapInfo << "\tBestRef\tscore\t1bIndel\t2bI"
                     "ndel\t>2bIndel\tlqMismatch\thqMismatch";
        }
        mapInfo << "\n";
        std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
        std::map<std::string, double> repTotals;
        for (const auto & ref : readCounts) {
          for (const auto & read : ref.reads_) {
            if ("" != pars.aSetRepName) {
              repTotals[pars.aSetRepName] += read.cnt_;
            } else {
              auto toks = tokenizeString(read.name_, ".");
              repTotals[njh::replaceString(toks[0], "_Comp", "")] +=
                  read.cnt_;
            }
          }
        }

        for (auto & ref : readCounts) {
          std::map<std::string, double> repCounts;
          for (const auto & read : ref.reads_) {
            auto toks = tokenizeString(read.name_, ".");
            repCounts[njh::replaceString(toks[0], "_Comp", "")] += read.cnt_;
          }
          mapInfo << seqName << "\t" << ref.seqBase_.name_ << "\t"
                  << ref.seqBase_.cnt_ << "\t" << ref.seqBase_.cnt_ / totalReadCnt;
          for (const auto & n : repNames) {
            if (repCounts.find(n) != repCounts.end()) {
              mapInfo << "\t" << repCounts[n] << "\t"
                      << repCounts[n] / repTotals[n];
            } else {
              mapInfo << "\t0\t0";
            }
          }
          if (setUp.pars_.refIoOptions_.firstName_ != "") {
            bool eventBased = true;
            mapInfo << "\t"
                    << profiler::compareToRefSingle(refSeqs, ref, alignerObj,
                                                    setUp.pars_.local_, eventBased).front();
          }
          mapInfo << "\n";
        }
        mapInfo << seqName << "\t" << "unmappable\t" << unmappableAmount << "\t"
                << unmappableAmount / totalReadCnt << "\n";
        mapInfo << seqName << "\t" << "indeterminate\t" << indeterminateAmount
                << "\t" << indeterminateAmount / totalReadCnt << "\n";
        mapInfo << seqName << "\t" << "ties\t" << tiesAmount
                << "\t" << tiesAmount / totalReadCnt << "\n";
      }
    }
  }
  if(pars.map){

    setUp.rLog_.logCurrentTime("Setting Aligner for re-calc");
    //alignPoolMap.destoryAligners();
    aligner alignerObjReCalc(maxReadLength, setUp.pars_.gapInfoRef_,
                             setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                             setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                             setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
    alignerObjReCalc.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
    if(pars.recalcConsensus){
      setUp.rLog_.logCurrentTime("Re-calculating consensus");
      if(setUp.pars_.verbose_){
        std::cout << "Re-calculating Consensus Sequences" << std::endl;
      }
    }else{
      setUp.rLog_.logCurrentTime("Setting Map Info");
    }

    std::vector<uint32_t> refReadPositions(readCounts.size());
    njh::iota<uint32_t>(refReadPositions, 0);
    njh::concurrent::LockableQueue<uint32_t> refReadPosQueue(refReadPositions);
    concurrent::AlignerPool alignPoolReCacl(alignerObjReCalc, pars.numThreads);
    alignPoolReCacl.initAligners();
    alignPoolReCacl.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
    auto setMappingInfo = [&refReadPosQueue,&consensusReads,
        &alignPoolReCacl,&readCounts,&pars,
        &clustersizeCutOff](){
      uint32_t readPos = 0;
      auto curThreadAligner = alignPoolReCacl.popAligner();
      while(refReadPosQueue.getVal(readPos)){

        auto foundReadPos = readVec::getReadIndexByName(consensusReads,
                                                        readCounts[readPos].seqBase_.name_);
//				std::cout << "foundReadPos: " << foundReadPos << std::endl;
//				std::cout << "readCounts[readPos].reads_.size(): " << readCounts[readPos].reads_.size() << std::endl;
        if (foundReadPos != std::string::npos) {
          if(readCounts[readPos].reads_.empty()){
            //if no reads were mapped to reference, remove it
            consensusReads[foundReadPos].remove = true;
            continue;
          }
          //update the the counts and fractions for the orignal clusters
          consensusReads[foundReadPos].seqBase_.cnt_ = readCounts[readPos].seqBase_.cnt_;
          consensusReads[foundReadPos].seqBase_.frac_ = readCounts[readPos].seqBase_.frac_;
          //add the reads from the refReads
          consensusReads[foundReadPos].reads_.clear();
          for(const auto & seq : readCounts[readPos].reads_){
            consensusReads[foundReadPos].reads_.emplace_back(std::make_shared<readObject>(seq));
          }
          consensusReads[foundReadPos].updateName();
          //re-calculate the consensus if indicated
          if(pars.recalcConsensus){
            readVecSorter::sortByTotalCountAE<seqInfo>(readCounts[readPos].reads_,true);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
            auto recalcConSeqInfo = getConsensusWithAReCalc(readCounts[readPos].reads_, *curThreadAligner, consensusReads[foundReadPos].seqBase_.name_);
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
            consensusReads[foundReadPos].seqBase_.seq_ = recalcConSeqInfo.seqBase_.seq_;
            consensusReads[foundReadPos].seqBase_.qual_ = recalcConSeqInfo.seqBase_.qual_;
            if(consensusReads[foundReadPos].seqBase_.cnt_ < clustersizeCutOff){
              consensusReads[foundReadPos].remove = true;
            }
          }
        }else{
          std::stringstream ss;
          ss << njh::bashCT::red << "Error, couldn't find: " << readCounts[readPos].seqBase_.name_ << njh::bashCT::reset << "\n";
          throw std::runtime_error{ss.str()};
        }
      }
    };

    std::vector<std::thread> mapInfoThreads;
    for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
      mapInfoThreads.emplace_back(std::thread(setMappingInfo));
    }

    for(auto & t : mapInfoThreads){
      t.join();
    }

    /*
    for (const auto & refRead : readCounts) {
      const auto foundReadPos = readVec::getReadIndexByName(consensusReads,
          refRead.seqBase_.name_);
      if (std::string::npos != foundReadPos) {
        if(refRead.reads_.empty()){
          //if no reads were mapped to reference, remove it
          consensusReads[foundReadPos].remove = true;
          continue;
        }
        std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        //update the the counts and fractions for the orignal clusters
        consensusReads[foundReadPos].seqBase_.cnt_ = refRead.seqBase_.cnt_;
        consensusReads[foundReadPos].seqBase_.frac_ = refRead.seqBase_.frac_;
        //add the reads from the refReads
        std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        consensusReads[foundReadPos].reads_.clear();
        std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        for(const auto & seq : refRead.reads_){
          consensusReads[foundReadPos].reads_.emplace_back(std::make_shared<readObject>(seq));
        }
        std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        consensusReads[foundReadPos].updateName();
        std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        //re-calculate the consensus if indicated
        if(pars.recalcConsensus){
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          std::cout << "refRead.reads_.size() : " << refRead.reads_.size() << std::endl;
          auto recalcConSeqInfo = getConsensusWithAReCalc(refRead.reads_, alignerObj, consensusReads[foundReadPos].seqBase_.name_);
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          consensusReads[foundReadPos].seqBase_.seq_ = recalcConSeqInfo.seqBase_.seq_;
          consensusReads[foundReadPos].seqBase_.qual_ = recalcConSeqInfo.seqBase_.qual_;
        }
        std::cout << __FILE__ << "  " << __LINE__ << std::endl;
      }else{
        std::stringstream ss;
        ss << njh::bashCT::red << "Error, couldn't find: " << refRead.seqBase_.name_ << njh::bashCT::reset << "\n";
        throw std::runtime_error{ss.str()};
      }
    }*/



    //get rid of the empty clusters
    consensusReads = readVecSplitter::splitVectorOnRemove(consensusReads).first;

    if(pars.recalcConsensus){
      setUp.rLog_.logCurrentTime("Collapsing after recalculating consensus sequences");
      if (setUp.pars_.verbose_) {
        setUp.rLog_.logCurrentTime("Comparing After Mapping");
        std::cout << njh::bashCT::bold << "Comparing sequences after mapping!" << njh::bashCT::reset
                  << std::endl;
      }
      //collapse clusters again if the consensus had to be calculated
      simpleCollapseWithPars(consensusReads, alignerObj, postCollapsePars);
      if (setUp.pars_.verbose_) {
        std::cout << njh::bashCT::bold << "There are " << consensusReads.size() << " after comparing" << njh::bashCT::reset
                  << std::endl;
      }
      if(pars.writeInitialClusters){
        //write out initial clusters if needed, this is mostly a debugging option
        std::string initialDir = njh::files::makeDir(njh::files::make_path(setUp.pars_.directoryName_, "initialClusters").string(), njh::files::MkdirPar("initialClustersCollapsed")).string();
        std::string initialClustersDir = njh::files::makeDir(initialDir, njh::files::MkdirPar("clusters")).string();
        SeqOutput::write(consensusReads, SeqIOOptions(initialDir + "allConsensus",setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
        for(const auto & clus : consensusReads){
          clus.writeClustersInDir(initialClustersDir, setUp.pars_.ioOptions_);
        }
      }
    }
  }


  if(setUp.pars_.debug_){
    std::string snpDir = njh::files::makeDir(setUp.pars_.directoryName_,
                                             njh::files::MkdirPar("debugInternalSnpInfo")).string();
    readVec::allSetFractionByTotalCount(consensusReads);
    for (const auto readPos : iter::range(consensusReads.size())) {
      //write out cluster info
      //log snp information
      /**@todo add indel information, less informative since pacbio have a lot*/
      std::unordered_map<uint32_t, std::unordered_map<char, VecStr>> mismatches;
      for (const auto subReadPos : iter::range(
          consensusReads[readPos].reads_.size())) {
        const auto & subRead = consensusReads[readPos].reads_[subReadPos];
        alignerObj.alignCache(consensusReads[readPos], subRead, false);
        //count gaps and mismatches and get identity
        alignerObj.profilePrimerAlignment(consensusReads[readPos], subRead);
        for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
          if(m.second.highQuality(setUp.pars_.qScorePars_)){
            mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(
                subRead->seqBase_.name_);
          }
        }
      }
      table misTab { VecStr { "refPos", "refBase", "seqBase", "freq", "fraction",
                              "seqs" } };
      for (const auto & m : mismatches) {
        for (const auto & seqM : m.second) {

          misTab.content_.emplace_back(
              toVecStr(m.first, consensusReads[readPos].seqBase_.seq_[m.first],
                       seqM.first, seqM.second.size(),
                       seqM.second.size() / consensusReads[readPos].seqBase_.cnt_,
                       vectorToString(seqM.second, ",")));
        }
      }
      misTab.sortTable("seqBase", false);
      misTab.sortTable("refPos", false);
      misTab.outPutContents(
          TableIOOpts(
              OutOptions(snpDir + consensusReads[readPos].seqBase_.name_,
                         ".tab.txt"), "\t", misTab.hasHeader_));
    }
  }


  setUp.rLog_.logCurrentTime("Determining Breakouts");
  bool breakout = false;
  uint32_t breakoutCount = 0;
  for (const auto readPos : iter::range(consensusReads.size())) {
    //log snp information
    auto breakouts = consensusReads[readPos].breakoutClustersBasedOnSnps(alignerObj, pars.breakoutPars);
    if(!breakouts.empty()){
      breakout = true;
      breakoutCount+= breakouts.size();
      addOtherVec(consensusReads, breakouts);
      if(setUp.pars_.debug_){
        std::cout << consensusReads[readPos].seqBase_.name_ << std::endl;
        for(const auto & breakoutPos : iter::range(breakouts.size())){
          std::cout << "\t" << breakouts[breakoutPos].seqBase_.name_ << std::endl;
        }
      }
    }
  }



  if(breakout){
    std::vector<cluster> aboveCutOffs;
    uint32_t belowCutOff = 0;
    for(const auto & seq : consensusReads){
      if(seq.seqBase_.cnt_ > clustersizeCutOff){
        aboveCutOffs.emplace_back(seq);
      }else{
        ++belowCutOff;
      }
    }
    if(belowCutOff > 0){
      consensusReads = aboveCutOffs;
    }

    simpleCollapseWithPars(consensusReads, alignerObj, postCollapsePars);
    renameReadNames(consensusReads,
                    bfs::basename(setUp.pars_.ioOptions_.firstName_), true, true,
                    true, "totalCount");
    if(pars.writeInitialClusters){
      //write out initial clusters if needed, this is mostly a debugging option
      std::string initialDir = njh::files::makeDir(njh::files::make_path(setUp.pars_.directoryName_, "initialClusters").string(), njh::files::MkdirPar("initialClustersAfrerBreakoutCollapsed")).string();
      std::string initialClustersDir = njh::files::makeDir(initialDir, njh::files::MkdirPar("clusters")).string();
      SeqOutput::write(consensusReads, SeqIOOptions(initialDir + "allConsensus",setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      for(const auto & clus : consensusReads){
        clus.writeClustersInDir(initialClustersDir, setUp.pars_.ioOptions_);
      }
    }
  }

  if(breakout){
    if (setUp.pars_.verbose_) {
      std::cout << njh::bashCT::bold << "Breakout!"
                << njh::bashCT::reset << std::endl;
    }

    if(setUp.pars_.debug_){
      auto breakOpts = SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, "debug_breakoutConsensus.fastq"));
      breakOpts.out_.overWriteFile_ = true;
      SeqOutput::write(consensusReads, breakOpts);
      OutputStream infoFile(njh::files::make_path(setUp.pars_.directoryName_, "debug_breakoutConsensus_info.tab.txt"));
      infoFile << "ClusterNumber\tClusterId\tTotalCount\tTotalClusters\tReadFraction\n";
      for (const auto readPos : iter::range(consensusReads.size())) {
        //write out cluster info
        infoFile << readPos
                 << '\t' << consensusReads[readPos].seqBase_.name_
                 << '\t' << consensusReads[readPos].seqBase_.cnt_
                 << '\t' << consensusReads[readPos].reads_.size()
                 << '\t' << consensusReads[readPos].seqBase_.frac_ << '\n';
      }
    }

    if (pars.map) {

      readCounts.clear();
      //njh::files::rmDirForce(mappingDir);
      njh::files::bfs::rename(mappingDir, njh::files::make_path(setUp.pars_.directoryName_, "initialMappingInfo"));

      unmappableAmount = 0.0;
      indeterminateAmount = 0.0;
      tiesAmount = 0.0;
      setUp.rLog_.logCurrentTime("Determining Variation ...again");
      if (setUp.pars_.verbose_) {
        std::cout << njh::bashCT::bold << "There was a " << njh::bashCT::red
                  << "breakout" << njh::bashCT::resetAdd(njh::bashCT::bold)
                  << " of " << breakoutCount << " clusters, now there are " << consensusReads.size() << " clusters after collapsing, re-mapping" << njh::bashCT::reset << std::endl;
      }

      //create aligner with gap info from ref so a different gap scoring can be used for mapping determination
      readVec::getMaxLength(consensusReads, maxReadLength);
      aligner alignerObjMapper(maxReadLength, setUp.pars_.gapInfoRef_,
                               setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                               setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                               setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);

      alignerObjMapper.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

      //set all reads' fraction
      njh::for_each(allInputReads, [&totalReadCnt](auto& seq) { getSeqBase(seq).setFractionByCount(totalReadCnt); });
      //create ref map containers to keep counts
      readCounts = refMapContainer<seqInfo>::createContainers<seqInfo, cluster>(consensusReads);
      //determine variation in final consensus sequences
      std::vector<refVariants> refVariationInfo;

      for (const auto refPos : iter::range(readCounts.size())) {
        refVariationInfo.emplace_back(readCounts[refPos].seqBase_);
      }

      for (const auto refPos : iter::range(readCounts.size())) {
        for (const auto refSubPos : iter::range(readCounts.size())) {
          if (refPos == refSubPos) {
            continue;
          }
          if(uAbsdiff(len(readCounts[refSubPos].seqBase_), len(readCounts[refPos].seqBase_)) > pars.readLengthMinDiff){
            continue;
          }
          refVariationInfo[refPos].addVariant(readCounts[refSubPos].seqBase_,
                                              alignerObjMapper, false);
        }
      }

      std::vector<kmerInfo> conKInfos(consensusReads.size());

      setUp.rLog_.logCurrentTime("Calculating Consensus Kmer info for comparison ...again");

      {
        std::vector<uint32_t> conReadPositions(consensusReads.size());
        njh::iota<uint32_t>(conReadPositions, 0);
        njh::concurrent::LockableQueue<uint32_t> conReadPosQueue(
            conReadPositions);
        auto fillInfos =
            [&consensusReads,&conKInfos,&conReadPosQueue,&kLenComp]() {
              uint32_t pos = std::numeric_limits<uint32_t>::max();
              while(conReadPosQueue.getVal(pos)) {
                conKInfos[pos] = kmerInfo(consensusReads[pos].seqBase_.seq_, kLenComp, false);
              }
            };
        std::vector<std::thread> threads;
        for (uint32_t tNum = 0; tNum < pars.numThreads; ++tNum) {
          threads.emplace_back(std::thread(fillInfos));
        }
        for (auto & t : threads) {
          t.join();
        }
      }


      setUp.rLog_.logCurrentTime("Mapping Sequences ...again");

      std::vector<std::pair<uint64_t, std::vector<uint64_t>>> ties;
      std::vector<seqInfo> tiesReads;
      std::mutex tiesMut;
      std::vector<seqInfo> unmappable;
      std::mutex unmappableMut;
      std::vector<seqInfo> indeterminate;
      std::mutex indeterminateMut;

      njh::ProgressBar pbar(allInputReads.size());
      pbar.progColors_ = pbar.RdYlGn_;
      //set the alignment best difference and function
      //the best alignment is considered within a certain range due such high error rates the best alignment is not
      //Necessary the one with the best score, it is the one that contains the segregating differences
      double difference = 0;
      std::function<double(const aligner &)> getScoreFunc;
      if (pars.byScore) {
        getScoreFunc = [](const aligner & alignerObject){
          return alignerObject.parts_.score_;
        };
        difference = 10;
      } else {
        getScoreFunc = [](const aligner & alignerObject){
          //return alignerObjMapper.comp_.distances_.eventBasedIdentity_;
          return alignerObject.comp_.distances_.eventBasedIdentityHq_;
        };
        difference = 0.001;
      }

      concurrent::AlignerPool alignPoolMap(alignerObjMapper, pars.numThreads);
      alignPoolMap.initAligners();
      alignPoolMap.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

      std::vector<uint32_t> readPositions(allInputReads.size());
      njh::iota<uint32_t>(readPositions, 0);
      njh::concurrent::LockableQueue<uint32_t> readPosQueue(readPositions);

      auto mapRead =
          [&pbar, &readPosQueue, &alignPoolMap, &pars, &getScoreFunc, &readCounts, &allInputReads,
              &tiesReads, &ties, &tiesMut, &unmappable, &unmappableMut, &indeterminate, &indeterminateMut,
              &difference, &refVariationInfo,&threadChunkNum, &setUp,
              &checkIndelsWhenMapping,&checkIndelsAgainstSNPsWhenMapping]() {
            std::vector<uint32_t> readPositionsChunk;
            auto curThreadAligner = alignPoolMap.popAligner();
            while(readPosQueue.getVals(readPositionsChunk, threadChunkNum)) {
              for(const auto readPos : readPositionsChunk){
                //best score to be used to get all the references that map within a certain distance to this latter
                double bestScore = 0;
                bool mappable = false;
                for (const auto refPos : iter::range(readCounts.size())) {
                  /*
                   * 							if(conKInfos[refPos].compareKmers(allInputReads[readPos]->kInfo_).second < kDistCutOff){
                    continue;
                  }
                   */
                  const auto & ref = readCounts[refPos];
                  if(uAbsdiff(len(getSeqBase(ref)), len(allInputReads[readPos]->seqBase_)) > pars.readLengthMinDiff){
                    continue;
                  }
                  curThreadAligner->alignCacheGlobal(ref.seqBase_,allInputReads[readPos]->seqBase_);
//                  curThreadAligner->alignCacheGlobalDiag(ref.seqBase_,allInputReads[readPos]->seqBase_);
                  curThreadAligner->profilePrimerAlignment(ref.seqBase_,
                                                           allInputReads[readPos]->seqBase_);
                  //check for very large gaps that wouldn't be expected no matter the technology
                  bool passLargeGaps = true;
                  for(const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                    if(g.second.size_ > pars.largeGapSizeCutOff) {
                      passLargeGaps = false;
                      break;
                    }
                  }
                  //check to see if the read is mappable by comparing to the percent identity cut off and has no large gaps
                  if (curThreadAligner->comp_.distances_.eventBasedIdentity_ >= pars.idCutOff && passLargeGaps) {
                    mappable = true;
                    double currentScore = getScoreFunc(*curThreadAligner);
                    if(currentScore > bestScore) {
                      bestScore = currentScore;
                    }
                  }
                }
                // if mappable get the best references
                if (mappable) {
                  std::vector<uint64_t> bestRefs;
                  for (const auto refPos : iter::range(readCounts.size())) {
                    const auto & ref = readCounts[refPos];
                    curThreadAligner->alignCacheGlobal(ref.seqBase_,allInputReads[readPos]->seqBase_);
//                    curThreadAligner->alignCacheGlobalDiag(ref.seqBase_,allInputReads[readPos]->seqBase_);
                    curThreadAligner->profilePrimerAlignment(ref.seqBase_,
                                                             allInputReads[readPos]->seqBase_);
                    if (curThreadAligner->comp_.distances_.eventBasedIdentity_ >= pars.idCutOff) {
                      mappable = true;
                      double currentScore = getScoreFunc(*curThreadAligner);
                      bool passLargeGaps = true;
                      for(const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                        if(g.second.size_ > pars.largeGapSizeCutOff) {
                          passLargeGaps = false;
                          break;
                        }
                      }
                      //if the current score is a certain amount within the best score and doesn't have large gaps then put it
                      //into the vector for possible match
                      if (std::abs(currentScore - bestScore) < difference && passLargeGaps) {
                        bestRefs.emplace_back(refPos);
                      }
                    }
                  }

                  if (bestRefs.size() == 1) {
                    //if only matching one, then place with that one
                    readCounts[bestRefs.front()].addRead(allInputReads[readPos]->seqBase_);
                  } else if (bestRefs.size() == 0) {
                    //if no mappable candidates found, put into un-mappable
                    std::lock_guard<std::mutex> unmLock(unmappableMut);
                    unmappable.emplace_back(allInputReads[readPos]->seqBase_);
                  } else {
                    VecStr matchingRefNames;
                    for(const auto & refPos : bestRefs) {
                      matchingRefNames.emplace_back(readCounts[refPos].seqBase_.name_);
                    }
                    std::vector<uint64_t> matchingRefs;
                    for (const auto & refPos : bestRefs) {
                      const auto & ref = readCounts[refPos];
                      //get the current snp info for the current ref to the other matching refs
                      //auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 1);
                      auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 0);
                      curThreadAligner->alignCacheGlobal(ref.seqBase_, allInputReads[readPos]->seqBase_);
//                      curThreadAligner->alignCacheGlobalDiag(ref.seqBase_, allInputReads[readPos]->seqBase_);
                      curThreadAligner->profilePrimerAlignment(ref.seqBase_, allInputReads[readPos]->seqBase_);
                      //determine the snps for this current read to this ref
                      std::unordered_map<uint32_t, char> currentSnps;
                      for (const auto & m : curThreadAligner->comp_.distances_.mismatches_) {
                        currentSnps[m.second.refBasePos] = m.second.seqBase;
                      }
                      bool pass = true;
                      //iterate over segregating snps locations
                      for(const auto & loc : seqSnpPosBases) {
                        //determine if current read has a snp at location
                        auto search = currentSnps.find(loc.first);
                        if(search != currentSnps.end()) {
                          //if it does have a snp here, check to see if it is a known variant, if it is the read doesn't pass
                          if(njh::in(search->second, loc.second)) {
                            pass = false;
                            break;
                          }
                        }
                      }
                      if(pass && checkIndelsWhenMapping){
                        //check indels
                        //get the current insertions info for the current ref to the other matching refs
                        auto deletionsInOtherMatching = refVariationInfo[refPos].getVariantDeletionLociMap(matchingRefNames, 2, true);
                        auto insertionsInOtherMatching = refVariationInfo[refPos].getVariantInsertionLociMap(matchingRefNames, 2, true);

                        //determine the insertions for this current read to this ref
                        std::unordered_map<uint32_t, gap> currentInsertions;
                        std::unordered_map<uint32_t, gap> currentDeletions;
                        for (const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                          if(g.second.ref_){
                            currentInsertions.emplace(g.second.refPos_, g.second);
                          }else{
                            currentDeletions.emplace(g.second.refPos_, g.second);
                          }
                        }

                        //check deletions
                        //iterate over segregating deletions locations
                        for(const auto & loc : deletionsInOtherMatching) {
                          //determine if current read has a snp at location
                          auto search = currentDeletions.find(loc.first);
                          if(search != currentDeletions.end()) {
                            //if it does have a deletion here, check to see if it is a known variant, if it is the read doesn't pass
                            if(njh::contains(loc.second, search->second,  [](const gap & gap1, const gap& gap2){
                              return gap1.gapedSequence_ == gap2.gapedSequence_;
                            })) {
                              pass = false;
                              break;
                            }
                          }
                        }
                        //check insertions
                        //iterate over segregating insertions locations
                        for(const auto & loc : insertionsInOtherMatching) {
                          //determine if current read has a snp at location
                          auto search = currentInsertions.find(loc.first);
                          if(search != currentInsertions.end()) {
                            //if it does have a insertion here, check to see if it is a known variant, if it is the read doesn't pass
                            if(njh::contains(loc.second, search->second,  [](const gap & gap1, const gap& gap2){
                              return gap1.gapedSequence_ == gap2.gapedSequence_;
                            })) {
                              pass = false;
                              break;
                            }
                          }
                        }
                      }
                      if(pass && checkIndelsAgainstSNPsWhenMapping) {
                        auto seqSnpPos = refVariationInfo[refPos].getVariantSnpLoci(matchingRefNames,0);
                        for(const auto & loc : seqSnpPos) {
                          //determine if current read has a snp at location
                          if(currentSnps.find(loc) == currentSnps.end()) {
                            //if there is a gap here, don't believe it either
                            if(curThreadAligner->alignObjectB_.seqBase_.seq_[curThreadAligner->getAlignPosForSeqAPos(loc)] == '-') {
                              pass = false;
                              break;
                            }
                          }
                        }
                      }
                      //if all the segregating check out, put this ref in the matchingRefs
                      if(pass) {
                        matchingRefs.emplace_back(refPos);
                      }
                    }
                    if(matchingRefs.size() == 1) {
                      //if only one matching ref, add to it
                      readCounts[matchingRefs.front()].addRead(allInputReads[readPos]->seqBase_);
                    } else if(matchingRefs.size() == 0) {
                      //if none of the ref work out, put it in indeterminate
                      std::lock_guard<std::mutex> indeterminateLock(indeterminateMut);
                      indeterminate.emplace_back(allInputReads[readPos]->seqBase_);
                    } else {
                      //if there is more than one matching ref, put it in ties
                      std::lock_guard<std::mutex> tiesLock(tiesMut);
                      tiesReads.emplace_back(allInputReads[readPos]->seqBase_);
                      ties.emplace_back(readPos, matchingRefs);
                      const seqInfo & info = allInputReads[readPos]->seqBase_;
                      seqInfo tempObj(info.name_, info.seq_, info.qual_,
                                      info.cnt_ / matchingRefs.size());
                      tempObj.updateName();
                      if(pars.tiesDuringMapping){
                        for (const auto & best : matchingRefs) {
                          readCounts[best].addRead(tempObj);
                        } 
                      }
                    }
                  }
                } else {
                  std::lock_guard<std::mutex> unmLock(unmappableMut);
                  unmappable.emplace_back(allInputReads[readPos]->seqBase_);
                }
              }
              if (setUp.pars_.verbose_) {
                pbar.outputProgAdd(std::cout, readPositionsChunk.size(), true);
              }
            }
          };


      std::vector<std::thread> threads;
      for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
        threads.emplace_back(std::thread(mapRead));
      }
      for(auto & t : threads){
        t.join();
      }


      mappingDir = njh::files::makeDir(setUp.pars_.directoryName_,
                                       njh::files::MkdirPar("mappingInfo")).string();
      SeqOutput::write(unmappable, SeqIOOptions(mappingDir + "unmappable", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      SeqOutput::write(indeterminate, SeqIOOptions(mappingDir + "indeterminate", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      std::vector<readObject> tiesObj;
      std::ofstream tiesInfo;
      openTextFile(tiesInfo, mappingDir + "tiesInfo", ".tab.txt",
                   setUp.pars_.ioOptions_.out_);
      tiesInfo << "readName\trefs\n";
      for (const auto & tie : ties) {
        tiesObj.emplace_back(allInputReads[tie.first]->seqBase_);
        tiesInfo << allInputReads[tie.first]->seqBase_.name_ << "\t";
        VecStr refNames;
        for (const auto & r : tie.second) {
          refNames.emplace_back(readCounts[r].seqBase_.name_);
        }
        tiesInfo << vectorToString(refNames, ",") << "\n";
      }
      SeqOutput::write(tiesObj, SeqIOOptions(mappingDir + "ties", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      unmappableAmount = std::accumulate(unmappable.begin(),
                                         unmappable.end(), 0.0,
                                         [](double init, const seqInfo & seq) {return init + seq.cnt_;});
      indeterminateAmount = std::accumulate(indeterminate.begin(),
                                            indeterminate.end(), 0.0,
                                            [](double init, const seqInfo & seq) {return init + seq.cnt_;});
      tiesAmount  = std::accumulate(tiesReads.begin(),
                                    tiesReads.end(), 0.0,
                                   [](double init, const seqInfo & seq) {return init + seq.cnt_;});

      setUp.rLog_.logCurrentTime("Setting Aligner for re-calc ...again" );
      //alignPoolMap.destoryAligners();
      aligner alignerObjReCalc(maxReadLength, setUp.pars_.gapInfoRef_,
                               setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                               setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                               setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
      alignerObjReCalc.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
      if(pars.recalcConsensus){
        setUp.rLog_.logCurrentTime("Re-calculating consensus ...again");
        if(setUp.pars_.verbose_){
          std::cout << "Re-calculating Consensus Sequences ...again" << std::endl;
        }
      }else{
        setUp.rLog_.logCurrentTime("Setting Map Info ...again");
      }

      std::vector<uint32_t> refReadPositions(readCounts.size());
      njh::iota<uint32_t>(refReadPositions, 0);
      njh::concurrent::LockableQueue<uint32_t> refReadPosQueue(refReadPositions);
      concurrent::AlignerPool alignPoolReCacl(alignerObjReCalc, pars.numThreads);
      alignPoolReCacl.initAligners();
      alignPoolReCacl.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
      auto setMappingInfo = [&refReadPosQueue,&consensusReads,
          &alignPoolReCacl,&readCounts,&pars,
          &clustersizeCutOff](){
        uint32_t readPos = 0;
        auto curThreadAligner = alignPoolReCacl.popAligner();
        while(refReadPosQueue.getVal(readPos)){

          auto foundReadPos = readVec::getReadIndexByName(consensusReads,
                                                          readCounts[readPos].seqBase_.name_);
          if (foundReadPos != std::string::npos) {
            if(readCounts[readPos].reads_.empty()){
              //if no reads were mapped to reference, remove it
              consensusReads[foundReadPos].remove = true;
              continue;
            }
            //update the the counts and fractions for the orignal clusters
            consensusReads[foundReadPos].seqBase_.cnt_ = readCounts[readPos].seqBase_.cnt_;
            consensusReads[foundReadPos].seqBase_.frac_ = readCounts[readPos].seqBase_.frac_;
            //add the reads from the refReads
            consensusReads[foundReadPos].reads_.clear();
            for(const auto & seq : readCounts[readPos].reads_){
              consensusReads[foundReadPos].reads_.emplace_back(std::make_shared<readObject>(seq));
            }
            consensusReads[foundReadPos].updateName();
            //re-calculate the consensus if indicated
            if(pars.recalcConsensus){
              readVecSorter::sortByTotalCountAE<seqInfo>(readCounts[readPos].reads_, true);
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
              auto recalcConSeqInfo = getConsensusWithAReCalc(readCounts[readPos].reads_, *curThreadAligner, consensusReads[foundReadPos].seqBase_.name_);
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
              consensusReads[foundReadPos].seqBase_.seq_ = recalcConSeqInfo.seqBase_.seq_;
              consensusReads[foundReadPos].seqBase_.qual_ = recalcConSeqInfo.seqBase_.qual_;
              if(consensusReads[foundReadPos].seqBase_.cnt_ < clustersizeCutOff){
                consensusReads[foundReadPos].remove = true;
              }
            }
          }else{
            std::stringstream ss;
            ss << njh::bashCT::red << "Error, couldn't find: " << readCounts[readPos].seqBase_.name_ << njh::bashCT::reset << "\n";
            throw std::runtime_error{ss.str()};
          }
        }
      };

      std::vector<std::thread> mapInfoThreads;
      for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
        mapInfoThreads.emplace_back(std::thread(setMappingInfo));
      }

      for(auto & t : mapInfoThreads){
        t.join();
      }

      /*
      for (const auto & refRead : readCounts) {
        const auto foundReadPos = readVec::getReadIndexByName(consensusReads,
            refRead.seqBase_.name_);
        if (std::string::npos != foundReadPos) {
          if(refRead.reads_.empty()){
            //if no reads were mapped to reference, remove it
            consensusReads[foundReadPos].remove = true;
            continue;
          }
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          //update the the counts and fractions for the orignal clusters
          consensusReads[foundReadPos].seqBase_.cnt_ = refRead.seqBase_.cnt_;
          consensusReads[foundReadPos].seqBase_.frac_ = refRead.seqBase_.frac_;
          //add the reads from the refReads
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          consensusReads[foundReadPos].reads_.clear();
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          for(const auto & seq : refRead.reads_){
            consensusReads[foundReadPos].reads_.emplace_back(std::make_shared<readObject>(seq));
          }
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          consensusReads[foundReadPos].updateName();
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          //re-calculate the consensus if indicated
          if(pars.recalcConsensus){
            std::cout << __FILE__ << "  " << __LINE__ << std::endl;
            std::cout << "refRead.reads_.size() : " << refRead.reads_.size() << std::endl;
            auto recalcConSeqInfo = getConsensusWithAReCalc(refRead.reads_, alignerObj, consensusReads[foundReadPos].seqBase_.name_);
            std::cout << __FILE__ << "  " << __LINE__ << std::endl;
            consensusReads[foundReadPos].seqBase_.seq_ = recalcConSeqInfo.seqBase_.seq_;
            consensusReads[foundReadPos].seqBase_.qual_ = recalcConSeqInfo.seqBase_.qual_;
          }
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        }else{
          std::stringstream ss;
          ss << njh::bashCT::red << "Error, couldn't find: " << refRead.seqBase_.name_ << njh::bashCT::reset << "\n";
          throw std::runtime_error{ss.str()};
        }
      }*/



      //get rid of the empty clusters
      consensusReads = readVecSplitter::splitVectorOnRemove(consensusReads).first;
      if(pars.recalcConsensus){
        if (setUp.pars_.verbose_) {
          std::cout << njh::bashCT::bold << "Collapsing after recalculating consensus sequences ...again" << njh::bashCT::reset
                    << std::endl;
        }
        setUp.rLog_.logCurrentTime("Collapsing after recalculating consensus sequences ...again");
        //collapse clusters again if the consensus had to be calculated
        simpleCollapseWithPars(consensusReads, alignerObj, postCollapsePars);
      }

      if (setUp.pars_.verbose_) {
        std::cout << njh::bashCT::bold << "Finished with "
                  << consensusReads.size() << " Consensus Sequences after Breakout"
                  << njh::bashCT::reset << std::endl;
      }
      setUp.rLog_ << "Finished with " << consensusReads.size()
                  << " Consensus Sequences after Breakout" << "\n";
      std::ofstream mapInfo;
      openTextFile(mapInfo, mappingDir + "mapInfo", ".tab.txt", setUp.pars_.ioOptions_.out_);
      std::set<std::string> repNames;
      for (const auto & read : allInputReads) {
        if("" != pars.aSetRepName){
          repNames.emplace(pars.aSetRepName);
        }else{
          auto toks = tokenizeString(read->seqBase_.name_, ".");
          repNames.emplace(njh::replaceString(toks[0], "_Comp", ""));
        }

      }
      mapInfo << "seqName\trefId\tftotalReadsMapped\ttotalMappedFraction";
      for (const auto & n : repNames) {
        mapInfo << "\t" << n << "_readsMapped\t" << n << "_mappedFraction";
      }
      if (setUp.pars_.refIoOptions_.firstName_ != "") {
        mapInfo << "\tBestRef\tscore\t1bIndel\t2bI"
                   "ndel\t>2bIndel\tlqMismatch\thqMismatch";
      }
      mapInfo << "\n";
      std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
      std::map<std::string, double> repTotals;
      for (const auto & ref : readCounts) {
        for (const auto & read : ref.reads_) {
          if ("" != pars.aSetRepName) {
            repTotals[pars.aSetRepName] += read.cnt_;
          } else {
            auto toks = tokenizeString(read.name_, ".");
            repTotals[njh::replaceString(toks[0], "_Comp", "")] +=
                read.cnt_;
          }
        }
      }

      for (auto & ref : readCounts) {
        std::map<std::string, double> repCounts;
        for (const auto & read : ref.reads_) {
          auto toks = tokenizeString(read.name_, ".");
          repCounts[njh::replaceString(toks[0], "_Comp", "")] += read.cnt_;
        }
        mapInfo << seqName << "\t" << ref.seqBase_.name_ << "\t"
                << ref.seqBase_.cnt_ << "\t" << ref.seqBase_.cnt_ / totalReadCnt;
        for (const auto & n : repNames) {
          if (repCounts.find(n) != repCounts.end()) {
            mapInfo << "\t" << repCounts[n] << "\t"
                    << repCounts[n] / repTotals[n];
          } else {
            mapInfo << "\t0\t0";
          }
        }
        if (setUp.pars_.refIoOptions_.firstName_ != "") {
          bool eventBased = true;
          mapInfo << "\t"
                  << profiler::compareToRefSingle(refSeqs, ref, alignerObj,
                                                  setUp.pars_.local_, eventBased).front();
        }
        mapInfo << "\n";
      }

      mapInfo << seqName << "\t" << "unmappable\t" << unmappableAmount << "\t"
              << unmappableAmount / totalReadCnt << "\n";
      mapInfo << seqName << "\t" << "indeterminate\t" << indeterminateAmount
              << "\t" << indeterminateAmount / totalReadCnt << "\n";
      mapInfo << seqName << "\t" << "ties\t" << tiesAmount
              << "\t" << tiesAmount / totalReadCnt << "\n";

    }


    renameReadNames(consensusReads,
                    bfs::basename(setUp.pars_.ioOptions_.firstName_), true, true,
                    true, "totalCount");


    if(pars.map){
      readCounts.clear();
      //njh::files::rmDirForce(mappingDir);
//      njh::files::bfs::rename(mappingDir, njh::files::make_path(setUp.pars_.directoryName_, "initialMappingInfo"));

      unmappableAmount = 0.0;
      indeterminateAmount = 0.0;
      tiesAmount = 0.0;
      setUp.rLog_.logCurrentTime("Determining Variation ...again final");
      if (setUp.pars_.verbose_) {
        std::cout << njh::bashCT::bold << "There was a " << njh::bashCT::red
                  << "breakout" << njh::bashCT::resetAdd(njh::bashCT::bold)
                  << " of " << breakoutCount << " clusters, now there are " << consensusReads.size() << " clusters after a remap and re-collapsing, re-mapping final time" << njh::bashCT::reset << std::endl;
      }

      //create aligner with gap info from ref so a different gap scoring can be used for mapping determination
      readVec::getMaxLength(consensusReads, maxReadLength);
      aligner alignerObjMapper(maxReadLength, setUp.pars_.gapInfoRef_,
                               setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                               setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                               setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);

      alignerObjMapper.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

      //set all reads' fraction
      njh::for_each(allInputReads, [&totalReadCnt](auto& seq) { getSeqBase(seq).setFractionByCount(totalReadCnt); });
      //create ref map containers to keep counts
      readCounts = refMapContainer<seqInfo>::createContainers<seqInfo, cluster>(consensusReads);
      //determine variation in final consensus sequences
      std::vector<refVariants> refVariationInfo;

      for (const auto refPos : iter::range(readCounts.size())) {
        refVariationInfo.emplace_back(readCounts[refPos].seqBase_);
      }

      for (const auto refPos : iter::range(readCounts.size())) {
        for (const auto refSubPos : iter::range(readCounts.size())) {
          if (refPos == refSubPos) {
            continue;
          }
          if(uAbsdiff(len(readCounts[refSubPos].seqBase_), len(readCounts[refPos].seqBase_)) > pars.readLengthMinDiff){
            continue;
          }
          refVariationInfo[refPos].addVariant(readCounts[refSubPos].seqBase_,
                                              alignerObjMapper, false);
        }
      }

      std::vector<kmerInfo> conKInfos(consensusReads.size());

      setUp.rLog_.logCurrentTime("Calculating Consensus Kmer info for comparison ...again final");

      {
        std::vector<uint32_t> conReadPositions(consensusReads.size());
        njh::iota<uint32_t>(conReadPositions, 0);
        njh::concurrent::LockableQueue<uint32_t> conReadPosQueue(
            conReadPositions);
        auto fillInfos =
            [&consensusReads,&conKInfos,&conReadPosQueue,&kLenComp]() {
              uint32_t pos = std::numeric_limits<uint32_t>::max();
              while(conReadPosQueue.getVal(pos)) {
                conKInfos[pos] = kmerInfo(consensusReads[pos].seqBase_.seq_, kLenComp, false);
              }
            };
        std::vector<std::thread> threads;
        for (uint32_t tNum = 0; tNum < pars.numThreads; ++tNum) {
          threads.emplace_back(std::thread(fillInfos));
        }
        for (auto & t : threads) {
          t.join();
        }
      }


      setUp.rLog_.logCurrentTime("Mapping Sequences ...again final");

      std::vector<std::pair<uint64_t, std::vector<uint64_t>>> ties;
      std::vector<seqInfo> tiesReads;
      std::mutex tiesMut;
      std::vector<seqInfo> unmappable;
      std::mutex unmappableMut;
      std::vector<seqInfo> indeterminate;
      std::mutex indeterminateMut;

      njh::ProgressBar pbar(allInputReads.size());
      pbar.progColors_ = pbar.RdYlGn_;
      //set the alignment best difference and function
      //the best alignment is considered within a certain range due such high error rates the best alignment is not
      //Necessary the one with the best score, it is the one that contains the segregating differences
      double difference = 0;
      std::function<double(const aligner &)> getScoreFunc;
      if (pars.byScore) {
        getScoreFunc = [](const aligner & alignerObject){
          return alignerObject.parts_.score_;
        };
        difference = 10;
      } else {
        getScoreFunc = [](const aligner & alignerObject){
          //return alignerObjMapper.comp_.distances_.eventBasedIdentity_;
          return alignerObject.comp_.distances_.eventBasedIdentityHq_;
        };
        difference = 0.001;
      }

      concurrent::AlignerPool alignPoolMap(alignerObjMapper, pars.numThreads);
      alignPoolMap.initAligners();
      alignPoolMap.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

      std::vector<uint32_t> readPositions(allInputReads.size());
      njh::iota<uint32_t>(readPositions, 0);
      njh::concurrent::LockableQueue<uint32_t> readPosQueue(readPositions);

      auto mapRead =
          [&pbar, &readPosQueue, &alignPoolMap, &pars, &getScoreFunc, &readCounts, &allInputReads,
              &tiesReads, &ties, &tiesMut, &unmappable, &unmappableMut, &indeterminate, &indeterminateMut,
              &difference, &refVariationInfo,&threadChunkNum, &setUp,
              &checkIndelsWhenMapping,&checkIndelsAgainstSNPsWhenMapping]() {
            std::vector<uint32_t> readPositionsChunk;
            auto curThreadAligner = alignPoolMap.popAligner();
            while(readPosQueue.getVals(readPositionsChunk, threadChunkNum)) {
              for(const auto readPos : readPositionsChunk){
                //best score to be used to get all the references that map within a certain distance to this latter
                double bestScore = 0;
                bool mappable = false;
                for (const auto refPos : iter::range(readCounts.size())) {
                  /*
                   * 							if(conKInfos[refPos].compareKmers(allInputReads[readPos]->kInfo_).second < kDistCutOff){
                    continue;
                  }
                   */
                  const auto & ref = readCounts[refPos];
                  if(uAbsdiff(len(getSeqBase(ref)), len(allInputReads[readPos]->seqBase_)) > pars.readLengthMinDiff){
                    continue;
                  }
                  curThreadAligner->alignCacheGlobal(ref.seqBase_,allInputReads[readPos]->seqBase_);
//                  curThreadAligner->alignCacheGlobalDiag(ref.seqBase_,allInputReads[readPos]->seqBase_);
                  curThreadAligner->profilePrimerAlignment(ref.seqBase_,
                                                           allInputReads[readPos]->seqBase_);
                  //check for very large gaps that wouldn't be expected no matter the technology
                  bool passLargeGaps = true;
                  for(const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                    if(g.second.size_ > pars.largeGapSizeCutOff) {
                      passLargeGaps = false;
                      break;
                    }
                  }
                  //check to see if the read is mappable by comparing to the percent identity cut off and has no large gaps
                  if (curThreadAligner->comp_.distances_.eventBasedIdentity_ >= pars.idCutOff && passLargeGaps) {
                    mappable = true;
                    double currentScore = getScoreFunc(*curThreadAligner);
                    if(currentScore > bestScore) {
                      bestScore = currentScore;
                    }
                  }
                }
                // if mappable get the best references
                if (mappable) {
                  std::vector<uint64_t> bestRefs;
                  for (const auto refPos : iter::range(readCounts.size())) {
                    const auto & ref = readCounts[refPos];
                    curThreadAligner->alignCacheGlobal(ref.seqBase_,allInputReads[readPos]->seqBase_);
//                    curThreadAligner->alignCacheGlobalDiag(ref.seqBase_,allInputReads[readPos]->seqBase_);
                    curThreadAligner->profilePrimerAlignment(ref.seqBase_,
                                                             allInputReads[readPos]->seqBase_);
                    if (curThreadAligner->comp_.distances_.eventBasedIdentity_ >= pars.idCutOff) {
                      mappable = true;
                      double currentScore = getScoreFunc(*curThreadAligner);
                      bool passLargeGaps = true;
                      for(const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                        if(g.second.size_ > pars.largeGapSizeCutOff) {
                          passLargeGaps = false;
                          break;
                        }
                      }
                      //if the current score is a certain amount within the best score and doesn't have large gaps then put it
                      //into the vector for possible match
                      if (std::abs(currentScore - bestScore) < difference && passLargeGaps) {
                        bestRefs.emplace_back(refPos);
                      }
                    }
                  }

                  if (bestRefs.size() == 1) {
                    //if only matching one, then place with that one
                    readCounts[bestRefs.front()].addRead(allInputReads[readPos]->seqBase_);
                  } else if (bestRefs.size() == 0) {
                    //if no mappable candidates found, put into un-mappable
                    std::lock_guard<std::mutex> unmLock(unmappableMut);
                    unmappable.emplace_back(allInputReads[readPos]->seqBase_);
                  } else {
                    VecStr matchingRefNames;
                    for(const auto & refPos : bestRefs) {
                      matchingRefNames.emplace_back(readCounts[refPos].seqBase_.name_);
                    }
                    std::vector<uint64_t> matchingRefs;
                    for (const auto & refPos : bestRefs) {
                      const auto & ref = readCounts[refPos];
                      //get the current snp info for the current ref to the other matching refs
                      //auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 1);
                      auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 0);
                      curThreadAligner->alignCacheGlobal(ref.seqBase_, allInputReads[readPos]->seqBase_);
//                      curThreadAligner->alignCacheGlobalDiag(ref.seqBase_, allInputReads[readPos]->seqBase_);
                      curThreadAligner->profilePrimerAlignment(ref.seqBase_, allInputReads[readPos]->seqBase_);
                      //determine the snps for this current read to this ref
                      std::unordered_map<uint32_t, char> currentSnps;
                      for (const auto & m : curThreadAligner->comp_.distances_.mismatches_) {
                        currentSnps[m.second.refBasePos] = m.second.seqBase;
                      }
                      bool pass = true;
                      //iterate over segregating snps locations
                      for(const auto & loc : seqSnpPosBases) {
                        //determine if current read has a snp at location
                        auto search = currentSnps.find(loc.first);
                        if(search != currentSnps.end()) {
                          //if it does have a snp here, check to see if it is a known variant, if it is the read doesn't pass
                          if(njh::in(search->second, loc.second)) {
                            pass = false;
                            break;
                          }
                        }
                      }
                      if(pass && checkIndelsWhenMapping){
                        //check indels
                        //get the current insertions info for the current ref to the other matching refs
                        auto deletionsInOtherMatching = refVariationInfo[refPos].getVariantDeletionLociMap(matchingRefNames, 2, true);
                        auto insertionsInOtherMatching = refVariationInfo[refPos].getVariantInsertionLociMap(matchingRefNames, 2, true);

                        //determine the insertions for this current read to this ref
                        std::unordered_map<uint32_t, gap> currentInsertions;
                        std::unordered_map<uint32_t, gap> currentDeletions;
                        for (const auto & g : curThreadAligner->comp_.distances_.alignmentGaps_) {
                          if(g.second.ref_){
                            currentInsertions.emplace(g.second.refPos_, g.second);
                          }else{
                            currentDeletions.emplace(g.second.refPos_, g.second);
                          }
                        }

                        //check deletions
                        //iterate over segregating deletions locations
                        for(const auto & loc : deletionsInOtherMatching) {
                          //determine if current read has a snp at location
                          auto search = currentDeletions.find(loc.first);
                          if(search != currentDeletions.end()) {
                            //if it does have a deletion here, check to see if it is a known variant, if it is the read doesn't pass
                            if(njh::contains(loc.second, search->second,  [](const gap & gap1, const gap& gap2){
                              return gap1.gapedSequence_ == gap2.gapedSequence_;
                            })) {
                              pass = false;
                              break;
                            }
                          }
                        }
                        //check insertions
                        //iterate over segregating insertions locations
                        for(const auto & loc : insertionsInOtherMatching) {
                          //determine if current read has a snp at location
                          auto search = currentInsertions.find(loc.first);
                          if(search != currentInsertions.end()) {
                            //if it does have a insertion here, check to see if it is a known variant, if it is the read doesn't pass
                            if(njh::contains(loc.second, search->second,  [](const gap & gap1, const gap& gap2){
                              return gap1.gapedSequence_ == gap2.gapedSequence_;
                            })) {
                              pass = false;
                              break;
                            }
                          }
                        }
                      }
                      if(pass && checkIndelsAgainstSNPsWhenMapping) {
                        auto seqSnpPos = refVariationInfo[refPos].getVariantSnpLoci(matchingRefNames,0);
                        for(const auto & loc : seqSnpPos) {
                          //determine if current read has a snp at location
                          if(currentSnps.find(loc) == currentSnps.end()) {
                            //if there is a gap here, don't believe it either
                            if(curThreadAligner->alignObjectB_.seqBase_.seq_[curThreadAligner->getAlignPosForSeqAPos(loc)] == '-') {
                              pass = false;
                              break;
                            }
                          }
                        }
                      }
                      //if all the segregating check out, put this ref in the matchingRefs
                      if(pass) {
                        matchingRefs.emplace_back(refPos);
                      }
                    }
                    if(matchingRefs.size() == 1) {
                      //if only one matching ref, add to it
                      readCounts[matchingRefs.front()].addRead(allInputReads[readPos]->seqBase_);
                    } else if(matchingRefs.size() == 0) {
                      //if none of the ref work out, put it in indeterminate
                      std::lock_guard<std::mutex> indeterminateLock(indeterminateMut);
                      indeterminate.emplace_back(allInputReads[readPos]->seqBase_);
                    } else {
                      //if there is more than one matching ref, put it in ties
                      std::lock_guard<std::mutex> tiesLock(tiesMut);
                      tiesReads.emplace_back(allInputReads[readPos]->seqBase_);
                      ties.emplace_back(readPos, matchingRefs);
                      const seqInfo & info = allInputReads[readPos]->seqBase_;
                      seqInfo tempObj(info.name_, info.seq_, info.qual_,
                                      info.cnt_ / matchingRefs.size());
                      tempObj.updateName();
                      if(pars.tiesDuringMapping){
                        for (const auto & best : matchingRefs) {
                          readCounts[best].addRead(tempObj);
                        }
                      }
                    }
                  }
                } else {
                  std::lock_guard<std::mutex> unmLock(unmappableMut);
                  unmappable.emplace_back(allInputReads[readPos]->seqBase_);
                }
              }
              if (setUp.pars_.verbose_) {
                pbar.outputProgAdd(std::cout, readPositionsChunk.size(), true);
              }
            }
          };


      std::vector<std::thread> threads;
      for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
        threads.emplace_back(std::thread(mapRead));
      }
      for(auto & t : threads){
        t.join();
      }


      mappingDir = njh::files::makeDir(setUp.pars_.directoryName_,
                                       njh::files::MkdirPar("finalMappingInfo")).string();
      SeqOutput::write(unmappable, SeqIOOptions(mappingDir + "unmappable", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      SeqOutput::write(indeterminate, SeqIOOptions(mappingDir + "indeterminate", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      std::vector<readObject> tiesObj;
      std::ofstream tiesInfo;
      openTextFile(tiesInfo, mappingDir + "tiesInfo", ".tab.txt",
                   setUp.pars_.ioOptions_.out_);
      tiesInfo << "readName\trefs\n";
      for (const auto & tie : ties) {
        tiesObj.emplace_back(allInputReads[tie.first]->seqBase_);
        tiesInfo << allInputReads[tie.first]->seqBase_.name_ << "\t";
        VecStr refNames;
        for (const auto & r : tie.second) {
          refNames.emplace_back(readCounts[r].seqBase_.name_);
        }
        tiesInfo << vectorToString(refNames, ",") << "\n";
      }
      SeqOutput::write(tiesObj, SeqIOOptions(mappingDir + "ties", setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      unmappableAmount = std::accumulate(unmappable.begin(),
                                         unmappable.end(), 0.0,
                                         [](double init, const seqInfo & seq) {return init + seq.cnt_;});
      indeterminateAmount = std::accumulate(indeterminate.begin(),
                                            indeterminate.end(), 0.0,
                                            [](double init, const seqInfo & seq) {return init + seq.cnt_;});
      tiesAmount  = std::accumulate(tiesReads.begin(),
                                    tiesReads.end(), 0.0,
                                    [](double init, const seqInfo & seq) {return init + seq.cnt_;});

      setUp.rLog_.logCurrentTime("Setting Aligner for re-calc ...again final" );
      //alignPoolMap.destoryAligners();
      aligner alignerObjReCalc(maxReadLength, setUp.pars_.gapInfoRef_,
                               setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                               setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                               setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
      alignerObjReCalc.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
      if(pars.recalcConsensus){
        setUp.rLog_.logCurrentTime("Re-calculating consensus ...again final");
        if(setUp.pars_.verbose_){
          std::cout << "Re-calculating Consensus Sequences ...again final" << std::endl;
        }
      }else{
        setUp.rLog_.logCurrentTime("Setting Map Info ...again final");
      }

      std::vector<uint32_t> refReadPositions(readCounts.size());
      njh::iota<uint32_t>(refReadPositions, 0);
      njh::concurrent::LockableQueue<uint32_t> refReadPosQueue(refReadPositions);
      concurrent::AlignerPool alignPoolReCacl(alignerObjReCalc, pars.numThreads);
      alignPoolReCacl.initAligners();
      alignPoolReCacl.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
      auto setMappingInfo = [&refReadPosQueue,&consensusReads,
          &alignPoolReCacl,&readCounts,&pars,
          &clustersizeCutOff](){
        uint32_t readPos = 0;
        auto curThreadAligner = alignPoolReCacl.popAligner();
        while(refReadPosQueue.getVal(readPos)){

          auto foundReadPos = readVec::getReadIndexByName(consensusReads,
                                                          readCounts[readPos].seqBase_.name_);
          if (foundReadPos != std::string::npos) {
            if(readCounts[readPos].reads_.empty()){
              //if no reads were mapped to reference, remove it
              consensusReads[foundReadPos].remove = true;
              continue;
            }
            //update the the counts and fractions for the orignal clusters
            consensusReads[foundReadPos].seqBase_.cnt_ = readCounts[readPos].seqBase_.cnt_;
            consensusReads[foundReadPos].seqBase_.frac_ = readCounts[readPos].seqBase_.frac_;
            //add the reads from the refReads
            consensusReads[foundReadPos].reads_.clear();
            for(const auto & seq : readCounts[readPos].reads_){
              consensusReads[foundReadPos].reads_.emplace_back(std::make_shared<readObject>(seq));
            }
            consensusReads[foundReadPos].updateName();
            //re-calculate the consensus if indicated
            if(pars.recalcConsensus){
              readVecSorter::sortByTotalCountAE<seqInfo>(readCounts[readPos].reads_, true);
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
              auto recalcConSeqInfo = getConsensusWithAReCalc(readCounts[readPos].reads_, *curThreadAligner, consensusReads[foundReadPos].seqBase_.name_);
//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
              consensusReads[foundReadPos].seqBase_.seq_ = recalcConSeqInfo.seqBase_.seq_;
              consensusReads[foundReadPos].seqBase_.qual_ = recalcConSeqInfo.seqBase_.qual_;
              if(consensusReads[foundReadPos].seqBase_.cnt_ < clustersizeCutOff){
                consensusReads[foundReadPos].remove = true;
              }
            }
          }else{
            std::stringstream ss;
            ss << njh::bashCT::red << "Error, couldn't find: " << readCounts[readPos].seqBase_.name_ << njh::bashCT::reset << "\n";
            throw std::runtime_error{ss.str()};
          }
        }
      };

      std::vector<std::thread> mapInfoThreads;
      for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
        mapInfoThreads.emplace_back(std::thread(setMappingInfo));
      }

      for(auto & t : mapInfoThreads){
        t.join();
      }

      /*
      for (const auto & refRead : readCounts) {
        const auto foundReadPos = readVec::getReadIndexByName(consensusReads,
            refRead.seqBase_.name_);
        if (std::string::npos != foundReadPos) {
          if(refRead.reads_.empty()){
            //if no reads were mapped to reference, remove it
            consensusReads[foundReadPos].remove = true;
            continue;
          }
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          //update the the counts and fractions for the orignal clusters
          consensusReads[foundReadPos].seqBase_.cnt_ = refRead.seqBase_.cnt_;
          consensusReads[foundReadPos].seqBase_.frac_ = refRead.seqBase_.frac_;
          //add the reads from the refReads
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          consensusReads[foundReadPos].reads_.clear();
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          for(const auto & seq : refRead.reads_){
            consensusReads[foundReadPos].reads_.emplace_back(std::make_shared<readObject>(seq));
          }
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          consensusReads[foundReadPos].updateName();
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
          //re-calculate the consensus if indicated
          if(pars.recalcConsensus){
            std::cout << __FILE__ << "  " << __LINE__ << std::endl;
            std::cout << "refRead.reads_.size() : " << refRead.reads_.size() << std::endl;
            auto recalcConSeqInfo = getConsensusWithAReCalc(refRead.reads_, alignerObj, consensusReads[foundReadPos].seqBase_.name_);
            std::cout << __FILE__ << "  " << __LINE__ << std::endl;
            consensusReads[foundReadPos].seqBase_.seq_ = recalcConSeqInfo.seqBase_.seq_;
            consensusReads[foundReadPos].seqBase_.qual_ = recalcConSeqInfo.seqBase_.qual_;
          }
          std::cout << __FILE__ << "  " << __LINE__ << std::endl;
        }else{
          std::stringstream ss;
          ss << njh::bashCT::red << "Error, couldn't find: " << refRead.seqBase_.name_ << njh::bashCT::reset << "\n";
          throw std::runtime_error{ss.str()};
        }
      }*/



      //get rid of the empty clusters
      consensusReads = readVecSplitter::splitVectorOnRemove(consensusReads).first;
      if(pars.recalcConsensus){
        if (setUp.pars_.verbose_) {
          std::cout << njh::bashCT::bold << "Collapsing after recalculating consensus sequences ...again final" << njh::bashCT::reset
                    << std::endl;
        }
        setUp.rLog_.logCurrentTime("Collapsing after recalculating consensus sequences ...again final");
        //collapse clusters again if the consensus had to be calculated
        simpleCollapseWithPars(consensusReads, alignerObj, postCollapsePars);
      }

      if (setUp.pars_.verbose_) {
        std::cout << njh::bashCT::bold << "Finished with "
                  << consensusReads.size() << " Consensus Sequences after Breakout remap collapsed and final remap"
                  << njh::bashCT::reset << std::endl;
      }
      setUp.rLog_ << "Finished with " << consensusReads.size()
                  << " Consensus Sequences after Breakout, final" << "\n";
      std::ofstream mapInfo;
      openTextFile(mapInfo, mappingDir + "mapInfo", ".tab.txt", setUp.pars_.ioOptions_.out_);
      std::set<std::string> repNames;
      for (const auto & read : allInputReads) {
        if("" != pars.aSetRepName){
          repNames.emplace(pars.aSetRepName);
        }else{
          auto toks = tokenizeString(read->seqBase_.name_, ".");
          repNames.emplace(njh::replaceString(toks[0], "_Comp", ""));
        }

      }
      mapInfo << "seqName\trefId\tftotalReadsMapped\ttotalMappedFraction";
      for (const auto & n : repNames) {
        mapInfo << "\t" << n << "_readsMapped\t" << n << "_mappedFraction";
      }
      if (setUp.pars_.refIoOptions_.firstName_ != "") {
        mapInfo << "\tBestRef\tscore\t1bIndel\t2bI"
                   "ndel\t>2bIndel\tlqMismatch\thqMismatch";
      }
      mapInfo << "\n";
      std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
      std::map<std::string, double> repTotals;
      for (const auto & ref : readCounts) {
        for (const auto & read : ref.reads_) {
          if ("" != pars.aSetRepName) {
            repTotals[pars.aSetRepName] += read.cnt_;
          } else {
            auto toks = tokenizeString(read.name_, ".");
            repTotals[njh::replaceString(toks[0], "_Comp", "")] +=
                read.cnt_;
          }
        }
      }

      for (auto & ref : readCounts) {
        std::map<std::string, double> repCounts;
        for (const auto & read : ref.reads_) {
          auto toks = tokenizeString(read.name_, ".");
          repCounts[njh::replaceString(toks[0], "_Comp", "")] += read.cnt_;
        }
        mapInfo << seqName << "\t" << ref.seqBase_.name_ << "\t"
                << ref.seqBase_.cnt_ << "\t" << ref.seqBase_.cnt_ / totalReadCnt;
        for (const auto & n : repNames) {
          if (repCounts.find(n) != repCounts.end()) {
            mapInfo << "\t" << repCounts[n] << "\t"
                    << repCounts[n] / repTotals[n];
          } else {
            mapInfo << "\t0\t0";
          }
        }
        if (setUp.pars_.refIoOptions_.firstName_ != "") {
          bool eventBased = true;
          mapInfo << "\t"
                  << profiler::compareToRefSingle(refSeqs, ref, alignerObj,
                                                  setUp.pars_.local_, eventBased).front();
        }
        mapInfo << "\n";
      }

      mapInfo << seqName << "\t" << "unmappable\t" << unmappableAmount << "\t"
              << unmappableAmount / totalReadCnt << "\n";
      mapInfo << seqName << "\t" << "indeterminate\t" << indeterminateAmount
              << "\t" << indeterminateAmount / totalReadCnt << "\n";
      mapInfo << seqName << "\t" << "ties\t" << tiesAmount
              << "\t" << tiesAmount / totalReadCnt << "\n";
    }
    renameReadNames(consensusReads,
                    bfs::basename(setUp.pars_.ioOptions_.firstName_), true, true,
                    true, "totalCount");
  }
  renameReadNames(consensusReads,
                  bfs::basename(setUp.pars_.ioOptions_.firstName_), true, true,
                  true, "totalCount");


  if(pars.visualize){
    //output visualization where clusters are written and colored either by a comparison to a reference sequence
    //or to the best matching consensus sequence
    /** @todo change the color to the read it gets mapped to if remapping*/
    setUp.rLog_.logCurrentTime("Visualizing");
    std::string visDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("vis")).string();
    std::unordered_map<std::string, std::string> nameToColor;
    std::mutex nameToColorMut;
    std::vector<uint32_t> positions(allInputReads.size());
    njh::iota<uint32_t>(positions,0);
    njh::concurrent::LockableQueue<uint32_t> positionsQueue(positions);

    //create aligner with gap info from ref so a different gap scoring can be used for mapping determination
    aligner alignerObjMapper(maxReadLength, setUp.pars_.gapInfoRef_,
                             setUp.pars_.scoring_, KmerMaps(), setUp.pars_.qScorePars_,
                             setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
                             setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
    alignerObjMapper.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
    concurrent::AlignerPool alignPool(alignerObjMapper, pars.numThreads);
    alignPool.initAligners();
    alignPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;


    table colorToName;
    std::string backgroundColor = "#000000";
    VecStr refNames;
    if(!refSeqs.empty()){
      refNames = readVec::getNames(refSeqs);
    }else{
      refNames = readVec::getNames(consensusReads);
    }
    auto popColors = getColorsForNames(refNames);
    if(colorFile != ""){
      VecStr adding;
      colorToName = table(colorFile,"\t", true);
      if(colorToName.columnNames_.size() != 2){
        std::stringstream ss;
        ss << "Wrong number of columns in " << colorFile << ", , should only have two, "
           << colorFile << " has " << colorToName.columnNames_.size()
           << std::endl;
        throw std::runtime_error{ss.str()};
      }
      for(const auto & row : colorToName.content_){
        if(stringToLowerReturn(row[0]) == "backgroundcolor"){
          backgroundColor = row[1];
          continue;
        }
        adding.emplace_back(row[0]);
        if(popColors.find(row[0]) == popColors.end()){
          std::cerr << "Warning, adding " << row[0] << " but it wasn't found in the reads" << std::endl;
        }
        popColors[row[0]] = njh::color(row[1]);
      }
      VecStr didNotAdd;
      for(const auto & n : refNames){
        if(!njh::in(n, adding)){
          didNotAdd.emplace_back(n);
          std::cerr << "Warning, add a file to specify colors for reads but didn't have a color set for " << n << std::endl;
        }
        auto popColorsNoColor = getColorsForNames(didNotAdd);
        for(const auto & color : popColorsNoColor){
          popColors[color.first] = color.second;
        }
      }
    }
    if(setUp.pars_.verbose_){
      std::cout << "SeqName\tcolor" << std::endl;
      for(const auto & popColor : popColors){
        std::cout << popColor.first <<"\t" << popColor.second.getHexStr() << std::endl;
      }
    }
    auto compareReads = [&alignPool, &refSeqs,&nameToColorMut,&nameToColor,&allInputReads,&consensusReads,&positionsQueue,&popColors](){
      uint32_t pos = std::numeric_limits<uint32_t>::max();
      std::unordered_map<std::string, std::string> nameToColorCurrent;
      auto currentAligner = alignPool.popAligner();
      while(positionsQueue.getVal(pos)){
        const auto & seq = allInputReads[pos];
        if (!refSeqs.empty()) {
          auto bestRef = profiler::compareToRefSingle(refSeqs, seq, *currentAligner,
                                                      false, false);
          auto bestRefToks = tokenizeString(bestRef.front(), ",");
          nameToColorCurrent[seq->seqBase_.name_] = popColors[bestRefToks[0]].getHexStr();
        } else {
          auto bestRef = profiler::compareToRefSingle(consensusReads, seq,
                                                      *currentAligner, false, false);
          auto bestRefToks = tokenizeString(bestRef.front(), ",");
          nameToColorCurrent[seq->seqBase_.name_] = popColors[bestRefToks[0]].getHexStr();
        }
      }

      {
        std::lock_guard<std::mutex> lock(nameToColorMut);
        for(const auto & name : nameToColorCurrent){
          nameToColor[name.first] = name.second;
        }
      }
    };

    std::vector<std::thread> threads;
    for(uint32_t threadNum = 0; threadNum < pars.numThreads; ++threadNum){
      threads.emplace_back(std::thread(compareReads));
    }
    for(auto & t : threads){
      t.join();
    }
    uint32_t subGroupingCount = 0;
    for(const auto & disGraph : disGraphs){
      auto subGroupVisDir = njh::files::make_path(visDir, njh::pasteAsStr("subgroup-", subGroupingCount));
      njh::files::makeDir(njh::files::MkdirPar{subGroupVisDir});
      ++subGroupingCount;
      auto treeJson = disGraph->toJson(clustersizeCutOff, nameToColor);
      auto & nodes = treeJson["nodes"];
      for(auto & node : nodes){
        node["size"] = disGraph->nodes_[disGraph->nameToNodePos_[node["name"].asString()]]->value_->cnt_ * pars.nodeSize;
      }
      std::ofstream treeJsonFile(visDir + "tree.json");
      treeJsonFile << treeJson;

      std::ofstream treeHtmlFile(visDir + "tree.html");
      genTreeHtml(treeHtmlFile, "tree.json", "tree.js");

      std::ofstream treeJsFile(visDir + "tree.js");
      genSimpleTreeJs(treeJsFile);
    }
  }

  std::string snpDir = njh::files::makeDir(setUp.pars_.directoryName_,
                                           njh::files::MkdirPar("internalSnpInfo")).string();
  //output info file on consensus
  std::ofstream infoFile;
  openTextFile(infoFile, setUp.pars_.directoryName_ + "outputInfo",
               ".tab.txt", false, false);
  readVec::allSetFractionByTotalCount(consensusReads);
  infoFile << "ClusterNumber\tClusterId\tClusterSize\tClusterFraction\n";
  for (const auto readPos : iter::range(consensusReads.size())) {
    //write out cluster info
    infoFile << readPos
             << '\t' << consensusReads[readPos].seqBase_.name_
             << '\t' << consensusReads[readPos].seqBase_.cnt_
             << '\t' << consensusReads[readPos].seqBase_.frac_ << '\n';
    SeqOutput::write(consensusReads[readPos].reads_,
                     SeqIOOptions(clusDir.string() + consensusReads[readPos].seqBase_.name_,
                                  setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
    //log snp information
    /**@todo add indel information, less informative since pacbio have a lot*/
    std::unordered_map<uint32_t, std::unordered_map<char, VecStr>> mismatches;
    for (const auto subReadPos : iter::range(
        consensusReads[readPos].reads_.size())) {
      const auto & subRead = consensusReads[readPos].reads_[subReadPos];
      alignerObj.alignCache(consensusReads[readPos], subRead, false);
      //count gaps and mismatches and get identity
      alignerObj.profilePrimerAlignment(consensusReads[readPos], subRead);
      for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
        if(m.second.highQuality(setUp.pars_.qScorePars_)){
          mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(
              subRead->seqBase_.name_);
        }
      }
    }
    table misTab { VecStr { "refPos", "refBase", "seqBase", "freq", "fraction",
                            "seqs" } };
    for (const auto & m : mismatches) {
      for (const auto & seqM : m.second) {

        misTab.content_.emplace_back(
            toVecStr(m.first, consensusReads[readPos].seqBase_.seq_[m.first],
                     seqM.first, seqM.second.size(),
                     seqM.second.size() / consensusReads[readPos].seqBase_.cnt_,
                     vectorToString(seqM.second, ",")));
      }
    }
    misTab.sortTable("seqBase", false);
    misTab.sortTable("refPos", false);
    misTab.outPutContents(
        TableIOOpts(
            OutOptions(snpDir + consensusReads[readPos].seqBase_.name_,
                       ".tab.txt"), "\t", misTab.hasHeader_));
  }


  if (setUp.pars_.chiOpts_.checkChimeras_) {
    setUp.rLog_.logCurrentTime("Checking Chimeras");

    setUp.pars_.chiOpts_.chiOverlap_.oneBaseIndel_ = 0.99;
    setUp.pars_.chiOpts_.chiOverlap_.twoBaseIndel_ = 0.99;
    setUp.pars_.chiOpts_.chiOverlap_.lqMismatches_ = pars.chiAllowableError;
    setUp.pars_.chiOpts_.chiOverlap_.hqMismatches_ = pars.chiAllowableError;
    collapser collapserObj = collapser(setUp.pars_.colOpts_);
    collapserObj.opts_.verboseOpts_.debug_ = setUp.pars_.debug_;
    auto chiInfoTab = collapserObj.markChimeras(consensusReads, alignerObj,
                                                setUp.pars_.chiOpts_);
    chiInfoTab.outPutContents(
        TableIOOpts(
            OutOptions(setUp.pars_.directoryName_ + "chiParentsInfo.txt",
                       ".txt"), "\t", true));
  }

  //write final consensus reads
  SeqOutput::write(consensusReads,
                   SeqIOOptions(setUp.pars_.directoryName_ + "output",
                                setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
  //write out which reads went into which cluster
  std::ofstream clusterNamesFile;
  openTextFile(clusterNamesFile, setUp.pars_.directoryName_ + "clusterNames",
               ".tab.txt", false, false);
  clusterNamesFile << "clusterName\treadName\n";
  for (const auto & con : consensusReads) {
    for (const auto & read : con.reads_) {
      clusterNamesFile << con.seqBase_.name_ << "\t" << read->seqBase_.name_
                       << std::endl;
    }
  }

  //check if writing out additional location file
  if (pars.additionalOutLocationFile != "") {
    auto fnp = setUp.pars_.ioOptions_.firstName_.filename().string();
    if(njh::endsWith(fnp, ".gz")){
      fnp = fnp.substr(0, fnp.rfind(".gz"));
    }
    std::string additionalOutDir = findAdditonalOutLocation(
        pars.additionalOutLocationFile, fnp);
    if (additionalOutDir == "") {
      std::cerr << njh::bashCT::red << njh::bashCT::bold;
      std::cerr << "No additional out directory found for: "
                << setUp.pars_.ioOptions_.firstName_ << std::endl;
      std::cerr << njh::bashCT::reset;
      std::cerr << processFileNameForID(setUp.pars_.ioOptions_.firstName_.string())<< std::endl;
      std::cerr << "Options:" << std::endl;
      table inTab(pars.additionalOutLocationFile, "\t");
      MapStrStr additionalOutNames;
      for (const auto& fIter : inTab.content_) {
        additionalOutNames[makeIDNameComparable(fIter[0])] = fIter[1];
      }
      std::cerr << njh::conToStr(njh::getVecOfMapKeys(additionalOutNames)) << std::endl;
    } else {
      SeqOutput::write(consensusReads, SeqIOOptions(additionalOutDir + setUp.pars_.ioOptions_.out_.outFilename_.string(),
                                                    setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
      std::ofstream metaDataFile;
      openTextFile(metaDataFile, additionalOutDir + "/" + "metaData", ".json",
                   setUp.pars_.ioOptions_.out_);
      metaDataFile << metaData;
    }
  }



  return 0;
}


}  // namespace njh

