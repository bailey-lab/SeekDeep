//
// Created by Nicholas Hathaway on 3/6/23.
//
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"

#include <njhseq/objects/kmer/SimpleKmerHash.hpp>

namespace njhseq {

int SeekDeepRunner::extractorByKmerMatching(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path idFnp;
  bfs::path lenCutOffsFnp;
  bfs::path uniqueKmersPerTargetFnp;
  bfs::path refSeqDir;
  CoreExtractorPars corePars;
  corePars.smallFragmentCutoff = 100;

  uint32_t numThreads = 1;

  std::string sampleName = "sample";
  bool rename = false;

  SeekDeepSetUp setUp(inputCommands);
  //id
  setUp.setOption(idFnp, "--ids,--id", "Primers file", true, "IDs");

  setUp.setOption(lenCutOffsFnp, "--lenCutOffsFnp,--lenCutOffs", "length Cut Offs per target", true, "IDs");
  setUp.setOption(uniqueKmersPerTargetFnp, "--uniqueKmersPerTarget", "unique Kmers Per Target", true, "IDs");
  setUp.setOption(refSeqDir, "--refSeqDir", "reference Seq Dir, should have a .fasta file named by each target", false, "IDs");

  //in and out
  setUp.processDefaultReader(VecStr{"--fastq", "--fastqgz", "--fasta", "--fastagz"});
  setUp.processDirectoryOutputName(true);
  if("sample" == sampleName){
    sampleName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
    if(njh::endsWith(setUp.pars_.ioOptions_.firstName_.string(), ".gz")){
      sampleName = bfs::path(sampleName).replace_extension("").string();
    }
  }
  setUp.setOption(sampleName, "--sampleName", "sample name");
  setUp.setOption(rename, "--rename", "rename input extracted reads");

  // filtering
  setUp.setOption(corePars.smallFragmentCutoff, "--minLenCutOff", "Hard cut off min length", false, "filtering");

  //running
  setUp.setOption(numThreads, "--numThreads", "number of threads");


  //primers
  corePars.pDetPars.useMotif_ = false;
  corePars.pDetPars.allowable_.lqMismatches_ = 5;
  corePars.pDetPars.allowable_.hqMismatches_ = 5;
  corePars.pDetPars.allowable_.distances_.query_.coverage_ = .50;
  corePars.pDetPars.allowable_.largeBaseIndel_ = .99;
  corePars.pDetPars.allowable_.oneBaseIndel_ = 4;
  corePars.pDetPars.allowable_.twoBaseIndel_ = 2;
  corePars.pDetPars.primerWithin_ = 50;
  bool primerToUpperCase = false;
  setUp.setOption(primerToUpperCase, "--primerUpper",
                  "Leave primers in upper case", false, "Primer");
  corePars.pDetPars.primerToLowerCase_ = !primerToUpperCase;
  setUp.setOption(corePars.pDetPars.allowable_.distances_.query_.coverage_, "--primerCoverage",
                  "Amount of primers found", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.hqMismatches_, "--primerNumOfMismatches",
                  "Number of Mismatches to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.oneBaseIndel_, "--primerOneBaseIndels",
                  "Number Of One base indels to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.twoBaseIndel_, "--primerTwoBaseIndels",
                  "Number Of Two base indels to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.largeBaseIndel_, "--primerMoreThan2BaseIndels",
                  "Number Of greater than two base indels to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.primerWithin_, "--primerWithinStart",
                  "By default the primers are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives",
                  false, "Primers");
  setUp.setOption(corePars.pDetPars.primerStart_, "--primerSearchStart",
                  "By default the primers are searched at the very beginning of seq, use this flag to start the search here",
                  false, "Primers");
  //copy over the regular determined primers to the back end primers
  corePars.backEndpDetPars.primerWithin_ = corePars.pDetPars.primerWithin_;
  corePars.backEndpDetPars.primerToLowerCase_ = corePars.pDetPars.primerToLowerCase_ ;

  corePars.backEndpDetPars.allowable_.hqMismatches_ = corePars.pDetPars.allowable_.hqMismatches_;
  corePars.backEndpDetPars.allowable_.lqMismatches_ = corePars.pDetPars.allowable_.lqMismatches_;
  corePars.backEndpDetPars.allowable_.distances_.query_.coverage_ = corePars.pDetPars.allowable_.distances_.query_.coverage_;
  corePars.backEndpDetPars.allowable_.largeBaseIndel_ = corePars.pDetPars.allowable_.largeBaseIndel_;
  corePars.backEndpDetPars.allowable_.oneBaseIndel_ = corePars.pDetPars.allowable_.oneBaseIndel_;
  corePars.backEndpDetPars.allowable_.twoBaseIndel_ = corePars.pDetPars.allowable_.twoBaseIndel_;
  corePars.backEndpDetPars.useMotif_ = corePars.pDetPars.useMotif_;

  corePars.qPars_.qualCheck_ = 30;
  corePars.qPars_.qualCheckCutOff_ = 0.15;
  setUp.setOption(corePars.qPars_.qualCheck_, "--qualCheckLevel",
                  "Bin qualities at this quality to do filtering on fraction above this", false, "Post Processing");
  setUp.setOption(corePars.qPars_.qualCheckCutOff_, "--qualCheckCutOff",
                  "The fractions of bases that have to be above the qualCheckLevel to be kept", false, "Post Processing");

  setUp.finishSetUp(std::cout);
  // run log
  setUp.startARunLog(setUp.pars_.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

  // create Primers and MIDs
  PrimersAndMids ids(idFnp);

  if(ids.getTargets().empty()){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "no targets read in from " << idFnp << "\n";
    throw std::runtime_error { ss.str() };
  }

  //init primer determinator
  ids.initPrimerDeterminator();

  //read in extra info
  ids.addLenCutOffs(lenCutOffsFnp);
  if(!refSeqDir.empty()){
    ids.addRefSeqs(refSeqDir);
  }
  uint32_t extractionKmer = ids.addUniqKmerCounts(uniqueKmersPerTargetFnp);


  // set up input
  seqInfo seq;
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();



  //qual check off;
  std::unique_ptr<ReadChecker> qualChecker = std::make_unique<ReadCheckerQualCheck>(corePars.qPars_.qualCheck_,
                                                                                    corePars.qPars_.qualCheckCutOff_, true);

  //set up extraction counts outputs
  ExtractionStator masterCounts;
  OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
  OutputStream outCounts(outCountsOpts);

  OutOptions outStatsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
  OutputStream outStats(outStatsOpts);

  std::unordered_map<std::string, uint32_t> readsPerSet;
  std::unordered_map<std::string, uint32_t> readsPerSetRevComp;

  VecStr names = ids.getTargets();
  MultiSeqIO seqOut;
  auto initialExtractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"extractedReads"});
  auto failedExtractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"failedExtractedReads"});
  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(initialExtractionDir, name) );
    seqOut.addReader(name, seqOutOpts);
  }

  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(name, sampleName) ) );
    seqOut.addReader(njh::pasteAsStr(name, "-passed"), seqOutOpts);
  }

  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(failedExtractionDir, name) );
    seqOut.addReader(njh::pasteAsStr(name, "-failed"), seqOutOpts);
  }

  seqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));

  //set up aligner
  // creating aligner
  // create aligner for primer identification
  auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
      setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);

  setUp.pars_.gapInfo_.gapLeftQueryOpen_ = 0;
  setUp.pars_.gapInfo_.gapLeftQueryExtend_ = 0;

  setUp.pars_.gapInfo_.gapRightQueryOpen_ = 0;
  setUp.pars_.gapInfo_.gapRightQueryExtend_ = 0;

  setUp.pars_.gapInfo_.gapLeftRefOpen_ = 0;
  setUp.pars_.gapInfo_.gapLeftRefExtend_ = 0;

  setUp.pars_.gapInfo_.gapRightRefOpen_ = 0;
  setUp.pars_.gapInfo_.gapRightRefExtend_ = 0;

  gapScoringParameters gapPars(setUp.pars_.gapInfo_);
  KmerMaps emptyMaps;
  bool countEndGaps = false;
  uint64_t maxPrimerSize = ids.pDeterminator_->getMaxPrimerSize();
  uint64_t maxReadSize =  maxPrimerSize * 4 + corePars.pDetPars.primerWithin_;
  aligner alignObj(maxReadSize, gapPars, scoreMatrix, emptyMaps, setUp.pars_.qScorePars_, countEndGaps, false);
  alignObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

  concurrent::AlignerPool alnPool(alignObj, numThreads);
  alnPool.initAligners();

  std::mutex mut;
  std::function<void()> readInComp = [&reader, &ids, &readsPerSet,&readsPerSetRevComp,&mut,
                                      &extractionKmer,&seqOut,&sampleName,&rename,
                                      &corePars,&alnPool,
                                      &qualChecker,
                                      &masterCounts]() {
    SimpleKmerHash hasher;
    seqInfo seq;
    std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;


    ExtractionStator masterCountsCurrent;

    auto currentAlignerObj = alnPool.popAligner();
    while(reader.readNextReadLock(seq)){
      ++masterCountsCurrent.totalReadCount_;
      if(len(seq) < corePars.smallFragmentCutoff){
        ++masterCountsCurrent.smallFrags_;
        continue;
      }
      std::unordered_map<uint64_t, uint64_t> hashedInputKmers;
      std::unordered_map<uint64_t, uint64_t> hashedInputKmersRevComp;
      if(len(seq.seq_) > extractionKmer){
        for(uint32_t pos = 0; pos < len(seq.seq_) - extractionKmer + 1; ++pos){
          auto hash = hasher.hash(seq.seq_.substr(pos, extractionKmer));
          ++hashedInputKmers[hash];
        }
      }
      if(len(seq.seq_) > extractionKmer){
        for(uint32_t pos = 0; pos < len(seq.seq_) - extractionKmer + 1; ++pos){
          auto hash = hasher.revCompHash(seq.seq_.substr(pos, extractionKmer));
          ++hashedInputKmersRevComp[hash];
        }
      }

      std::unordered_map<std::string, uint32_t> foundPerSet;
      std::unordered_map<std::string, uint32_t> foundPerSetRevComp;
      for(const auto & setName  : njh::getVecOfMapKeys(ids.uniqueKmersPerTarget_)){
        foundPerSet[setName] = 0;
        foundPerSetRevComp[setName] = 0;
      }
      for(const auto & hashedKmer : hashedInputKmers){
        for(const auto & uniqueKmers : ids.uniqueKmersPerTarget_){
          if(njh::in(hashedKmer.first, uniqueKmers.second)){
            ++foundPerSet[uniqueKmers.first];
          }
        }
      }
      for(const auto & hashedKmer : hashedInputKmersRevComp){
        for(const auto & uniqueKmers : ids.uniqueKmersPerTarget_){
          if(njh::in(hashedKmer.first, uniqueKmers.second)){
            ++foundPerSetRevComp[uniqueKmers.first];
          }
        }
      }
      std::string winnerSet = "undetermined";
      double bestFrac = 0;
      bool winnerRevComp = false;

      for(const auto & setName  : njh::getVecOfMapKeys(ids.uniqueKmersPerTarget_)){
        if(static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
          bestFrac = static_cast<double>(foundPerSet[setName])/static_cast<double>(hashedInputKmers.size());
          winnerSet = setName;
        }
        if(static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size()) > bestFrac){
          bestFrac = static_cast<double>(foundPerSetRevComp[setName])/static_cast<double>(hashedInputKmers.size());
          winnerSet = setName;
          winnerRevComp = true;
        }
      }
      if(winnerRevComp){
        ++readsPerSetRevCompCurrent[winnerSet];
        seq.reverseComplementRead(true, true);
      }else{
        ++readsPerSetCurrent[winnerSet];
      }
      if(rename){
        auto threadId = estd::to_string(std::this_thread::get_id());
        seq.name_ = njh::pasteAsStr(sampleName, ".",winnerSet, ".", readsPerSetCurrent[winnerSet] + readsPerSetRevCompCurrent[winnerSet], ".", threadId);
        if(winnerRevComp){
          seq.name_.append("_Comp");
        }
      }
      seqOut.openWrite(winnerSet, seq);
      if("undetermined" == winnerSet){
        ++masterCountsCurrent.readsUnrecBarcode_;
        continue;
      }
      auto extractorCase = ExtractionStator::extractCase::GOOD;
      //min len
      ids.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(seq);
      if(!seq.on_){
        extractorCase = ExtractionStator::extractCase::MINLENBAD;
      }
      if(seq.on_){

        std::string frontPrimerName = ids.pDeterminator_->determineForwardPrimer(seq, corePars.pDetPars, *currentAlignerObj, VecStr{winnerSet});
//        if("unrecognized" == frontPrimerName){
//          currentAlignerObj->alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//          currentAlignerObj->alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//          std::cout << njh::bashCT::cyan ;
//          std::cout << "currentAlignerObj->comp_.distances_.query_.coverage_: " << currentAlignerObj->comp_.distances_.query_.coverage_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.oneBaseIndel_:      " << currentAlignerObj->comp_.oneBaseIndel_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.twoBaseIndel_:      " << currentAlignerObj->comp_.twoBaseIndel_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.largeBaseIndel_:    " << currentAlignerObj->comp_.largeBaseIndel_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.hqMismatches_:      " << currentAlignerObj->comp_.hqMismatches_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.lqMismatches_:      " << currentAlignerObj->comp_.lqMismatches_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.lowKmerMismatches_: " << currentAlignerObj->comp_.lowKmerMismatches_ << std::endl;
//          std::cout << njh::bashCT::reset ;
//        } else {
//          currentAlignerObj->alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//          currentAlignerObj->alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//          std::cout << njh::bashCT::green ;
//          std::cout << "currentAlignerObj->comp_.distances_.query_.coverage_: " << currentAlignerObj->comp_.distances_.query_.coverage_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.oneBaseIndel_:      " << currentAlignerObj->comp_.oneBaseIndel_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.twoBaseIndel_:      " << currentAlignerObj->comp_.twoBaseIndel_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.largeBaseIndel_:    " << currentAlignerObj->comp_.largeBaseIndel_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.hqMismatches_:      " << currentAlignerObj->comp_.hqMismatches_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.lqMismatches_:      " << currentAlignerObj->comp_.lqMismatches_ << std::endl;
//          std::cout << "currentAlignerObj->comp_.lowKmerMismatches_: " << currentAlignerObj->comp_.lowKmerMismatches_ << std::endl;
//          std::cout << njh::bashCT::reset ;
//        }
        seq.reverseComplementRead(false, true);
        std::string backPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq, corePars.backEndpDetPars, *currentAlignerObj, VecStr{winnerSet});
        seq.reverseComplementRead(false, true);
        bool passesFront = frontPrimerName == winnerSet;
        bool passesBack = backPrimerName == winnerSet;

        if (!passesBack && !passesFront) {
					extractorCase = ExtractionStator::extractCase::FAILEDBOTHPRIMERS;
        } else if(!passesFront){
					extractorCase = ExtractionStator::extractCase::BADFORWARD;
        } else if(!passesBack){
					extractorCase = ExtractionStator::extractCase::BADREVERSE;
        }
        if (passesBack && passesFront) {
          //min len
          ids.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(seq);
          if(!seq.on_){
            extractorCase = ExtractionStator::extractCase::MINLENBAD;
          }

          //max len
          if(seq.on_){
            ids.targets_.at(winnerSet).lenCuts_->maxLenChecker_.checkRead(seq);
            if(!seq.on_){
              extractorCase = ExtractionStator::extractCase::MAXLENBAD;
            }
          }

          //quality
          if(seq.on_){
            qualChecker->checkRead(seq);
            if(!seq.on_){
              extractorCase = ExtractionStator::extractCase::QUALITYFAILED;
            }
          }
          if (seq.on_) {
            //passed all filters
            seqOut.openWrite(njh::pasteAsStr(winnerSet, "-passed"), seq);
          } else {
            //failed
            seqOut.openWrite(njh::pasteAsStr(winnerSet, "-failed"), seq);
          }
        } else {
          seq.name_.append(njh::pasteAsStr("[", "frontPrimerName=", frontPrimerName, ";", "backPrimerName=", backPrimerName, ";","]"));
          //failed primer check
          seqOut.openWrite(njh::pasteAsStr(winnerSet, "-failed"), seq);
        }
      }else{
        //failed first min length pass
        seqOut.openWrite(njh::pasteAsStr(winnerSet, "-failed"), seq);
      }
      masterCountsCurrent.increaseCounts(winnerSet, seq.name_, extractorCase);
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);

      masterCounts.addOtherExtractorCounts(masterCountsCurrent);

      for(const auto & readsPerSetCount : readsPerSetCurrent){
        readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
      }
      for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
        readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
      }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);


  uint64_t totalReadsProcessed = 0;
  for(const auto & readsPerSetCount : readsPerSet){
    totalReadsProcessed += readsPerSetCount.second;
  }
  for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
    totalReadsProcessed += readsPerSetRevCompCount.second;
  }

  uint64_t totalReadsPassedAllFilters = 0;


  outCounts << "sample\ttotalReadsProcessed\ttarget\tcount\tfrac\tforwardCount\tfracForward";
  outCounts << "\tminLenFailed\tminLenFailedFrac";
  outCounts << "\tmaxLenFailed\tmaxLenFailedFrac";
  outCounts << "\tqualityFailed\tqualityFailedFrac";
  outCounts << "\tforwardPrimerFailed\tforwardPrimerFailedFrac";
  outCounts << "\treversePrimerFailed\treversePrimerFailedFrac";
  outCounts << "\tbothForRevPrimerFailed\tbothForRevPrimerFailedFrac";
  outCounts << "\ttotalFailed\ttotalFailedFrac";
  outCounts << "\tpasssed\tpassedFrac";
  outCounts << std::endl;


  uint64_t totalExtractedAllTargets = 0;
  uint64_t totalExtractedAllTargetsForward = 0;
  uint64_t totalExtractedUndetermined = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];


  for(const auto & setName : ids.getTargets()){
    uint64_t totalExtracted = readsPerSet[setName] + readsPerSetRevComp[setName];
    totalExtractedAllTargets += totalExtracted;
		totalExtractedAllTargetsForward += readsPerSet[setName];

    uint64_t totalBad = masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_ +
        masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_+
        masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_ +
        masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_ +
        masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_ +
        masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_;

    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << setName
              << "\t" << totalExtracted
              << "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed)
              << "\t" << readsPerSet[setName]
              << "\t" << static_cast<double>(readsPerSet[setName]) / static_cast<double>(totalExtracted);
    //minlen
    outCounts<< "\t" << masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_
    << "\t" << static_cast<double>(masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_)/static_cast<double>(totalBad);
    //maxlen
    outCounts << "\t" << masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_
    << "\t" << static_cast<double>(masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_)/static_cast<double>(totalBad);
    //quality
    outCounts << "\t" << masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_
              << "\t" << static_cast<double>(masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_)/static_cast<double>(totalBad);
    //failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_
              << "\t" << static_cast<double>(masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_)/static_cast<double>(totalBad);
    //bad reverse
    outCounts << "\t" << masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_
              << "\t" << static_cast<double>(masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_)/static_cast<double>(totalBad);
    //both bad reverse and failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_
              << "\t" << static_cast<double>(masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_)/static_cast<double>(totalBad);
    //bad
    outCounts << "\t" << totalBad
              << "\t" << static_cast<double>(totalBad) / static_cast<double>(totalExtracted);
    //good
    outCounts << "\t" << masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_
              << "\t" << static_cast<double>(masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_) / static_cast<double>(totalExtracted)
              << std::endl;
    totalReadsPassedAllFilters += masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_;
  }
  {
    uint64_t totalExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << "undetermined"
              << "\t" << totalExtracted
              << "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed)
              << "\t" << readsPerSet["undetermined"]
              << "\t" << static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalExtracted)
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
              << std::endl;
  }

  outStats << "sampleName\ttotalReadsProcessed\tfailedMinLen_" << corePars.smallFragmentCutoff << "\tfailedMinLenFrac\tundetermined\tundeterminedFrac\textracted\textractedFrac\textractedForward\textractedForwardFrac\tpassed\tpassedFrac" << std::endl;
  outStats << sampleName
					 << "\t" << totalReadsProcessed + masterCounts.smallFrags_
					 << "\t" << masterCounts.smallFrags_
					 << "\t" << static_cast<double>(masterCounts.smallFrags_) / static_cast<double>(totalReadsProcessed + masterCounts.smallFrags_)
					 << "\t" << totalExtractedUndetermined
					 << "\t" << static_cast<double>(totalExtractedUndetermined) / static_cast<double>(totalExtractedUndetermined + totalExtractedAllTargets)
					 << "\t" << totalExtractedAllTargets
					 << "\t" << static_cast<double>(totalExtractedAllTargets) / static_cast<double>(totalReadsProcessed + masterCounts.smallFrags_)
					 << "\t" << totalExtractedAllTargetsForward
					 << "\t" << static_cast<double>(totalExtractedAllTargetsForward) / static_cast<double>(totalExtractedAllTargets)
           << "\t" << totalReadsPassedAllFilters
           << "\t" << static_cast<double>(totalReadsPassedAllFilters) / static_cast<double>(totalReadsProcessed + masterCounts.smallFrags_)
           << std::endl;

  return 0;
}

} // namespace njhseq
