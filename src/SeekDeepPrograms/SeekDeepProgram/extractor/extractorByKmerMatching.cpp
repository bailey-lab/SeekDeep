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
  corePars.primIdsPars.mPars_.allowableErrors_ = 1;
  setUp.setOption(corePars.primIdsPars.mPars_.allowableErrors_, "--MIDAllowableErrors", "Number of errors to allow in MIDs, more errors allowed can lead to mis-sorting of reads by MIDs", false, "MID");
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

  //init determinators
  ids.initPrimerDeterminator();

  if(ids.containsMids()) {
    ids.initMidDeterminator(corePars.primIdsPars.mPars_);
  }
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
    if(ids.containsMids()) {
      for(const auto & midname : getVectorOfMapKeys(ids.mids_)) {
        auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(name, midname) ) );
        seqOut.addReader(njh::pasteAsStr(name, midname, "-passed"), seqOutOpts);

        auto failedSeqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(failedExtractionDir, name + midname) );
        seqOut.addReader(njh::pasteAsStr(name, midname, "-failed"), failedSeqOutOpts);
      }
      seqOut.addReader(njh::pasteAsStr(name, "-MIDUnrecognized", "-failed"),
                       SeqIOOptions::genFastqOutGz(
                         njh::files::make_path(failedExtractionDir, name+ "-MIDUnrecognized")));
      seqOut.addReader(njh::pasteAsStr(name, "-MIDunknown", "-failed"),
                 SeqIOOptions::genFastqOutGz(
                   njh::files::make_path(failedExtractionDir, name + "-MIDunknown")));
    } else {
      auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(name, sampleName) ) );
      seqOut.addReader(njh::pasteAsStr(name, "-passed"), seqOutOpts);

      auto failedSeqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(failedExtractionDir, name) );
      seqOut.addReader(njh::pasteAsStr(name, "-failed"), failedSeqOutOpts);

    }
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
    auto maxMidSize = ids.getMaxMIDSize() + 2;
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

      std::string MIDunknownName = ids.containsMids() ? "-MIDunknown" : "";
      if(seq.on_){
        std::string frontPrimerName = "unrecognized";
        std::string backPrimerName = "unrecognized";
        std::string midName;
        bool passesMID = false;
        if(ids.containsMids()) {
          //if also contains mids check for those, the previous extraction by kmer step would have re-oriented everything into forward primer direction
          auto positionForwardResults = ids.pDeterminator_->determineBestForwardPrimerPosFront(seq, corePars.pDetPars, *currentAlignerObj, VecStr{winnerSet});
          seq.reverseComplementRead(false, true);
          auto positionReverseResults = ids.pDeterminator_->determineBestReversePrimerPosFront(seq, corePars.backEndpDetPars, *currentAlignerObj, VecStr{winnerSet});
          seq.reverseComplementRead(false, true);
          frontPrimerName = positionForwardResults.primerName_;
          backPrimerName = positionReverseResults.primerName_;
          //check for mid only if the read passed for both primers
          if(winnerSet == positionForwardResults.primerName_ && winnerSet == positionReverseResults.primerName_) {
            MidDeterminator::MidDeterminePars midpars;
            midpars.searchStart_ = positionForwardResults.start_ > maxMidSize ? positionForwardResults.start_ - maxMidSize : 0;
            midpars.searchStop_ = midpars.searchStart_ + maxMidSize;
            midpars.allowableErrors_ = corePars.primIdsPars.mPars_.allowableErrors_;
            auto midResults = ids.mDeterminator_->searchRead(seq, midpars, midpars);
            auto midResultsProcessed = ids.mDeterminator_->processSearchRead(seq, midResults);
            if(midResultsProcessed.case_ == MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH) {
              //process the primers locations
              //this will both trim and determine the primers, redundant to do this twice for the forward primer with the trimming of the MID it's easier to just recheck
              //the modification below only makes sense if there's only 1 MID and it's attached to the forward primer which is the most common scernario
              auto modifiedSearchPars = corePars.pDetPars;
              modifiedSearchPars.primerStart_ = 0;
              modifiedSearchPars.primerWithin_ = ids.pDeterminator_->getMaxPrimerSize();
              frontPrimerName = ids.pDeterminator_->determineForwardPrimer(seq, modifiedSearchPars, *currentAlignerObj, VecStr{winnerSet});
              seq.reverseComplementRead(false, true);
              backPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq, corePars.backEndpDetPars, *currentAlignerObj, VecStr{winnerSet});
              seq.reverseComplementRead(false, true);
              midName = midResultsProcessed.midName_;
              passesMID = true;
            }
          }
        } else {
          //this will both trim and determine the primers
          frontPrimerName = ids.pDeterminator_->determineForwardPrimer(seq, corePars.pDetPars, *currentAlignerObj, VecStr{winnerSet});
          seq.reverseComplementRead(false, true);
          backPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq, corePars.backEndpDetPars, *currentAlignerObj, VecStr{winnerSet});
          seq.reverseComplementRead(false, true);
          passesMID = true;
        }
        //passes if the front and end primer both match the expected primer based on the winder set from the kmer matching
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
          //check for MID, if no mids the pass MIDs is always true
          if(!passesMID) {
            //failed first min length pass
            seqOut.openWrite(njh::pasteAsStr(winnerSet, "-MIDUnrecognized", "-failed"), seq);
            extractorCase = ExtractionStator::extractCase::BADMID;
          } else {
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
              seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, midName), "-passed"), seq);
            } else {
              //failed
              seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, midName), "-failed"), seq);
            }
          }
        } else {
          seq.name_.append(njh::pasteAsStr("[", "frontPrimerName=", frontPrimerName, ";", "backPrimerName=", backPrimerName, ";","]"));
          //failed primer check
          seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, MIDunknownName), "-failed"), seq);
        }
      } else {
        //failed first min length pass
        seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, MIDunknownName), "-failed"), seq);
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
  outCounts << "\tMIDFailed\tMIDFailedFrac";
  outCounts << "\ttotalFailed\ttotalFailedFrac";
  outCounts << "\tpassed\tpassedFrac";
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
        masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_+
        masterCounts.counts_[setName][false].badmid_ + masterCounts.counts_[setName][true].badmid_;

    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << setName
              << "\t" << totalExtracted
              << "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed)
              << "\t" << readsPerSet[setName]
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(readsPerSet[setName]) / static_cast<double>(totalExtracted));
    //minlen
    outCounts<< "\t" << masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_)/static_cast<double>(totalBad));
    //maxlen
    outCounts << "\t" << masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_)/static_cast<double>(totalBad));
    //quality
    outCounts << "\t" << masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_)/static_cast<double>(totalBad));
    //failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_)/static_cast<double>(totalBad));
    //bad reverse
    outCounts << "\t" << masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_)/static_cast<double>(totalBad));
    //both bad reverse and failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_)/static_cast<double>(totalBad));
    //failed MID
    outCounts << "\t" << masterCounts.counts_[setName][false].badmid_ + masterCounts.counts_[setName][true].badmid_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badmid_ + masterCounts.counts_[setName][true].badmid_)/static_cast<double>(totalBad));

    //bad
    outCounts << "\t" << totalBad
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(totalBad) / static_cast<double>(totalExtracted));
    //good
    outCounts << "\t" << masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_) / static_cast<double>(totalExtracted))
              << std::endl;
    totalReadsPassedAllFilters += masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_;
  }
  {
    uint64_t totalExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << "undetermined"
              << "\t" << totalExtracted
              << "\t" << (totalReadsProcessed == 0 ? 0 : static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed))
              << "\t" << readsPerSet["undetermined"]
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalExtracted))
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
