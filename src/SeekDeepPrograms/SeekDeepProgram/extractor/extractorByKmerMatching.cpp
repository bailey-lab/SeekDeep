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

  setUp.setOption(lenCutOffsFnp, "--lenCutOffsFnp", "length Cut Offs per target", true, "IDs");
  setUp.setOption(uniqueKmersPerTargetFnp, "--uniqueKmersPerTarget", "unique Kmers Per Target", true, "IDs");
  setUp.setOption(refSeqDir, "--refSeqDir", "reference Seq Dir, should have a .fasta file named by each target", true, "IDs");

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
  corePars.pDetPars.allowable_.hqMismatches_ = 3;
  corePars.pDetPars.allowable_.distances_.query_.coverage_ = .60;
  corePars.pDetPars.allowable_.largeBaseIndel_ = .99;
  corePars.pDetPars.allowable_.oneBaseIndel_ = 2;
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
  ids.addRefSeqs(refSeqDir);
  uint32_t extractionKmer = ids.addUniqKmerCounts(uniqueKmersPerTargetFnp);


  // set up input
  seqInfo seq;
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();

  uint32_t failedMinLength = 0;
  std::unordered_map<std::string, uint32_t> totalPassedByTarget;

  //set up outputs
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
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, name) );
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
  std::function<void()> readInComp = [&reader, &ids, &readsPerSet,&readsPerSetRevComp,&mut,&extractionKmer,&seqOut,&sampleName,&rename,&failedMinLength,&corePars,&alnPool,&totalPassedByTarget]() {
    SimpleKmerHash hasher;
    seqInfo seq;
    std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
    std::unordered_map<std::string, uint32_t> totalPassedByTargetCurrent;
    uint32_t failedMinLengthCurrent = 0;
    auto currentAlignerObj = alnPool.popAligner();
    while(reader.readNextReadLock(seq)){
      if(len(seq) < corePars.smallFragmentCutoff){
        ++failedMinLengthCurrent;
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
      }
      seqOut.openWrite(winnerSet, seq);
      if("undetermined" == winnerSet){
        continue;
      }
      std::string frontPrimerName = ids.pDeterminator_->determineForwardPrimer(seq, corePars.pDetPars, *currentAlignerObj, VecStr{winnerSet});
      seq.reverseComplementRead(false, true);
      std::string backPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq, corePars.backEndpDetPars, *currentAlignerObj, VecStr{winnerSet});
      seq.reverseComplementRead(false, true);
      bool passesFront = frontPrimerName == winnerSet;
      bool passesBack = backPrimerName == winnerSet;
      if(passesBack && passesFront){
        //min len
        ids.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(seq);
        if(seq.on_){
          //max len
          ids.targets_.at(winnerSet).lenCuts_->maxLenChecker_.checkRead(seq);
          if(seq.on_){
            seqOut.openWrite(njh::pasteAsStr(winnerSet, "-passed"), seq);
            ++totalPassedByTargetCurrent[winnerSet];
          }else{
            //failed max length
            seqOut.openWrite(njh::pasteAsStr(winnerSet, "-failed"), seq);
          }
        }else{
          //failed min length
          seqOut.openWrite(njh::pasteAsStr(winnerSet, "-failed"), seq);
        }
      } else {
        seq.name_.append(njh::pasteAsStr("[", "frontPrimerName=", frontPrimerName, ";", "backPrimerName=", backPrimerName, ";", "]"));
        //failed primer check
        seqOut.openWrite(njh::pasteAsStr(winnerSet, "-failed"), seq);
      }
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);
      for(const auto & readsPerSetCount : readsPerSetCurrent){
        readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
      }
      for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
        readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
      }
      for(const auto & totalPassedByTargetCount : totalPassedByTargetCurrent){
        totalPassedByTarget[totalPassedByTargetCount.first] += totalPassedByTargetCount.second;
      }
      failedMinLength += failedMinLengthCurrent;
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

  for(const auto & totalPassedByTargetCount : totalPassedByTarget){
    totalReadsPassedAllFilters += totalPassedByTargetCount.second;
  }

  outCounts << "sample\ttotalReadsProcessed\ttarget\tcount\tfrac\tforwardCount\tfracForward";
  outCounts << "\tpasssed\tpassedFrac";
  outCounts << std::endl;
  uint64_t totalExtractedAllTargets = 0;
  uint64_t totalExtractedAllTargetsForard = 0;
  uint64_t totalExtractedUndetermined = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
  for(const auto & setName : ids.getTargets()){
    uint64_t totalExtracted = readsPerSet[setName] + readsPerSetRevComp[setName];
    totalExtractedAllTargets += totalExtracted;
    totalExtractedAllTargetsForard += readsPerSet[setName];
    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << setName
              << "\t" << totalExtracted
              << "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed)
              << "\t" << readsPerSet[setName]
              << "\t" << static_cast<double>(readsPerSet[setName]) / static_cast<double>(totalExtracted)
              << "\t" << totalPassedByTarget[setName]
              << "\t" << static_cast<double>(totalPassedByTarget[setName]) / static_cast<double>(totalExtracted)
              << std::endl;
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
              << std::endl;
  }

  outStats << "sampleName\ttotalReadsProcessed\tfailedMinLen_" << corePars.smallFragmentCutoff << "\tfailedMinLenFrac\tundetermined\tundeterminedFrac\textracted\textractedFrac\textractedForward\textractedForwardFrac\tpassed\tpassedFrac" << std::endl;
  outStats << sampleName
           << "\t" << totalReadsProcessed + failedMinLength
           << "\t" << failedMinLength
           << "\t" << static_cast<double>(failedMinLength) / static_cast<double>(totalReadsProcessed + failedMinLength)
           << "\t" << totalExtractedUndetermined
           << "\t" << static_cast<double>(totalExtractedUndetermined) / static_cast<double>(totalExtractedUndetermined + totalExtractedAllTargets)
           << "\t" << totalExtractedAllTargets
           << "\t" << static_cast<double>(totalExtractedAllTargets) / static_cast<double>(totalReadsProcessed + failedMinLength)
           << "\t" << totalExtractedAllTargetsForard
           << "\t" << static_cast<double>(totalExtractedAllTargetsForard) / static_cast<double>(totalExtractedAllTargets)
           << "\t" << totalReadsPassedAllFilters
           << "\t" << static_cast<double>(totalReadsPassedAllFilters) / static_cast<double>(totalReadsProcessed + failedMinLength)
           << std::endl;

  return 0;
}

} // namespace njhseq
