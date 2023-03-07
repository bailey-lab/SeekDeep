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

  uint32_t minLenCutOff = 100;
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
  setUp.setOption(minLenCutOff, "--minLenCutOff", "Hard cut off min length", false, "filtering");

  //running
  setUp.setOption(numThreads, "--numThreads", "number of threads");


  setUp.finishSetUp(std::cout);
  // run log
  setUp.startARunLog(setUp.pars_.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

  // create Primers and MIDs
  PrimersAndMids ids(idFnp);

  //read in extra info
  ids.addLenCutOffs(lenCutOffsFnp);
  ids.addRefSeqs(refSeqDir);
  uint32_t extractionKmer = ids.addUniqKmerCounts(uniqueKmersPerTargetFnp);

  // set up input
  seqInfo seq;
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();

  uint32_t failedMinLength = 0;

  //set up outputs
  OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
  OutputStream outCounts(outCountsOpts);

  OutOptions outStatsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
  OutputStream outStats(outStatsOpts);

  std::unordered_map<std::string, uint32_t> readsPerSet;
  std::unordered_map<std::string, uint32_t> readsPerSetRevComp;

  VecStr names = ids.getTargets();
  MultiSeqIO seqOut;
  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, name) );
    seqOut.addReader(name, seqOutOpts);
  }
  seqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));

  std::mutex mut;
  std::function<void()> readInComp = [&reader, &ids, &readsPerSet,&readsPerSetRevComp,&mut,&extractionKmer,&seqOut,&sampleName,&rename,&failedMinLength,&minLenCutOff]() {
    SimpleKmerHash hasher;
    seqInfo seq;
    std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
    uint32_t failedMinLengthCurrent = 0;
    while(reader.readNextReadLock(seq)){
      if(len(seq) < minLenCutOff){
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
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);
      for(const auto & readsPerSetCount : readsPerSetCurrent){
        readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
      }
      for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
        readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
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

  outCounts << "sample\ttotalReadsProcessed\ttarget\tcount\tfrac\tforwardCount\tfracForward";
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
              << std::endl;
  }

  outStats << "sampleName\ttotalReadsProcessed\tfailedMinLen_" << minLenCutOff << "\tfailedMinLenFrac\tundetermined\tundeterminedFrac\textracted\textractedFrac\textractedForward\textractedForwardFrac" << std::endl;
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
           << std::endl;

  return 0;
}

} // namespace njhseq
