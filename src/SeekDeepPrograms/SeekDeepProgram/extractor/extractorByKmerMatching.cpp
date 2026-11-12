//
// Created by Nicholas Hathaway on 3/6/23.
//
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"

#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <njhseq/objects/kmer/KmerSetInvertedIndex.hpp>
#include <njhseq/objects/kmer/KmerCollectionReadMatcher.hpp>

namespace njhseq {

namespace {

/**
 * @brief Where each target sits along a read that hit several targets, from the positions of k-mers unique to one target
 *
 * Reads that join several amplicons end to end (e.g. amplicons ligated to each other before barcoding) hold each amplicon's
 * unique k-mers in one stretch. Hits of the same target and orientation are chained while they are within maxGap_ of each
 * other and the chain stays no longer than the target can be, so two copies of one target in a row stay separate; chains with
 * fewer than minKmers_ hits are dropped as noise.
 */
struct MultihitSplitter {
	struct Pars {
		bool split_ = true;
		uint32_t minKmers_ = 20;      ///< unique k-mer hits a stretch needs to count as an amplicon
		uint32_t maxGap_ = 100;       ///< largest gap between unique k-mer hits within one stretch
		uint32_t primerSlop_ = 20;    ///< sequence kept past the longest primer on each side of a stretch
		std::string oneEnd_ = "none"; ///< reads demultiplexed on one end: split none of them, keep the amplicon nearest that barcode, or keep all
	};

	struct Segment {
		uint32_t setIdx_;
		bool revComp_;
		uint32_t start_;  ///< first base of the first k-mer
		uint32_t end_;    ///< one past the last base of the last k-mer
		uint32_t kmers_;
		uint32_t lastPos_;
	};

	/**
	 * @param allowed by index order, the targets a stretch may come from
	 * @param maxSpan by index order, the longest a stretch of that target can be
	 */
	static std::vector<Segment> findSegments(const std::string & seq, const KmerSetInvertedIndex & index, const SimpleKmerHash & hasher,
			uint32_t klen, const std::vector<char> & allowed, const std::vector<uint32_t> & maxSpan, const Pars & pars) {
		std::vector<Segment> done;
		if (seq.size() <= klen) {
			return done;
		}
		//the stretch currently open for each target and orientation, index*2 + revComp
		std::vector<int64_t> open(index.setNames_.size() * 2, -1);
		std::vector<Segment> all;
		for (uint32_t pos = 0; pos < seq.size() - klen + 1; ++pos) {
			for (const bool revComp : {false, true}) {
				const uint64_t hash = revComp ? hasher.revCompHash(seq, pos, klen) : hasher.hash(seq, pos, klen);
				const auto search = index.kmerToSet_.find(hash);
				if (search == index.kmerToSet_.end() || (search->second & KmerSetInvertedIndex::SHARED_FLAG)) {
					continue;
				}
				const uint32_t setIdx = search->second;
				if (!allowed[setIdx]) {
					continue;
				}
				auto & openIdx = open[setIdx * 2 + revComp];
				if (openIdx >= 0) {
					auto & seg = all[openIdx];
					if (pos <= seg.lastPos_ + pars.maxGap_ && pos + klen - seg.start_ <= maxSpan[setIdx]) {
						seg.lastPos_ = pos;
						seg.end_ = pos + klen;
						++seg.kmers_;
						continue;
					}
				}
				openIdx = static_cast<int64_t>(all.size());
				all.emplace_back(Segment{setIdx, revComp, pos, pos + klen, 1, pos});
			}
		}
		for (const auto & seg : all) {
			if (seg.kmers_ >= pars.minKmers_) {
				done.emplace_back(seg);
			}
		}
		std::sort(done.begin(), done.end(), [](const Segment & a, const Segment & b) { return a.start_ < b.start_; });
		return done;
	}
};

}  // namespace

int SeekDeepRunner::extractorByKmerMatching(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path idFnp;
  bfs::path lenCutOffsFnp;
  bfs::path uniqueKmersPerTargetFnp;
  bfs::path refSeqDir;
  CoreExtractorPars corePars;
  corePars.smallFragmentCutoff = 100;
  corePars.largeFragmentCutOff = 10000;
  uint32_t numThreads = 1;

  std::string sampleName = "sample";
  bool rename = false;
  bool ignoreMissingKmerSets = false;
  SeekDeepSetUp setUp(inputCommands);
  setUp.setOption(ignoreMissingKmerSets, "--ignoreMissingKmerSets", "ignore Missing Kmer Sets");

  //id
  setUp.setOption(idFnp, "--ids,--id", "Primers file", true, "IDs");

  setUp.setOption(lenCutOffsFnp, "--lenCutOffsFnp,--lenCutOffs", "length Cut Offs per target", true, "IDs");
  setUp.setOption(uniqueKmersPerTargetFnp, "--uniqueKmersPerTarget", "unique Kmers Per Target", true, "IDs");
  setUp.setOption(refSeqDir, "--refSeqDir", "reference Seq Dir, should have a .fasta file named by each target", false, "IDs");

  //in and out
  setUp.processDefaultReader(VecStr{"--fastq", "--fastqgz", "--fasta", "--fastagz", "--bam"});
  setUp.processDirectoryOutputName(true);
  if("sample" == sampleName){
    sampleName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
    if(njh::endsWith(setUp.pars_.ioOptions_.firstName_.string(), ".gz")){
      sampleName = bfs::path(sampleName).replace_extension("").string();
    }
    //for unaligned bams (e.g. from nanopore dorado demux), default to the SM of the read group(s) in the header
    if(SeqIOOptions::inFormats::BAM == setUp.pars_.ioOptions_.inFormat_){
      setSampleNameFromBamHeader(setUp.pars_.ioOptions_.firstName_, sampleName, "--sampleName");
    }
  }
  setUp.setOption(sampleName, "--sampleName", "sample name, if input is a bam this defaults to the SM of the read group in the header otherwise the input file's basename");
  setUp.setOption(rename, "--rename", "rename input extracted reads");

  // filtering
  setUp.setOption(corePars.smallFragmentCutoff, "--minLenCutOff", "Hard cut off min length", false, "filtering");
  setUp.setOption(corePars.largeFragmentCutOff, "--maxLenCutOff", "Hard cut off max length", false, "filtering");

  uint32_t hardKmersPerTarget = 5;
  double kmerFracPerTarget = 0.10;
  setUp.setOption(hardKmersPerTarget, "--hardKmersPerTarget", "hard occurencet cut off Kmers Per Target", false, "filtering");
  setUp.setOption(kmerFracPerTarget, "--kmerFracPerTarget", "at minimum have to have this fraction of kmers in a Target", false, "filtering");

  //running
  setUp.setOption(numThreads, "--numThreads", "number of threads");
  bool noInvertedKmerIndex = false;
  setUp.setOption(noInvertedKmerIndex, "--noInvertedKmerIndex",
      "Probe every target's k-mer set for every read k-mer instead of using one k-mer -> target(s) index. The index gives identical counts, counting a k-mer shared between targets toward each of them exactly as the probe does; this exists to fall back if that is ever in doubt");
  const bool useInvertedKmerIndex = !noInvertedKmerIndex;

  MultihitSplitter::Pars splitPars;
  bool noMultihitSplit = false;
  setUp.setOption(noMultihitSplit, "--noMultihitSplit",
      "Don't split reads that hit several targets; by default each target's stretch of the read (found by k-mers unique to that target) is extracted as its own read, e.g. amplicons ligated end to end", false, "Multihit Split");
  splitPars.split_ = !noMultihitSplit;
  setUp.setOption(splitPars.minKmers_, "--multihitSplitMinKmers", "Unique k-mer hits a stretch of a multihit read needs to count as a target's amplicon", false, "Multihit Split");
  setUp.setOption(splitPars.maxGap_, "--multihitSplitMaxGap", "Largest gap between a target's unique k-mer hits within one stretch of a multihit read", false, "Multihit Split");
  setUp.setOption(splitPars.primerSlop_, "--multihitSplitPrimerSlop", "Sequence kept on each side of a stretch past the longest primer, so the primers are still found", false, "Multihit Split");
  setUp.setOption(splitPars.oneEnd_, "--multihitSplitOneEnd",
      "Multihit reads demultiplexed on one end's barcode (barcodeEnds=1 in the name): none leaves them as multihit, nearest keeps only the amplicon nearest that barcode (needs assignedOn and direction in the name), all keeps every piece. "
      "One end multihit reads are often from other samples: in ligation Pool 1 they put 268-444 reads on 5 targets in the no template controls even keeping only the nearest piece, and in run 2 pieces past the nearest were ~5% from the other panel. Reads without barcodeEnds are always split", false, "Multihit Split");
  if(!njh::in(splitPars.oneEnd_, VecStr{"none", "nearest", "all"})){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "--multihitSplitOneEnd should be none, nearest or all, not " << splitPars.oneEnd_ << "\n";
    throw std::runtime_error { ss.str() };
  }


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
  uint32_t extractionKmer = ids.addKmerCounts(uniqueKmersPerTargetFnp, ignoreMissingKmerSets);
  //the table can hold all kmers per target (overlapping targets share kmers) or only unique kmers. All kmers pick the
  //winning target and kmers unique to a target count how many targets a read hits. A shared kmer counts toward every
  //target holding it.
  const njhseq::KmerCollectionReadMatcher kmerMatcher(ids.allKmersPerTarget_, extractionKmer, useInvertedKmerIndex);
  const auto & kmerTargetNames = kmerMatcher.names();
  if(useInvertedKmerIndex){
    auto invIndex = njhseq::KmerSetInvertedIndex::get(ids.allKmersPerTarget_);
    std::cout << "kmerSetExclusivity\tsets\t" << invIndex->nSets_
        << "\ttotalKmers\t" << invIndex->totalKmers_
        << "\tsharedKmers\t" << invIndex->nSharedKmers_
        << "\texclusive\t" << njh::boolToStr(0 == invIndex->nSharedKmers_) << std::endl;
  }
  //for splitting multihit reads, positions of k-mers unique to one target are looked up in the inverted index
  std::shared_ptr<const KmerSetInvertedIndex> splitIndex;
  std::vector<uint32_t> matcherToSplitIdx(kmerTargetNames.size(), 0);
  std::vector<uint32_t> splitMaxSpan;
  const SimpleKmerHash splitHasher;
  if(splitPars.split_){
    splitIndex = KmerSetInvertedIndex::get(ids.allKmersPerTarget_);
    std::unordered_map<std::string, uint32_t> splitIdxByName;
    for(uint32_t idx = 0; idx < splitIndex->setNames_.size(); ++idx){
      splitIdxByName[splitIndex->setNames_[idx]] = idx;
    }
    for(uint32_t idx = 0; idx < kmerTargetNames.size(); ++idx){
      matcherToSplitIdx[idx] = splitIdxByName.at(kmerTargetNames[idx]);
    }
    //a stretch can't be longer than its target's longest allowed amplicon, which keeps two copies of a target in a row apart
    splitMaxSpan.assign(splitIndex->setNames_.size(), corePars.largeFragmentCutOff);
    for(uint32_t idx = 0; idx < splitIndex->setNames_.size(); ++idx){
      const auto & name = splitIndex->setNames_[idx];
      if(njh::in(name, ids.targets_) && nullptr != ids.targets_.at(name).lenCuts_){
        splitMaxSpan[idx] = ids.targets_.at(name).lenCuts_->maxLenChecker_.maxLen_;
      }
    }
  }


  // set up input
  seqInfo seq;
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();



  //qual check off;
  std::unique_ptr<ReadChecker> qualChecker = std::make_unique<ReadCheckerQualCheck>(corePars.qPars_.qualCheck_,
                                                                                    corePars.qPars_.qualCheckCutOff_, true);

  //set up extraction counts outputs
  ExtractionStator masterCounts;
  //counts per target and MID, keyed by the name of the extracted file, <target><MID>, or <target>-MIDunknown/-MIDUnrecognized when the MID wasn't determined
  ExtractionStator perMidCounts;
  OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
  OutputStream outCounts(outCountsOpts);

  OutOptions outStatsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
  OutputStream outStats(outStatsOpts);

  std::unordered_map<std::string, uint32_t> readsPerSet;
  std::unordered_map<std::string, uint32_t> readsPerSetRevComp;
  //pieces cut from multihit reads, per target; the reads they came from are counted in readsPerSet under "multihitSplit"
  std::unordered_map<std::string, uint32_t> piecesPerSet;
  std::unordered_map<std::string, uint32_t> piecesPerSetRevComp;
  uint64_t multihitSplitFoldback = 0;

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
  seqOut.addReader("multihit", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "multihit")));
  if(splitPars.split_){
    seqOut.addReader("multihitSplit", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "multihitSplit")));
  }

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
                                      &seqOut,&sampleName,&rename,
                                      &corePars,&alnPool,
                                      &qualChecker,
                                      &masterCounts,
                                      &perMidCounts,
                                      &hardKmersPerTarget,
                                      &kmerFracPerTarget,
                                      &kmerMatcher,
                                      &kmerTargetNames,
                                      &splitPars,
                                      &splitIndex,
                                      &splitHasher,
                                      &splitMaxSpan,
                                      &matcherToSplitIdx,
                                      &extractionKmer,
                                      &maxPrimerSize,
                                      &piecesPerSet,
                                      &piecesPerSetRevComp,
                                      &multihitSplitFoldback
                                      ]() {
    auto maxMidSize = ids.getMaxMIDSize() + 2;
    seqInfo seq;
    std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
    std::unordered_map<std::string, uint32_t> piecesPerSetCurrent;
    std::unordered_map<std::string, uint32_t> piecesPerSetRevCompCurrent;
    uint64_t multihitSplitFoldbackCurrent = 0;


    ExtractionStator masterCountsCurrent;
    ExtractionStator perMidCountsCurrent;

    auto kmerMatch = kmerMatcher.genMatch();

    auto currentAlignerObj = alnPool.popAligner();
    //primer, MID, length and quality checks for a read assigned to winnerSet, writing it to that target's passed or failed output
    auto extractTarget = [&](seqInfo & seq, const std::string & winnerSet, ExtractionStator::extractCase extractorCase){
        //only check min len if INVERSECHIMERA wasn't true
        if(seq.on_) {
          //min len
          ids.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(seq);
          if(!seq.on_){
            extractorCase = ExtractionStator::extractCase::MINLENBAD;
          }
        }

        std::string MIDunknownName = ids.containsMids() ? "-MIDunknown" : "";
        std::string perMidName = njh::pasteAsStr(winnerSet, MIDunknownName);
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
                perMidName = njh::pasteAsStr(winnerSet, midName);
              } else {
                perMidName = njh::pasteAsStr(winnerSet, "-MIDUnrecognized");
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
        if(ids.containsMids()){
          perMidCountsCurrent.increaseCounts(perMidName, seq.name_, extractorCase);
        }
    };

    while(reader.readNextReadLock(seq)){
      ++masterCountsCurrent.totalReadCount_;
      if(len(seq) < corePars.smallFragmentCutoff){
        ++masterCountsCurrent.smallFrags_;
        continue;
      }
      if (len(seq) > corePars.largeFragmentCutOff) {
        ++masterCountsCurrent.largeFrags_;
        continue;
      }
      //all kmers per target pick the winner, kmers unique to a target count the targets hit
      kmerMatch.reset();
      kmerMatcher.addKmers(seq.seq_, true, kmerMatch);
      kmerMatcher.countFound(kmerMatch);

      //filter
      const double nInputKmers = static_cast<double>(kmerMatch.kmers.size());
      for(auto * found : {&kmerMatch.found, &kmerMatch.foundRevComp, &kmerMatch.uniqueFound, &kmerMatch.uniqueFoundRevComp}){
        for (auto & count : *found) {
          if(count < hardKmersPerTarget || static_cast<double>(count)/nInputKmers <= kmerFracPerTarget) {
            count = 0;
          }
        }
      }
      kmerMatcher.determineWinner(kmerMatch);
      std::string winnerSet = kmerMatch.hasWinner() ? kmerTargetNames[kmerMatch.winnerIdx] : "undetermined";
      const bool winnerRevComp = kmerMatch.winnerRevComp;
      // doesn't take into account inverse hits, 1 target with large hits in forward and reverse
      uint32_t target_hit_count = njhseq::KmerCollectionReadMatcher::countCollectionsHit(kmerMatch, hardKmersPerTarget);
      bool multi_hit = target_hit_count > 1;
      if (multi_hit) {
        winnerSet = "multihit";
      }
      if(multi_hit && splitPars.split_){
        //stretches of the targets that made this a multihit
        std::vector<char> allowed(splitIndex->setNames_.size(), 0);
        for(uint32_t idx = 0; idx < kmerTargetNames.size(); ++idx){
          if(kmerMatch.uniqueFound[idx] >= hardKmersPerTarget || kmerMatch.uniqueFoundRevComp[idx] >= hardKmersPerTarget){
            allowed[matcherToSplitIdx[idx]] = 1;
          }
        }
        const auto segments = MultihitSplitter::findSegments(seq.seq_, *splitIndex, splitHasher, extractionKmer, allowed, splitMaxSpan, splitPars);
        std::set<uint32_t> segmentTargets;
        for(const auto & seg : segments){
          segmentTargets.emplace(seg.setIdx_);
        }
        MetaDataInName inMeta;
        const bool hasMeta = MetaDataInName::nameHasMetaData(seq.name_);
        if(hasMeta){
          inMeta = MetaDataInName(seq.name_);
        }
        const bool oneEndRead = hasMeta && inMeta.containsMeta("barcodeEnds") && "1" == inMeta.getMeta("barcodeEnds");
        if(segmentTargets.size() > 1 && !(oneEndRead && "none" == splitPars.oneEnd_)){
          //a target seen again in the other orientation means the molecule was read through and back (a foldback),
          //so only the stretches before it turns are separate molecules
          std::vector<MultihitSplitter::Segment> kept;
          std::unordered_map<uint32_t, bool> orientationSeen;
          bool foldback = false;
          for(const auto & seg : segments){
            auto seen = orientationSeen.find(seg.setIdx_);
            if(seen != orientationSeen.end() && seen->second != seg.revComp_){
              foldback = true;
              break;
            }
            orientationSeen[seg.setIdx_] = seg.revComp_;
            kept.emplace_back(seg);
          }
          if(foldback){
            ++multihitSplitFoldbackCurrent;
          }
          //demultiplexed on one end's barcode, the amplicon nearest that barcode is the one it labels best
          if("nearest" == splitPars.oneEnd_ && oneEndRead && inMeta.containsMeta("assignedOn") && inMeta.containsMeta("direction") && kept.size() > 1){
            const bool barcodeAtFront = foldback || (("back" != inMeta.getMeta("assignedOn")) == ("rev" != inMeta.getMeta("direction")));
            const auto nearest = barcodeAtFront ? kept.front() : kept.back();
            kept = {nearest};
          }
          //a piece has the read's barcode at one of its own ends at most, so it doesn't count as support from both ends
          //(processClusters --minBothEndsBarcodeReads), otherwise a few joined molecules from other samples let bleed through, e.g. into no template controls
          std::string keptMeta;
          if(hasMeta && inMeta.containsMeta("barcodeEnds")){
            MetaDataInName barcodeEndsMeta;
            barcodeEndsMeta.addMeta("barcodeEnds", 1);
            keptMeta = barcodeEndsMeta.createMetaName();
          }
          ++readsPerSetCurrent["multihitSplit"];
          seqOut.openWrite("multihitSplit", seq);
          const uint64_t pad = maxPrimerSize + splitPars.primerSlop_;
          const auto threadId = estd::to_string(std::this_thread::get_id());
          for(const auto & seg : kept){
            const uint32_t start = seg.start_ > pad ? seg.start_ - pad : 0;
            const uint32_t end = std::min<uint64_t>(len(seq), seg.end_ + pad);
            auto piece = seq.getSubRead(start, end - start);
            const std::string & pieceSet = splitIndex->setNames_[seg.setIdx_];
            if(seg.revComp_){
              ++piecesPerSetRevCompCurrent[pieceSet];
            } else {
              ++piecesPerSetCurrent[pieceSet];
            }
            if(rename){
              piece.name_ = njh::pasteAsStr(sampleName, ".", pieceSet, ".mhSplit", piecesPerSetCurrent[pieceSet] + piecesPerSetRevCompCurrent[pieceSet], ".", threadId, keptMeta);
              if(seg.revComp_){
                piece.reverseComplementRead(false, true);
                piece.name_.append("_Comp");
              }
            } else {
              piece.name_ = njh::pasteAsStr(seq.name_, "_mhSplit", start, "-", end);
              if(seg.revComp_){
                piece.reverseComplementRead(true, true);
              }
            }
            seqOut.openWrite(pieceSet, piece);
            extractTarget(piece, pieceSet, ExtractionStator::extractCase::GOOD);
          }
          continue;
        }
      }
      if(winnerRevComp){
        ++readsPerSetRevCompCurrent[winnerSet];
        seq.reverseComplementRead(true, true);
      } else {
        ++readsPerSetCurrent[winnerSet];
      }

      if(rename){
        //keep how many of the read's ends had the sample's barcode (nanopore demultiplexing), for clustering to count
        std::string keptMeta;
        if(MetaDataInName::nameHasMetaData(seq.name_)){
          MetaDataInName inMeta(seq.name_);
          if(inMeta.containsMeta("barcodeEnds")){
            MetaDataInName barcodeEndsMeta;
            barcodeEndsMeta.addMeta("barcodeEnds", inMeta.getMeta("barcodeEnds"));
            keptMeta = barcodeEndsMeta.createMetaName();
          }
        }
        auto threadId = estd::to_string(std::this_thread::get_id());
        seq.name_ = njh::pasteAsStr(sampleName, ".",winnerSet, ".", readsPerSetCurrent[winnerSet] + readsPerSetRevCompCurrent[winnerSet], ".", threadId, keptMeta);
        if(winnerRevComp){
          seq.name_.append("_Comp");
        }
      }

      seqOut.openWrite(winnerSet, seq);
      if("undetermined" == winnerSet){
        ++masterCountsCurrent.readsUnrecBarcode_;
        continue;
      }
      if("multihit" == winnerSet){
        ++masterCountsCurrent.multihit_;
        continue;
      }

      auto extractorCase = ExtractionStator::extractCase::GOOD;

      //check for inverse hit (has hit of the same target in both the forward and reverse directions)
      if (kmerMatch.found[kmerMatch.winnerIdx] >= hardKmersPerTarget && kmerMatch.foundRevComp[kmerMatch.winnerIdx]  >= hardKmersPerTarget) {
        extractorCase = ExtractionStator::extractCase::INVERSECHIMERA;
        seq.on_ = false;
      }
      extractTarget(seq, winnerSet, extractorCase);
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);

      masterCounts.addOtherExtractorCounts(masterCountsCurrent);
      perMidCounts.addOtherExtractorCounts(perMidCountsCurrent);

      for(const auto & readsPerSetCount : readsPerSetCurrent){
        readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
      }
      for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
        readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
      }
      for(const auto & count : piecesPerSetCurrent){
        piecesPerSet[count.first] += count.second;
      }
      for(const auto & count : piecesPerSetRevCompCurrent){
        piecesPerSetRevComp[count.first] += count.second;
      }
      multihitSplitFoldback += multihitSplitFoldbackCurrent;
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
  outCounts << "\tinverseChimeraFailed\tinverseChimeraFailedFrac";
  outCounts << "\tmaxLenFailed\tmaxLenFailedFrac";
  outCounts << "\tqualityFailed\tqualityFailedFrac";
  outCounts << "\tforwardPrimerFailed\tforwardPrimerFailedFrac";
  outCounts << "\treversePrimerFailed\treversePrimerFailedFrac";
  outCounts << "\tbothForRevPrimerFailed\tbothForRevPrimerFailedFrac";
  outCounts << "\tMIDFailed\tMIDFailedFrac";
  outCounts << "\ttotalFailed\ttotalFailedFrac";
  outCounts << "\tpassed\tpassedFrac";
  outCounts << "\tfromMultihitSplit";
  outCounts << std::endl;


  uint64_t totalExtractedAllTargets = 0;
  uint64_t totalExtractedAllTargetsForward = 0;
  uint64_t totalExtractedUndetermined = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
  uint64_t totalExtractedMultihit = readsPerSet["multihit"] + readsPerSetRevComp["multihit"];


  for(const auto & setName : ids.getTargets()){
    const uint64_t fromMultihitSplit = piecesPerSet[setName] + piecesPerSetRevComp[setName];
    uint64_t totalExtracted = readsPerSet[setName] + readsPerSetRevComp[setName] + fromMultihitSplit;
    totalExtractedAllTargets += totalExtracted;
		totalExtractedAllTargetsForward += readsPerSet[setName] + piecesPerSet[setName];

    uint64_t totalBad = masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_ +
          masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_+
          masterCounts.counts_[setName][false].inverse_chimera_ + masterCounts.counts_[setName][true].inverse_chimera_+
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
              << "\t" << readsPerSet[setName] + piecesPerSet[setName]
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(readsPerSet[setName] + piecesPerSet[setName]) / static_cast<double>(totalExtracted));
    //minlen
    outCounts<< "\t" << masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_)/static_cast<double>(totalBad));
    //inverse chimera
    outCounts<< "\t" << masterCounts.counts_[setName][false].inverse_chimera_ + masterCounts.counts_[setName][true].inverse_chimera_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].inverse_chimera_ + masterCounts.counts_[setName][true].inverse_chimera_)/static_cast<double>(totalBad));
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
              << "\t" << fromMultihitSplit
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
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
              << std::endl;
  }

  outStats << "sampleName\ttotalReadsProcessed\tfailedMinLen_" << corePars.smallFragmentCutoff << "\tfailedMinLenFrac\tfailedMaxLen_" << corePars.largeFragmentCutOff << "\tfailedMaxLenFrac\tundetermined\tundeterminedFrac\tmultihit\tmultihitFrac\textracted\textractedFrac\textractedForward\textractedForwardFrac\tpassed\tpassedFrac\tmultihitSplit\tmultihitSplitFrac\tpiecesFromMultihitSplit\tmultihitSplitFoldback" << std::endl;
  //reads that hit several targets and were split into pieces, the pieces are counted in extracted and passed
  const uint64_t totalMultihitSplit = readsPerSet["multihitSplit"];
  uint64_t totalPiecesFromMultihitSplit = 0;
  for(const auto & count : piecesPerSet){
    totalPiecesFromMultihitSplit += count.second;
  }
  for(const auto & count : piecesPerSetRevComp){
    totalPiecesFromMultihitSplit += count.second;
  }
  uint32_t grand_total_reads = totalReadsProcessed + masterCounts.smallFrags_ + masterCounts.largeFrags_;
  outStats << sampleName
					 << "\t" << grand_total_reads
					 << "\t" << masterCounts.smallFrags_
					 << "\t" << static_cast<double>(masterCounts.smallFrags_) / static_cast<double>(grand_total_reads)
           << "\t" << masterCounts.largeFrags_
           << "\t" << static_cast<double>(masterCounts.largeFrags_) / static_cast<double>(grand_total_reads)
					 << "\t" << totalExtractedUndetermined
					 << "\t" << static_cast<double>(totalExtractedUndetermined) / static_cast<double>(totalReadsProcessed)
           << "\t" << totalExtractedMultihit
           << "\t" << static_cast<double>(totalExtractedMultihit) / static_cast<double>(totalReadsProcessed)
					 << "\t" << totalExtractedAllTargets
					 << "\t" << static_cast<double>(totalExtractedAllTargets) / static_cast<double>(grand_total_reads)
					 << "\t" << totalExtractedAllTargetsForward
					 << "\t" << static_cast<double>(totalExtractedAllTargetsForward) / static_cast<double>(totalExtractedAllTargets)
           << "\t" << totalReadsPassedAllFilters
           << "\t" << static_cast<double>(totalReadsPassedAllFilters) / static_cast<double>(grand_total_reads)
           << "\t" << totalMultihitSplit
           << "\t" << static_cast<double>(totalMultihitSplit) / static_cast<double>(totalReadsProcessed)
           << "\t" << totalPiecesFromMultihitSplit
           << "\t" << multihitSplitFoldback
           << std::endl;

  if(ids.containsMids()){
    //the counts per target per MID, name matches the extracted file <target><MID> so it can be linked to the sample the MID belongs to
    //reads failing before the MID is looked for (inverse chimera, min length, either primer) are counted under MIDunknown, reads whose MID wasn't found under MIDUnrecognized
    OutputStream outPerMid(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfilePerMID.tab.txt"));
    outPerMid << "sample\ttarget\tmid\tname\ttotal"
              << "\tpassed\tpassedFrac"
              << "\tminLenFailed\tmaxLenFailed\tqualityFailed\tinverseChimeraFailed"
              << "\tforwardPrimerFailed\treversePrimerFailed\tbothForRevPrimerFailed\tMIDFailed"
              << "\ttotalFailed" << std::endl;
    auto midNames = getVectorOfMapKeys(ids.mids_);
    njh::sort(midNames);
    midNames.emplace_back("MIDUnrecognized");
    midNames.emplace_back("MIDunknown");
    for(const auto & setName : ids.getTargets()){
      for(const auto & midName : midNames){
        const bool undeterminedMid = "MIDUnrecognized" == midName || "MIDunknown" == midName;
        const std::string name = undeterminedMid ? njh::pasteAsStr(setName, "-", midName) : njh::pasteAsStr(setName, midName);
        //both read directions
        ExtractionStator::extractCounts counts;
        if(njh::in(name, perMidCounts.counts_)){
          for(const auto & dir : perMidCounts.counts_.at(name)){
            counts.good_ += dir.second.good_;
            counts.bad_ += dir.second.bad_;
            counts.minLenBad_ += dir.second.minLenBad_;
            counts.maxLenBad_ += dir.second.maxLenBad_;
            counts.qualityFailed_ += dir.second.qualityFailed_;
            counts.inverse_chimera_ += dir.second.inverse_chimera_;
            counts.badForward_ += dir.second.badForward_;
            counts.badReverse_ += dir.second.badReverse_;
            counts.failedBothPrimers_ += dir.second.failedBothPrimers_;
            counts.badmid_ += dir.second.badmid_;
          }
        }
        const uint32_t total = counts.good_ + counts.bad_;
        outPerMid << sampleName
                  << "\t" << setName
                  << "\t" << midName
                  << "\t" << name
                  << "\t" << total
                  << "\t" << counts.good_
                  << "\t" << (0 == total ? 0 : static_cast<double>(counts.good_) / static_cast<double>(total))
                  << "\t" << counts.minLenBad_
                  << "\t" << counts.maxLenBad_
                  << "\t" << counts.qualityFailed_
                  << "\t" << counts.inverse_chimera_
                  << "\t" << counts.badForward_
                  << "\t" << counts.badReverse_
                  << "\t" << counts.failedBothPrimers_
                  << "\t" << counts.badmid_
                  << "\t" << counts.bad_
                  << std::endl;
      }
    }
  }

  return 0;
}

} // namespace njhseq
