/*
 * testForFoldbacks.cpp
 *
 *  Created on: Nov 5, 2026
 *      Author: nicholashathaway
 */

#include "SeekDeepUtilsRunner.hpp"
#include <njhseq/seqToolsUtils/FoldbackDetector.hpp>
#include <njhseq/IO/OutputStream.hpp>
#include <atomic>

namespace njhseq {

int SeekDeepUtilsRunner::testForFoldbacks(const njh::progutils::CmdArgs & inputCommands){
	FoldbackDetector::Pars pars;
	std::string sampleName;
	uint32_t numThreads = 1;
	bool writeAllCandidates = false;

	seqSetUp setUp(inputCommands);
	setUp.description_ = "Find reads that are a foldback (inverse chimera), a stretch of sequence and then its reverse complement, and where on the read the fold is, using only the read so it works on reads from anywhere and doesn't need a reference or target. Reads should have their adapters and barcodes trimmed, the same adapter on both ends of a molecule reads as a short inverted repeat and gets called";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastqgz sample.fastq.gz --sampleName sample --numThreads 8 --dout foldbackOut");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--fastq", "--fastqgz", "--fasta", "--fastagz"}, true);
	setUp.setOption(sampleName, "--sampleName", "Name of the sample for the output, by default the name of the input file");
	setUp.setOption(numThreads, "--numThreads", "Number of threads, the output is the same no matter the number of threads");
	setUp.setOption(writeAllCandidates, "--writeAllCandidates", "Write every read with any reverse complement kmers in it, with the column called saying which passed, for trying out other cut offs");
	const std::string detectGroup = "Detection";
	setUp.setOption(pars.kmerLen_, "--kmerLen", "The size of the kmers that vote for a fold", false, detectGroup);
	setUp.setOption(pars.seedBand_, "--seedBand", "The width of the range of kmer position sums that the strongest fold is looked for within", false, detectGroup);
	setUp.setOption(pars.maxKmerOccurrences_, "--maxKmerOccurrences", "Kmers seen more than this many times in a read don't vote", false, detectGroup);
	setUp.setOption(pars.minArmLen_, "--minArmLen", "The shorter arm of a fold has to be at least this long", false, detectGroup);
	setUp.setOption(pars.minAnchors_, "--minAnchors", "The fewest kmers that have to pair between the arms", false, detectGroup);
	setUp.setOption(pars.minIdentity_, "--minIdentity", "The lowest identity between the arms", false, detectGroup);
	setUp.setOption(pars.minAnchorsPerBase_, "--minAnchorsPerBase", "The fewest kmers that have to pair per base of the anchored part of the arm", false, detectGroup);
	setUp.setOption(pars.maxChainGap_, "--maxChainGap", "The furthest apart two paired kmers can be and still be on the same arms", false, detectGroup);
	setUp.setOption(pars.maxLoopPerArm_, "--maxLoopPerArm", "The loop between the arms can be at most this times the longer arm plus --maxLoopFlank", false, detectGroup);
	setUp.setOption(pars.maxLoopFlank_, "--maxLoopFlank", "See --maxLoopPerArm", false, detectGroup);
	setUp.setOption(pars.extensionBand_, "--extensionBand", "How far apart the arms can get from insertions and deletions when finding the arms' ends", false, detectGroup);
	setUp.setOption(pars.extensionXDrop_, "--extensionXDrop", "Stop extending an arm when the score has fallen this far below its best", false, detectGroup);
	setUp.processDirectoryOutputName("testForFoldbacks_TODAY", true);
	if(pars.kmerLen_ < 5 || pars.kmerLen_ > 31){
		setUp.failed_ = true;
		setUp.addWarning("--kmerLen has to be between 5 and 31");
	}
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);
	if(sampleName.empty()){
		sampleName = bfs::path(setUp.pars_.ioOptions_.firstName_).filename().string();
		for(const auto & extension : VecStr{".gz", ".fastq", ".fq", ".fasta", ".fa"}){
			if(njh::endsWith(sampleName, extension)){
				sampleName.erase(sampleName.size() - extension.size());
			}
		}
	}

	const FoldbackDetector detector(pars);
	OutputStream readsOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "foldbackReads.tsv.gz")));
	readsOut << "sample\tread\tlength\tcalled\tleftStart\tleftEnd\trightStart\trightEnd\tleftArmLen\trightArmLen\tloopLen\tanchors\tanchorsPerBase\totherVotes\tarmIdentity\n";

	uint64_t totalReads = 0;
	uint64_t totalBases = 0;
	uint64_t readsWithCandidates = 0;
	uint64_t readsCalled = 0;
	uint64_t basesCalled = 0;
	std::vector<uint32_t> shorterArmLens;
	std::vector<uint32_t> loopLens;
	std::vector<double> identities;

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	//the reads are worked on a batch at a time and written serially so that the order doesn't depend on the threads
	const size_t batchSize = 20000;
	std::vector<seqInfo> batch;
	std::vector<FoldbackDetector::Fold> folds;
	bool moreReads = true;
	while(moreReads){
		batch.clear();
		seqInfo seq;
		while(batch.size() < batchSize){
			if(!reader.readNextRead(seq)){
				moreReads = false;
				break;
			}
			batch.emplace_back(seq);
		}
		if(batch.empty()){
			break;
		}
		folds.assign(batch.size(), FoldbackDetector::Fold{});
		std::atomic<size_t> nextRead{0};
		std::function<void()> detectFolds = [&](){
			size_t pos = 0;
			while((pos = nextRead++) < batch.size()){
				folds[pos] = detector.detect(batch[pos].seq_);
			}
		};
		njh::concurrent::runVoidFunctionThreaded(detectFolds, std::max<uint32_t>(numThreads, 1));
		for(size_t pos = 0; pos < batch.size(); ++pos){
			const auto & fold = folds[pos];
			++totalReads;
			totalBases += len(batch[pos]);
			if(fold.hasCandidate()){
				++readsWithCandidates;
			}
			if(fold.called_){
				++readsCalled;
				basesCalled += len(batch[pos]);
				shorterArmLens.emplace_back(fold.shorterArmLen());
				loopLens.emplace_back(fold.loopLen());
				identities.emplace_back(fold.identity_);
			}
			if(fold.called_ || (writeAllCandidates && fold.hasCandidate())){
				//the demultiplexing programs add information to the name, leave only the read's name
				const auto & name = batch[pos].name_;
				readsOut << sampleName
						<< "\t" << name.substr(0, name.find_first_of(" \t["))
						<< "\t" << len(batch[pos])
						<< "\t" << (fold.called_ ? 1 : 0)
						<< "\t" << fold.leftStart_
						<< "\t" << fold.leftEnd_
						<< "\t" << fold.rightStart_
						<< "\t" << fold.rightEnd_
						<< "\t" << fold.leftArmLen()
						<< "\t" << fold.rightArmLen()
						<< "\t" << fold.loopLen()
						<< "\t" << fold.anchors_
						<< "\t" << fold.anchorsPerBase()
						<< "\t" << fold.otherVotes_
						<< "\t" << fold.identity_ << "\n";
			}
		}
	}

	OutputStream countsOut(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, "foldbackCounts.tsv")));
	countsOut << "sample\ttotalReads\treadsCalled\tpercReadsCalled\tpercBasesCalled\treadsWithCandidates\tmedianShorterArmLen\tmedianLoopLen\tmedianArmIdentity\n";
	countsOut << sampleName
			<< "\t" << totalReads
			<< "\t" << readsCalled
			<< "\t" << (0 == totalReads ? 0.0 : 100.0 * readsCalled / totalReads)
			<< "\t" << (0 == totalBases ? 0.0 : 100.0 * basesCalled / totalBases)
			<< "\t" << readsWithCandidates
			<< "\t" << (shorterArmLens.empty() ? 0.0 : vectorMedianRef(shorterArmLens))
			<< "\t" << (loopLens.empty() ? 0.0 : vectorMedianRef(loopLens))
			<< "\t" << (identities.empty() ? 0.0 : vectorMedianRef(identities)) << std::endl;
	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

}  // namespace njhseq
