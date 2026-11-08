/*
 * nanoporeAmpliconReadQC.cpp
 *
 *  Created on: Nov 6, 2026
 *      Author: nicholashathaway
 */

#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/TarAmpSetupUtils/NanoporeAmpliconReadQC.hpp"
#include "SeekDeep/objects/TarAmpSetupUtils/NanoporeReadPlacer.hpp"
#include <njhseq/IO/OutputStream.hpp>
#include <atomic>
#include <regex>

namespace njhseq {

namespace {

struct SampleInput{
	std::string name_;
	SeqIOOptions opts_;
	bool capped_{false}; /**< not a sample, only look at some of its reads */
};

std::string sampleNameFromFile(const bfs::path & fnp){
	auto name = fnp.filename().string();
	for(const auto & extension : VecStr{".gz", ".fastq", ".fq", ".fasta", ".fa"}){
		if(njh::endsWith(name, extension)){
			name.erase(name.size() - extension.size());
		}
	}
	return name;
}

/**@brief The fastq files of a directory, a demultiplexing's Undetermined file doesn't belong to a sample so it's skipped unless asked for
 */
std::vector<SampleInput> readInSamplesFromDir(const bfs::path & fastqDir, bool includeUndetermined){
	std::vector<SampleInput> ret;
	for(const auto & file : njh::files::listAllFiles(fastqDir, false, {std::regex{R"(.*\.(fastq|fq)(\.gz)?$)"}})){
		const auto name = sampleNameFromFile(file.first);
		if(!includeUndetermined && njh::beginsWith(njh::strToLowerRet(name), "undetermined")){
			continue;
		}
		ret.emplace_back(SampleInput{name, njh::endsWith(file.first.string(), ".gz") ? SeqIOOptions::genFastqInGz(file.first) : SeqIOOptions::genFastqIn(file.first),
			njh::beginsWith(njh::strToLowerRet(name), "undetermined")});
	}
	njh::sort(ret, [](const SampleInput & s1, const SampleInput & s2){
		return s1.name_ < s2.name_;
	});
	return ret;
}

/**@brief The same hash everywhere, picks which reads of a capped input are looked at
 */
uint64_t hashReadName(const std::string & name){
	uint64_t ret = 14695981039346656037ULL;
	for(const auto c : name){
		ret ^= static_cast<unsigned char>(c);
		ret *= 1099511628211ULL;
	}
	return ret;
}

std::string readNameOnly(const std::string & name){
	return name.substr(0, name.find_first_of(" \t["));
}

/**@brief A value the demultiplexing added to a read's name, name[tag=value;tag=value]
 */
std::string getNameTag(const std::string & name, const std::string & tag){
	const auto pat = tag + "=";
	auto pos = name.find(pat);
	while(pos != std::string::npos){
		if(pos > 0 && ('[' == name[pos - 1] || ';' == name[pos - 1])){
			const auto start = pos + pat.size();
			return name.substr(start, name.find_first_of(";]", start) - start);
		}
		pos = name.find(pat, pos + 1);
	}
	return "";
}

/**@brief The barcode numbers of each sample, a file with the columns barcode,sample with a number in the barcode's name, a header is fine
 */
std::map<std::string, std::set<int32_t>> readBarcodeKey(const bfs::path & keyFnp){
	std::map<std::string, std::set<int32_t>> ret;
	InputStream key{InOptions{keyFnp}};
	std::string line;
	while(njh::files::crossPlatGetline(key, line)){
		const auto toks = tokenizeString(line, "\t");
		if(toks.size() < 2){
			continue;
		}
		std::string digits;
		for(const auto c : toks[0]){
			if(std::isdigit(static_cast<unsigned char>(c))){
				digits.push_back(c);
			}
		}
		if(!digits.empty()){
			ret[toks[1]].emplace(std::stoi(digits));
		}
	}
	return ret;
}

} // namespace

int SeekDeepUtilsRunner::nanoporeAmpliconReadQC(const njh::progutils::CmdArgs & inputCommands){
	NanoporeAmpliconReadQC::Pars pars;
	bfs::path primersFnp;
	bfs::path targetInfoDir;
	bfs::path refSeqsDir;
	bfs::path fastqDir;
	std::string sampleName;
	uint32_t numThreads = 1;
	bool includeUndetermined = false;
	bool noFoldbackReads = false;
	bool noSalvageReads = false;
	bool writeReadCategories = false;
	NanoporeReadPlacer::Pars placerPars;
	bfs::path hostGenomeFnp;
	uint64_t maxCappedReads = 500000;
	bfs::path mismatchedSamplesFastq;
	bfs::path barcodesFnp;
	uint32_t barcodeStart = 0;
	uint32_t barcodeLen = 0;
	bfs::path barcodeKeyFnp;
	uint32_t minOffTargetLocusReads = 5;
	uint32_t offTargetLocusBinSize = 1000;

	seqSetUp setUp(inputCommands);
	setUp.description_ = "Report what the nanopore reads of targeted amplicon samples are, complete amplicons, joined amplicons, folds, and off target products, directly from the reads without extraction or clustering. Reads should have adapters and barcodes trimmed";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastqgz sample.fastq.gz --sampleName sample --primers primers.tsv --targetInfoDir genome_extraction --numThreads 8 --dout qcOut");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastqDir demux/fastq/ --primers primers.tsv --refSeqsDir refSeqs/ --numThreads 8 --dout qcOut");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(VecStr{"--fastq", "--fastqgz", "--fasta", "--fastagz"}, false);
	setUp.setOption(sampleName, "--sampleName", "Name of the sample when giving just one, by default the name of the input file", false, "Input");
	setUp.setOption(fastqDir, "--fastqDir", "A directory of samples, a file for each, the file's name is the sample's name", false, "Input");
	setUp.setOption(includeUndetermined, "--includeUndetermined", "Include files in --fastqDir that start with Undetermined, by default they're skipped since they aren't a sample, only --maxCappedReads of their reads are looked at", false, "Input");
	setUp.setOption(mismatchedSamplesFastq, "--mismatchedSamplesFastq", "The file of reads the demultiplexing found with different samples at each end, only --maxCappedReads of them are looked at", false, "Input");
	setUp.setOption(maxCappedReads, "--maxCappedReads", "At most about this many reads of inputs that aren't a sample, chosen by a hash of the read's name so the same reads are chosen every time", false, "Input");
	setUp.setOption(primersFnp, "--primers", "A file with the columns target,forward primer,reverse primer", true, "Input");
	setUp.setOption(targetInfoDir, "--targetInfoDir", "The output directory of SeekDeep genTargetInfoFromGenomes, for the targets' reference sequences", false, "Input");
	setUp.setOption(refSeqsDir, "--refSeqsDir", "A directory of TARGET.fasta reference sequences without their primers, instead of --targetInfoDir", false, "Input");
	setUp.setOption(numThreads, "--numThreads", "Number of threads, the output is the same no matter the number of threads");

	const std::string ampliconGroup = "Complete Amplicons";
	setUp.setOption(pars.maxPrimerError_, "--maxPrimerError", "The fraction of a primer that can be edits, reads with joined amplicons are also counted at 0.15 and 0.10 to see how much depends on this", false, ampliconGroup);
	setUp.setOption(pars.minInnerKmerFrac_, "--minInnerKmerFrac", "The fraction of kmers between the primers that have to be in the target's reference sequences", false, ampliconGroup);
	setUp.setOption(pars.minLenTolerance_, "--minLenTolerance", "How far off the expected length a complete amplicon can be, at least", false, ampliconGroup);
	setUp.setOption(pars.lenToleranceFrac_, "--lenToleranceFrac", "How far off the expected length a complete amplicon can be as a fraction of it, whichever of this and --minLenTolerance is larger", false, ampliconGroup);
	setUp.setOption(pars.otherBasesCutOff_, "--otherBasesCutOff", "The most bases outside a read's only complete amplicon before it's ampliconPlusOther", false, ampliconGroup);
	setUp.setOption(pars.decoyEvery_, "--decoyEvery", "Also search for shuffled primers on 1 in this many reads to see how often the search finds something by chance, 0 to skip", false, ampliconGroup);

	const std::string otherGroup = "Other Categories";
	setUp.setOption(pars.minTargetKmerFrac_, "--minTargetKmerFrac", "The fraction of a read or fold arm's kmers that have to be a target's to say it's that target", false, otherGroup);
	setUp.setOption(pars.otherTargetKmerFrac_, "--otherTargetKmerFrac", "The fraction of a read's kmers that have to be a second target's to say it holds both", false, otherGroup);
	setUp.setOption(pars.maxDimerInsert_, "--maxDimerInsert", "The most bases between two primers on a read for it to be a dimer", false, otherGroup);

	const std::string foldGroup = "Folds";
	setUp.setOption(pars.foldPars_.kmerLen_, "--kmerLen", "The size of the kmers that vote for a fold", false, foldGroup);
	setUp.setOption(pars.foldPars_.minArmLen_, "--minArmLen", "The shorter arm of a fold has to be at least this long", false, foldGroup);
	setUp.setOption(pars.foldPars_.minAnchors_, "--minAnchors", "The fewest kmers that have to pair between the arms", false, foldGroup);
	setUp.setOption(pars.foldPars_.minIdentity_, "--minIdentity", "The lowest identity between the arms", false, foldGroup);
	setUp.setOption(pars.foldPars_.minAnchorsPerBase_, "--minAnchorsPerBase", "The fewest kmers that have to pair per base of the anchored part of the arm", false, foldGroup);
	setUp.setOption(pars.turnBinSize_, "--turnBinSize", "The bin size for where folds turn along an amplicon", false, foldGroup);
	setUp.setOption(pars.primerAtStartWithin_, "--primerAtStartWithin", "A read's first primer has to start within this many bases to place its fold on the amplicon", false, foldGroup);

	const std::string barcodeGroup = "Barcodes";
	setUp.setOption(barcodesFnp, "--barcodes", "A fasta of the sample barcodes to look for at the junctions of joined amplicons, another sample's barcode there means the amplicons can't be said to belong to this sample", false, barcodeGroup);
	setUp.setOption(barcodeStart, "--barcodeStart", "Where the barcode starts in the sequences in --barcodes, they can have flanking sequence", false, barcodeGroup);
	setUp.setOption(barcodeLen, "--barcodeLen", "The length of the barcode in the sequences in --barcodes", false, barcodeGroup);
	setUp.setOption(barcodeKeyFnp, "--barcodeKey", "A file with the columns barcode,sample, the number in the barcode's name is the barcode's number", false, barcodeGroup);

	const std::string placementGroup = "Placing Reads";
	setUp.setOption(placerPars.mappingIndex_, "--mappingIndex", "A minimap2 index of the target genome and the host, to say where reads and folds are, reads are only placed if this is given", false, placementGroup);
	setUp.setOption(placerPars.targetGenomeFnp_, "--targetGenome", "A fasta of the target genome, says which contigs of --mappingIndex are the target's, the rest are host", false, placementGroup);
	setUp.setOption(placerPars.targetBedFnp_, "--targetBed", "Where the panel is on the target genome, the columns chrom,start,end,name with the primers in the coordinates", false, placementGroup);
	setUp.setOption(minOffTargetLocusReads, "--minOffTargetLocusReads", "Only write the places twoPrimerOffTarget reads map to with at least this many reads over all samples, for offTargetPrimerPairLoci.tsv.gz", false, placementGroup);
	setUp.setOption(offTargetLocusBinSize, "--offTargetLocusBinSize", "The bin size for where twoPrimerOffTarget reads map to, by the start of the read's longest alignment", false, placementGroup);

	const std::string outputGroup = "Output";
	setUp.setOption(noFoldbackReads, "--noFoldbackReads", "Don't write foldbackReads.tsv.gz, every read with a fold, can be most of a sample", false, outputGroup);
	setUp.setOption(noSalvageReads, "--noSalvageReads", "Don't write salvageReads.tsv.gz, where the complete amplicons are on the reads with 2 or more", false, outputGroup);
	setUp.setOption(writeReadCategories, "--writeReadCategories", "Write readCategories.tsv.gz, a line for every read", false, outputGroup);
	setUp.processDirectoryOutputName("nanoporeAmpliconReadQC_TODAY", true);

	const uint32_t inputsGiven = (setUp.pars_.ioOptions_.firstName_.empty() ? 0 : 1) + (fastqDir.empty() ? 0 : 1);
	if(1 != inputsGiven){
		setUp.failed_ = true;
		setUp.addWarning("Need to give exactly one of --fastq, --fastqgz, --fasta, --fastagz or --fastqDir");
	}
	if(targetInfoDir.empty() == refSeqsDir.empty()){
		setUp.failed_ = true;
		setUp.addWarning("Need to give exactly one of --targetInfoDir or --refSeqsDir, the reference sequences say a target's expected length and sequence");
	}
	if(pars.foldPars_.kmerLen_ < 5 || pars.foldPars_.kmerLen_ > 31){
		setUp.failed_ = true;
		setUp.addWarning("--kmerLen has to be between 5 and 31");
	}
	if(!placerPars.mappingIndex_.empty() && (placerPars.targetGenomeFnp_.empty() || placerPars.targetBedFnp_.empty())){
		setUp.failed_ = true;
		setUp.addWarning("Need to give --targetGenome and --targetBed with --mappingIndex");
	}
	if(!barcodesFnp.empty() && (0 == barcodeLen || barcodeKeyFnp.empty())){
		setUp.failed_ = true;
		setUp.addWarning("Need to give --barcodeLen and --barcodeKey with --barcodes");
	}
	if(0 == maxCappedReads){
		setUp.failed_ = true;
		setUp.addWarning("--maxCappedReads has to be above 0");
	}
	if(0 == pars.turnBinSize_){
		setUp.failed_ = true;
		setUp.addWarning("--turnBinSize has to be above 0");
	}
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

	std::vector<SampleInput> samples;
	if(!fastqDir.empty()){
		samples = readInSamplesFromDir(fastqDir, includeUndetermined);
	}else{
		samples.emplace_back(SampleInput{sampleName.empty() ? sampleNameFromFile(setUp.pars_.ioOptions_.firstName_) : sampleName, setUp.pars_.ioOptions_});
	}
	if(samples.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no samples found in input" << "\n";
		throw std::runtime_error{ss.str()};
	}

	PrimersAndMids ids(primersFnp);
	if(ids.getTargets().empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << primersFnp << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.addRefSeqs(refSeqsDir.empty() ? njh::files::make_path(targetInfoDir, "forSeekDeep", "refSeqs") : refSeqsDir);
	VecStr targetsWithoutRefs;
	for(const auto & target : ids.targets_){
		if(target.second.refs_.empty()){
			targetsWithoutRefs.emplace_back(target.first);
		}
	}
	if(!targetsWithoutRefs.empty()){
		std::cerr << "warning, no reference sequences for " << njh::conToStr(targetsWithoutRefs, ", ") << ", no complete amplicons can be called for them" << std::endl;
	}

	if(!mismatchedSamplesFastq.empty()){
		samples.emplace_back(SampleInput{sampleNameFromFile(mismatchedSamplesFastq), njh::endsWith(mismatchedSamplesFastq.string(), ".gz") ? SeqIOOptions::genFastqInGz(mismatchedSamplesFastq) : SeqIOOptions::genFastqIn(mismatchedSamplesFastq), true});
	}

	NanoporeAmpliconReadQC qc(ids, pars);
	std::map<std::string, std::set<int32_t>> barcodeKey;
	if(!barcodesFnp.empty()){
		qc.setBarcodes(barcodesFnp, barcodeStart, barcodeLen);
		barcodeKey = readBarcodeKey(barcodeKeyFnp);
	}
	const auto & targetNames = qc.targetNames();
	std::unique_ptr<NanoporeReadPlacer> placer;
	bfs::path mappingDir;
	if(!placerPars.mappingIndex_.empty()){
		placerPars.numThreads_ = numThreads;
		placer = std::make_unique<NanoporeReadPlacer>(placerPars);
		mappingDir = njh::files::make_path(setUp.pars_.directoryName_, "mapping");
		njh::files::makeDir(njh::files::MkdirPar{mappingDir});
	}

	const auto outFnp = [&setUp](const std::string & name){
		return OutOptions(njh::files::make_path(setUp.pars_.directoryName_, name));
	};
	OutputStream categoryOut(outFnp("categoryCounts.tsv.gz"));
	categoryOut << "sample\tcategory\treads\tpercReads\tbases\tpercBases\n";
	OutputStream targetCategoryOut(outFnp("targetCategoryCounts.tsv.gz"));
	targetCategoryOut << "sample\ttarget\tcategory\treads\n";
	OutputStream ampliconsOut(outFnp("targetAmplicons.tsv.gz"));
	ampliconsOut << "sample\ttarget\texpectedLen\tsingleAmpliconReads\tunitsInJoinedReads\tunitsCountedOnce\tusableAmplicons\n";
	OutputStream salvageOut(outFnp("salvageCounts.tsv.gz"));
	salvageOut << "sample\treads\treadsWith2+\tpercReadsWith2+\twith2\twith3\twith4\twith5+";
	for(const auto & primerError : qc.stricterPrimerErrors_){
		salvageOut << "\treadsWith2+AtPrimerError" << primerError;
	}
	salvageOut << "\treadsWithForeignBarcode\tsingleAmpliconReads\tunitsInJoinedReads\tunitsCountedOnce\n";
	OutputStream junctionsOut(outFnp("junctions.tsv.gz"));
	junctionsOut << "sample\ttargets\torientation\tgapBin\tbarcode\tjunctions\tmedianGap\n";
	OutputStream pairStatusOut(outFnp("pairStatusCounts.tsv.gz"));
	pairStatusOut << "sample\tpairStatus\tcategory\treads\n";
	OutputStream cappedOut(outFnp("cappedInputs.tsv.gz"));
	cappedOut << "sample\ttotalReads\tlookedAtOneIn\treadsLookedAt\n";
	OutputStream ratesOut(outFnp("foldbackTargetRates.tsv.gz"));
	ratesOut << "sample\ttarget\ttargetReads\tfolds\tpercFolds\tfoldsPlaced\tfoldsStartedOnForwardPrimer\tpercMainStrand\tmedianLoop\tmedianShorterArm\n";
	OutputStream pileupOut(outFnp("foldbackTurnPileup.tsv.gz"));
	pileupOut << "sample\ttarget\treadStartedWith\tampliconPosBinStart\treads\n";
	OutputStream dimersOut(outFnp("primerDimers.tsv.gz"));
	dimersOut << "sample\tfirstPrimer\tsecondPrimer\treads\n";
	OutputStream offTargetPairsOut(outFnp("offTargetPrimerPairs.tsv.gz"));
	offTargetPairsOut << "sample\tprimer1\tprimer2\tpairType\treads\tpercOfTwoPrimerOffTarget\tpercOfSampleReads\tbases\tmedianReadLen\tmedianInsert\treadsBothPrimersAtEnds\n";
	OutputStream offTargetPrimersOut(outFnp("offTargetPrimers.tsv.gz"));
	offTargetPrimersOut << "sample\tprimer\ttwoPrimerOffTargetReads\tpercOfTwoPrimerOffTarget\tonePrimerOffTargetReads\tpercOfOnePrimerOffTarget\n";
	OutputStream calibrationOut(outFnp("primerCalibration.tsv.gz"));
	calibrationOut << "sample\tprimers\tprimerError\treadsChecked\tprimerHitsPerRead\treadsWith1+Units\treadsWith2+Units\n";
	std::unique_ptr<OutputStream> foldReadsOut;
	if(!noFoldbackReads){
		foldReadsOut = std::make_unique<OutputStream>(outFnp("foldbackReads.tsv.gz"));
		*foldReadsOut << "sample\tread\tlength\tcategory\ttarget\tleftStart\tleftEnd\trightStart\trightEnd\tleftArmLen\trightArmLen\tloopLen\tanchors\tanchorsPerBase\totherVotes\tarmIdentity\treadStartedWith\tturnAmpliconPos\n";
	}
	std::unique_ptr<OutputStream> salvageReadsOut;
	if(!noSalvageReads){
		salvageReadsOut = std::make_unique<OutputStream>(outFnp("salvageReads.tsv.gz"));
		*salvageReadsOut << "sample\tread\tlength\tcategory\tnUnits\tunits\tjunctions\tbefore\tafter\n";
	}

	std::unique_ptr<OutputStream> readCategoriesOut;
	if(writeReadCategories){
		readCategoriesOut = std::make_unique<OutputStream>(outFnp("readCategories.tsv.gz"));
		*readCategoriesOut << "sample\tread\tlength\tcategory\ttarget\tnUnits\tprimerHits\tfoldCalled\n";
	}

	std::unique_ptr<OutputStream> placementOut;
	std::unique_ptr<OutputStream> foldPlacementCountsOut;
	std::unique_ptr<OutputStream> targetTurnsOut;
	std::unique_ptr<OutputStream> foldPlacementsOut;
	std::unique_ptr<OutputStream> offTargetPairPlacementOut;
	if(placer){
		offTargetPairPlacementOut = std::make_unique<OutputStream>(outFnp("offTargetPrimerPairPlacement.tsv.gz"));
		*offTargetPairPlacementOut << "sample\tprimer1\tprimer2\tpairType\tplacement\treads\treadsTwoLoci\n";
		placementOut = std::make_unique<OutputStream>(outFnp("placementCounts.tsv.gz"));
		*placementOut << "sample\tcategory\tplacement\treads\tbases\treadsTwoLoci\treadsWith2+WholeUnits\n";
		foldPlacementCountsOut = std::make_unique<OutputStream>(outFnp("foldPlacementCounts.tsv.gz"));
		*foldPlacementCountsOut << "sample\tstatus\torganism\tturnClass\tloopSide\tfolds\tmedianReach\n";
		targetTurnsOut = std::make_unique<OutputStream>(outFnp("foldbackGenomeTurnPileup.tsv.gz"));
		*targetTurnsOut << "sample\ttarget\textendStrand\tgenomeBinStart\treads\n";
		if(!noFoldbackReads){
			foldPlacementsOut = std::make_unique<OutputStream>(outFnp("foldPlacements.tsv.gz"));
			*foldPlacementsOut << "sample\tread\tlength\tleftStart\tleftEnd\trightStart\trightEnd\tloopLen\tarmIdentity\tstatus\torganism\tchrom\tmapqL\tmapqR\textendStrand\tfar\tentry\treach\tentryL\tentryR\tqgap\tloopSide\tspanLo\tspanHi\ttarget\ttargetStart\ttargetEnd\tfarToEnd\tentryToEnd\tturnClass\n";
		}
	}
	//all the samples go through minimap2 together at the end since loading the index takes a while with a host in it
	struct ReadForPlacing{
		uint32_t sample_{0};
		NanoporeAmpliconReadQC::ReadCategory category_{NanoporeAmpliconReadQC::ReadCategory::NO_PRIMER};
		uint32_t len_{0};
		bool seen_{false};
		uint32_t offTargetPair_{std::numeric_limits<uint32_t>::max()}; /**< primer1 * number of primers + primer2 for twoPrimerOffTarget */
	};
	const uint32_t numPrimers = targetNames.size() * 2;
	struct FoldForPlacing{
		uint32_t sample_{0};
		NanoporeReadPlacer::FoldInput input_;
	};
	std::unordered_map<std::string, ReadForPlacing> readsForPlacing;
	std::unordered_map<std::string, FoldForPlacing> foldsForPlacing;
	std::vector<bfs::path> seqsToMap;
	std::unique_ptr<OutputStream> segmentsOut;
	std::unique_ptr<OutputStream> cappedSeqsOut;
	const auto segmentsFnp = placer ? njh::files::make_path(mappingDir, "foldSides.fasta") : bfs::path{};
	const auto cappedSeqsFnp = placer ? njh::files::make_path(mappingDir, "readsOfCappedInputs.fasta") : bfs::path{};
	if(placer){
		segmentsOut = std::make_unique<OutputStream>(OutOptions(segmentsFnp));
	}
	uint64_t readsInMoreThanOneInput = 0;
	//off panel folds are spread over the host genome, it's only worth saying which places they pile up at over all the samples
	std::map<std::tuple<std::string, std::string, int64_t>, std::pair<uint64_t, std::set<std::string>>> offPanelTurns;

	for(uint32_t samplePos = 0; samplePos < samples.size(); ++samplePos){
		const auto & sample = samples[samplePos];
		if(setUp.pars_.verbose_){
			std::cout << sample.name_ << std::endl;
		}
		if(placer){
			if(!sample.capped_){
				seqsToMap.emplace_back(sample.opts_.firstName_);
			}else if(!cappedSeqsOut){
				cappedSeqsOut = std::make_unique<OutputStream>(OutOptions(cappedSeqsFnp));
				seqsToMap.emplace_back(cappedSeqsFnp);
			}
		}
		uint64_t keepEvery = 1;
		uint64_t totalInInput = 0;
		if(sample.capped_){
			SeqInput counter(sample.opts_);
			counter.openIn();
			seqInfo seq;
			while(counter.readNextRead(seq)){
				++totalInInput;
			}
			keepEvery = (totalInInput + maxCappedReads - 1) / maxCappedReads;
			keepEvery = std::max<uint64_t>(keepEvery, 1);
		}
		const auto sampleBarcodesFound = barcodeKey.find(sample.name_);
		const std::set<int32_t> sampleBarcodes = sampleBarcodesFound == barcodeKey.end() ? std::set<int32_t>{} : sampleBarcodesFound->second;
		//a read that isn't in a sample's file says which samples it could be
		const auto ownBarcodesForRead = [&](const std::string & name){
			if(!sample.capped_){
				return sampleBarcodes;
			}
			std::set<int32_t> ret;
			for(const auto & tag : VecStr{"fwdSample", "revSample", "frontSample", "backSample"}){
				const auto found = barcodeKey.find(getNameTag(name, tag));
				if(found != barcodeKey.end()){
					ret.insert(found->second.begin(), found->second.end());
				}
			}
			return ret;
		};
		std::map<std::pair<std::string, NanoporeAmpliconReadQC::ReadCategory>, uint64_t> pairStatusCounts;
		NanoporeAmpliconReadQC::SampleTallies tallies;
		tallies.sampleName_ = sample.name_;
		SeqInput reader(sample.opts_);
		reader.openIn();
		//the reads are worked on a batch at a time and written serially so that the order doesn't depend on the threads
		const size_t batchSize = 20000;
		std::vector<seqInfo> batch;
		std::vector<NanoporeAmpliconReadQC::ReadResult> results;
		bool moreReads = true;
		uint64_t readsSeen = 0;
		while(moreReads){
			batch.clear();
			seqInfo seq;
			while(batch.size() < batchSize){
				if(!reader.readNextRead(seq)){
					moreReads = false;
					break;
				}
				if(keepEvery > 1 && 0 != hashReadName(readNameOnly(seq.name_)) % keepEvery){
					continue;
				}
				batch.emplace_back(seq);
			}
			if(batch.empty()){
				break;
			}
			results.assign(batch.size(), NanoporeAmpliconReadQC::ReadResult{});
			std::atomic<size_t> nextRead{0};
			std::function<void()> processReads = [&](){
				size_t pos = 0;
				while((pos = nextRead++) < batch.size()){
					results[pos] = qc.processRead(batch[pos].seq_, pars.decoyEvery_ > 0 && 0 == (readsSeen + pos) % pars.decoyEvery_);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(processReads, std::max<uint32_t>(numThreads, 1));
			for(size_t pos = 0; pos < batch.size(); ++pos){
				const auto & result = results[pos];
				//the demultiplexing programs add information to the name, leave only the read's name
				const auto & name = batch[pos].name_;
				const auto readName = readNameOnly(name);
				const auto ownBarcodes = ownBarcodesForRead(name);
				qc.addToTallies(result, tallies, ownBarcodes);
				{
					auto pairStatus = getNameTag(name, "pairStatus");
					if(pairStatus.empty()){
						pairStatus = getNameTag(name, "sample").empty() ? "noTag" : "bothEndsCalled";
					}
					++pairStatusCounts[{pairStatus, result.category_}];
				}
				if(placer){
					ReadForPlacing forPlacing{samplePos, result.category_, static_cast<uint32_t>(result.len_), false};
					if(NanoporeAmpliconReadQC::ReadCategory::TWO_PRIMER_OFF_TARGET == result.category_){
						forPlacing.offTargetPair_ = result.offTargetPrimer1_ * numPrimers + result.offTargetPrimer2_;
					}
					if(!readsForPlacing.emplace(readName, forPlacing).second){
						++readsInMoreThanOneInput;
					}
					if(sample.capped_){
						*cappedSeqsOut << ">" << readName << "\n" << batch[pos].seq_ << "\n";
					}
					if(result.fold_.called_){
						foldsForPlacing[readName] = FoldForPlacing{samplePos, NanoporeReadPlacer::FoldInput{static_cast<uint32_t>(result.len_), result.fold_}};
						*segmentsOut << ">" << readName << "|L\n" << batch[pos].seq_.substr(0, result.fold_.rightStart_) << "\n";
						*segmentsOut << ">" << readName << "|R\n" << batch[pos].seq_.substr(result.fold_.leftEnd_) << "\n";
					}
				}
				const auto category = NanoporeAmpliconReadQC::getCategoryStr(result.category_);
				if(readCategoriesOut){
					*readCategoriesOut << sample.name_
							<< "\t" << readName
							<< "\t" << result.len_
							<< "\t" << category
							<< "\t" << (result.target_ >= 0 ? targetNames[result.target_] : std::string{})
							<< "\t" << result.units_.size()
							<< "\t" << result.primerHits_
							<< "\t" << (result.fold_.called_ ? 1 : 0) << "\n";
				}
				if(result.units_.size() >= 2 && salvageReadsOut){
					*salvageReadsOut << sample.name_
							<< "\t" << readName
							<< "\t" << result.len_
							<< "\t" << category
							<< "\t" << result.units_.size()
							<< "\t" << qc.unitsStr(result)
							<< "\t" << qc.junctionsStr(result, ownBarcodes)
							<< "\t" << result.units_.front().start_
							<< "\t" << result.len_ - result.units_.back().end_ << "\n";
				}
				if(result.fold_.called_ && foldReadsOut){
					uint32_t turnPos = 0;
					const bool placed = qc.foldTurnPos(result, turnPos);
					const auto & fold = result.fold_;
					*foldReadsOut << sample.name_
							<< "\t" << readName
							<< "\t" << result.len_
							<< "\t" << category
							<< "\t" << (result.target_ >= 0 ? targetNames[result.target_] : std::string{})
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
							<< "\t" << fold.identity_
							<< "\t" << (placed ? (qc.probe(result.startProbe_).forwardPrimer_ ? "F" : "R") : "")
							<< "\t";
					if(placed){
						*foldReadsOut << turnPos;
					}
					*foldReadsOut << "\n";
				}
			}
			readsSeen += batch.size();
		}

		if(sample.capped_){
			cappedOut << sample.name_ << "\t" << totalInInput << "\t" << keepEvery << "\t" << readsSeen << "\n";
		}
		for(const auto & count : pairStatusCounts){
			pairStatusOut << sample.name_ << "\t" << count.first.first << "\t" << NanoporeAmpliconReadQC::getCategoryStr(count.first.second) << "\t" << count.second << "\n";
		}
		const auto perc = [](uint64_t num, uint64_t den){
			return 0 == den ? 0.0 : 100.0 * num / den;
		};
		for(const auto & category : NanoporeAmpliconReadQC::allCategories_){
			categoryOut << sample.name_
					<< "\t" << NanoporeAmpliconReadQC::getCategoryStr(category)
					<< "\t" << tallies.categoryReads_[category]
					<< "\t" << perc(tallies.categoryReads_[category], tallies.reads_)
					<< "\t" << tallies.categoryBases_[category]
					<< "\t" << perc(tallies.categoryBases_[category], tallies.bases_) << "\n";
		}
		for(const auto & count : tallies.targetCategoryReads_){
			targetCategoryOut << sample.name_
					<< "\t" << targetNames[count.first.first]
					<< "\t" << NanoporeAmpliconReadQC::getCategoryStr(count.first.second)
					<< "\t" << count.second << "\n";
		}
		uint64_t singles = 0;
		uint64_t unitsAll = 0;
		uint64_t unitsOnce = 0;
		for(uint32_t target = 0; target < targetNames.size(); ++target){
			const auto single = tallies.singleAmpliconReads_[target];
			const auto all = tallies.unitsInJoins_[target];
			const auto once = tallies.unitsCountedOnce_[target];
			singles += single;
			unitsAll += all;
			unitsOnce += once;
			ampliconsOut << sample.name_
					<< "\t" << targetNames[target]
					<< "\t" << qc.expectedLen(target)
					<< "\t" << single
					<< "\t" << all
					<< "\t" << once
					<< "\t" << single + once << "\n";
			const auto folds = tallies.targetFolds_[target];
			const auto forward = tallies.targetFoldStrand_[{target, true}];
			const auto reverse = tallies.targetFoldStrand_[{target, false}];
			auto & loops = tallies.targetFoldLoops_[target];
			auto & arms = tallies.targetFoldArms_[target];
			ratesOut << sample.name_
					<< "\t" << targetNames[target]
					<< "\t" << tallies.targetReads_[target]
					<< "\t" << folds
					<< "\t" << perc(folds, tallies.targetReads_[target])
					<< "\t" << forward + reverse
					<< "\t" << forward
					<< "\t" << perc(std::max(forward, reverse), forward + reverse)
					<< "\t" << (loops.empty() ? 0.0 : vectorMedianRef(loops))
					<< "\t" << (arms.empty() ? 0.0 : vectorMedianRef(arms)) << "\n";
		}
		uint64_t joinedReads = 0;
		for(const auto & count : tallies.readsWithUnits_){
			joinedReads += count.second;
		}
		salvageOut << sample.name_
				<< "\t" << tallies.reads_
				<< "\t" << joinedReads
				<< "\t" << perc(joinedReads, tallies.reads_)
				<< "\t" << tallies.readsWithUnits_[2]
				<< "\t" << tallies.readsWithUnits_[3]
				<< "\t" << tallies.readsWithUnits_[4]
				<< "\t" << tallies.readsWithUnits_[5];
		for(size_t pos = 0; pos < qc.stricterPrimerErrors_.size(); ++pos){
			salvageOut << "\t" << (tallies.joinedReadsAtStricterErrors_.empty() ? 0 : tallies.joinedReadsAtStricterErrors_[pos]);
		}
		salvageOut << "\t" << tallies.readsWithForeignBarcode_ << "\t" << singles << "\t" << unitsAll << "\t" << unitsOnce << "\n";
		for(auto & junction : tallies.junctions_){
			auto & gaps = tallies.junctionGaps_[junction.first];
			junctionsOut << sample.name_
					<< "\t" << (std::get<0>(junction.first) ? "sameTarget" : "diffTarget")
					<< "\t" << (std::get<1>(junction.first) ? "inverted" : "headToTail")
					<< "\t" << std::get<2>(junction.first)
					<< "\t" << std::get<3>(junction.first)
					<< "\t" << junction.second
					<< "\t" << vectorMedianRef(gaps) << "\n";
		}
		for(const auto & bin : tallies.turnPileup_){
			pileupOut << sample.name_
					<< "\t" << targetNames[std::get<0>(bin.first)]
					<< "\t" << (std::get<1>(bin.first) ? "F" : "R")
					<< "\t" << std::get<2>(bin.first)
					<< "\t" << bin.second << "\n";
		}
		for(const auto & dimer : tallies.dimers_){
			dimersOut << sample.name_ << "\t" << dimer.first.first << "\t" << dimer.first.second << "\t" << dimer.second << "\n";
		}
		{
			const auto twoPrimerReads = tallies.categoryReads_[NanoporeAmpliconReadQC::ReadCategory::TWO_PRIMER_OFF_TARGET];
			const auto onePrimerReads = tallies.categoryReads_[NanoporeAmpliconReadQC::ReadCategory::ONE_PRIMER_OFF_TARGET];
			//most reads first so the pairs that matter are at the top
			std::vector<std::pair<std::pair<uint32_t, uint32_t>, NanoporeAmpliconReadQC::SampleTallies::OffTargetPair *>> pairs;
			for(auto & pair : tallies.offTargetPairs_){
				pairs.emplace_back(pair.first, &pair.second);
			}
			std::stable_sort(pairs.begin(), pairs.end(), [](const auto & p1, const auto & p2){
				return p1.second->reads_ > p2.second->reads_;
			});
			for(const auto & pair : pairs){
				offTargetPairsOut << sample.name_
						<< "\t" << qc.primerName(pair.first.first)
						<< "\t" << qc.primerName(pair.first.second)
						<< "\t" << NanoporeAmpliconReadQC::primerPairType(pair.first.first, pair.first.second)
						<< "\t" << pair.second->reads_
						<< "\t" << perc(pair.second->reads_, twoPrimerReads)
						<< "\t" << perc(pair.second->reads_, tallies.reads_)
						<< "\t" << pair.second->bases_
						<< "\t" << vectorMedianRef(pair.second->lens_)
						<< "\t" << vectorMedianRef(pair.second->inserts_)
						<< "\t" << pair.second->atEnds_ << "\n";
			}
			for(const auto & primer : tallies.offTargetPrimerReads_){
				offTargetPrimersOut << sample.name_
						<< "\t" << qc.primerName(primer.first)
						<< "\t" << primer.second[0]
						<< "\t" << perc(primer.second[0], twoPrimerReads)
						<< "\t" << primer.second[1]
						<< "\t" << perc(primer.second[1], onePrimerReads) << "\n";
			}
		}
		for(const auto & calibration : tallies.calibration_){
			calibrationOut << sample.name_
					<< "\t" << (calibration.first.first ? "shuffled" : "real")
					<< "\t" << calibration.first.second
					<< "\t" << tallies.decoyReads_
					<< "\t" << (0 == tallies.decoyReads_ ? 0.0 : static_cast<double>(calibration.first.first ? tallies.decoyHits_ : tallies.realHitsOnDecoyReads_) / tallies.decoyReads_)
					<< "\t" << calibration.second[0]
					<< "\t" << calibration.second[1] << "\n";
		}
	}
	if(placer){
		if(readsInMoreThanOneInput > 0){
			std::cerr << "warning, " << readsInMoreThanOneInput << " reads have a name that's in more than one input, they're placed as the first one" << std::endl;
		}
		segmentsOut.reset();
		cappedSeqsOut.reset();
		struct PlacementCount{
			uint64_t reads_{0};
			uint64_t bases_{0};
			uint64_t twoLoci_{0};
			uint64_t wholeUnits_{0};
		};
		std::map<std::tuple<uint32_t, NanoporeAmpliconReadQC::ReadCategory, std::string>, PlacementCount> placementCounts;
		std::map<std::tuple<uint32_t, uint32_t, std::string>, std::pair<uint64_t, uint64_t>> pairPlacements;
		//where off target products map, over all samples, the reads of a sample come together so the samples are only added when they change
		std::unordered_map<std::string, uint32_t> chromIds;
		VecStr chromNames;
		std::map<std::tuple<uint32_t, uint32_t, int64_t>, std::pair<uint64_t, std::vector<uint32_t>>> pairLoci;
		placer->mapSeqs(seqsToMap, njh::files::make_path(mappingDir, "reads.log"), false,
				[&](const std::string & name, const std::vector<NanoporeReadPlacer::Alignment> & alns){
			const auto found = readsForPlacing.find(name);
			if(found == readsForPlacing.end() || found->second.seen_){
				return;
			}
			found->second.seen_ = true;
			const auto placement = placer->placeRead(alns);
			auto & count = placementCounts[std::make_tuple(found->second.sample_, found->second.category_, placement.placement_)];
			++count.reads_;
			count.bases_ += found->second.len_;
			count.twoLoci_ += placement.twoLoci_ ? 1 : 0;
			count.wholeUnits_ += placement.wholeUnits_ >= 2 ? 1 : 0;
			if(std::numeric_limits<uint32_t>::max() != found->second.offTargetPair_){
				auto & pairCount = pairPlacements[std::make_tuple(found->second.sample_, found->second.offTargetPair_, placement.placement_)];
				++pairCount.first;
				pairCount.second += placement.twoLoci_ ? 1 : 0;
				const NanoporeReadPlacer::Alignment * longest = nullptr;
				for(const auto & aln : alns){
					if(nullptr == longest || aln.matches_ > longest->matches_){
						longest = &aln;
					}
				}
				if(nullptr != longest && longest->matches_ >= placer->pars_.minAlnMatches_){
					auto chromId = chromIds.find(longest->chrom_);
					if(chromId == chromIds.end()){
						chromId = chromIds.emplace(longest->chrom_, chromNames.size()).first;
						chromNames.emplace_back(longest->chrom_);
					}
					auto & locus = pairLoci[std::make_tuple(found->second.offTargetPair_, chromId->second, static_cast<int64_t>(longest->rStart_ / offTargetLocusBinSize * offTargetLocusBinSize))];
					++locus.first;
					if(locus.second.empty() || locus.second.back() != found->second.sample_){
						locus.second.emplace_back(found->second.sample_);
					}
				}
			}
		});
		for(const auto & read : readsForPlacing){
			if(!read.second.seen_){
				auto & count = placementCounts[std::make_tuple(read.second.sample_, read.second.category_, std::string{"unplaced"})];
				++count.reads_;
				count.bases_ += read.second.len_;
				if(std::numeric_limits<uint32_t>::max() != read.second.offTargetPair_){
					++pairPlacements[std::make_tuple(read.second.sample_, read.second.offTargetPair_, std::string{"unplaced"})].first;
				}
			}
		}
		for(const auto & count : pairPlacements){
			const uint32_t primer1 = std::get<1>(count.first) / numPrimers;
			const uint32_t primer2 = std::get<1>(count.first) % numPrimers;
			*offTargetPairPlacementOut << samples[std::get<0>(count.first)].name_
					<< "\t" << qc.primerName(primer1)
					<< "\t" << qc.primerName(primer2)
					<< "\t" << NanoporeAmpliconReadQC::primerPairType(primer1, primer2)
					<< "\t" << std::get<2>(count.first)
					<< "\t" << count.second.first
					<< "\t" << count.second.second << "\n";
		}
		offTargetPairPlacementOut.reset();
		{
			OutputStream pairLociOut(outFnp("offTargetPrimerPairLoci.tsv.gz"));
			pairLociOut << "primer1\tprimer2\tpairType\tchrom\tbinStart\treads\tsamples\n";
			std::vector<std::pair<std::tuple<uint32_t, uint32_t, int64_t>, std::pair<uint64_t, std::vector<uint32_t>> *>> loci;
			for(auto & locus : pairLoci){
				if(locus.second.first >= minOffTargetLocusReads){
					loci.emplace_back(locus.first, &locus.second);
				}
			}
			std::stable_sort(loci.begin(), loci.end(), [](const auto & l1, const auto & l2){
				return l1.second->first > l2.second->first;
			});
			for(const auto & locus : loci){
				const uint32_t primer1 = std::get<0>(locus.first) / numPrimers;
				const uint32_t primer2 = std::get<0>(locus.first) % numPrimers;
				auto & sampleIds = locus.second->second;
				std::sort(sampleIds.begin(), sampleIds.end());
				sampleIds.erase(std::unique(sampleIds.begin(), sampleIds.end()), sampleIds.end());
				pairLociOut << qc.primerName(primer1)
						<< "\t" << qc.primerName(primer2)
						<< "\t" << NanoporeAmpliconReadQC::primerPairType(primer1, primer2)
						<< "\t" << chromNames[std::get<1>(locus.first)]
						<< "\t" << std::get<2>(locus.first)
						<< "\t" << locus.second->first
						<< "\t" << sampleIds.size() << "\n";
			}
		}
		pairLoci.clear();
		for(const auto & count : placementCounts){
			*placementOut << samples[std::get<0>(count.first)].name_
					<< "\t" << NanoporeAmpliconReadQC::getCategoryStr(std::get<1>(count.first))
					<< "\t" << std::get<2>(count.first)
					<< "\t" << count.second.reads_
					<< "\t" << count.second.bases_
					<< "\t" << count.second.twoLoci_
					<< "\t" << count.second.wholeUnits_ << "\n";
		}
		readsForPlacing.clear();

		std::map<std::tuple<uint32_t, std::string, std::string, std::string, std::string>, std::vector<int64_t>> foldCounts;
		std::map<std::tuple<uint32_t, std::string, char, int64_t>, uint64_t> targetTurns;
		std::unordered_map<std::string, std::vector<NanoporeReadPlacer::Alignment>> pendingLeft;
		const auto writeFold = [&](const std::string & name, const std::vector<NanoporeReadPlacer::Alignment> & leftAlns, const std::vector<NanoporeReadPlacer::Alignment> & rightAlns){
			const auto & forPlacing = foldsForPlacing.at(name);
			const auto & input = forPlacing.input_;
			const auto & sampleName = samples[forPlacing.sample_].name_;
			const auto place = placer->placeFold(input, leftAlns, rightAlns);
			foldCounts[std::make_tuple(forPlacing.sample_, place.status_, place.organism_, place.turnClass_, place.loopSide_)].emplace_back(place.reach_);
			if(place.placed_ && place.reach_ > 10 && !place.target_.empty() && "internal" == place.turnClass_){
				//the same size bins as along the amplicon so they can be compared
				++targetTurns[std::make_tuple(forPlacing.sample_, place.target_, place.extendStrand_, place.far_ / pars.turnBinSize_ * pars.turnBinSize_)];
			}
			if(place.placed_ && "offPanel" == place.turnClass_ && place.mapqLeft_ >= 1 && place.mapqRight_ >= 1){
				auto & bin = offPanelTurns[std::make_tuple(place.chrom_, std::string(1, place.extendStrand_), place.far_ / 100 * 100)];
				++bin.first;
				bin.second.emplace(sampleName);
			}
			if(foldPlacementsOut){
				const auto & fold = input.fold_;
				*foldPlacementsOut << sampleName
						<< "\t" << name
						<< "\t" << input.len_
						<< "\t" << fold.leftStart_ << "\t" << fold.leftEnd_ << "\t" << fold.rightStart_ << "\t" << fold.rightEnd_
						<< "\t" << fold.loopLen()
						<< "\t" << fold.identity_
						<< "\t" << place.status_
						<< "\t" << place.organism_
						<< "\t" << place.chrom_;
				if(place.placed_){
					*foldPlacementsOut << "\t" << place.mapqLeft_ << "\t" << place.mapqRight_
							<< "\t" << place.extendStrand_ << "\t" << place.far_ << "\t" << place.entry_ << "\t" << place.reach_
							<< "\t" << place.entryLeft_ << "\t" << place.entryRight_ << "\t" << place.queryGap_ << "\t" << place.loopSide_
							<< "\t" << place.spanLow_ << "\t" << place.spanHigh_ << "\t" << place.target_;
					if(!place.target_.empty()){
						*foldPlacementsOut << "\t" << place.targetStart_ << "\t" << place.targetEnd_ << "\t" << place.farToEnd_ << "\t" << place.entryToEnd_;
					}else{
						*foldPlacementsOut << "\t\t\t\t";
					}
					*foldPlacementsOut << "\t" << place.turnClass_;
				}else{
					*foldPlacementsOut << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
				}
				*foldPlacementsOut << "\n";
			}
		};
		//a fold's two sides come out one after the other, the same order they went in
		std::unordered_set<std::string> foldsPlaced;
		placer->mapSeqs(std::vector<bfs::path>{segmentsFnp}, njh::files::make_path(mappingDir, "foldSides.log"), true,
				[&](const std::string & name, const std::vector<NanoporeReadPlacer::Alignment> & alns){
			const auto read = name.substr(0, name.size() - 2);
			if('L' == name.back()){
				pendingLeft[read] = alns;
			}else{
				const auto left = pendingLeft.find(read);
				writeFold(read, left == pendingLeft.end() ? std::vector<NanoporeReadPlacer::Alignment>{} : left->second, alns);
				foldsPlaced.emplace(read);
				if(left != pendingLeft.end()){
					pendingLeft.erase(left);
				}
			}
		});
		//the sides with only a left alignment or with neither
		for(const auto & left : pendingLeft){
			writeFold(left.first, left.second, std::vector<NanoporeReadPlacer::Alignment>{});
			foldsPlaced.emplace(left.first);
		}
		for(const auto & fold : foldsForPlacing){
			if(!njh::in(fold.first, foldsPlaced)){
				writeFold(fold.first, std::vector<NanoporeReadPlacer::Alignment>{}, std::vector<NanoporeReadPlacer::Alignment>{});
			}
		}
		for(auto & count : foldCounts){
			*foldPlacementCountsOut << samples[std::get<0>(count.first)].name_
					<< "\t" << std::get<1>(count.first)
					<< "\t" << std::get<2>(count.first)
					<< "\t" << std::get<3>(count.first)
					<< "\t" << std::get<4>(count.first)
					<< "\t" << count.second.size()
					<< "\t" << vectorMedianRef(count.second) << "\n";
		}
		for(const auto & turn : targetTurns){
			*targetTurnsOut << samples[std::get<0>(turn.first)].name_ << "\t" << std::get<1>(turn.first) << "\t" << std::get<2>(turn.first) << "\t" << std::get<3>(turn.first) << "\t" << turn.second << "\n";
		}
		foldsForPlacing.clear();
		bfs::remove(segmentsFnp);
		if(bfs::exists(cappedSeqsFnp)){
			bfs::remove(cappedSeqsFnp);
		}
		OutputStream offPanelOut(outFnp("foldbackOffPanelPileup.tsv.gz"));
		offPanelOut << "chrom\textendStrand\tgenomeBinStart\treads\tsamples\n";
		for(const auto & bin : offPanelTurns){
			offPanelOut << std::get<0>(bin.first) << "\t" << std::get<1>(bin.first) << "\t" << std::get<2>(bin.first) << "\t" << bin.second.first << "\t" << bin.second.second.size() << "\n";
		}
	}
	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

}  // namespace njhseq
