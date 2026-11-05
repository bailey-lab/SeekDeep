/*
 * ampliconReadQC.cpp
 *
 *  Created on: Nov 4, 2026
 *      Author: nicholashathaway
 */

#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/TarAmpSetupUtils/AmpliconReadQC.hpp"
#include <njhseq/IO/OutputStream.hpp>
#include <atomic>
#include <regex>

namespace njhseq {

namespace {

struct SampleInput{
	std::string name_;
	bfs::path r1_;
	bfs::path r2_;

	[[nodiscard]] SeqIOOptions genOpts() const{
		return njh::endsWith(r1_.string(), ".gz") ? SeqIOOptions::genPairedInGz(r1_, r2_) : SeqIOOptions::genPairedIn(r1_, r2_);
	}
};

/**@brief Pair up the files in a directory, a file is a first mate if it ends in _R1_001, _R1, or _1 and the second mate
 * is named the same with a 2, what comes before that is the sample name
 */
std::vector<SampleInput> readInSamplesFromDir(const bfs::path & fastqDir){
	const std::regex firstMatePat{"^(.+?)_(R?)1(_001)?(\\.(fastq|fq)(\\.gz)?)$"};
	std::vector<SampleInput> ret;
	std::stringstream errors;
	for(const auto & file : njh::files::listAllFiles(fastqDir, false, VecStr{})){
		const auto fileName = file.first.filename().string();
		std::smatch match;
		if(std::regex_match(fileName, match, firstMatePat)){
			const auto mateFnp = njh::files::make_path(fastqDir, njh::pasteAsStr(match[1], "_", match[2], "2", match[3], match[4]));
			if(bfs::exists(mateFnp)){
				ret.emplace_back(SampleInput{match[1], file.first, mateFnp});
			}else{
				errors << "no mate file " << mateFnp << " for " << file.first << "\n";
			}
		}
	}
	if(!errors.str().empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << errors.str();
		throw std::runtime_error{ss.str()};
	}
	njh::sort(ret, [](const SampleInput & s1, const SampleInput & s2){
		return s1.name_ < s2.name_;
	});
	return ret;
}

std::vector<SampleInput> readInSamplesFromTable(const bfs::path & sampleTableFnp){
	table samplesTab(sampleTableFnp, "\t", true);
	samplesTab.checkForColumnsThrow(VecStr{"sample", "r1", "r2"}, __PRETTY_FUNCTION__);
	std::vector<SampleInput> ret;
	for(const auto & row : samplesTab.content_){
		ret.emplace_back(SampleInput{row[samplesTab.getColPos("sample")], row[samplesTab.getColPos("r1")], row[samplesTab.getColPos("r2")]});
	}
	return ret;
}

} // namespace

int SeekDeepUtilsRunner::ampliconReadQC(const njh::progutils::CmdArgs & inputCommands){
	AmpliconReadQC::QcPars pars;
	bfs::path primersFnp;
	bfs::path targetInfoDir;
	bfs::path overlapStatusFnp;
	bfs::path fastqDir;
	bfs::path sampleTableFnp;
	std::string sampleName;
	uint32_t numThreads = 1;
	bool noTrimLowQualWindows = false;
	bfs::path genomeDir;
	std::string primaryGenome;
	std::string selectedGenomesStr;
	bfs::path gffDir;
	uint32_t pairedEndLength = 0;
	uint32_t targetInfoErrors = 2;
	bfs::path refSeqsDir;
	bool noKmerComparison = false;
	bfs::path adaptersFnp;
	AmpliconProductPlacer::PlacerPars placerPars;

	seqSetUp setUp(inputCommands);
	setUp.description_ = "Report what the read pairs of targeted amplicon samples are, on their expected target or primer dimers and other off target products, directly from the reads without extraction or clustering";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq1gz sample_R1.fastq.gz --fastq2gz sample_R2.fastq.gz --sampleName sample --primers primers.tsv --targetInfoDir genome_extraction --dout qcOut");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastqDir fastqs/ --primers primers.tsv --overlapStatusFnp overlapStatuses.txt --numThreads 8 --dout qcOut");
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processReadInNames(seqSetUp::pairedReadInFormatsAvailable_, false);
	setUp.setOption(sampleName, "--sampleName", "Name of the sample when giving just one with the fastq options", false, "Input");
	setUp.setOption(fastqDir, "--fastqDir", "A directory of samples, first mates end in _R1_001, _R1 or _1 and the second mates are the same with a 2, what comes before that is the sample name", false, "Input");
	setUp.setOption(sampleTableFnp, "--sampleTable", "A tab delimited table of samples with the columns sample,r1,r2", false, "Input");
	setUp.setOption(primersFnp, "--primers", "A file with the columns target,forward primer,reverse primer", true, "Input");
	setUp.setOption(targetInfoDir, "--targetInfoDir", "The output directory of SeekDeep genTargetInfoFromGenomes, for how each target is expected to overlap, its reference sequences, and where it is, if this isn't given it gets created from --genomeDir", false, "Input");
	setUp.setOption(overlapStatusFnp, "--overlapStatusFnp", "A file with the columns target,status for how each target is expected to overlap, instead of --targetInfoDir", false, "Input");
	setUp.setOption(numThreads, "--numThreads", "Number of samples to process at the same time and threads for BLAST, the output is the same no matter the number of threads");

	const std::string genomeGroup = "Genomes";
	setUp.setOption(genomeDir, "--genomeDir", "A directory of the target genomes, named GENOME.fasta, to create --targetInfoDir if it's not given and to BLAST products that aren't what was expected against", false, genomeGroup);
	setUp.setOption(primaryGenome, "--primaryGenome", "The primary reference genome in --genomeDir", false, genomeGroup);
	setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Other genomes in --genomeDir to use, separated by commas, leave blank to use all of them", false, genomeGroup);
	setUp.setOption(gffDir, "--gffDir", "A directory with a gff for the genomes in --genomeDir, only used when creating --targetInfoDir", false, genomeGroup);
	setUp.setOption(pairedEndLength, "--pairedEndLength", "The length of the reads, only used when creating --targetInfoDir, by default it's taken from the reads", false, genomeGroup);
	setUp.setOption(targetInfoErrors, "--errors", "Number of errors to allow in primers when finding them in the genomes, only used when creating --targetInfoDir", false, genomeGroup);
	setUp.setOption(placerPars.hostGenomeFnp_, "--hostGenome", "A fasta of the host genome to BLAST products that aren't what was expected against", false, genomeGroup);

	const std::string kmerGroup = "Reference Kmer Comparison";
	setUp.setOption(refSeqsDir, "--refSeqsDir", "A directory of TARGET.fasta reference sequences to compare reads to, by default the ones in --targetInfoDir", false, kmerGroup);
	setUp.setOption(noKmerComparison, "--noKmerComparison", "Don't compare reads to the reference sequences of their target", false, kmerGroup);
	setUp.setOption(pars.compKmerLen_, "--compKmerLen", "The size of the kmers for comparing to the reference sequences", false, kmerGroup);
	setUp.setOption(pars.compKmerSimCutOff_, "--compKmerSimCutOff", "The fraction of kmers that have to be shared with a reference sequence", false, kmerGroup);

	const std::string placementGroup = "Placing Products";
	setUp.setOption(adaptersFnp, "--adapters", "A fasta of the adapter sequences to check products for instead of the TruSeq ones, a hit to adapter says nothing about where a sequence came from since many database records have some left in them", false, placementGroup);
	setUp.setOption(placerPars.adapterKmerLen_, "--adapterKmerLen", "How many bases in a row have to match an adapter", false, placementGroup);
	setUp.setOption(placerPars.minBlastLen_, "--minBlastLen", "Sequences shorter than this are too short to place with any confidence", false, placementGroup);
	setUp.setOption(placerPars.minPercentId_, "--blastMinPercentId", "The lowest percent identity for a BLAST hit to count", false, placementGroup);
	setUp.setOption(placerPars.minCoverage_, "--blastMinCoverage", "The lowest fraction of the sequence a BLAST hit has to cover to count", false, placementGroup);
	setUp.setOption(placerPars.targetRegionFlank_, "--targetRegionFlank", "How far from a target a BLAST hit can be and still count as being in that target", false, placementGroup);
	setUp.setOption(placerPars.longSeqDust_, "--blastDust", "The low complexity masking for BLASTing longer sequences, the blast default of \"20 64 1\" hides too much of an AT rich genome and turning it off with \"no\" is orders of magnitude slower", false, placementGroup);
	setUp.setOption(placerPars.maxSeqsToBlast_, "--maxSeqsToBlast", "Only the most common sequences across all the samples are BLASTed if there's more than this, the rest are reported as notBlasted", false, placementGroup);

	const std::string stitchGroup = "Stitching";
	setUp.setOption(pars.stitchPars_.errorAllowed_, "--overLapErrorAllowed", "The fraction of high quality mismatches to allow in the overlap when stitching, rounds down so an overlap under 1/this bases is allowed none, should be the same as what the extraction will use", false, stitchGroup);
	setUp.setOption(pars.dimerOverLapErrorAllowed_, "--dimerOverLapErrorAllowed", "The same as --overLapErrorAllowed but for a second attempt at stitching the pairs that failed, only used to find primer dimers", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.minOverlap_, "--minOverlap", "The minimal amount of overlap when stitching", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.hardMismatchCutOff_, "--hardMismatchCutOff", "A hard cut off for number of mismatches in the overlap when stitching", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.lqMismatchCutOff, "--lqMismatchCutOff", "A hard cut off for number of low quality mismatches in the overlap when stitching", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.hqMismatchCutOff, "--hqMismatchCutOff", "A hard cut off for number of high quality mismatches in the overlap when stitching", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.matchScore_, "--stitchMatchScore", "The score for a match when aligning the mates to stitch them", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.mismatchPenalty_, "--stitchMismatchPenalty", "The penalty for a mismatch when aligning the mates to stitch them, harsher than the match score so a long low identity placement of the mates doesn't out score a short true overlap", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.r1Trim_, "--r1Trim", "Remove this many bases off of the end of r1 reads", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.r2Trim_, "--r2Trim", "Remove this many bases off of the end of r2 reads", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.qualWindowPar_.avgQualCutOff_, "--qWindowTrimAvgQualCutOff", "Quality Window Trim Avg Qual Cut Off", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.qualWindowPar_.windowSize_, "--qWindowSize", "Quality Window Trim Size", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.qualWindowPar_.windowStep_, "--qWindowStep", "Quality Window Trim Step", false, stitchGroup);
	setUp.setOption(noTrimLowQualWindows, "--noTrimLowQualWindows", "Don't Trim Low Qual Windows", false, stitchGroup);
	setUp.setOption(pars.stitchPars_.stitchCacheMaxEntries_, "--stitchCacheMaxEntries", "Max number of mate pairs per sample to remember the alignment of so that they only have to be aligned once, 0 to not remember any, the output is the same either way", false, stitchGroup);
	pars.stitchPars_.trimLowQaulWindows_ = !noTrimLowQualWindows;

	const std::string dimerGroup = "Primer Dimers";
	setUp.setOption(pars.dimerPars_.maxGapBetweenPrimers_, "--maxGapBetweenPrimers", "The most bases allowed between the primers at each end of a read for it to be a dimer, the primers can also overlap or abut", false, dimerGroup);
	setUp.setOption(pars.dimerPars_.minPrimerMatchBases_, "--minPrimerMatchBases", "The fewest bases of a primer that have to match at a read's end for it to be trusted, 0 to work it out from the number of primers since the more primers there are the more likely a short match is by chance", false, dimerGroup);
	setUp.setOption(pars.dimerPars_.minCandidateMatchBases_, "--minCandidateMatchBases", "The fewest bases of a primer that have to match at a read's end for it to be kept as a possible partner of a trusted primer at the other end, which makes the read an unresolved dimer", false, dimerGroup);

	const std::string primerGroup = "Primers";
	setUp.setOption(pars.primerPars_.primerWithin_, "--primerWithinStart", "By default the primers are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives", false, primerGroup);
	setUp.setOption(pars.primerPars_.primerStart_, "--primerSearchStart", "By default the primers are searched at the very beginning of seq, use this flag to start the search here", false, primerGroup);
	setUp.setOption(pars.primerPars_.allowable_.hqMismatches_, "--primerNumOfMismatches", "The number of mismatches to allow in primers", false, primerGroup);
	setUp.setOption(pars.primerCallCacheMaxEntries_, "--primerCallCacheMaxEntries", "Max number of distinct read starts per sample to remember the primer found for so it's not searched for again, 0 to not remember any, the output is the same either way", false, primerGroup);
	bool useAlnPrimerSearch = false;
	setUp.setOption(useAlnPrimerSearch, "--useAlnPrimerSearch", "Use alignment to search for primers rather than motif search, motif search is faster but alignment search allows indels", false, primerGroup);
	pars.primerPars_.useMotif_ = !useAlnPrimerSearch;

	setUp.processDirectoryOutputName("ampliconReadQC_TODAY", true);

	const uint32_t inputsGiven = (setUp.pars_.ioOptions_.firstName_.empty() ? 0 : 1) + (fastqDir.empty() ? 0 : 1) + (sampleTableFnp.empty() ? 0 : 1);
	if(1 != inputsGiven){
		setUp.failed_ = true;
		setUp.addWarning("Need to give exactly one of --fastq1/--fastq2, --fastq1gz/--fastq2gz, --fastqDir, or --sampleTable");
	}
	if(!setUp.pars_.ioOptions_.firstName_.empty() && sampleName.empty()){
		setUp.failed_ = true;
		setUp.addWarning("Need to give --sampleName when giving a single sample");
	}
	if(!targetInfoDir.empty() && !overlapStatusFnp.empty()){
		setUp.failed_ = true;
		setUp.addWarning("Can't give both --targetInfoDir and --overlapStatusFnp");
	}
	if(targetInfoDir.empty() && overlapStatusFnp.empty() && genomeDir.empty()){
		setUp.failed_ = true;
		setUp.addWarning("Need to give one of --targetInfoDir, --overlapStatusFnp, or --genomeDir, it's needed to know how each target is expected to overlap");
	}
	if(!genomeDir.empty() && primaryGenome.empty()){
		setUp.failed_ = true;
		setUp.addWarning("Need to give --primaryGenome with --genomeDir");
	}
	if(!genomeDir.empty() && !overlapStatusFnp.empty()){
		setUp.failed_ = true;
		setUp.addWarning("Can't give --genomeDir with --overlapStatusFnp, BLAST hits need the target locations from --targetInfoDir");
	}
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

	std::vector<SampleInput> samples;
	if(!fastqDir.empty()){
		samples = readInSamplesFromDir(fastqDir);
	}else if(!sampleTableFnp.empty()){
		samples = readInSamplesFromTable(sampleTableFnp);
	}else{
		samples.emplace_back(SampleInput{sampleName, setUp.pars_.ioOptions_.firstName_, setUp.pars_.ioOptions_.secondName_});
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
	if(ids.containsMids()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << primersFnp << " has MIDs in it but MIDs aren't handled yet, the barcodes would stop the primers from being found" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();

	if(targetInfoDir.empty() && overlapStatusFnp.empty()){
		//nothing given about the targets, get it from the genomes
		targetInfoDir = njh::files::make_path(setUp.pars_.directoryName_, "genome_extraction");
		if(0 == pairedEndLength){
			uint64_t maxReadLen = 0;
			PairedRead pair;
			SeqInput reader(samples.front().genOpts());
			reader.openIn();
			uint32_t readCount = 0;
			while(readCount < 1000 && reader.readNextRead(pair)){
				readVec::getMaxLength(pair.seqBase_, maxReadLen);
				++readCount;
			}
			pairedEndLength = maxReadLen;
		}
		std::map<std::string, std::string> targetInfoArgs{
			{"--primers", primersFnp.string()},
			{"--genomedir", genomeDir.string()},
			{"--primarygenome", primaryGenome},
			{"--dout", targetInfoDir.string()},
			{"--pairedendlength", njh::pasteAsStr(pairedEndLength)},
			{"--useblast", ""},
			{"--errors", njh::pasteAsStr(targetInfoErrors)},
			{"--numthreads", njh::pasteAsStr(numThreads)}
		};
		if(!selectedGenomesStr.empty()){
			targetInfoArgs["--selectedgenomes"] = selectedGenomesStr;
		}
		if(!gffDir.empty()){
			targetInfoArgs["--gffdir"] = gffDir.string();
		}
		const njh::progutils::CmdArgs targetInfoCommands(inputCommands.masterProgram_, "genTargetInfoFromGenomes", targetInfoArgs, inputCommands.commandLine_, inputCommands.workingDir_);
		SeekDeepUtilsRunner::genTargetInfoFromGenomes(targetInfoCommands);
	}
	ids.addOverLapStatuses(overlapStatusFnp.empty() ?
			njh::files::make_path(targetInfoDir, "forSeekDeep", "overlapStatuses.txt") : overlapStatusFnp);
	if(refSeqsDir.empty() && !targetInfoDir.empty()){
		refSeqsDir = njh::files::make_path(targetInfoDir, "forSeekDeep", "refSeqs");
	}
	if(!noKmerComparison && !refSeqsDir.empty()){
		ids.addRefSeqs(refSeqsDir);
		ids.setRefSeqsKInfos(pars.compKmerLen_, true);
	}

	if(!adaptersFnp.empty()){
		placerPars.adapters_.clear();
		SeqInput adapterReader(SeqIOOptions::genFastaIn(adaptersFnp));
		for(const auto & adapter : adapterReader.readAllReads<seqInfo>()){
			placerPars.adapters_.emplace_back(adapter.seq_);
		}
	}
	placerPars.numThreads_ = numThreads;
	if(!genomeDir.empty()){
		//the same as genTargetInfoFromGenomes, no selected genomes means all of them, hits outside the primary genome can be what places strain specific sequence
		VecStr genomes{primaryGenome};
		VecStr otherGenomes;
		if(selectedGenomesStr.empty()){
			for(const auto & fastaFile : njh::files::listAllFiles(genomeDir, false, {std::regex{R"(.*\.fasta$)"}})){
				otherGenomes.emplace_back(bfs::basename(fastaFile.first));
			}
			njh::sort(otherGenomes);
		}else{
			otherGenomes = tokenizeString(selectedGenomesStr, ",");
		}
		for(const auto & genome : otherGenomes){
			if(!njh::in(genome, genomes)){
				genomes.emplace_back(genome);
			}
		}
		for(const auto & genome : genomes){
			//could be a --targetInfoDir that was made from fewer genomes
			if(!bfs::exists(njh::files::make_path(targetInfoDir, "locationsByGenome", genome + ".bed"))){
				std::cerr << __PRETTY_FUNCTION__ << ", warning no target locations for " << genome << " in " << targetInfoDir << ", not BLASTing against it" << "\n";
				continue;
			}
			placerPars.targetGenomes_.emplace_back(AmpliconProductPlacer::TargetGenome{genome,
				njh::files::make_path(genomeDir, genome + ".fasta"),
				njh::files::make_path(targetInfoDir, "locationsByGenome", genome + ".bed")});
		}
	}

	const AmpliconReadQC qc(ids, pars);
	const AmpliconProductPlacer placer(placerPars);

	std::vector<AmpliconReadQC::SampleResults> results(samples.size());
	std::atomic<size_t> nextSample{0};
	std::function<void()> runSamples = [&](){
		size_t sampleIdx = 0;
		while((sampleIdx = nextSample++) < samples.size()){
			if(setUp.pars_.verbose_){
				std::cout << "Processing " << samples[sampleIdx].name_ << std::endl;
			}
			results[sampleIdx] = qc.processSample(samples[sampleIdx].genOpts(), samples[sampleIdx].name_);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(runSamples, std::min<uint64_t>(std::max<uint32_t>(numThreads, 1), samples.size()));

	//each sequence only needs to be placed once no matter how many samples it's in
	std::unordered_map<std::string, uint64_t> seqsToPlace;
	for(const auto & res : results){
		for(const auto & product : res.shortProducts_){
			if(AmpliconReadQC::ReadPairCategory::SHORT_PRODUCT == product.second.category_ && !product.second.gapSeq_.empty()){
				seqsToPlace[product.second.gapSeq_] += product.second.count_;
			}
		}
		for(const auto & seqCount : res.kmerMismatchSeqs_){
			seqsToPlace[seqCount.first.second] += seqCount.second;
		}
		for(const auto & seqCount : res.mismatchedPrimerSeqs_){
			seqsToPlace[std::get<2>(seqCount.first)] += seqCount.second;
		}
	}
	const auto placementDir = njh::files::make_path(setUp.pars_.directoryName_, "placingProducts");
	njh::files::makeDirP(njh::files::MkdirPar{placementDir});
	const auto placements = placer.placeSeqs(seqsToPlace, placementDir);
	const auto writePlacement = [](std::ostream & out, const AmpliconProductPlacer::Placement & placement, const std::string & label){
		out << "\t" << label
				<< "\t" << placement.genome_
				<< "\t" << placement.chrom_
				<< "\t" << (placement.chrom_.empty() ? std::string("") : estd::to_string(placement.start_))
				<< "\t" << (placement.chrom_.empty() ? std::string("") : estd::to_string(placement.end_))
				<< "\t" << (placement.chrom_.empty() ? std::string("") : estd::to_string(placement.percentId_))
				<< "\t" << (placement.chrom_.empty() ? std::string("") : estd::to_string(placement.coverage_))
				<< "\t" << njh::conToStr(placement.overlappingTargets_, ";");
	};
	const std::string placementHeader = "\tplacement\thitGenome\thitChrom\thitStart\thitEnd\thitPercentId\thitCoverage\toverlappingTargets";
	const auto targetOfPrimer = [](const std::string & primerName){
		return primerName.substr(0, primerName.rfind('_'));
	};

	const auto outPath = [&setUp](const std::string & fileName){
		return OutOptions(njh::files::make_path(setUp.pars_.directoryName_, fileName));
	};
	OutputStream categoryOut(outPath("readCategoryCounts.tsv"));
	categoryOut << "sample\ttotal";
	for(const auto & category : AmpliconReadQC::allCategories_){
		categoryOut << "\t" << AmpliconReadQC::getCategoryStr(category);
	}
	for(const auto & category : AmpliconReadQC::allCategories_){
		categoryOut << "\t" << AmpliconReadQC::getCategoryStr(category) << "Perc";
	}
	categoryOut << "\tprimersMatched\tprimersMismatched\tprimersFailed\tdimersFromRetryStitch\toverLapErrorAllowed\tdimerOverLapErrorAllowed\n";

	OutputStream dimerOut(outPath("primerDimerCounts.tsv"));
	PrimerDimerTallies::writeHeader(dimerOut);

	OutputStream targetOut(outPath("targetCounts.tsv"));
	targetOut << "sample\ttarget\texpectedStatus\tstitchStatus\tcategory\tcount\tmedianLength\n";

	OutputStream otherOut(outPath("otherPrimerPairCounts.tsv"));
	otherOut << "sample\tforwardPrimer\treversePrimer\tstitchStatus\tcategory\tcount\tmedianLength\n";

	OutputStream shortOut(outPath("shortProducts.tsv"));
	shortOut << "sample\tcategory\tfrontPrimer\tendPrimer\tgapBetweenPrimers\tcount\tlength\tseq\tgapSeq" << placementHeader << "\n";

	OutputStream kmerMismatchOut(outPath("kmerMismatchProducts.tsv"));
	kmerMismatchOut << "sample\ttarget\tcount\tlength\tseq" << placementHeader << "\n";

	OutputStream mismatchedPrimerOut(outPath("mismatchedPrimerProducts.tsv"));
	mismatchedPrimerOut << "sample\tforwardPrimer\treversePrimer\tcount\tlength\tseq" << placementHeader << "\n";

	OutputStream placementCountsOut(outPath("placementCounts.tsv"));
	placementCountsOut << "sample\tcategory\tplacement\tcount\tperc\n";

	const auto medianLen = [](AmpliconReadQC::LengthTally & tally) -> std::string {
		return tally.lengths_.empty() ? std::string("") : estd::to_string(vectorMedianRef(tally.lengths_));
	};
	for(size_t sampleIdx = 0; sampleIdx < samples.size(); ++sampleIdx){
		auto & res = results[sampleIdx];
		uint64_t categorized = 0;
		categoryOut << res.sampleName_ << "\t" << res.total_;
		for(const auto & category : AmpliconReadQC::allCategories_){
			categorized += res.categoryCounts_[category];
			categoryOut << "\t" << res.categoryCounts_[category];
		}
		for(const auto & category : AmpliconReadQC::allCategories_){
			categoryOut << "\t" << (0 == res.total_ ? 0.0 : 100.0 * res.categoryCounts_[category] / res.total_);
		}
		categoryOut << "\t" << res.primersMatched_
				<< "\t" << res.primersMismatched_
				<< "\t" << res.primersFailed_
				<< "\t" << res.dimersFromRetryStitch_
				<< "\t" << pars.stitchPars_.errorAllowed_
				<< "\t" << pars.dimerOverLapErrorAllowed_ << std::endl;
		if(categorized != res.total_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error for " << res.sampleName_ << " categorized " << categorized << " read pairs out of " << res.total_ << "\n";
			throw std::runtime_error{ss.str()};
		}

		res.dimers_.write(dimerOut, samples[sampleIdx].r1_.filename().string(), res.sampleName_,
				[&res](const std::string &){ return res.total_; }, pars.dimerOverLapErrorAllowed_, qc.detector_.pars_);

		for(auto & target : res.targetCounts_){
			VecStr expected;
			for(const auto & status : ids.targets_.at(std::get<0>(target.first)).overlapStatuses_){
				expected.emplace_back(PairedReadProcessor::getOverlapStatusStr(status));
			}
			targetOut << res.sampleName_
					<< "\t" << std::get<0>(target.first)
					<< "\t" << njh::conToStr(expected, ",")
					<< "\t" << std::get<1>(target.first)
					<< "\t" << AmpliconReadQC::getCategoryStr(std::get<2>(target.first))
					<< "\t" << target.second.count_
					<< "\t" << medianLen(target.second) << std::endl;
		}
		for(auto & other : res.otherPrimerCounts_){
			otherOut << res.sampleName_
					<< "\t" << std::get<0>(other.first)
					<< "\t" << std::get<1>(other.first)
					<< "\t" << std::get<2>(other.first)
					<< "\t" << AmpliconReadQC::getCategoryStr(std::get<3>(other.first))
					<< "\t" << other.second.count_
					<< "\t" << medianLen(other.second) << std::endl;
		}
		std::map<std::pair<std::string, std::string>, uint64_t> placementCounts;
		for(const auto & product : res.shortProducts_){
			shortOut << res.sampleName_
					<< "\t" << AmpliconReadQC::getCategoryStr(product.second.category_)
					<< "\t" << product.second.frontPrimer_
					<< "\t" << product.second.endPrimer_
					<< "\t" << product.second.gap_
					<< "\t" << product.second.count_
					<< "\t" << product.first.size()
					<< "\t" << product.first
					<< "\t" << product.second.gapSeq_;
			const auto placement = placements.find(product.second.gapSeq_);
			if(AmpliconReadQC::ReadPairCategory::SHORT_PRODUCT == product.second.category_ && placements.end() != placement){
				const auto label = placement->second.getLabel(VecStr{targetOfPrimer(product.second.frontPrimer_), targetOfPrimer(product.second.endPrimer_)});
				writePlacement(shortOut, placement->second, label);
				placementCounts[std::make_pair(AmpliconReadQC::getCategoryStr(product.second.category_), label)] += product.second.count_;
			}else{
				writePlacement(shortOut, AmpliconProductPlacer::Placement{}, "");
			}
			shortOut << std::endl;
		}
		for(const auto & seqCount : res.kmerMismatchSeqs_){
			const auto & placement = placements.at(seqCount.first.second);
			const auto label = placement.getLabel(VecStr{seqCount.first.first});
			kmerMismatchOut << res.sampleName_
					<< "\t" << seqCount.first.first
					<< "\t" << seqCount.second
					<< "\t" << seqCount.first.second.size()
					<< "\t" << seqCount.first.second;
			writePlacement(kmerMismatchOut, placement, label);
			kmerMismatchOut << std::endl;
			placementCounts[std::make_pair(AmpliconReadQC::getCategoryStr(AmpliconReadQC::ReadPairCategory::KMER_MISMATCH), label)] += seqCount.second;
		}
		for(const auto & seqCount : res.mismatchedPrimerSeqs_){
			const auto & seq = std::get<2>(seqCount.first);
			const auto & placement = placements.at(seq);
			const auto label = placement.getLabel(VecStr{std::get<0>(seqCount.first), std::get<1>(seqCount.first)});
			mismatchedPrimerOut << res.sampleName_
					<< "\t" << std::get<0>(seqCount.first)
					<< "\t" << std::get<1>(seqCount.first)
					<< "\t" << seqCount.second
					<< "\t" << seq.size()
					<< "\t" << seq;
			writePlacement(mismatchedPrimerOut, placement, label);
			mismatchedPrimerOut << std::endl;
			placementCounts[std::make_pair(AmpliconReadQC::getCategoryStr(AmpliconReadQC::ReadPairCategory::MISMATCHED_PRIMERS), label)] += seqCount.second;
		}
		for(const auto & placementCount : placementCounts){
			placementCountsOut << res.sampleName_
					<< "\t" << placementCount.first.first
					<< "\t" << placementCount.first.second
					<< "\t" << placementCount.second
					<< "\t" << (0 == res.total_ ? 0.0 : 100.0 * placementCount.second / res.total_) << std::endl;
		}
	}
	return 0;
}

}  // namespace njhseq
