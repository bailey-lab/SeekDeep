//  SeekDeepUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/06/24.
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of SeekDeep.
//
// SeekDeep is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SeekDeep is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SeekDeep.  If not, see <http://www.gnu.org/licenses/>.
//

#include "SeekDeepUtilsRunner.hpp"

#include <bibseq/ProgramRunners.h>

namespace bibseq {

SeekDeepUtilsRunner::SeekDeepUtilsRunner() :
		bib::progutils::ProgramRunner(
				{ addFunc("dryRunQualityFiltering", dryRunQualityFiltering, false),
					addFunc("runMultipleCommands",    runMultipleCommands, false),
					addFunc("setupTarAmpAnalysis", setupTarAmpAnalysis, false),
					addFunc("replaceUnderscores", replaceUnderscores, false),
				  addFunc("rBind", ManipulateTableRunner::rBind, false),
					addFunc("genTargetInfoFromGenomes", genTargetInfoFromGenomes, false)}, //
				"SeekDeepUtils") {
}

int SeekDeepUtilsRunner::genTargetInfoFromGenomes(const bib::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	MultiGenomeMapper::inputParameters pars;
	std::string gffExtraAttributesStr = "descriptions";
	bfs::path primersFile = "";
	std::string forwardPrimer = "";
	std::string reversePrimer = "";
	std::string targetName = "";
	uint32_t errors = 0;
	uint32_t sizeLimit = 1000;
	uint32_t lenCutOffSizeExpand = 20;
	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();
	uint32_t barcodeSize = 0;
	std::string selectedGenomesStr = "";
	bool removeRefAlignments = false;
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(lenCutOffSizeExpand, "--lenCutOffSizeExpand", "When creating length cut off file how much to expand the length of the found targets");
	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Read Length", true);
	setUp.setOption(barcodeSize, "--barcodeSize", "Barcode Size, if on both primer, the sum of the two barcodes");
	setUp.setOption(errors, "--errors", "Number of errors to allow in primers");
	setUp.setOption(sizeLimit, "--sizeLimit", "Output target extractions for only targets below this size");
	setUp.setOption(pars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(removeRefAlignments, "--removeRefAlignments", "Remove Ref Alignments");
  setUp.setOption(pars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
  setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome");
  setUp.setOption(pars.gffDir_, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
  setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");

  setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	setUp.setOption(primersFile, "--primers", "A file that contains three columns, target,forwardPrimer,reversePrimer 5` to 3` directions, same file as the input to SeekDeep", true);
	setUp.processDirectoryOutputName("extractedRegions_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_ );

	PrimersAndMids ids(primersFile);
	if(0 == 	ids.getTargets().size() ){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << primersFile << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();
	bib::concurrent::LockableQueue<std::string> targetsQueue(getVectorOfMapKeys(ids.targets_));

	std::unique_ptr<MultiGenomeMapper> gMapper;
	//check for existence of genome directory
	bib::files::checkExistenceThrow(pars.genomeDir_, __PRETTY_FUNCTION__);

	//set up genome mapper;
	gMapper = std::make_unique<MultiGenomeMapper>(pars);

	//set up selected genomes
	gMapper->setSelectedGenomes(selectedGenomesStr);
	//init is threaded
	gMapper->pars_.numThreads_ = pars.numThreads_;
	gMapper->init();
	gMapper->pars_.numThreads_ = 1;
	//set threads;
	if(pars.numThreads_ >= 4){
		gMapper->pars_.numThreads_ = 2;
		pars.numThreads_ = pars.numThreads_/2;
	}
	if(1 == ids.targets_.size()){
		gMapper->pars_.numThreads_ = pars.numThreads_;
		pars.numThreads_ = 1;
	}
	bfs::path outputDir(setUp.pars_.directoryName_);

	auto alignToGenome = [](bib::concurrent::LockableQueue<bfs::path> & genomesQuque,
			const bfs::path & forwardFnp,
			const bfs::path & reverseFnp,
			const bfs::path & alignDir){
		bfs::path genomeFnp = "";
		while(genomesQuque.getVal(genomeFnp)){
			{
				BioCmdsUtils bioRunner;
				auto seqOpts = SeqIOOptions::genFastaIn(forwardFnp, false);
				seqOpts.out_.outFilename_ = bib::files::make_path(alignDir, bfs::basename(genomeFnp) + "_" + bfs::basename(forwardFnp) + ".sorted.bam");
				bioRunner.bowtie2Align(seqOpts, genomeFnp	, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -a --end-to-end");
			}
			{
				BioCmdsUtils bioRunner;
				auto seqOpts = SeqIOOptions::genFastaIn(reverseFnp, false);
				seqOpts.out_.outFilename_ = bib::files::make_path(alignDir, bfs::basename(genomeFnp) + "_" + bfs::basename(reverseFnp) + ".sorted.bam");
				bioRunner.bowtie2Align(seqOpts, genomeFnp, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -a --end-to-end"	);
			}
		}
	};

	comparison allowableErrors;
	allowableErrors.hqMismatches_ = errors;


	auto extractPathway =
			[&targetsQueue,&removeRefAlignments,&outputDir,&ids,&alignToGenome,&allowableErrors,&sizeLimit](
					const std::unique_ptr<MultiGenomeMapper> & gMapper) {
				std::string target;
				while(targetsQueue.getVal(target)) {
					const auto & primerInfo = ids.pDeterminator_->primers_.at(target);
					auto primerDirectory = bib::files::makeDir(outputDir, bib::files::MkdirPar(primerInfo.primerPairName_));
					auto bedDirectory = bib::files::makeDir(primerDirectory, bib::files::MkdirPar("genomeLocations"));
					auto forwardOpts = SeqIOOptions::genFastaOut(bib::files::make_path(primerDirectory, "forwardPrimer"));
					auto reverseOpts = SeqIOOptions::genFastaOut(bib::files::make_path(primerDirectory, "reversePrimer"));
					auto forDegens = createDegenStrs(primerInfo.forwardPrimerInfo_.seq_);
					auto revDegens = createDegenStrs(primerInfo.reversePrimerInfoForDir_.seq_);
					//
					std::vector<seqInfo> forSeqs;
					if(forDegens.size() == 1){
						forSeqs.emplace_back(primerInfo.forwardPrimerInfo_);
					}else{
						forSeqs = vecStrToReadObjs<seqInfo>(forDegens, primerInfo.forwardPrimerInfo_.name_);
					}
					std::vector<seqInfo> revSeqs;
					if(revDegens.size() == 1){
						revSeqs.emplace_back(primerInfo.reversePrimerInfoForDir_);
					}else{
						revSeqs = vecStrToReadObjs<seqInfo>(revDegens, primerInfo.reversePrimerInfoForDir_.name_);
					}
					SeqOutput::write(forSeqs, forwardOpts);
					SeqOutput::write(revSeqs, reverseOpts);
					auto refAlignmentDir = bib::files::makeDir(primerDirectory,bib::files::MkdirPar{ "refAlignments"});
					std::vector<std::thread> threads;
					bib::concurrent::LockableQueue<bfs::path> genomesQueue(gMapper->getGenomeFnps());
					auto forOutFnp = forwardOpts.out_.outName();
					auto revOutFnp = reverseOpts.out_.outName() ;
					for(uint32_t t = 0; t < gMapper->pars_.numThreads_; ++t){
						threads.emplace_back(std::thread(alignToGenome,
								std::ref(genomesQueue),
								std::cref(forOutFnp),
								std::cref(revOutFnp),
								std::cref(refAlignmentDir)));
					}
					for(auto & t : threads){
						t.join();
					}
					struct GenExtracRes{
						uint32_t forwardHits_{0};
						uint32_t reverseHits_{0};
						uint32_t extractCounts_{0};
					};
					std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
					for(const auto & genome : gMapper->genomes_){
						genomeExtractionsResults[genome.first] = GenExtracRes{};
					}
					std::unordered_map<std::string, std::vector<GenomeExtractResult>> genomeExtracts;
					for(const auto & genome : gMapper->genomes_) {
						auto forBamFnp = bib::files::make_path(refAlignmentDir,
								bfs::basename(genome.second->fnp_) + "_" + bfs::basename(forwardOpts.out_.outName()) + ".sorted.bam");
						auto revBamFnp = bib::files::make_path(refAlignmentDir,
														bfs::basename(genome.second->fnp_) + "_" + bfs::basename(reverseOpts.out_.outName()) + ".sorted.bam");
						auto forResults = gatherMapResults(forBamFnp, genome.second->fnpTwoBit_, allowableErrors);
						auto revResults = gatherMapResults(revBamFnp, genome.second->fnpTwoBit_, allowableErrors);
						genomeExtractionsResults[genome.first].forwardHits_ = forResults.size();
						genomeExtractionsResults[genome.first].reverseHits_ = revResults.size();
						if(!forResults.empty() && !revResults.empty()){
							auto uniForRes = getUniqueLocationResults(forResults);
							auto uniRevRes = getUniqueLocationResults(revResults);
							genomeExtracts[genome.first] = getPossibleGenomeExtracts(uniForRes, uniRevRes, sizeLimit);
						}
					}
					std::vector<seqInfo> refSeqs;
					std::vector<seqInfo> refTrimmedSeqs;
					for(auto & genome : genomeExtracts){
						if(genome.second.empty()){
							continue;
						}
						genomeExtractionsResults[genome.first].extractCounts_ = genome.second.size();
						OutputStream bedOut{OutOptions(bib::files::make_path(bedDirectory, genome.first + ".bed"))};
						uint32_t extractionCount = 0;
						for(auto & extract : genome.second){
							extract.setRegion();
							auto bedRegion = extract.gRegion_->genBedRecordCore();
							bedRegion.name_ = genome.first + "-" + target;
							auto name = genome.first;
							if(0 != extractionCount){
								name.append("." + estd::to_string(extractionCount));
							}
							MetaDataInName meta;
							meta.addMeta("genome", genome.first);
							meta.addMeta("target", target);
							meta.addMeta("extractionCount", extractionCount);
							name += " " + meta.createMetaName();
							bedRegion.name_ = name;
							bedOut << bedRegion.toDelimStr() << std::endl;
							TwoBit::TwoBitFile tReader(gMapper->genomes_.at(genome.first)->fnpTwoBit_);
							auto eSeq = extract.gRegion_->extractSeq(tReader);
							eSeq.name_ = name;
							bool refFound = false;
							for(auto & rSeq : refSeqs){
								if(rSeq.seq_ == eSeq.seq_){
									refFound = true;
									rSeq.name_.append("-" + eSeq.name_);
								}
							}
							if(!refFound){
								refSeqs.emplace_back(eSeq);
							}
							bool trimmed_refFound = false;

							auto innerSeq = extract.gRegionInner_->extractSeq(tReader);
							innerSeq.name_ = name;
							for(auto & rSeq : refTrimmedSeqs){
								if(rSeq.seq_ == innerSeq.seq_){
									trimmed_refFound = true;
									rSeq.name_.append("-" + innerSeq.name_);
								}
							}
							if(!trimmed_refFound){
								refTrimmedSeqs.emplace_back(innerSeq);
							}
							++extractionCount;
						}
					}

					if(!refSeqs.empty()){
						auto fullSeqOpts = SeqIOOptions::genFastaOut(bib::files::make_path(primerDirectory, primerInfo.primerPairName_ +".fasta"));
						SeqOutput::write(refSeqs, fullSeqOpts);
					}
					if(!refTrimmedSeqs.empty()){
						auto innerSeqOpts = SeqIOOptions::genFastaOut(bib::files::make_path(primerDirectory, primerInfo.primerPairName_ +"_primersRemoved.fasta"));
						SeqOutput::write(refTrimmedSeqs, innerSeqOpts);
					}

					table performanceTab(VecStr{"genome", "forwardPrimerHits", "reversePrimerHits", "extractionCounts"});
					auto genomeKeys = getVectorOfMapKeys(genomeExtractionsResults);
					bib::sort(genomeKeys);
					for(const auto & genomeKey : genomeKeys){
						performanceTab.addRow(genomeKey,
								genomeExtractionsResults[genomeKey].forwardHits_,
								genomeExtractionsResults[genomeKey].reverseHits_,
								genomeExtractionsResults[genomeKey].extractCounts_);
					}
					auto perTabOpts = TableIOOpts::genTabFileOut(bib::files::make_path(primerDirectory, "extractionCounts"),true);
					performanceTab.outPutContents(perTabOpts);
					if(removeRefAlignments){
						bib::files::rmDirForce(refAlignmentDir);
					}
				}
			};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars.numThreads_; ++t){
		threads.emplace_back(std::thread(extractPathway,
				std::cref(gMapper)));
	}
	for(auto & t : threads){
		t.join();
	}
	auto locationsCombined = bib::files::make_path(setUp.pars_.directoryName_, "locationsByGenome");
	bib::files::makeDir(bib::files::MkdirPar{locationsCombined});

	for(const auto & genome : gMapper->genomes_){
		auto genomeBedOpts = bib::files::make_path(locationsCombined, genome.first + ".bed");
		std::vector<std::shared_ptr<Bed3RecordCore>> allRegions;
		for(const auto & tar : ids.targets_){
			auto bedForTarFnp = bib::files::make_path(setUp.pars_.directoryName_, tar.first, "genomeLocations", genome.first + ".bed");
			if(bfs::exists(bedForTarFnp)){
				auto locs = getBed3s(bedForTarFnp);
				addOtherVec(allRegions, locs);
			}
		}
		if(!allRegions.empty()){
			if("" != genome.second->gffFnp_){
				intersectBedLocsWtihGffRecordsPars pars(genome.second->gffFnp_);
				pars.selectFeatures_ = VecStr{"gene"};
				pars.extraAttributes_ = VecStr{"description"};
				intersectBedLocsWtihGffRecords(allRegions, pars);
			}
			OutputStream genomeBedOut(genomeBedOpts);
			bib::sort(allRegions, [](
					const std::shared_ptr<Bed3RecordCore> & bed1,
					const std::shared_ptr<Bed3RecordCore> & bed2
					){
				if(bed1->chrom_ == bed2->chrom_){
					return bed1->chromStart_ < bed2->chromStart_;
				}else{
					return bed1->chrom_ < bed2->chrom_;
				}
			});
			for(const auto & loc : allRegions){
				genomeBedOut << loc->toDelimStrWithExtra() << std::endl;
			}
		}
	}


	auto forSeekDeepDir = bib::files::make_path(setUp.pars_.directoryName_, "forSeekDeep");

	bib::files::makeDir(bib::files::MkdirPar{forSeekDeepDir});
	auto refSeqsDir = bib::files::make_path(forSeekDeepDir, "refSeqs");
	bib::files::makeDir(bib::files::MkdirPar{refSeqsDir});
	OutOptions lenCutOffsOpts(bib::files::make_path(forSeekDeepDir, "lenCutOffs.txt"));
	OutputStream lenCutOffsOut(lenCutOffsOpts);
	lenCutOffsOut << "target\tminlen\tmaxlen" << "\n";

	OutOptions overlapStatusOpts(bib::files::make_path(forSeekDeepDir, "overlapStatuses.txt"));
	OutputStream overlapStatusOut(overlapStatusOpts);
	overlapStatusOut << "target\tstatus" << "\n";


	for(const auto & tar : ids.getTargets()){
		auto primersRemovedFnp = bib::files::make_path(setUp.pars_.directoryName_, tar, tar + "_primersRemoved.fasta");
		auto extractedSeqsFnp = bib::files::make_path(setUp.pars_.directoryName_, tar, tar + ".fasta");
		auto primersRemovedFinalFnp = bib::files::make_path(refSeqsDir, tar + ".fasta");
		if(bfs::exists(primersRemovedFnp)){
			{
				SeqIO reader(SeqIOOptions::genFastaInOut(primersRemovedFnp, primersRemovedFinalFnp));
				reader.openIn();
				reader.openOut();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					reader.write(seq);
				}
			}
			{
				std::vector<uint32_t> readLengths;
				SeqInput reader(SeqIOOptions::genFastaIn(extractedSeqsFnp));
				reader.openIn();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					readLengths.emplace_back(len(seq));
				}
				auto minlen = vectorMinimum(readLengths);
				auto maxlen = vectorMaximum(readLengths);
				lenCutOffsOut << tar
						<< "\t" << (minlen > lenCutOffSizeExpand ? minlen - lenCutOffSizeExpand : 0)
						<< "\t" << maxlen + lenCutOffSizeExpand << std::endl;
				uint32_t finalSize = maxlen + barcodeSize;
				if(finalSize > pairedEndLength && finalSize < (2* pairedEndLength - 10)){
					overlapStatusOut << tar
							<< "\t" << "R1EndsInR2" << std::endl;
				} else if(finalSize < pairedEndLength){
					overlapStatusOut << tar
							<< "\t" << "R1BeginsInR2" << std::endl;
				} else {
					overlapStatusOut << tar
							<< "\t" << "NoOverLap" << std::endl;
				}
			}
		}else{
			std::cerr << "Warning, no sequences extracted for " << tar << std::endl;
		}
	}


	return 0;
}


int SeekDeepUtilsRunner::dryRunQualityFiltering(
		const bib::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	std::string qualWindow;
	uint32_t qualityWindowLength;
	uint32_t qualityWindowStep;
	uint32_t qualityWindowThres;
	if (setUp.setOption(qualWindow, "--qualWindow", "SlidingQualityWindow")) {
		seqUtil::processQualityWindowString(qualWindow, qualityWindowLength,
				qualityWindowStep, qualityWindowThres);
	} else {
		qualityWindowLength = 50;
		qualityWindowStep = 5;
		qualityWindowThres = 25;
	}
	uint32_t qualCheck = 30;
	setUp.setOption(qualCheck, "--qualCheck", "Qual Check Level");
	double qualCheckCutOff = 0.90;
	setUp.setOption(qualCheckCutOff, "--qualCheckCutOff",
			"Cut Off for fraction of bases above qual check");
	bool plot = false;
	setUp.setOption(plot, "--plot", "Plot");
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::unordered_map<uint32_t, uint32_t> qualWindowCounts;
	uint32_t count = 0;
	std::vector<double> qualChecks;
	uint32_t failedQualCheck = 0;
	uint32_t failedQualWindow = 0;
	readObject read;
	while (reader.readNextRead(read)) {
		read.setBaseCountOnQualCheck(qualCheck);
		std::cout << "Currently on " << count << "\r";
		std::cout.flush();
		qualChecks.emplace_back(read.fractionAboveQualCheck_);
		if (read.fractionAboveQualCheck_ < qualCheckCutOff) {
			++failedQualCheck;
		}
		if (!seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
				qualityWindowStep, read.seqBase_.qual_)) {
			++failedQualWindow;
		}
		++count;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << bib::bashCT::boldBlack("Sliding Quality Window") << std::endl;
	std::cout << "Window Size: " << qualityWindowLength << ", Window Step: "
			<< qualityWindowStep << ", Window Thresdhold: " << qualityWindowThres
			<< std::endl;
	std::cout << "FailedWindow: " << getPercentageString(failedQualWindow, count)
			<< std::endl;
	std::cout << std::endl;
	std::cout << bib::bashCT::boldBlack("Quality Fraction above cut off")
			<< std::endl;
	std::cout << "Q" << qualCheck << ">" << qualCheckCutOff << std::endl;
	std::cout << "FailedCheck: " << getPercentageString(failedQualCheck, count)
			<< std::endl;
	auto qualCheckStats = getStatsOnVec(qualChecks);
	table qualChecksTab(qualCheckStats, VecStr { "stat", "value" });
	qualChecksTab.outPutContentOrganized(std::cout);
	return 0;
}

int SeekDeepUtilsRunner::runMultipleCommands(
		const bib::progutils::CmdArgs & inputCommands) {
	std::string filename = "";
	std::string logFile = "";
	uint32_t numThreads = 1;
	bool raw = false;
	bool noFilesInReplacementToks = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(logFile, "--logFile",
			"Name of a file to log the output of the commands");
	setUp.setOption(filename, "--cmdFile",
			"Name of the file, first line is command with REPLACETHIS, the next lines are the cmd to run with that line replacing REPLACETHIS",
			true);
	setUp.setOption(raw, "--raw",
			"If if the file is simply just a list of commands to run");
	setUp.setOption(noFilesInReplacementToks, "--noFilesInReplacementToks",
			"The replacement tokens are by default checked to see if there are files and more replacement strings are read in where each line is a replacement, this turns off that behavior");
	setUp.processVerbose();
	setUp.processWritingOptions();
	setUp.processDebug();
	if (setUp.needsHelp()) {
		std::cout
				<< "Input cmdFile should start with CMD: and then command and "
						"each line after that should be some sort of replacement to run the commnd"
				<< std::endl;
		std::cout << "Example" << std::endl;
		std::cout << "CMD:echo hello REPLACETHIS1 from REPLACETHIS2" << std::endl;
		std::cout << "REPLACETHIS1:nick,jon,mike" << std::endl;
		std::cout << "REPLACETHIS2:world,everyone" << std::endl;
		setUp.printFlags(std::cout);
		exit(1);
	}
	if (setUp.commands_.hasFlagCaseInsenNoDash("--gen")) {
		std::cout << "CMD:echo hello REPLACETHIS1 from REPLACETHIS2" << std::endl;
		std::cout << "REPLACETHIS1:nick,jon,mike" << std::endl;
		std::cout << "REPLACETHIS2:world,everyone" << std::endl;
		exit(1);
	}
	setUp.finishSetUp(std::cout);
	if (logFile == "") {
		logFile = bfs::path(filename).filename().replace_extension("").string()
				+ "Log.json";
	}
	std::ofstream outFile;
	if (!setUp.pars_.debug_) {
		openTextFile(outFile, logFile, ".json", setUp.pars_.ioOptions_.out_);
	}

	std::ifstream inFile(filename);
	if (!inFile) {
		std::stringstream ss;
		ss << bib::bashCT::bold << "Error in opening "
				<< bib::bashCT::boldRed(filename) << std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::string cmd;
	VecStr cmds;
	std::string line;
	if (raw) {
		while (bib::files::crossPlatGetline(inFile, line)) {
			if ("" == line || allWhiteSpaceStr(line) || bib::beginsWith(line, "#")) {
				continue;
			}
			cmds.emplace_back(line);
		}
	} else {
		std::map<std::string, VecStr> replacements;
		while (bib::files::crossPlatGetline(inFile, line)) {
			if ("" == line || allWhiteSpaceStr(line) || bib::beginsWith(line, "#")) {
				continue;
			}
			auto colonPos = line.find(":");

			if (std::string::npos == colonPos) {
				std::stringstream ss;
				ss << "Error in processing line: " << bib::bashCT::boldRed(line)
						<< std::endl;
				ss << "Need at least one colon" << std::endl;
				throw std::runtime_error { ss.str() };
			}
			VecStr toks { line.substr(0, colonPos), line.substr(colonPos + 1) };
			if (toks.front() == "CMD") {
				cmd = toks.back();
			} else {
				if(noFilesInReplacementToks){
					replacements[toks.front()] = tokenizeString(toks.back(), ",");
				}else{
					auto repToks = tokenizeString(toks.back(), ",");
					VecStr finalReplacements;
					for(const auto & tok : repToks){
						addOtherVec(finalReplacements, getInputValues(tok, ","));
					}
					replacements[toks.front()] = finalReplacements;
				}
			}
		}
		for (const auto & r : replacements) {
			if (cmds.empty()) {
				for (const auto & subR : r.second) {
					cmds.emplace_back(bib::replaceString(cmd, r.first, subR));
				}
			} else {
				VecStr newCmds;
				for (const auto & c : cmds) {
					for (const auto & subR : r.second) {
						newCmds.emplace_back(bib::replaceString(c, r.first, subR));
					}
				}
				cmds = newCmds;
			}
		}
	}

	if (setUp.pars_.verbose_) {
		printVector(cmds, "\n");
	}

	auto allRunOutputs = bib::sys::runCmdsThreaded(cmds, numThreads,
			setUp.pars_.verbose_, setUp.pars_.debug_);
	Json::Value allLog;

	allLog["totalTime"] = setUp.timer_.totalTime();
	allLog["cmdsfile"] = bib::json::toJson(bfs::absolute(filename));
	auto & cmdsLog = allLog["cmdsLog"];
	for (const auto & out : allRunOutputs) {
		cmdsLog.append(out.toJson());
	}
	if (!setUp.pars_.debug_) {
		outFile << allLog << std::endl;
	}

	return 0;
}

int SeekDeepUtilsRunner::replaceUnderscores(
		const bib::progutils::CmdArgs & inputCommands) {
	/**@todo add ability to search up to a pattern to replace
	 *  on rather than just number of underscores*/

	std::string pat = "";
	bfs::path dir = "./";
	uint32_t upTo = std::numeric_limits<uint32_t>::max();
	std::string replacement = "-";
	bool keepOriginals = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(dir, "--dir", "Directory with the files to rename", true);
	setUp.setOption(pat, "--pat,-p", "Rename files that match this pattern",
			true);
	setUp.setOption(upTo, "--upTo",
			"Replace up this many underscores from the begging of the filename");
	setUp.setOption(keepOriginals, "--keepOriginals",
			"Keep a copy of the original files");
	setUp.setOption(replacement, "--replacement",
			"What to replace the underscores with");
	setUp.finishSetUp(std::cout);

	auto files = bib::files::listAllFiles(dir.string(), false, {
			std::regex { pat } });

	std::map<bfs::path, bfs::path> renameKey;

	for (const auto & f : files) {
		if (bfs::is_regular_file(f.first)) {
			auto parDir = f.first.parent_path();
			auto bName = f.first.filename().string();
			uint32_t count = 0;
			while (count < upTo) {
				auto underPos = bName.find("_");
				if (std::string::npos == underPos) {
					break;
				}
				bName.erase(bName.begin() + underPos);
				bName.insert(bName.begin() + underPos, replacement.begin(),
						replacement.end());
				++count;
			}
			auto outPath = bib::files::make_path(parDir, bName);
			renameKey[f.first] = outPath;
		}
	}
	if (setUp.pars_.debug_) {
		for (const auto & rk : renameKey) {
			std::cout << bib::bashCT::green << rk.first << bib::bashCT::reset
					<< " -> " << bib::bashCT::blue << rk.second << bib::bashCT::reset
					<< std::endl;
		}
		return 0;
	}
	for (const auto & rk : renameKey) {
		if (keepOriginals) {
			bfs::copy(rk.first, rk.second);
		} else {
			bfs::rename(rk.first, rk.second);
		}
	}

	return 0;
}

table getStitchReport(const Json::Value & logs) {

	//report for stitching
	table flashOutTab { VecStr { "SampleName", "TotalPairs", "CombinedPairs",
			"UncombinedPairs", "PercentCombined" } };

	for (const auto & samp : logs) {
		std::stringstream ss;
		ss << samp["stitch"]["stdOut_"].asString();
		std::unordered_map<std::string, std::string> output;
		bool log = false;
		std::string line = "";
		while (bib::files::crossPlatGetline(ss, line)) {
			if (std::string::npos != line.find("Writing histogram files") && log) {
				break;
			}
			if (log) {
				if (std::string::npos != line.find(":")) {
					line = bib::replaceString(line, "[FLASH]", "");
					line.erase(std::remove_if(line.begin(), line.end(), isspace),
							line.end());
					auto toks = tokenizeString(line, ":");
					if (toks.size() != 2) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__
								<< ": Error in processing stitch output for sample: "
								<< samp["name"].asString() << " expected two values on line "
								<< " line " << " but got " << toks.size() << "\n";
						throw std::runtime_error { ss.str() };
					}
					output[toks[0]] = toks[1];
				}
			}
			if (std::string::npos != line.find("Read combination statistics:")) {
				log = true;
			}
		}
		VecStr headers = { "Totalpairs", "Combinedpairs", "Uncombinedpairs",
				"Percentcombined" };
		for (const auto & head : headers) {
			if (!bib::in(head, output)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__
						<< ": Error in processing stitch output for sample: "
						<< samp["name"].asString() << "\n";
				ss << "Coulnd't find: " << head << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		flashOutTab.addRow(samp["name"].asString(), output["Totalpairs"],
				output["Combinedpairs"], output["Uncombinedpairs"],
				bib::replaceString(output["Percentcombined"], "%", ""));
	}
	return flashOutTab;
}




//

}// namespace bibseq
