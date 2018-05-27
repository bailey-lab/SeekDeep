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
  //setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome", true);
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
	gMapper = std::make_unique<MultiGenomeMapper>(pars.genomeDir_, pars.primaryGenome_);

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
								std::cref( forOutFnp),
								std::cref( revOutFnp),
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
					for( auto & genome : genomeExtracts){
						genomeExtractionsResults[genome.first].extractCounts_ = genome.second.size();
						OutputStream bedOut{OutOptions(bib::files::make_path(bedDirectory, genome.first + ".bed"))};
						uint32_t extractionCount = 0;
						for(auto & extract : genome.second){
							extract.setRegion();
							auto bedRegion = extract.gRegion_->genBedRecordCore();
							bedRegion.name_ = genome.first + "-" + target;
							bedOut << bedRegion.toDelimStr() << std::endl;
							auto name = genome.first;
							if(0 != extractionCount){
								name.append("." + estd::to_string(extractionCount));
							}
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
					auto fullSeqOpts = SeqIOOptions::genFastaOut(bib::files::make_path(primerDirectory, primerInfo.primerPairName_ +".fasta"));
					auto innerSeqOpts = SeqIOOptions::genFastaOut(bib::files::make_path(primerDirectory, primerInfo.primerPairName_ +"_primersRemoved.fasta"));
					SeqOutput::write(refSeqs, fullSeqOpts);
					SeqOutput::write(refTrimmedSeqs, innerSeqOpts);
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



int SeekDeepUtilsRunner::setupTarAmpAnalysis(
		const bib::progutils::CmdArgs & inputCommands) {
	TarAmpAnalysisSetup::TarAmpPars pars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.debug = setUp.pars_.debug_;
	setUp.setOption(pars.technology, "--technology",
			"Sequencing Technology (should be 454,IonTorrent, or Illumina",
			false, "Technology");
	stringToLower(pars.technology);
	if (pars.technology != "454" && pars.technology != "iontorrent"
			&& pars.technology != "illumina") {
		setUp.failed_ = true;
		std::stringstream ss;
		ss
				<< "Error in setting technology, should be 454, IonTorrent, or Illumina not "
				<< pars.technology << "\n";
		setUp.addWarning(ss.str());
	}
	setUp.setOption(pars.samplesNamesFnp, "--samples",
			"A file containing the samples names, columns should go target,sample,pcr1,pcr2 (optional) etc",
			true, "ID Files");
	setUp.setOption(pars.outDir, "--outDir", "Directory to setup for analysis",
			true, "Output");
	setUp.setOption(pars.inputDir, "--inputDir",
			"Input Directory of raw data files", true, "Input");
	setUp.setOption(pars.idFile, "--idFile",
			"ID file containing primer and possible additional MIDs", true, "ID Files");
	setUp.setOption(pars.byIndex, "--byIndex",
			"If the input sample names are by index and not targets", false, "ID Files");
	setUp.setOption(pars.targetsToIndexFnp, "--targetsToIndex",
			"A tsv file with two columns named targets and index where targets in a comma separated value with the targets for the index in index",
			false, "ID Files");

	setUp.setOption(pars.numThreads, "--numThreads", "Number of CPUs to use");

	setUp.setOption(pars.extraExtractorCmds, "--extraExtractorCmds",
			"Extra extractor cmds to add to the defaults", false, "Extra Commands");
	setUp.setOption(pars.extraQlusterCmds, "--extraQlusterCmds",
			"Extra qluster cmds to add to the defaults", false, "Extra Commands");
	setUp.setOption(pars.extraProcessClusterCmds, "--extraProcessClusterCmds",
			"Extra process clusters cmds to add to the defaults", false, "Extra Commands");

	setUp.setOption(pars.noQualTrim, "--noQualTrim", "No Quality Trim", false, "Pre Processing");

	setUp.setOption(pars.r1Trim, "--r1Trim", "Trim R1 Reads to this length", false, "Pre Processing");
	setUp.setOption(pars.r2Trim, "--r2Trim", "Trim R2 Reads to this length", false, "Pre Processing");

	setUp.setOption(pars.groupMeta, "--groupMeta", "Group Metadata", false, "Meta");

	setUp.setOption(pars.lenCutOffsFnp, "--lenCutOffs",
			"A file with 3 columns, target,minlen,maxlen to supply length cut off specifically for each target", false, "Extractor");
	setUp.setOption(pars.refSeqsDir, "--refSeqsDir",
			"A directory of fasta files where each file is named with the input target names", false, "Extractor");
	if (pars.techIs454() || pars.techIsIonTorrent()) {
		pars.inputFilePat = ".*.fastq";
	}
	setUp.setOption(pars.overlapStatusFnp, "--overlapStatusFnp",
			"A file with two columns, target,status; status column should contain 1 of 3 values (capitalization doesn't matter): r1BegOverR2End,r1EndOverR2Beg,NoOverlap. r1BegOverR2End=target size < read length (causes read through),r1EndOverR2Beg= target size > read length less than 2 x read length, NoOverlap=target size > 2 x read length",
			pars.techIsIllumina(), "Illumina Stitching");

	setUp.setOption(pars.inputFilePat, "--inputFilePat",
			"The input file pattern in the input directory to work on", false, "Input");

	setUp.finishSetUp(std::cout);

	VecStr warnings;

	pars.allChecks(warnings);

	if (!warnings.empty()) {
		std::stringstream ss;
		if (1 == warnings.size()) {
			ss << __PRETTY_FUNCTION__ << ": There is " << warnings.size() << " error."
					<< "\n";
		} else {
			ss << __PRETTY_FUNCTION__ << ": There are " << warnings.size()
					<< " errors." << "\n";
		}
		for (const auto & warn : warnings) {
			ss << warn << std::endl;
		}
		throw std::runtime_error { ss.str() };
	}

	bool foundErrors = false;
	std::stringstream errorOutput;

	TarAmpAnalysisSetup analysisSetup(pars);
	setUp.startARunLog(bib::appendAsNeededRet(analysisSetup.dir_.string(), "/"));
	auto targets = analysisSetup.idsMids_->getTargets();
	auto sampNamesTargets = analysisSetup.getTargets();
	bib::sort(targets);
	bib::sort(sampNamesTargets);
	VecStr idMissingTargets;
	VecStr sampNamesTargetsTarMissing;
	for (const auto & tar : sampNamesTargets) {
		if (!bib::in(tar, targets)) {
			idMissingTargets.emplace_back(tar);
		}
	}
	for (const auto & tar : targets) {
		if (!bib::in(tar, sampNamesTargets)) {
			sampNamesTargetsTarMissing.emplace_back(tar);
		}
	}

	if (!idMissingTargets.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": error, missing the following targets from the id file "
				<< pars.idFile << "\n";
		ss << "Targets: " << bib::conToStr(idMissingTargets, ", ") << "\n";
		ss << "AvailableTargets: " << bib::conToStr(targets, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}

	if (!sampNamesTargetsTarMissing.empty()) {
		foundErrors = true;
		errorOutput << __PRETTY_FUNCTION__
				<< ": warning, missing the following targets from the sample name file"
				<< pars.samplesNamesFnp << "\n";
		errorOutput << "Targets: "
				<< bib::conToStr(sampNamesTargetsTarMissing, ", ") << "\n";
		errorOutput << "AvailableTargets: " << bib::conToStr(sampNamesTargets, ", ")
				<< "\n";
	}

	if ("" != pars.groupMeta) {
		if (!analysisSetup.groupMetaData_->missingSamples_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, the following samples were in the meta data file, "
					<< pars.groupMeta
					<< " but were missing from the input sample names file, "
					<< pars.samplesNamesFnp << std::endl;
			for (const auto & samp : analysisSetup.groupMetaData_->missingSamples_) {
				errorOutput << "\tSample: " << samp << "\n";
			}
		}
		if (!analysisSetup.groupMetaData_->missingMetaForSamples_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, the following samples were in the input samples name file, "
					<< pars.samplesNamesFnp
					<< " but were missing from the meta data file, " << pars.groupMeta
					<< std::endl;
			for (const auto & samp : analysisSetup.groupMetaData_->missingMetaForSamples_) {
				errorOutput << "\tSample: " << samp << "\n";
			}
		}
	}

	if ("" != pars.refSeqsDir) {
		if (!analysisSetup.forRefSeqs_.missing_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, reference sequences were supplied but reference files for the following targets couldn't be found"
					<< "\n";
			for (const auto & tar : analysisSetup.forRefSeqs_.missing_) {
				errorOutput << "\tTarget: " << tar << "\n";
			}
		}
		if (!analysisSetup.forRefSeqs_.notMatching_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, reference sequences were supplied but there are no matching reference targets in project for the following"
					<< "\n";
			for (const auto & tar : analysisSetup.forRefSeqs_.notMatching_) {
				errorOutput << "\tTarget: " << tar << "\n";
			}
			errorOutput << "\tOptions are: "
					<< bib::conToStr(analysisSetup.getTargets(), ", ") << "\n";
		}
	}

	if ("" != pars.lenCutOffsFnp) {
		if (!analysisSetup.forLenCutOffs_.missing_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, length cut offs were supplied but cut offs for the following targets couldn't be found"
					<< "\n";
			for (const auto & tar : analysisSetup.forLenCutOffs_.missing_) {
				errorOutput << "\tTarget: " << tar << "\n";
			}
		}
		if (!analysisSetup.forLenCutOffs_.notMatching_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, length cut offs were supplied but cut offs for targets didn't match targets in project for the following targets"
					<< "\n";
			for (const auto & tar : analysisSetup.forRefSeqs_.notMatching_) {
				errorOutput << "\tTarget: " << tar << "\n";
			}
			errorOutput << "\tOptions are: "
					<< bib::conToStr(analysisSetup.getTargets(), ", ") << "\n";
		}
	}

	//now write id files
	analysisSetup.writeOutIdFiles();

	auto files = bib::files::listAllFiles(pars.inputDir.string(), false, {
			std::regex { analysisSetup.pars_.inputFilePat } });
	if (setUp.pars_.debug_) {
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}
	auto expectedSamples = analysisSetup.getExpectantInputNames();
	if(setUp.pars_.debug_){
		std::cout << "Expected input files: " << std::endl;
		std::cout << bib::conToStr(expectedSamples, "\n") << std::endl;
	}
	VecStr unrecognizedInput;
	VecStr sampleFilesFound;
	VecStr sampleFilesNotFound;
	//VecStr samplesExtracted;
	//VecStr samplesEmpty;
	Json::Value logs;
	std::mutex logsMut;
	std::unordered_map<std::string, std::pair<VecStr, VecStr>> readsByPairs ;
	std::unordered_map<std::string, bfs::path> filesByPossibleName;
	ReadPairsOrganizer rpOrganizer(expectedSamples);

	if (analysisSetup.pars_.techIsIllumina()) {

		rpOrganizer.processFiles(files);
		readsByPairs = rpOrganizer.processReadPairs();
		auto keys = getVectorOfMapKeys(readsByPairs);
		bib::sort(keys);
		sampleFilesFound = keys;
		if(setUp.pars_.debug_){
			std::cout << "Samples found: " << std::endl;
			std::cout << bib::conToStr(keys, "\n") << std::endl;
		}
		for(const auto & expected : expectedSamples){
			if(!bib::in(expected, sampleFilesFound)){
				sampleFilesNotFound.emplace_back(expected);
			}
		}
		if (!rpOrganizer.readPairsUnrecognized_.empty()) {
			foundErrors = true;
			errorOutput
					<< "The following files were found but didn't match any input sample names in "
					<< analysisSetup.pars_.samplesNamesFnp << std::endl;
			for (const auto & samp : rpOrganizer.readPairsUnrecognized_) {
				errorOutput << "\tPossible Sample Name: " << samp.first << std::endl;
				for (const auto & sampFiles : samp.second) {
					errorOutput << "\t\t" << sampFiles << std::endl;
				}
			}
		}
		if (!sampleFilesNotFound.empty()) {
			foundErrors = true;
			errorOutput << "The following files were expected but were not found " << std::endl;
			for (const auto & samp : sampleFilesNotFound) {
				errorOutput << "\tSample: " << samp << std::endl;
			}
		}
//		if (!samplesEmpty.empty()) {
//			foundErrors = true;
//			errorOutput << "The following files were found to be empty " << std::endl;
//			for (const auto & samp : samplesEmpty) {
//				errorOutput << "\tSample: " << samp << std::endl;
//				for (const auto & sampFiles : rpOrganizer.readPairs_.at(samp)) {
//					errorOutput << "\t\t" << sampFiles << std::endl;
//				}
//			}
//		}
	} else {
		// ion torrent and 454 extraction and moving goes here
		// need to add to extractSamples if found
		// need to find samples that are empty
		// add to inputPassed
		// just checking for possible compression with file ending, might consider changing to libmagic or something to make sure
		bool compressed = false;
		if(bib::endsWith(analysisSetup.pars_.inputFilePat, ".gz")){
			compressed = true;
		}

		for(const auto & file : files){
			auto fNameNoExt = file.first.filename().replace_extension("");
			if(compressed){
				fNameNoExt.replace_extension("");
			}
			if(bib::in(fNameNoExt.string(), expectedSamples)){
				filesByPossibleName[fNameNoExt.string()] = file.first;
			}else{
				unrecognizedInput.emplace_back(file.first.string());
			}
		}
		auto keys = bib::getVecOfMapKeys(filesByPossibleName);
		bib::sort(keys);
	}

//	VecStr samplesNotFound;
//	for (const auto & expected : expectedSamples) {
//		if (!bib::in(expected, samplesExtracted)
//				&& !bib::in(bib::replaceString(expected, "MID", ""),
//						samplesExtracted)) {
//			samplesNotFound.emplace_back(expected);
//		}
//	}
//
//	if (!samplesNotFound.empty()) {
//		foundErrors = true;
//		errorOutput << "The following input files were not found" << std::endl;
//		for (const auto & samp : samplesNotFound) {
//			errorOutput << "\tSample: " << samp << std::endl;
//		}
//	}

	OutOptions wrningsOpts(
			bib::files::make_path(analysisSetup.dir_, "WARNINGS_PLEASE_READ.txt"));
	if (foundErrors) {
		std::ofstream outWarnings;
		openTextFile(outWarnings, wrningsOpts);
		outWarnings << errorOutput.str() << std::endl;
	}

	//make population clustering directory
	analysisSetup.setUpPopClusteringDirs(setUp.pars_.verbose_);

	VecStr extractorCmds;
	VecStr qlusterCmds;
	if (analysisSetup.pars_.byIndex) {
		if(setUp.pars_.debug_){
			std::cout << "Samples:" << std::endl;
			std::cout << bib::conToStr(getVectorOfMapKeys(analysisSetup.samples_), "\n") << std::endl;

		}


		//extractor cmds
		std::string extractorCmdTemplate;
		if(analysisSetup.pars_.techIsIllumina()){
			 extractorCmdTemplate =
							setUp.commands_.masterProgram_
									+ " extractorPairedEnd --dout {INDEX}_extraction --overWriteDir  ";
		}else{
			 extractorCmdTemplate =
							setUp.commands_.masterProgram_
									+ " extractor --dout {INDEX}_extraction --overWriteDir  ";
		}
		if (analysisSetup.pars_.techIsIllumina()) {
			extractorCmdTemplate += "--fastq1 {INDEX_R1} --fastq2 /{INDEX_R2} ";
		} else {
			extractorCmdTemplate += "--fastq {INDEX_InFile} ";
		}
		std::string lenCutOffsTemplate ="--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt ";
		std::string idTemplate = "--id info/ids/{TARS}.id.txt ";
		std::string overLapStatusTemplate = "--overlapStatusFnp info/ids/{TARS}_overlapStatus.tab.txt ";
		std::string refSeqsDir = "--compareSeq info/refs/ ";
		//qluster cmds;
		auto extractionDirs = bib::files::make_path(
				bfs::absolute(analysisSetup.dir_), "{INDEX}_extraction");
		std::string qlusterCmdTemplate =
				"cd \"" + extractionDirs.string() + "\" && "
						+ "if [ -f {TARGET}{MIDREP}.fastq  ]; then "
						+ setUp.commands_.masterProgram_
						+ " qluster "
								"--fastq {TARGET}{MIDREP}.fastq "
								"--alnInfoDir {TARGET}{MIDREP}_alnCache --overWriteDir "
								"--additionalOut \"../popClustering/{TARGET}/locationByIndex/{INDEX}.tab.txt\" "
								"--overWrite --dout {TARGET}{MIDREP}_qlusterOut ";
		if (analysisSetup.pars_.techIsIllumina()) {
			qlusterCmdTemplate += "--illumina --qualThres 25,20";
		} else if (analysisSetup.pars_.techIs454()) {
			qlusterCmdTemplate += "--454";
		} else if (analysisSetup.pars_.techIsIonTorrent()) {
			qlusterCmdTemplate += "--ionTorrent";
		}
		auto indexes = analysisSetup.getIndexes();
		if(setUp.pars_.verbose_){
			std::cout << "indexes" << std::endl;
			printVector(indexes);
		}
		for (const auto & index : indexes) {
			//if (bib::in(index, inputPassed)) {
			if (true) {

				auto tarsNames = bib::conToStr(analysisSetup.indexToTars_[index], "_");
				auto cmds = VecStr { extractorCmdTemplate, idTemplate };
				if (bfs::exists(
						bib::files::make_path(analysisSetup.idsDir_,
								tarsNames + "_lenCutOffs.tab.txt"))) {
					cmds.emplace_back(lenCutOffsTemplate);
				}
				bool anyRefs = false;
				for (const auto & tar : analysisSetup.indexToTars_[index]) {
					if (bfs::exists(
							bib::files::make_path(analysisSetup.refsDir_, tar + ".fasta"))) {
						anyRefs = true;
						break;
					}
				}
				if(anyRefs){
					cmds.emplace_back(refSeqsDir);
				}

				if ("" != analysisSetup.pars_.extraExtractorCmds) {
					cmds.emplace_back(analysisSetup.pars_.extraExtractorCmds);
				}
				if(pars.techIsIllumina()){
					cmds.emplace_back(overLapStatusTemplate);
				}
				auto currentExtractCmd = bib::conToStr(cmds, " ");
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{INDEX}",
						index);
				if(pars.techIsIllumina()){
					//INDEX_R1, INDEX_R2
					auto searchForPairs = readsByPairs.find(index);
					if(searchForPairs != readsByPairs.end()){
						currentExtractCmd = bib::replaceString(currentExtractCmd, "{INDEX_R1}",
								bfs::absolute(searchForPairs->second.first.front()).string() );
						currentExtractCmd = bib::replaceString(currentExtractCmd, "{INDEX_R2}",
								bfs::absolute(searchForPairs->second.second.front()).string() );
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error index: " << index << " was not found in readsByPairs" << "\n";
						ss << "Options: " << bib::conToStr(bib::getVecOfMapKeys(readsByPairs), ", ");
					}
				}else{
					//INDEX_InFile
					auto searchForFile = filesByPossibleName.find(index);
					if(searchForFile != filesByPossibleName.end()){
						currentExtractCmd = bib::replaceString(currentExtractCmd, "{INDEX_InFile}",
								bfs::absolute(searchForFile->second).string());
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error index: " << index << " was not found in filesByPossibleName" << "\n";
						ss << "Options: " << bib::conToStr(bib::getVecOfMapKeys(filesByPossibleName), ", ");
					}
				}
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{TARS}",
						tarsNames);
				extractorCmds.emplace_back(currentExtractCmd);
				for (const auto & mid : analysisSetup.idsMids_->getMids()) {
					for (const auto & tar : analysisSetup.indexToTars_[index]) {
						auto currentQlusterCmdTemplate = qlusterCmdTemplate;
						if ("" != analysisSetup.pars_.extraQlusterCmds) {
							currentQlusterCmdTemplate += " "
									+ analysisSetup.pars_.extraQlusterCmds;
						}
						currentQlusterCmdTemplate += "; fi";

						currentQlusterCmdTemplate = bib::replaceString(
								currentQlusterCmdTemplate, "{INDEX}", index);
						currentQlusterCmdTemplate = bib::replaceString(
								currentQlusterCmdTemplate, "{TARGET}", tar);
						currentQlusterCmdTemplate = bib::replaceString(
								currentQlusterCmdTemplate, "{MIDREP}", mid);
						qlusterCmds.emplace_back(currentQlusterCmdTemplate);
					}
				}
			}
		}
	} else {
		//extractor cmds
		std::string extractorCmdTemplate;
		if(analysisSetup.pars_.techIsIllumina()){
			 extractorCmdTemplate =
							setUp.commands_.masterProgram_
									+ " extractorPairedEnd --dout {REP}_extraction --overWriteDir  ";
		}else{
			 extractorCmdTemplate =
							setUp.commands_.masterProgram_
									+ " extractor --dout {INDEX}_extraction --overWriteDir  ";
		}
		if (analysisSetup.pars_.techIsIllumina()) {
			extractorCmdTemplate += "--fastq1 {REP_R1} --fastq2 /{REP_R2} ";
		} else {
			extractorCmdTemplate += "--fastq {REP_InFile} ";
		}

		std::string lenCutOffsTemplate ="--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt ";
		std::string idTemplate = "--id info/ids/{TARS}.id.txt ";
		std::string overLapStatusTemplate = "--overlapStatusFnp info/ids/{TARS}_overlapStatus.tab.txt ";
		std::string refSeqsDir = "--compareSeq info/refs/ ";

		std::string sampleNameTemplate = "--sampleName {REP}";
		//qluster cmds;
		auto extractionDirs = bib::files::make_path(
				bfs::absolute(analysisSetup.dir_), "{REP}_extraction");

		std::string qlusterCmdTemplate = "cd " + extractionDirs.string() + " && "
				+ " if [ -f {TARGET}{MIDREP}.fastq  ]; then "
										+ setUp.commands_.masterProgram_
										+ " qluster "
						"--fastq {TARGET}{MIDREP}.fastq "
						"--alnInfoDir {TARGET}{MIDREP}_alnCache --overWriteDir "
						"--additionalOut ../popClustering/locationByIndex/{TARGET}.tab.txt "
						"--overWrite --dout {TARGET}{MIDREP}_qlusterOut ";
		if (analysisSetup.pars_.techIsIllumina()) {
			qlusterCmdTemplate += "--illumina --qualThres 25,20";
		} else if (analysisSetup.pars_.techIs454()) {
			qlusterCmdTemplate += "--454";
		} else if (analysisSetup.pars_.techIsIonTorrent()) {
			qlusterCmdTemplate += "--ionTorrent";
		}
		if(setUp.pars_.debug_){
			std::cout << "Samples:" << std::endl;
			std::cout << bib::conToStr(getVectorOfMapKeys(analysisSetup.samples_), "\n") << std::endl;

		}

		std::unordered_map<std::string, VecStr> targetsForReps;
		for (const auto & tars : analysisSetup.samples_) {
			for (const auto & rep : tars.second.getReps()) {
				targetsForReps[rep].emplace_back(tars.first);
			}
		}
		for (auto & rep : targetsForReps) {
			bib::sort(rep.second);
		}

		for (const auto & rep : targetsForReps) {
			if (bib::in(rep.first, sampleFilesFound)
					|| bib::in(bib::replaceString(rep.first, "MID", ""), sampleFilesFound)) {
				std::string fName = rep.first;
//				if (!bib::in(rep.first, samplesExtracted)) {
//					fName = bib::replaceString(rep.first, "MID", "");
//				}

				std::string sampName = rep.first;
				if (!bib::beginsWith(rep.first, "MID")) {
					sampName = "MID" + rep.first;
				}

				auto tarsNames = bib::conToStr(rep.second, "_");
				auto currentSampTemp = bib::replaceString(sampleNameTemplate, "{REP}",
						sampName);
				auto cmds = VecStr { extractorCmdTemplate, idTemplate, currentSampTemp };
				if (bfs::exists(
						bib::files::make_path(analysisSetup.idsDir_,
								tarsNames + "_lenCutOffs.tab.txt"))) {
					cmds.emplace_back(lenCutOffsTemplate);
				}
				if ("" != analysisSetup.pars_.extraExtractorCmds) {
					cmds.emplace_back(analysisSetup.pars_.extraExtractorCmds);
				}
				//add ref files
				bool anyRefs = false;
				for (const auto & tar : rep.second) {
					if (bfs::exists(
							bib::files::make_path(analysisSetup.refsDir_, tar + ".fasta"))) {
						anyRefs = true;
						break;
					}
				}
				if(anyRefs){
					cmds.emplace_back(refSeqsDir);
				}
				// add overlap status
				if(pars.techIsIllumina()){
					cmds.emplace_back(overLapStatusTemplate);
				}


				auto currentExtractCmd = bib::conToStr(cmds, " ");
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{REP}",
						fName);
				if(pars.techIsIllumina()){
					//INDEX_R1, INDEX_R2
					auto searchForPairs = readsByPairs.find(fName);
					if(searchForPairs != readsByPairs.end()){
						currentExtractCmd = bib::replaceString(currentExtractCmd, "{REP_R1}",
								bfs::absolute(searchForPairs->second.first.front()).string() );
						currentExtractCmd = bib::replaceString(currentExtractCmd, "{REP_R2}",
								bfs::absolute(searchForPairs->second.second.front()).string() );
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error REP: " << fName << " was not found in readsByPairs" << "\n";
						ss << "Options: " << bib::conToStr(bib::getVecOfMapKeys(readsByPairs), ", ");
					}
				}else{
					//INDEX_InFile
					auto searchForFile = filesByPossibleName.find(fName);
					if(searchForFile != filesByPossibleName.end()){
						currentExtractCmd = bib::replaceString(currentExtractCmd, "{REP_InFile}",
								bfs::absolute(searchForFile->second).string());
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error REP: " << fName << " was not found in filesByPossibleName" << "\n";
						ss << "Options: " << bib::conToStr(bib::getVecOfMapKeys(filesByPossibleName), ", ");
					}
				}
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{TARS}",
						tarsNames);
				extractorCmds.emplace_back(currentExtractCmd);
				for (const auto & tar : rep.second) {
					std::string currentQlusterCmdTemplate = qlusterCmdTemplate;
					if ("" != analysisSetup.pars_.extraQlusterCmds) {
						currentQlusterCmdTemplate += " "
								+ analysisSetup.pars_.extraQlusterCmds;
					}
					currentQlusterCmdTemplate += "; fi";
					currentQlusterCmdTemplate = bib::replaceString(
							currentQlusterCmdTemplate, "{REP}", fName);
					currentQlusterCmdTemplate = bib::replaceString(
							currentQlusterCmdTemplate, "{MIDREP}", sampName);
					currentQlusterCmdTemplate = bib::replaceString(
							currentQlusterCmdTemplate, "{TARGET}", tar);
					qlusterCmds.emplace_back(currentQlusterCmdTemplate);
				}
			}
		}
	}

	OutOptions extractorCmdsOpts(
			bib::files::make_path(analysisSetup.dir_, "extractorCmds.txt"));
	std::ofstream extractorCmdsFile;
	openTextFile(extractorCmdsFile, extractorCmdsOpts);
	printVector(extractorCmds, "\n", extractorCmdsFile);

	OutOptions qlusterCmdsOpts(
			bib::files::make_path(analysisSetup.dir_, "qlusterCmds.txt"));
	std::ofstream qlusterCmdsFile;
	openTextFile(qlusterCmdsFile, qlusterCmdsOpts);
	printVector(qlusterCmds, "\n", qlusterCmdsFile);

	//process cluster cmds
	std::string processClusterTemplate =
			setUp.commands_.masterProgram_
					+ " processClusters "
							"--alnInfoDir alnCache --strictErrors --dout analysis --fastq output.fastq --overWriteDir";
	if ("" != analysisSetup.pars_.extraProcessClusterCmds) {
		processClusterTemplate += " " + analysisSetup.pars_.extraProcessClusterCmds;
	}
	VecStr processClusterCmds;
	if (nullptr != analysisSetup.groupMetaData_) {
		processClusterTemplate += " --groupingsFile "
				+ bib::files::make_path(bfs::absolute(analysisSetup.infoDir_),
						"groupMeta.tab.txt").string();
	}

	auto popDir = bib::files::make_path(bfs::absolute(analysisSetup.dir_),
			"popClustering");

	for (const auto & tar : targets) {
		std::stringstream processClustersCmdsStream;
		processClustersCmdsStream
				<< "cd " + bib::files::make_path(popDir, tar).string() + " && "
						+ processClusterTemplate + " --experimentName " + tar;
		auto refSeqFnp = bib::files::make_path(pars.outDir, "info/refs/" + tar + ".fasta");
		if(bfs::exists(refSeqFnp)){
			processClustersCmdsStream << " --ref " << bfs::absolute(refSeqFnp);
		}
		processClusterCmds.emplace_back(processClustersCmdsStream.str());

	}
	OutOptions processClusterCmdsOpts(
			bib::files::make_path(analysisSetup.dir_, "processClusterCmds.txt"));
	std::ofstream processClusterCmdsFile;
	openTextFile(processClusterCmdsFile, processClusterCmdsOpts);
	printVector(processClusterCmds, "\n", processClusterCmdsFile);

	//gen analysis configs
	std::string genConfigTemplate =
			setUp.commands_.masterProgram_
					+ " genProjectConfig --projectName {TARGET} --out serverConfigs/{TARGET}.config";

	VecStr genConfigCmds;
	for (const auto & tar : targets) {
		genConfigCmds.emplace_back(
				bib::replaceString(
						genConfigTemplate + " --mainDir popClustering/{TARGET}/analysis",
						"{TARGET}", tar));
	}
	OutOptions genConfigCmdsOpts(
			bib::files::make_path(analysisSetup.dir_, "genConfigCmds.txt"));
	std::ofstream genConfigCmdsFile;
	openTextFile(genConfigCmdsFile, genConfigCmdsOpts);
	printVector(genConfigCmds, "\n", genConfigCmdsFile);

	//start server config
	OutOptions startServerCmdOpts(
			bib::files::make_path(analysisSetup.dir_, "startServerCmd.sh"));
	std::ofstream startServerCmdFile;
	openTextFile(startServerCmdFile, startServerCmdOpts);
	startServerCmdFile << "#!/usr/bin/env bash" << std::endl;
	startServerCmdFile
			<< "# Will automatically run the server in the background and with nohup so it will keep running"
			<< std::endl;

	startServerCmdFile << "if [[ $# -ne 2 ]] && [[ $# -ne 0 ]]; then"
			<< std::endl;
	startServerCmdFile
			<< "	echo \"Illegal number of parameters, needs either 0 or 2 arguments, if 2 args 1) port number to server on 2) the name to serve on\""
			<< std::endl;
	startServerCmdFile << "	echo \"Examples\"" << std::endl;
	startServerCmdFile << "	echo \"./startServerCmd.sh\"" << std::endl;
	startServerCmdFile << "	echo \"./startServerCmd.sh 9882 pcv2\"	"
			<< std::endl;
	startServerCmdFile << "	exit	" << std::endl;
	startServerCmdFile << "fi" << std::endl;
	startServerCmdFile << "" << std::endl;
	startServerCmdFile << "if [[ $# -eq 2 ]]; then" << std::endl;
	startServerCmdFile << "	nohup " << setUp.commands_.masterProgram_
			<< " popClusteringViewer --verbose --configDir $(pwd)/serverConfigs --port $1 --name $2 &"
			<< std::endl;
	startServerCmdFile << "else" << std::endl;
	startServerCmdFile << "	nohup " << setUp.commands_.masterProgram_
			<< " popClusteringViewer --verbose --configDir $(pwd)/serverConfigs & "
			<< std::endl;
	startServerCmdFile << "fi" << std::endl;
	//make file executable
	chmod(startServerCmdOpts.outFilename_.c_str(),
			S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	//combining extraction counts
	OutOptions combineExtractionCmdOpts(
			bib::files::make_path(analysisSetup.dir_, "combineExtractionCountsCmd.sh"));
	OutputStream combineExtractionCmdOut(combineExtractionCmdOpts);
	combineExtractionCmdOut << "#!/usr/bin/env bash" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains allPrimerCounts.tab.txt --delim tab --header --out reports/allAllPrimerCounts.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains extractionProfile.tab.txt --delim tab --header --out reports/allExtractionProfile.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains extractionStats.tab.txt --delim tab --header --out reports/allExtractionStats.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains midCounts.tab.txt --delim tab --header --out reports/allMidCounts.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains processPairsCounts.tab.txt --delim tab --header --out reports/allProcessPairsCounts.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains top_mostCommonR1Starts_for_unrecognizedBarcodes.tab.txt --delim tab --header --out reports/allTop_mostCommonR1Starts_for_unrecognizedBarcodes.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << "SeekDeep rBind --recursive --depth 1 --contains top_mostCommonR2Starts_for_unrecognizedBarcodes.tab.txt --delim tab --header --out reports/allTop_mostCommonR2Starts_for_unrecognizedBarcodes.tab.txt  --overWrite" << std::endl;
	chmod(combineExtractionCmdOpts.outFilename_.c_str(),
			S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	//run file analysis file
	OutOptions runAnalysisOpts(
			bib::files::make_path(analysisSetup.dir_, "runAnalysis.sh"));
	std::ofstream runAnalysisFile;
	openTextFile(runAnalysisFile, runAnalysisOpts);
	runAnalysisFile << "#!/usr/bin/env bash" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "##run all parts of the pipeline" << std::endl;
	runAnalysisFile
			<< "##except start the server which is done by running ./startServerCmd.sh"
			<< std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "numThreads=1" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "if [[ $# -eq 1 ]]; then" << std::endl;
	runAnalysisFile << "	numThreads=$1" << std::endl;
	runAnalysisFile << "fi" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile extractorCmds.txt      --numThreads $numThreads --raw"
			<< std::endl;
	runAnalysisFile << "./combineExtractionCountsCmd.sh"<< std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile qlusterCmds.txt        --numThreads $numThreads --raw"
			<< std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile processClusterCmds.txt --numThreads $numThreads --raw"
			<< std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile genConfigCmds.txt      --numThreads $numThreads --raw"
			<< std::endl;
	runAnalysisFile << "" << std::endl;
	//make file executable
	chmod(runAnalysisOpts.outFilename_.c_str(),
			S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	if (foundErrors) {
		std::cerr << bib::bashCT::flashing << bib::bashCT::red << bib::bashCT::bold
				<< "ERRORS FOUND!!" << bib::bashCT::reset << std::endl;
		std::cerr << "Read Warnings in "
				<< bib::bashCT::boldRed(wrningsOpts.outFilename_.string()) << std::endl;
	}

	return 0;
}
//

}// namespace bibseq
