//  SeekDeepUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/06/24.
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2019 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "SeekDeep/utils.h"

#include <njhseq/ProgramRunners.h>

namespace njhseq {

SeekDeepUtilsRunner::SeekDeepUtilsRunner() :
		njh::progutils::ProgramRunner(
				{ addFunc("dryRunQualityFiltering", dryRunQualityFiltering, false),
					addFunc("runMultipleCommands",    runMultipleCommands, false),
					addFunc("setupTarAmpAnalysis", setupTarAmpAnalysis, false),
					addFunc("replaceUnderscores", replaceUnderscores, false),
				  addFunc("rBind", ManipulateTableRunner::rBind, false),
					addFunc("genTargetInfoFromGenomes", genTargetInfoFromGenomes, false),
					addFunc("benchmarkControlMixtures", benchmarkControlMixtures, false),
					addFunc("gatherInfoOnTargetedAmpliconSeqFile", gatherInfoOnTargetedAmpliconSeqFile, false),
					addFunc("getPossibleSampleNamesFromRawInput", getPossibleSampleNamesFromRawInput, false),
					addFunc("SampleBarcodeFileToSeekDeepInput", SampleBarcodeFileToSeekDeepInput, false),
					addFunc("primersToFasta", primersToFasta, false),
				}, //
				"SeekDeepUtils") {
}

//
int SeekDeepUtilsRunner::primersToFasta(const njh::progutils::CmdArgs & inputCommands) {
	TarAmpAnalysisSetup::TarAmpPars pars;

	auto outOptions = SeqIOOptions::genFastaOut(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.idFile, "--primers", "Primers file", true);
	setUp.processWritingOptions(outOptions.out_);
	setUp.finishSetUp(std::cout);


	PrimersAndMids primers(pars.idFile);
	primers.initPrimerDeterminator();
	SeqOutput writer(outOptions);
	writer.openOut();
	for(const auto & primer : primers.pDeterminator_->primers_){
		{
			uint32_t fwdCount = 0;
			for(const auto & fwd : primer.second.fwds_){
				std::string name = primer.first + "-fwd";

				if(primer.second.fwds_.size() > 1){
					name += njh::leftPadNumStr<uint32_t>(fwdCount, primer.second.fwds_.size());
				}
				seqInfo out(name, fwd.primer_);
				writer.write(out);
				++fwdCount;
			}
		}

		{
			uint32_t revCount = 0;
			for(const auto & rev : primer.second.revs_){
				std::string name = primer.first + "-rev";

				if(primer.second.fwds_.size() > 1){
					name += njh::leftPadNumStr<uint32_t>(revCount, primer.second.revs_.size());
				}
				seqInfo out(name, rev.primer_);
				writer.write(out);
				++revCount;
			}
		}
	}

	return 0;
}


int SeekDeepUtilsRunner::getPossibleSampleNamesFromRawInput(const njh::progutils::CmdArgs & inputCommands) {
	TarAmpAnalysisSetup::TarAmpPars pars;
	auto tabOutOpts = TableIOOpts::genTabFileOut("", true);
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(pars.inputDir, "--inputDir", "Input Directory", true);
	setUp.setOption(pars.inputFilePat, "--inputFilePat", "Input File Pat");
	setUp.setOption(pars.technology, "--technology", "Technology");
	setUp.setOption(pars.ignoreSamples, "--ignoreSamples", "Ignore these Samples if found");
	setUp.setOption(pars.replicatePattern, "--replicatePattern", "Replicate Pattern to match on to indicate replicates when guessing sample names, should have two groups e.g. (.*)(-rep.*)");
	setUp.processWritingOptions(tabOutOpts.out_);
	setUp.finishSetUp(std::cout);

	auto sampleNames = GuessPossibleSamps(pars);
	sampleNames.outPutContents(std::cout, "\t");
	return 0;
}

int SeekDeepUtilsRunner::genTargetInfoFromGenomes(const njh::progutils::CmdArgs & inputCommands) {
	extractBetweenSeqsPars pars;
	uint32_t minOverlap  = 10;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.primersFile, "--primers", "A file that contains three columns, target,forwardPrimer,reversePrimer 5` to 3` directions, same file as the input to SeekDeep", true);
	pars.verbose_ = setUp.pars_.verbose_;
	pars.debug_ = setUp.pars_.debug_;
	pars.setUpCoreOptions(setUp, true);
	setUp.setOption(minOverlap, "--minOverlap", "Minimum overlap for stitching");
	setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("bowtie2");
	njh::sys::requireExternalProgramThrow("samtools");

	PrimersAndMids ids(pars.primersFile);

	if(0 == ids.getTargets().size() ){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << pars.primersFile << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();



	extractBetweenSeqs(ids, pars);

	setUp.startARunLog(pars.outputDirPars.dirName_.string());


	auto forSeekDeepDir = njh::files::make_path(pars.outputDirPars.dirName_, "forSeekDeep");

	njh::files::makeDir(njh::files::MkdirPar{forSeekDeepDir});
	auto refSeqsDir = njh::files::make_path(forSeekDeepDir, "refSeqs");
	njh::files::makeDir(njh::files::MkdirPar{refSeqsDir});
	OutOptions lenCutOffsOpts(njh::files::make_path(forSeekDeepDir, "lenCutOffs.txt"));
	OutputStream lenCutOffsOut(lenCutOffsOpts);
	lenCutOffsOut << "target\tminlen\tmaxlen" << "\n";

	OutOptions overlapStatusOpts(njh::files::make_path(forSeekDeepDir, "overlapStatuses.txt"));
	OutputStream overlapStatusOut(overlapStatusOpts);
	overlapStatusOut << "target\tstatus" << "\n";


	for(const auto & tar : ids.getTargets()){
		auto primersRemovedFnp = njh::files::make_path(pars.outputDirPars.dirName_, tar, tar + "_primersRemoved.fasta");
		auto extractedSeqsFnp = njh::files::make_path(pars.outputDirPars.dirName_, tar, tar + ".fasta");
		auto primersRemovedFinalFnp = njh::files::make_path(refSeqsDir, tar + ".fasta");
		if(bfs::exists(primersRemovedFnp)){
			{
				auto finalSeqOpts = SeqIOOptions::genFastaInOut(primersRemovedFnp, primersRemovedFinalFnp);
				SeqIO reader(finalSeqOpts);
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
						<< "\t" << (minlen > pars.lenCutOffSizeExpand ? minlen - pars.lenCutOffSizeExpand : 0)
						<< "\t" << maxlen + pars.lenCutOffSizeExpand << std::endl;
				uint32_t finalMaxSize = maxlen + pars.barcodeSize;
				uint32_t finalMinSize = minlen + pars.barcodeSize;

				uint32_t maxInsertSize = 2* pars.pairedEndLength - minOverlap;

				std::string status = "";
				VecStr statuses;
				if(finalMaxSize > maxInsertSize){
					statuses.emplace_back("NoOverLap");
				}else{
					if(finalMaxSize >= pars.pairedEndLength){
						statuses.emplace_back("R1EndsInR2");
					}
					if(finalMinSize < pars.pairedEndLength){
						if(!statuses.empty()){
							statuses.emplace_back("PerfectOverlap");
						}
						statuses.emplace_back("R1BeginsInR2");
					}
				}
				overlapStatusOut << tar << "\t" << njh::conToStr(statuses, ",") << std::endl;
			}
		}else{
			std::cerr << "Warning, no sequences extracted for " << tar << std::endl;
		}
	}


	return 0;
}


int SeekDeepUtilsRunner::dryRunQualityFiltering(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	std::string qualWindow;
	uint32_t qualityWindowLength;
	uint32_t qualityWindowStep;
	uint8_t qualityWindowThres;
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
	std::cout << njh::bashCT::boldBlack("Sliding Quality Window") << std::endl;
	std::cout << "Window Size: " << qualityWindowLength << ", Window Step: "
			<< qualityWindowStep << ", Window Thresdhold: " << qualityWindowThres
			<< std::endl;
	std::cout << "FailedWindow: " << getPercentageString(failedQualWindow, count)
			<< std::endl;
	std::cout << std::endl;
	std::cout << njh::bashCT::boldBlack("Quality Fraction above cut off")
			<< std::endl;
	std::cout << "Q" << qualCheck << ">" << qualCheckCutOff << std::endl;
	std::cout << "FailedCheck: " << getPercentageString(failedQualCheck, count)
			<< std::endl;
	auto qualCheckStats = getStatsOnVec(qualChecks);
	table qualChecksTab(qualCheckStats, VecStr { "stat", "value" });
	qualChecksTab.outPutContentOrganized(std::cout);
	return 0;
}


inline std::vector<njh::sys::RunOutput> runCmdsThreadedQueue(
		const std::vector<std::string> & cmds, uint32_t numThreads, bool verbose,
		bool debug) {
	if (debug) {
		for (const auto & cmd : cmds) {
			std::cout << cmd << std::endl;
		}
		exit(1);
	}

	std::mutex allOutputsMut;
	std::mutex coutMut;
	std::vector<njh::sys::RunOutput> ret;
	njh::concurrent::LockableQueue<std::string> cmdsQueue(cmds);
	std::function<void()> runCmds = [&coutMut, &allOutputsMut, &ret, &verbose,&cmdsQueue]() {
		std::string cmd = "";
		uint32_t cmdNum = 0;
		std::vector<njh::sys::RunOutput> currentRets;
		bool run = cmdsQueue.getVal(cmd);
		while(run) {
//		while(cmdsQueue.getVal(cmd)) {
//			auto toks = tokenizeString(cmd, " ");
//			OutOptions opts(bfs::path{toks[0]});
//			opts.overWriteFile_ = true;
//			OutputStream out(opts);
//			sleep(njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]));
//			auto currentOut = njh::sys::run(std::vector<std::string> {cmd});
//			if(0 == cmdNum) {
//				auto currentOut = njh::sys::run(std::vector<std::string> {"sleep 1"});
//			}
			if(verbose) {
				std::lock_guard<std::mutex> lock(coutMut);
				auto tId = std::this_thread::get_id();
				std::cout << "Thread: " << tId << std::endl;
				std::cout << "\tRunning: " << cmd << std::endl;
			}
			auto currentOut = njh::sys::run(std::vector<std::string> {cmd});
			if(verbose) {
				auto tId = std::this_thread::get_id();
				std::cout << "\tThread: " << tId << std::endl;
				std::cout << "\tDone running: " << cmd << std::endl;
			}
			{
				if(verbose) {
					auto tId = std::this_thread::get_id();
					std::cout << "\tThread: " << tId << std::endl;
					std::cout << "\tInserting Results from: " << cmd << std::endl;
				}
				currentRets.emplace_back(currentOut);
			}
			run = cmdsQueue.getVal(cmd);
			++cmdNum;
		}
		{
			std::lock_guard<std::mutex> lock(allOutputsMut);
			for(const auto & runOut : currentRets){
				ret.emplace_back(runOut);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(runCmds, numThreads);
	return ret;
}


int SeekDeepUtilsRunner::runMultipleCommands(
		const njh::progutils::CmdArgs & inputCommands) {
	std::string filename = "";
	std::string logFile = "";
	bfs::path logDir = "./";
	uint32_t numThreads = 1;
	bool raw = false;
	bool noFilesInReplacementToks = false;
	std::string additionalFields = "";
	bool replaceFields = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(logDir, "--logDir", "Directory to create log file");
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(logFile, "--logFile", "Name of a file to log the output of the commands, this will cause --logDir to be ignored");
	setUp.setOption(filename, "--cmdFile",
			"Name of the file, first line is command with REPLACETHIS, the next lines are the cmd to run with that line replacing REPLACETHIS",
			true);
	setUp.setOption(raw, "--raw",
			"If if the file is simply just a list of commands to run");
	setUp.setOption(noFilesInReplacementToks, "--noFilesInReplacementToks",
			"The replacement tokens are by default checked to see if there are files and more replacement strings are read in where each line is a replacement, this turns off that behavior");
	setUp.setOption(additionalFields, "--additionalFields", "Additional fields to be replaced in the CMD: line, need to be given in format Key1:values;Key2:values");
	setUp.setOption(replaceFields, "--replaceFields", "Replace Fields");

	setUp.processWritingOptions();

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
		logFile = njh::files::make_path(logDir, bfs::path(filename).filename().replace_extension("").string() + "_TODAY_Log.json").string();
	}
	std::ofstream outFile;

	if (!setUp.pars_.debug_) {
		logFile = njh::replaceString(logFile, "TODAY", njh::getCurrentDate());
		if(!setUp.pars_.ioOptions_.out_.overWriteFile_){
			uint32_t number = 1;
			while(bfs::exists(logFile)){
				logFile = njh::files::nameAppendBeforeExt(logFile, estd::to_string(number)).string();
				++number;
			}
		}
		openTextFile(outFile, logFile, ".json", setUp.pars_.ioOptions_.out_);
	}

	std::ifstream inFile(filename);
	if (!inFile) {
		std::stringstream ss;
		ss << njh::bashCT::bold << "Error in opening "
				<< njh::bashCT::boldRed(filename) << std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::string cmd;
	VecStr rawPrepCmds;
	VecStr finalPrepCmds;
	VecStr rawPostCmds;
	VecStr finalPostCmds;

	VecStr cmds;

	if (raw) {
		std::string line;
		while (njh::files::crossPlatGetline(inFile, line)) {
			if ("" == line || allWhiteSpaceStr(line) || njh::beginsWith(line, "#")) {
				continue;
			}
			cmds.emplace_back(line);
		}
	} else {
		std::string line;
		std::map<std::string, VecStr> replacements;
		while (njh::files::crossPlatGetline(inFile, line)) {
			if ("" == line || allWhiteSpaceStr(line) || njh::beginsWith(line, "#")) {
				continue;
			}

			auto colonPos = line.find(":");

			if (std::string::npos == colonPos) {
				std::stringstream ss;
				ss << "Error in processing line: " << njh::bashCT::boldRed(line)
						<< std::endl;
				ss << "Need at least one colon" << std::endl;
				throw std::runtime_error { ss.str() };
			}
			VecStr toks { line.substr(0, colonPos), line.substr(colonPos + 1) };
			if (toks.front() == "CMD") {
				cmd = toks.back();
			}else if (toks.front() == "PREP") {
				rawPrepCmds.emplace_back(toks.back());
			} else if (toks.front() == "POST") {
				rawPostCmds.emplace_back(toks.back());
			} else {
				if(noFilesInReplacementToks){
					replacements[toks.front()] = tokenizeString(toks.back(), ",");
				}else{
					auto repToks = tokenizeString(toks.back(), ",");
					VecStr finalReplacements;
					for(const auto & tok : repToks){
						addOtherVec(finalReplacements, getInputValues(tok, ","));
					}
					if(njh::in(toks.front(), replacements) && !replaceFields){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "already have key " << toks.front()<< "\n";
						throw std::runtime_error{ss.str()};
					}
					replacements[toks.front()] = finalReplacements;
				}
			}
		}
		if("" != additionalFields){
			auto addFieldToks = tokenizeString(additionalFields, ";");
			for(const auto & addFieldTok : addFieldToks){
				auto colonPos = addFieldTok.find(":");

				if (std::string::npos == colonPos) {
					std::stringstream ss;
					ss << "Error in processing argument: " << njh::bashCT::boldRed(addFieldTok)
							<< std::endl;
					ss << "Need at least one colon" << std::endl;
					throw std::runtime_error { ss.str() };
				}
				VecStr toks { addFieldTok.substr(0, colonPos), addFieldTok.substr(colonPos + 1) };
				if(noFilesInReplacementToks){
					replacements[toks.front()] = tokenizeString(toks.back(), ",");
				}else{
					auto repToks = tokenizeString(toks.back(), ",");
					VecStr finalReplacements;
					for(const auto & tok : repToks){
						addOtherVec(finalReplacements, getInputValues(tok, ","));
					}
					if(njh::in(toks.front(), replacements) && !replaceFields){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "already have key " << toks.front()<< "\n";
						throw std::runtime_error{ss.str()};
					}
					replacements[toks.front()] = finalReplacements;
				}
			}
		}
		for (const auto & r : replacements) {
			if (cmds.empty()) {
				for (const auto & subR : r.second) {
					cmds.emplace_back(njh::replaceString(cmd, r.first, subR));
				}
			} else {
				VecStr newCmds;
				for (const auto & c : cmds) {
					for (const auto & subR : r.second) {
						newCmds.emplace_back(njh::replaceString(c, r.first, subR));
					}
				}
				cmds = newCmds;
			}
		}
		//prep cmds
		if(!rawPrepCmds.empty()){
			finalPrepCmds = rawPrepCmds;
			for(const auto & r : replacements){
				VecStr newPreps;
				for(const auto & prep : finalPrepCmds){
					if(std::string::npos != prep.find(r.first)){
						for(const auto & subR : r.second){
							newPreps.emplace_back(njh::replaceString(prep, r.first, subR));
						}
					}else{
						newPreps.emplace_back(prep);
					}
				}
				finalPrepCmds = newPreps;
			}
		}
		//post cmds
		if(!rawPostCmds.empty()){
			finalPostCmds = rawPostCmds;
			for(const auto & r : replacements){
				VecStr newPreps;
				for(const auto & prep : finalPostCmds){
					if(std::string::npos != prep.find(r.first)){
						for(const auto & subR : r.second){
							newPreps.emplace_back(njh::replaceString(prep, r.first, subR));
						}
					}else{
						newPreps.emplace_back(prep);
					}
				}
				finalPostCmds = newPreps;
			}
		}
	}

	if (setUp.pars_.verbose_) {
		printVector(cmds, "\n");
	}
	Json::Value allLog;

	if(!finalPrepCmds.empty()){
		uint32_t prepCmdCount = 1;
		for(const auto & prepCmd : finalPrepCmds){
			if(setUp.pars_.debug_){
				std::cout << prepCmd << std::endl;
			}else{
				auto prepOut = njh::sys::run(VecStr{prepCmd});
				if(!prepOut.success_){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "in running prep cmd: " << prepCmd << "\n";
					ss << prepOut.stdOut_ << "\n";
					ss << prepOut.stdErr_ << "\n";
					throw std::runtime_error{ss.str()};
				}
				allLog[njh::pasteAsStr("00_PREP_CMD_", prepCmdCount)] = prepOut.toJson();
				++prepCmdCount;
			}
		}
	}

	auto allRunOutputs = runCmdsThreadedQueue(cmds, numThreads, setUp.pars_.verbose_, setUp.pars_.debug_);

	allLog["totalTime"] = setUp.timer_.totalTime();
	allLog["cmdsfile"] = njh::json::toJson(bfs::absolute(filename));
	allLog["numThreads"] = njh::json::toJson(numThreads);

	auto & cmdsLog = allLog["cmdsLog"];
	for (const auto & out : allRunOutputs) {
		cmdsLog.append(out.toJson());
	}

	if(!finalPostCmds.empty()){
		uint32_t postCmdCount = 1;
		for(const auto & postCmd : finalPostCmds){
			if(setUp.pars_.debug_){
				std::cout << postCmd << std::endl;
			}else{
				auto postOut = njh::sys::run(VecStr{postCmd});
				if(!postOut.success_){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "in running post cmd: " << postCmd << "\n";
					ss << postOut.stdOut_ << "\n";
					ss << postOut.stdErr_ << "\n";
					throw std::runtime_error{ss.str()};
				}
				allLog[njh::pasteAsStr("ZZ_POST_CMD_", postCmdCount)] = postOut.toJson();
				++postCmdCount;
			}
		}
	}

	if (!setUp.pars_.debug_) {
		outFile << allLog << std::endl;
	}

	return 0;
}

int SeekDeepUtilsRunner::replaceUnderscores(
		const njh::progutils::CmdArgs & inputCommands) {
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

	auto files = njh::files::listAllFiles(dir.string(), false, {
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
			auto outPath = njh::files::make_path(parDir, bName);
			renameKey[f.first] = outPath;
		}
	}
	if (setUp.pars_.debug_) {
		for (const auto & rk : renameKey) {
			std::cout << njh::bashCT::green << rk.first << njh::bashCT::reset
					<< " -> " << njh::bashCT::blue << rk.second << njh::bashCT::reset
					<< std::endl;
		}
		return 0;
	}
	for (const auto & rk : renameKey) {
		if (keepOriginals) {
			bfs::copy_file(rk.first, rk.second);
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
		while (njh::files::crossPlatGetline(ss, line)) {
			if (std::string::npos != line.find("Writing histogram files") && log) {
				break;
			}
			if (log) {
				if (std::string::npos != line.find(":")) {
					line = njh::replaceString(line, "[FLASH]", "");
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
			if (!njh::in(head, output)) {
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
				njh::replaceString(output["Percentcombined"], "%", ""));
	}
	return flashOutTab;
}




//

}// namespace njhseq
