
//  SeekDeepUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/06/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
    
#include "SeekDeepUtilsRunner.hpp"

    
namespace bibseq {

SeekDeepUtilsRunner::SeekDeepUtilsRunner()
    : bib::progutils::programRunner(
    		{
				addFunc("dryRunQaulityFiltering", dryRunQaulityFiltering, false),
				addFunc("runMultipleCommands", runMultipleCommands, false),
				addFunc("setupPairedEndDir", setupPairedEndDir, false)
    		},//
                    "SeekDeepUtils") {}
                    
int SeekDeepUtilsRunner::dryRunQaulityFiltering(const bib::progutils::CmdArgs & inputCommands){
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
  setUp.setOption(qualCheckCutOff, "--qualCheckCutOff", "Cut Off for fraction of bases above qual check"	);
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
	while(reader.readNextRead(read)){
		read.setBaseCountOnQualCheck(qualCheck);
		std::cout << "Currently on " << count << "\r";
		std::cout.flush();
		qualChecks.emplace_back(read.fractionAboveQualCheck_);
		if(read.fractionAboveQualCheck_ < qualCheckCutOff){
			++failedQualCheck;
		}
		if(!seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
        qualityWindowStep, read.seqBase_.qual_)){
			++failedQualWindow;
		}
		++count;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << bib::bashCT::boldBlack("Sliding Quality Window") << std::endl;
	std::cout << "Window Size: "  << qualityWindowLength << ", Window Step: " << qualityWindowStep << ", Window Thresdhold: " << qualityWindowThres << std::endl;
	std::cout << "FailedWindow: " << getPercentageString(failedQualWindow, count) << std::endl;
	std::cout << std::endl;
	std::cout << bib::bashCT::boldBlack("Quality Fraction above cut off") << std::endl;
	std::cout << "Q"<< qualCheck << ">" << qualCheckCutOff << std::endl;
	std::cout << "FailedCheck: " << getPercentageString(failedQualCheck, count) << std::endl;
	auto qualCheckStats = getStatsOnVec(qualChecks);
	table qualChecksTab(qualCheckStats, VecStr{"stat", "value"});
	qualChecksTab.outPutContentOrganized(std::cout);
	return 0;
}



int SeekDeepUtilsRunner::runMultipleCommands(
		const bib::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	std::string logFile = "";
	uint32_t numThreads = 1;
	bool raw = false;
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(logFile, "--logFile",
			"Name of a file to log the output of the commands");
	setUp.setOption(filename, "--cmdFile",
			"Name of the file, first line is command with REPLACETHIS, the next lines are the cmd to run with that line replacing REPLACETHIS",
			true);
	setUp.setOption(raw, "--raw,-r",
			"If if the file is simply just a list of commands to run");
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
		logFile = bfs::path(filename).filename().replace_extension("").string() + "Log.json";
	}
	std::ofstream outFile;
	openTextFile(outFile, logFile, ".json", setUp.pars_.ioOptions_.out_);
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
				ss << "Error in processing line: " << bib::bashCT::boldRed(line) << std::endl;
				ss << "Need at least one colon" << std::endl;
				throw std::runtime_error{ss.str()};
			}
			VecStr toks{line.substr(0,colonPos), line.substr(colonPos + 1)};
			if (toks.front() == "CMD") {
				cmd = toks.back();
			} else {
				replacements[toks.front()] = tokenizeString(toks.back(), ",");
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
	outFile << allLog << std::endl;
	return 0;
}







std::map<std::string, PrimersAndMids::Target::lenCutOffs> readInLenCutOffs(const bfs::path & lenCutOffsFnp){
	std::map<std::string, PrimersAndMids::Target::lenCutOffs> ret;
	table lenCutTab = table(lenCutOffsFnp.string(), "whitespace", true);
	bib::for_each(lenCutTab.columnNames_,
			[](std::string & str) {stringToLower(str);});
	lenCutTab.setColNamePositions();
	if (!bib::in(std::string("target"), lenCutTab.columnNames_)
			|| !bib::in(std::string("minlen"), lenCutTab.columnNames_)
			|| !bib::in(std::string("maxlen"), lenCutTab.columnNames_)) {
		std::stringstream ss;
		ss << "need to have columns " << "target, minlen, and maxlen"
				<< " when reading in a table for multiple cut off lengths"
				<< std::endl;
		ss << "only have " << vectorToString(lenCutTab.columnNames_, ",")
				<< std::endl;
		throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
	}
	for (const auto & row : lenCutTab.content_) {
		ret.emplace(row[lenCutTab.getColPos("target")],
				PrimersAndMids::Target::lenCutOffs { bib::lexical_cast<uint32_t>(
						row[lenCutTab.getColPos("minlen")]), bib::lexical_cast<uint32_t>(
						row[lenCutTab.getColPos("maxlen")]) });
	}
	return ret;
}




int SeekDeepUtilsRunner::setupPairedEndDir(const bib::progutils::CmdArgs & inputCommands){
	bfs::path samplesNamesFnp = "";
	bfs::path outDir = "";
	bfs::path inputDir = "";
	bfs::path groupMeta = "";
	bfs::path idFile = "";
	bfs::path lenCutOffsFnp = "";
	bfs::path refSeqsDir = "";
	uint32_t maxOverlap = 250;
	uint32_t numThreads = 1;
	uint32_t r1Trim = std::numeric_limits<uint32_t>::max();
	uint32_t r2Trim = std::numeric_limits<uint32_t>::max();
	bool noQualTrim = false;

	std::string extraExtractorCmds = "";
	std::string extraQlusterCmds = "";
	std::string extraProcessClusterCmds = "";


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(samplesNamesFnp, "--samples", "A file containing the samples names, columns should go target,sample,pcr1,(optional)pcr2 etc", true);
	setUp.setOption(outDir, "--outDir", "Directory to setup for analysis", true);
	setUp.setOption(inputDir, "--inputDir", "Input Directory of raw data files", true);
	setUp.setOption(maxOverlap, "--maxOverlap", "Max overlap allowed in stitcher");
	setUp.setOption(numThreads, "--numThreads", "Number of CPUs to use");

	setUp.setOption(extraExtractorCmds, "--extraExtractorCmds", "Extra extractor cmds to add to the defaults");
	setUp.setOption(extraQlusterCmds, "--extraQlusterCmds", "Extra qluster cmds to add to the defaults");
	setUp.setOption(extraProcessClusterCmds, "--extraProcessClusterCmds", "Extra process clusters cmds to add to the defaults");

	setUp.setOption(noQualTrim, "--noQualTrim", "No Quality Trim");

	setUp.setOption(r1Trim, "--r1Trim", "Trim R1 Reads to this length");
	setUp.setOption(r2Trim, "--r2Trim", "Trim R2 Reads to this length");

	setUp.setOption(groupMeta, "--groupMeta", "Group Metadata");
	setUp.setOption(idFile, "--idFile", "ID file containing primer and possible additional MIDs", true);
	setUp.setOption(lenCutOffsFnp, "--lenCutOffs", "A file with 3 columns, target,minlen,maxlen to supply length cut off specifically for each target");
	setUp.setOption(refSeqsDir, "--refSeqsDir", "A directory of fasta files where each file is named with the input target names");

	setUp.finishSetUp(std::cout);

	std::string stitcher = "flash";
#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
	//apple zcat is stupid and requires files end with .Z because why not
	//using brew install gnutls instead
	std::string zcatCmd = "gzcat";
#else
	std::string zcatCmd = "zcat";
#endif
	VecStr warnings;
	if (!bib::sys::hasSysCommand(stitcher)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, need to have " << stitcher
				<< " in path" << "\n";
		warnings.emplace_back(ss.str());
	}
	if (!bib::sys::hasSysCommand(zcatCmd)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, need to have " << zcatCmd
				<< " in path" << "\n";
		warnings.emplace_back(ss.str());
	}
	if (bfs::exists(outDir)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, out directory " << "\""
				<< bib::bashCT::boldRed(outDir.string()) << "\"" << " already exists"
				<< ", must remove by hand will not overwrite\n";
		warnings.emplace_back(ss.str());
	}
	if (!bfs::exists(inputDir)) {
		std::stringstream ss;
		ss << bib::bashCT::boldRed(inputDir.string()) << " doesn't exist" << "\n";
		warnings.emplace_back(ss.str());
	}
	if (!bfs::exists(idFile)) {
		std::stringstream ss;
		ss << bib::bashCT::boldRed(idFile.string()) << " doesn't exist" << "\n";
		warnings.emplace_back(ss.str());
	}

	if(!warnings.empty()){
		std::stringstream ss;
		if(1 == warnings.size()){
			ss << __PRETTY_FUNCTION__ << ": There is " << warnings.size() << " error." << "\n";
		}else{
			ss << __PRETTY_FUNCTION__ << ": There are " << warnings.size() << " errors." << "\n";
		}
		for(const auto & warn : warnings){
			ss << warn << std::endl;
		}
		throw std::runtime_error{ss.str()};
	}

	bool foundErrors = false;
	std::stringstream errorOutput;

	TarAmpPEAnalysisSetup analysisSetup(outDir);
	setUp.startARunLog(bib::appendAsNeededRet(analysisSetup.dir_.string(), "/"));

	analysisSetup.addSamplesNames(samplesNamesFnp);

	analysisSetup.writeSampleNamesFile();

	/**@todo consider putting a check for correct primers and mids in the input files*/

	PrimersAndMids ids(idFile);

	auto targets = ids.getTargets();
	auto sampNamesTargets = analysisSetup.getTargets();
	bib::sort(targets);
	bib::sort(sampNamesTargets);

	VecStr idMissingTargets;
	VecStr sampNamesTargetsTarMissing;
	for(const auto & tar : sampNamesTargets){
		if(!bib::in(tar, targets)){
			idMissingTargets.emplace_back(tar);
		}
	}
	for(const auto & tar : targets){
		if(!bib::in(tar, sampNamesTargets)){
			sampNamesTargetsTarMissing.emplace_back(tar);
		}
	}
	if(!idMissingTargets.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, missing the following targets from the id file" << idFile << "\n";
		ss << "Targets: " << bib::conToStr(idMissingTargets, ", ") << "\n";
		ss << "AvailableTargs: " << bib::conToStr(targets, ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}

	if(!sampNamesTargetsTarMissing.empty()){
		foundErrors = true;
		errorOutput << __PRETTY_FUNCTION__ << ": warning, missing the following targets from the sample name file" << samplesNamesFnp << "\n";
		errorOutput << "Targets: " << bib::conToStr(sampNamesTargetsTarMissing, ", ") << "\n";
		errorOutput << "AvailableTargs: " << bib::conToStr(sampNamesTargets, ", ") << "\n";
	}

	if("" != groupMeta){
		analysisSetup.addGroupingMetaData(groupMeta);

		if (!analysisSetup.groupMetaData_->missingSamples_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, the following samples were in the meta data file, "
					<< groupMeta << " but were missing from the input sample names file, "
					<< samplesNamesFnp << std::endl;
			for (const auto & samp : analysisSetup.groupMetaData_->missingSamples_) {
				errorOutput << "\tSample: " << samp << "\n";
			}
		}
		if (!analysisSetup.groupMetaData_->missingMetaForSamples_.empty()) {
			foundErrors = true;
			errorOutput
								<< "Warning, the following samples were in the input samples name file, "
								<< samplesNamesFnp << " but were missing from the meta data file, "
								<< groupMeta << std::endl;
			for(const auto & samp : analysisSetup.groupMetaData_->missingMetaForSamples_){
				errorOutput << "\tSample: " << samp << "\n";
			}
		}
		bfs::copy(analysisSetup.groupMetaData_->groupingsFile_,
				bib::files::make_path(analysisSetup.infoDir_, "groupMeta.tab.txt"));
	}

	if("" != refSeqsDir){
		/**@todo add ref seq functions and checks*/
		if (!bfs::exists(refSeqsDir)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, directory "
					<< bib::bashCT::boldRed(refSeqsDir.string()) << " doesn't exists"
					<< "\n";
			;
			throw std::runtime_error { ss.str() };
		}
		auto fastaFiles = bib::files::listAllFiles(refSeqsDir.string(), false, {
				std::regex { R"(.*\.fasta$)" } });
		VecStr found;
		VecStr missing;
		VecStr noMatching;
		auto tars = ids.getTargets();
		for (const auto & ff : fastaFiles) {
			auto tarName = bfs::basename(ff.first);
			found.emplace_back(tarName);
			if (!bib::in(tarName, tars)) {
				noMatching.emplace_back(tarName);
			} else {
				auto inOpts = SeqIOOptions::genFastaIn(ff.first.string(), false);
				SeqInput reader(inOpts);
				auto seqs = reader.readAllReads<seqInfo>();
				if (seqs.size() == 1) {
					ids.targets_.at(tarName).addSingleRef(seqs.front());
				} else {
					ids.targets_.at(tarName).addMultileRef(seqs);
				}
			}
		}
		for (const auto & tar : tars) {
			if (!bib::in(tar, found)) {
				missing.emplace_back(tar);
			}
		}
		if(!missing.empty()){
			foundErrors = true;
			errorOutput
								<< "Warning, reference sequences were supplied but reference files for the following targets couldn't be found" << "\n";
			for(const auto & tar : missing){
				errorOutput << "\tTarget: " << tar << "\n";
			}
			errorOutput << "\tOptions are: " << bib::conToStr(found, ", ") << "\n";
		}
		if(!noMatching.empty()){
			foundErrors = true;
			errorOutput
								<< "Warning, reference sequences were supplied but there are no matching reference targets in project for the following" << "\n";
			for(const auto & tar : noMatching){
				errorOutput << "\tTarget: " << tar << "\n";
			}
			errorOutput << "\tOptions are: " << bib::conToStr(tars, ", ") << "\n";
		}
	}

	if("" != lenCutOffsFnp){
		auto multipleLenCutOffs = readInLenCutOffs(lenCutOffsFnp);
		for(const auto & cutOff : multipleLenCutOffs){
			if(!bib::in(cutOff.first, analysisSetup.samplesForTargets_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, in processing " << lenCutOffsFnp << "\n";
				ss << "Supplying len cut off for a target not in input" << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		VecStr missingLensForTargets;
		for (const auto & tar : ids.getTargets()) {
			if (!bib::in(tar, multipleLenCutOffs)) {
				missingLensForTargets.emplace_back(tar);
			} else {
				ids.targets_.at(tar).addLenCutOff(
						multipleLenCutOffs.at(tar).minLenChecker_.minLen_,
						multipleLenCutOffs.at(tar).maxLenChecker_.maxLen_);
			}
		}
		if (!missingLensForTargets.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, supplying length cut offs for targets but didn't supply any cuts for the following targets\n";
			for (const auto & tar : missingLensForTargets) {
				errorOutput << "\tTarget: " << tar << "\n";
			}
		}
	}

	//now write id files

	std::unordered_map<std::string, VecStr> targetsForReps;
	for (const auto & tars : analysisSetup.samplesForTargets_) {
		for (const auto & rep : tars.second.getReps()) {
			targetsForReps[rep].emplace_back(tars.first);
		}
	}
	std::vector<VecStr> tarCombos;
	for (auto & sampTars : targetsForReps) {
		bib::sort(sampTars.second);
		if (!bib::in(sampTars.second, tarCombos)) {
			tarCombos.emplace_back(sampTars.second);
		}
	}

	for (const auto & tarCombo : tarCombos) {
		auto collapse = bib::conToStr(tarCombo, "_");
		auto refs = ids.getRefSeqs(tarCombo);
		auto lens = ids.genLenCutOffs(tarCombo);
		if (!lens.empty()) {
			auto lensOutOpts = TableIOOpts::genTabFileOut(
					bib::files::make_path(analysisSetup.idsDir_,
							collapse + "_lenCutOffs.tab.txt"));
			lens.outPutContents(lensOutOpts);
		}
		if (!refs.empty()) {
			SeqOutput::write(refs,
					SeqIOOptions::genFastaOut(
							bib::files::make_path(analysisSetup.idsDir_, collapse).string()));
		}
		ids.writeIdFile(
				OutOptions(
						bib::files::make_path(analysisSetup.idsDir_, collapse + ".id.txt").string()),
				tarCombo);
	}



	auto files = bib::files::listAllFiles(inputDir.string(), false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;

	auto reps = analysisSetup.getReps();

	for(const auto & f : files){
		auto filename = f.first.filename().string();
		auto underPos = filename.find("_");
		if(0 == underPos){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
					<< f.first <<  ", shouldn't start with an _, can't determine sample name\n";
			throw std::runtime_error{ss.str()};
		}
		if(std::string::npos == underPos){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
					<< f.first <<  ", should contain an _ to separate sample name from mate pair designations\n";
			ss << "Examples: Samp1_R1.fastq Samp1_R2.fastq Samp2_R1.fastq ..." << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::string sampName = filename.substr(0, underPos);

		if(bib::in(sampName, reps) || bib::in("MID" + sampName, reps)){
			readPairs[sampName].emplace_back(f.first.string());
		}else{
			readPairsUnrecognized[sampName].emplace_back(f.first.string());
		}
	}

	std::unordered_map<std::string, std::pair<VecStr, VecStr>> readsByPairs;


	for (const auto & reads : readPairs) {
		for (const auto & read : reads.second) {
			auto filename = bfs::path(read).filename().string();
			auto lastUnderPos = filename.rfind("_");
			auto periodPos = filename.find(".", lastUnderPos);
			if(std::string::npos == lastUnderPos){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
						<< read <<  ", should contain an _ before the read mate designation\n";
				throw std::runtime_error{ss.str()};
			}
			if(std::string::npos == periodPos){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
						<< read <<  ", should contain an . after the read designation, normally for the file extension\n";
				throw std::runtime_error{ss.str()};
			}

			std::string pairNum = filename.substr(lastUnderPos + 1, periodPos - lastUnderPos - 1);
			while("R1" != pairNum
					&& "1" != pairNum
					&& "R2" != pairNum
					&& "2" != pairNum){
				//find next _ to try to determine read pairs designation
				periodPos = lastUnderPos;
				lastUnderPos = filename.rfind("_", lastUnderPos - 1);
				if(std::string::npos == lastUnderPos){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
							<< read <<  ", couldn't find the mate designation \n";
					ss << "Examples: Samp1_R1.fastq Samp1_R2.fastq Samp2_R1.fastq ..." << "\n";
					ss << "Examples: Samp1_R1_001.fastq Samp1_R2_001.fastq Samp2_R1_001.fastq ..." << "\n";
					throw std::runtime_error{ss.str()};
				}
				pairNum = filename.substr(lastUnderPos + 1, periodPos - lastUnderPos - 1);
			}
			if("R1" == pairNum || "1" == pairNum){
				readsByPairs[reads.first].first.emplace_back(read);
			}else if("R2" == pairNum || "2" == pairNum){
				readsByPairs[reads.first].second.emplace_back(read);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, cannot determine pair number for "
						<< filename << ", found designation of " << pairNum << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
	}
	//file checks
	for(auto & reads : readsByPairs){
		//check size
		if(reads.second.first.size() != reads.second.second.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, read paired file numbers for " << reads.first << " don't match up \n";
			ss << "R1 file number: " << reads.second.first.size() << ", files: " << bib::conToStr(reads.second.first, ", ") << "\n";
			ss << "R2 file number: " << reads.second.second.size() << ", files: " << bib::conToStr(reads.second.second, ", ") << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::vector<size_t> needToErrase;
		//sort
		bib::sort(reads.second.first);
		bib::sort(reads.second.second);
		for(const auto pos : iter::range(reads.second.first.size())){

			auto lastUnderPos1 = reads.second.first[pos].rfind("_");
			{
				auto periodPos = reads.second.first[pos].find(".", lastUnderPos1);
				std::string pairNum = reads.second.first[pos].substr(lastUnderPos1 + 1,
						periodPos - lastUnderPos1 - 1);
				while ("R1" != pairNum && "1" != pairNum) {
					//find next _ to try to determine read pairs designation
					periodPos = lastUnderPos1;
					lastUnderPos1 = reads.second.first[pos].rfind("_", lastUnderPos1 - 1);
					pairNum = reads.second.first[pos].substr(lastUnderPos1 + 1,
							periodPos - lastUnderPos1 - 1);
				}
			}
			auto lastUnderPos2 = reads.second.second[pos].rfind("_");
			{
				auto periodPos = reads.second.second[pos].find(".", lastUnderPos2);
				std::string pairNum = reads.second.second[pos].substr(lastUnderPos2 + 1,
						periodPos - lastUnderPos2 - 1);
				while ("R2" != pairNum && "2" != pairNum) {
					//find next _ to try to determine read pairs designation
					periodPos = lastUnderPos2;
					lastUnderPos2 = reads.second.second[pos].rfind("_",
							lastUnderPos2 - 1);
					pairNum = reads.second.second[pos].substr(lastUnderPos2 + 1,
							periodPos - lastUnderPos2 - 1);
				}
			}
			auto nameStub1 = reads.second.first[pos].substr(0, lastUnderPos1);
			auto nameStub2 = reads.second.second[pos].substr(0, lastUnderPos2);
			//check name
			if(nameStub1 != nameStub2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, read paired file names don't match up for " << reads.first << "\n";
				ss << "for " << reads.second.first[pos] << " and " << reads.second.second[pos] << "\n";
				ss << "stub1: " << nameStub1 << "\n";
				ss << "stub2: " << nameStub2 << "\n";
				throw std::runtime_error{ss.str()};
			}
			//check to see if they are empty
			if(0 == bfs::file_size(reads.second.first[pos]) ||
					0 == bfs::file_size(reads.second.second[pos])){
				needToErrase.emplace_back(pos);
			}
		}
		if(!needToErrase.empty()){
			//sort the position so the last positions come first so positions don't get invalidated after erasing (back -> front errasing)
			std::sort(needToErrase.rbegin(), needToErrase.rend());
			for(const auto pos : needToErrase){
				reads.second.first.erase(reads.second.first.begin() + pos);
				reads.second.second.erase(reads.second.second.begin() + pos);
			}
		}
	}


	auto keys = getVectorOfMapKeys(readsByPairs);
	bib::sort(keys);
	VecStr samplesExtracted;
	VecStr samplesEmpty;
	bib::concurrent::LockableQueue<std::string> filesKeys(keys);
	Json::Value logs;
	std::mutex logsMut;
	const std::string zcatTempCmdR1 = zcatCmd + " " + "{FILES} > "
			+  "{OUTPUT}_R1.fastq";
	const std::string zcatTempCmdR2 = zcatCmd + " " + "{FILES} > "
			+  "{OUTPUT}_R2.fastq";
	const std::string stitchCmd = stitcher + " " +
			"{OUTPUT}_R1.fastq " +
			"{OUTPUT}_R2.fastq " +
			"-o {OUTPUT} " +
			"--max-overlap " + estd::to_string(maxOverlap);

	auto extractStitchFiles =
			[&analysisSetup ,&samplesExtracted, &samplesEmpty,
			 &zcatTempCmdR1, &zcatTempCmdR2, &stitchCmd,
			 &filesKeys, &readsByPairs,
			 &logs, &logsMut,
			 &r1Trim, &r2Trim]() {
				std::string key = "";
				VecStr currentEmptySamps;
				VecStr currentSamplesExtracted;
				std::unordered_map<std::string, Json::Value> currentLogs;

				while(filesKeys.getVal(key)){
					//check to see if all the files were erased and therefore sample was just empty files
					if(readsByPairs.at(key).first.empty()){
						currentEmptySamps.emplace_back(key);
						continue;
					}
					auto zcatR1 = bib::replaceString(zcatTempCmdR1, "{FILES}", bib::conToStr(readsByPairs.at(key).first, " "));
					zcatR1 = bib::replaceString(zcatR1, "{OUTPUT}", bib::files::make_path(analysisSetup.dir_, key).string());
					auto zcatR2 = bib::replaceString(zcatTempCmdR2, "{FILES}", bib::conToStr(readsByPairs.at(key).second, " "));
					zcatR2 = bib::replaceString(zcatR2, "{OUTPUT}", bib::files::make_path(analysisSetup.dir_, key).string());
					auto curStitchCmd =  bib::replaceString(stitchCmd, "{OUTPUT}",  bib::files::make_path(analysisSetup.dir_, key).string());

					auto zcatR1Out = bib::sys::run({zcatR1});
					auto zcatR2Out = bib::sys::run({zcatR2});
					if(std::numeric_limits<uint32_t>::max() != r1Trim){
						auto finalFnp = bib::files::make_path(analysisSetup.dir_, key).string() + "_R1.fastq";
						auto tempFnp = bib::files::make_path(analysisSetup.dir_, "tempTrim_" + key).string() + "_R1.fastq";
						auto opts = SeqIOOptions::genFastqInOut(finalFnp,tempFnp,false);
						SeqIO reader(opts);
						reader.openIn();
						reader.openOut();
						seqInfo r1seq;
						while(reader.readNextRead(r1seq)){
							readVecTrimmer::trimToMaxLength(r1seq, r1Trim);
							reader.write(r1seq);
						}
						reader.closeIn();
						reader.closeOut();
						bfs::remove(finalFnp);
						bfs::rename(tempFnp, finalFnp);
					}

					if(std::numeric_limits<uint32_t>::max() != r2Trim){
						auto finalFnp = bib::files::make_path(analysisSetup.dir_, key).string() + "_R2.fastq";
						auto tempFnp =  bib::files::make_path(analysisSetup.dir_, "tempTrim_" + key).string() + "_R2.fastq";
						auto opts = SeqIOOptions::genFastqInOut(finalFnp,tempFnp,false);
						SeqIO reader(opts);
						reader.openIn();
						reader.openOut();
						seqInfo r2seq;
						while(reader.readNextRead(r2seq)){
							readVecTrimmer::trimToMaxLength(r2seq, r2Trim);
							reader.write(r2seq);
						}
						reader.closeIn();
						reader.closeOut();
						bfs::remove(finalFnp);
						bfs::rename(tempFnp, finalFnp);
					}
					auto curStitchCmdOut = bib::sys::run({curStitchCmd});
					currentLogs[key]["zcat1"] = zcatR1Out.toJson();
					currentLogs[key]["zcat2"] = zcatR2Out.toJson();
					currentLogs[key]["stitch"] = curStitchCmdOut.toJson();
					currentLogs[key]["name"] = key;
					currentSamplesExtracted.emplace_back(key);
				}
				{
					std::lock_guard<std::mutex> lock(logsMut);
					addOtherVec(samplesEmpty, currentEmptySamps);
					addOtherVec(samplesExtracted, currentSamplesExtracted);
					for(const auto & keyLog : currentLogs){
						logs[keyLog.first] = keyLog.second;

					}
				}
			};
	std::vector<std::thread> threads;
	for(uint32_t tNum = 0; tNum < numThreads; ++tNum){
		threads.emplace_back(std::thread(extractStitchFiles));
	}

	for(auto & th : threads){
		th.join();
	}
	std::ofstream outLog;
	openTextFile(outLog, OutOptions(bib::files::make_path(analysisSetup.logsDir_, "extractLog").string(), ".json"));
	outLog << logs << std::endl;

	//report for stitching
	table flashOutTab{VecStr{"SampleName", "TotalPairs", "CombinedPairs", "UncombinedPairs", "PercentCombined"}};
	VecStr noReadsStitched;
	for(const auto & samp : logs){
		if(setUp.pars_.debug_){
			std::cout << samp["name"].asString() << std::endl;;
		}
		std::stringstream ss;
		ss << samp["stitch"]["stdOut_"].asString();
		std::unordered_map<std::string, std::string> output;
		bool log = false;
		std::string line = "";
		while(bib::files::crossPlatGetline(ss, line)){
			if(std::string::npos != line.find("Writing histogram files") && log){
				break;
			}
			if(log){
				if(std::string::npos!= line.find(":")){
					line = bib::replaceString(line, "[FLASH]", "");
					line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
					auto toks = tokenizeString(line, ":");
					if (toks.size() != 2) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__
								<< ": Error in processing stitch output for sample: "
								<< samp["name"].asString() << " expected two values on line " << " line "
								<< " but got " << toks.size() << "\n";
						throw std::runtime_error { ss.str() };
					}
					output[toks[0]] = toks[1];
				}
			}
			if(std::string::npos != line.find("Read combination statistics:")){
				log = true;
			}
		}
		VecStr headers = {"Totalpairs", "Combinedpairs", "Uncombinedpairs", "Percentcombined"};
		for(const auto & head : headers){
			if(!bib::in(head, output)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__
						<< ": Error in processing stitch output for sample: "
						<< samp["name"].asString() << "\n";
				ss << "Coulnd't find: " << head << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
		auto combinedPairsCnt = bib::lexical_cast<uint32_t>(output["Combinedpairs"]);
		if(combinedPairsCnt == 0){
			noReadsStitched.emplace_back(samp["name"].asString());
		}
		flashOutTab.addRow(samp["name"].asString(),
				output["Totalpairs"],
				output["Combinedpairs"],
				output["Uncombinedpairs"],
				bib::replaceString(output["Percentcombined"], "%", ""));
	}
	flashOutTab.outPutContents(
			TableIOOpts::genTabFileOut(
					bib::files::make_path(analysisSetup.reportsDir_,
							"StitchStats.tab.txt"), true));

	if(!readPairsUnrecognized.empty()){
		foundErrors = true;
		errorOutput << "The following files were found but didn't match any input sample names in " << samplesNamesFnp << std::endl;
		for(const auto & samp : readPairsUnrecognized){
			errorOutput << "\tPossible Sample Name: " << samp.first << std::endl;
			for(const auto & sampFiles : samp.second){
				errorOutput << "\t\t" << sampFiles << std::endl;
			}
		}
	}

	if(!samplesEmpty.empty()){
		foundErrors = true;
		errorOutput << "The following files were found to be empty "<< std::endl;
		for(const auto & samp : samplesEmpty){
			errorOutput << "\tSample: " << samp << std::endl;
			for(const auto & sampFiles : readPairs.at(samp)){
				errorOutput << "\t\t" << sampFiles << std::endl;
			}
		}
	}
	if(!noReadsStitched.empty()){
		foundErrors = true;
		errorOutput << "The following samples had no reads stitched "<< std::endl;
		for(const auto & samp : noReadsStitched){
			errorOutput << "\tSample: " << samp << std::endl;
			for(const auto & sampFiles : readPairs.at(samp)){
				errorOutput << "\t\t" << sampFiles << std::endl;
			}
		}
	}




	VecStr samplesNotFound;
	for(const auto & rep : reps){

		if(!bib::in(rep, samplesExtracted) &&
				!bib::in(bib::replaceString(rep, "MID", ""), samplesExtracted)){
			samplesNotFound.emplace_back(rep);
		}
	}

	if(!samplesNotFound.empty()){
		foundErrors = true;
		errorOutput << "The following input files were not found"<< std::endl;
		for(const auto & samp : samplesNotFound){
			errorOutput << "\tSample: " << samp << std::endl;
		}
	}


	OutOptions wrningsOpts(bib::files::make_path(analysisSetup.dir_, "WARNINGS_PLEASE_READ.txt").string());
	if(foundErrors){
		std::ofstream outWarnings;
		openTextFile(outWarnings,wrningsOpts );
		outWarnings << errorOutput.str() << std::endl;;
	}

	//make population clustering directory
	bool multipleTargets = analysisSetup.samplesForTargets_.size() > 1;
	setUpSampleDirs(
			bib::files::make_path(analysisSetup.infoDir_, "sampNames.tab.txt").string(),
			bib::files::make_path(analysisSetup.dir_, "popClustering").string(),
			multipleTargets,
			setUp.pars_.verbose_);

	//extractor cmds
	std::string extractorCmdTemplate = setUp.commands_.masterProgram_ + " extractor --fastq {REP}.extendedFrags.fastq --dout {REP}_extraction --overWriteDir --illumina";
	std::string qualTrimTemplate = "--trimAtQual 2";
	std::string lenCutOffsTemplate = "--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt";
	std::string idTemplate = "--id info/ids/{TARS}.id.txt";
	std::string sampleNameTemplate = "--sampleName {REP}";



	//qluster cmds;
	std::string qlusterCmdTemplate = "cd {REP}_extraction && " + setUp.commands_.masterProgram_ + " qluster "
			"--fastq {TARGET}{MIDREP}.fastq  --illumina --qualThres 25,20 "
			"--alnInfoDir alnCache_{TARGET}{MIDREP} --markChimeras --overWriteDir "
			"--additionalOut ../popClustering/locationByIndex/{TARGET}.tab.txt "
			"--overWrite --dout {TARGET}{MIDREP}_qlusterOut";

	VecStr extractorCmds;
	VecStr qlusterCmds;
	for(const auto & rep : targetsForReps){
		if((bib::in(rep.first, samplesExtracted) ||
				bib::in(bib::replaceString(rep.first, "MID", "") , samplesExtracted) )
				&& !(bib::in(rep.first, noReadsStitched) ||
						bib::in(bib::replaceString(rep.first, "MID", "") , noReadsStitched))){
			std::string fName = rep.first;
			if(!bib::in(rep.first, samplesExtracted)){
				fName = bib::replaceString(rep.first, "MID", "");
			}
			std::string sampName = rep.first;
			if(!bib::beginsWith(rep.first, "MID")){
				sampName = "MID" + rep.first;
			}
			auto tarsNames = bib::conToStr(rep.second, "_");
			auto currentSampTemp = bib::replaceString(sampleNameTemplate, "{REP}", sampName);
			auto cmds = VecStr{extractorCmdTemplate, idTemplate, currentSampTemp};
			if(!noQualTrim){
				cmds.emplace_back(qualTrimTemplate);
			}

			if(bfs::exists(bib::files::make_path(analysisSetup.idsDir_, tarsNames + "_lenCutOffs.tab.txt"))){
				cmds.emplace_back(lenCutOffsTemplate);
			}
			if(ids.containsMids()){
				cmds.emplace_back("--multiplex");
			}
			if(rep.second.size() > 1){
				cmds.emplace_back("--multipleTargets");
			}
			if("" != extraExtractorCmds){
				cmds.emplace_back(extraExtractorCmds);
			}
			auto currentExtractCmd = bib::conToStr(cmds, " ");
			currentExtractCmd = bib::replaceString(currentExtractCmd, "{REP}", fName);
			currentExtractCmd = bib::replaceString(currentExtractCmd, "{TARS}", tarsNames);
			extractorCmds.emplace_back(currentExtractCmd);
			for(const auto & tar : rep.second){
				std::string currentQlusterCmdTemplate = qlusterCmdTemplate;
				if("" != extraQlusterCmds){
					currentQlusterCmdTemplate += " " + extraQlusterCmds;
				}
				currentQlusterCmdTemplate = bib::replaceString(currentQlusterCmdTemplate,"{REP}", fName);
				currentQlusterCmdTemplate = bib::replaceString(currentQlusterCmdTemplate,"{MIDREP}", sampName);
				currentQlusterCmdTemplate = bib::replaceString(currentQlusterCmdTemplate,"{TARGET}", tar);
				qlusterCmds.emplace_back(currentQlusterCmdTemplate);
			}
		}
	}
	OutOptions extractorCmdsOpts(bib::files::make_path(analysisSetup.dir_, "extractorCmds.txt").string());
	std::ofstream extractorCmdsFile;
	openTextFile(extractorCmdsFile, extractorCmdsOpts);
	printVector(extractorCmds, "\n", extractorCmdsFile);

	OutOptions qlusterCmdsOpts(bib::files::make_path(analysisSetup.dir_, "qlusterCmds.txt").string());
	std::ofstream qlusterCmdsFile;
	openTextFile(qlusterCmdsFile, qlusterCmdsOpts);
	printVector(qlusterCmds, "\n", qlusterCmdsFile);

	//process cluster cmds
	std::string processClusterTemplate = setUp.commands_.masterProgram_ + " processClusters "
			"--alnInfoDir alnCache --strictErrors --dout analysis --fastq output.fastq --overWriteDir";
	if("" != extraProcessClusterCmds){
		processClusterTemplate += " " + extraProcessClusterCmds;
	}
	VecStr processClusterCmds;
	if(nullptr != analysisSetup.groupMetaData_){
		processClusterTemplate += " --groupingsFile " + bib::files::make_path(bfs::absolute(analysisSetup.infoDir_), "groupMeta.tab.txt").string();
	}

	if(analysisSetup.samplesForTargets_.size() > 1){
		for(const auto & tar : analysisSetup.samplesForTargets_){
			processClusterCmds.emplace_back("cd popClustering/" + tar.first + " && " + processClusterTemplate);
		}
	}else{
		processClusterCmds.emplace_back("cd popClustering && " + processClusterTemplate);
	}
	OutOptions processClusterCmdsOpts(bib::files::make_path(analysisSetup.dir_, "processClusterCmds.txt").string());
	std::ofstream processClusterCmdsFile;
	openTextFile(processClusterCmdsFile, processClusterCmdsOpts);
	printVector(processClusterCmds, "\n", processClusterCmdsFile);

	//gen analysis configs
	std::string genConfigTemplate = setUp.commands_.masterProgram_ + " genProjectConfig --projectName {TARGET} --out serverConfigs/{TARGET}.config";


	VecStr genConfigCmds;
	if (analysisSetup.samplesForTargets_.size() > 1) {
		for (const auto & tar : analysisSetup.samplesForTargets_) {
			genConfigCmds.emplace_back(
					bib::replaceString(
							genConfigTemplate + " --mainDir popClustering/{TARGET}/analysis",
							"{TARGET}", tar.first));
		}
	} else {
		genConfigCmds.emplace_back(
				bib::replaceString(genConfigTemplate, "{TARGET}",
						analysisSetup.samplesForTargets_.begin()->first)
						+ " --mainDir popClustering/analysis");
	}
	OutOptions genConfigCmdsOpts(bib::files::make_path(analysisSetup.dir_, "genConfigCmds.txt").string());
	std::ofstream genConfigCmdsFile;
	openTextFile(genConfigCmdsFile, genConfigCmdsOpts);
	printVector(genConfigCmds, "\n", genConfigCmdsFile);

	//start server config
	OutOptions startServerCmdOpts(
			bib::files::make_path(analysisSetup.dir_, "startServerCmd.sh").string());
	std::ofstream startServerCmdFile;
	openTextFile(startServerCmdFile, startServerCmdOpts);
	startServerCmdFile << "#!/usr/bin/env bash" << std::endl;
	startServerCmdFile << "# Will automatically run the server in the background and with nohup so it will keep running" << std::endl;

	startServerCmdFile << "if [[ $# -ne 2 ]] && [[ $# -ne 0 ]]; then" << std::endl;
	startServerCmdFile
			<< "	echo \"Illegal number of parameters, needs either 0 or 2 argument, if 2 args 1) port number to server on 2) the name to serve on\""
			<< std::endl;
	startServerCmdFile << "	echo \"Examples\"" << std::endl;
	startServerCmdFile << "	echo \"./startServerCmd.txt\"" << std::endl;
	startServerCmdFile << "	echo \"./startServerCmd.txt 9882 pcv2\"	"
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
	chmod(startServerCmdOpts.outFilename_.c_str(), S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	//start server config
	OutOptions runAnalysisOpts(
			bib::files::make_path(analysisSetup.dir_, "runAnalysis.sh").string());
	std::ofstream runAnalysisFile;
	openTextFile(runAnalysisFile, runAnalysisOpts);
	runAnalysisFile << "#!/usr/bin/env bash" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "##run all parts of the pipeline" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "numThreads=1" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "if [[ $# -eq 1 ]]; then" << std::endl;
	runAnalysisFile << "	numThreads=$1" << std::endl;
	runAnalysisFile << "fi" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "sequenceTools runMultipleCommands --cmdFile extractorCmds.txt      --numThreads $numThreads --raw" << std::endl;
	runAnalysisFile << "sequenceTools runMultipleCommands --cmdFile qlusterCmds.txt        --numThreads $numThreads --raw" << std::endl;
	runAnalysisFile << "sequenceTools runMultipleCommands --cmdFile processClusterCmds.txt --numThreads $numThreads --raw" << std::endl;
	runAnalysisFile << "sequenceTools runMultipleCommands --cmdFile genConfigCmds.txt      --numThreads $numThreads --raw" << std::endl;
	runAnalysisFile << "" << std::endl;
	//make file executable
	chmod(runAnalysisOpts.outFilename_.c_str(), S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	if (foundErrors) {
		std::cerr << bib::bashCT::flashing << bib::bashCT::red << bib::bashCT::bold
				<< "ERRORS FOUND!!" << bib::bashCT::reset << std::endl;
		std::cerr << "Read Warnings in " << bib::bashCT::boldRed(wrningsOpts.outFilename_) << std::endl;
	}

	return 0;
}
//


} // namespace bibseq
