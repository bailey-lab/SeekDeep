
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
				addFunc("setupPairedEndDir", setupPairedEndDir, false),
				addFunc("replaceUnderscores", replaceUnderscores, false)
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
	if(!setUp.pars_.debug_){
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
	if(!setUp.pars_.debug_){
		outFile << allLog << std::endl;
	}

	return 0;
}


int SeekDeepUtilsRunner::replaceUnderscores(const bib::progutils::CmdArgs & inputCommands){

	std::string pat = "";
	bfs::path dir = "./";
	uint32_t upTo = std::numeric_limits<uint32_t>::max();
	std::string replacement = "-";
	bool keepOriginals = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(dir, "--dir", "Directory with the files to rename", true);
	setUp.setOption(pat, "--pat,-p", "Rename files that match this pattern", true);
	setUp.setOption(upTo, "--upTo", "Replace up this many underscores from the begging of the filename");
	setUp.setOption(keepOriginals, "--keepOriginals", "Keep a copy of the original files");
	setUp.setOption(replacement, "--replacement", "What to replace the underscores with");
	setUp.finishSetUp(std::cout);


	auto files = bib::files::listAllFiles(dir.string(), false, {std::regex{pat}});


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
					<< " -> " << bib::bashCT::blue << rk.second << bib::bashCT::reset << std::endl;
		}
		return 0;
	}
	for (const auto & rk : renameKey) {
		if(keepOriginals){
			bfs::copy(rk.first, rk.second);
		}else{
			bfs::rename(rk.first, rk.second);
		}
	}


	return 0;
}








int SeekDeepUtilsRunner::setupPairedEndDir(const bib::progutils::CmdArgs & inputCommands){
	TarAmpAnalysisSetup::TarAmpPars pars;


	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.samplesNamesFnp, "--samples", "A file containing the samples names, columns should go target,sample,pcr1,(optional)pcr2 etc", true);
	setUp.setOption(pars.outDir, "--outDir", "Directory to setup for analysis", true);
	setUp.setOption(pars.inputDir, "--inputDir", "Input Directory of raw data files", true);
	setUp.setOption(pars.idFile, "--idFile", "ID file containing primer and possible additional MIDs", true);
	setUp.setOption(pars.byIndex, "--byIndex", "If the input sample names are by index and not targets");
	setUp.setOption(pars.targetsToIndexFnp, "--targetsToIndex", "A tsv file with two columns named targets and index where targets in a comma separated value with the targets for the index in index",pars.byIndex);


	setUp.setOption(pars.maxOverlap, "--maxOverlap", "Max overlap allowed in stitcher");
	setUp.setOption(pars.numThreads, "--numThreads", "Number of CPUs to use");


	setUp.setOption(pars.extraExtractorCmds, "--extraExtractorCmds", "Extra extractor cmds to add to the defaults");
	setUp.setOption(pars.extraQlusterCmds, "--extraQlusterCmds", "Extra qluster cmds to add to the defaults");
	setUp.setOption(pars.extraProcessClusterCmds, "--extraProcessClusterCmds", "Extra process clusters cmds to add to the defaults");

	setUp.setOption(pars.noQualTrim, "--noQualTrim", "No Quality Trim");

	setUp.setOption(pars.r1Trim, "--r1Trim", "Trim R1 Reads to this length");
	setUp.setOption(pars.r2Trim, "--r2Trim", "Trim R2 Reads to this length");

	setUp.setOption(pars.groupMeta, "--groupMeta", "Group Metadata");

	setUp.setOption(pars.lenCutOffsFnp, "--lenCutOffs", "A file with 3 columns, target,minlen,maxlen to supply length cut off specifically for each target");
	setUp.setOption(pars.refSeqsDir, "--refSeqsDir", "A directory of fasta files where each file is named with the input target names");

	setUp.finishSetUp(std::cout);


	VecStr warnings;

	pars.allChecks(warnings);

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

	TarAmpAnalysisSetup analysisSetup(pars);
	setUp.startARunLog(bib::appendAsNeededRet(analysisSetup.dir_.string(), "/"));
	auto targets = analysisSetup.idsMids_->getTargets();
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
		ss << __PRETTY_FUNCTION__ << ": error, missing the following targets from the id file " << pars.idFile << "\n";
		ss << "Targets: " << bib::conToStr(idMissingTargets, ", ") << "\n";
		ss << "AvailableTargets: " << bib::conToStr(targets, ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}

	if(!sampNamesTargetsTarMissing.empty()){
		foundErrors = true;
		errorOutput << __PRETTY_FUNCTION__ << ": warning, missing the following targets from the sample name file" << pars.samplesNamesFnp << "\n";
		errorOutput << "Targets: " << bib::conToStr(sampNamesTargetsTarMissing, ", ") << "\n";
		errorOutput << "AvailableTargets: " << bib::conToStr(sampNamesTargets, ", ") << "\n";
	}

	if("" != pars.groupMeta){
		if (!analysisSetup.groupMetaData_->missingSamples_.empty()) {
			foundErrors = true;
			errorOutput
					<< "Warning, the following samples were in the meta data file, "
					<< pars.groupMeta << " but were missing from the input sample names file, "
					<< pars.samplesNamesFnp << std::endl;
			for (const auto & samp : analysisSetup.groupMetaData_->missingSamples_) {
				errorOutput << "\tSample: " << samp << "\n";
			}
		}
		if (!analysisSetup.groupMetaData_->missingMetaForSamples_.empty()) {
			foundErrors = true;
			errorOutput
								<< "Warning, the following samples were in the input samples name file, "
								<< pars.samplesNamesFnp << " but were missing from the meta data file, "
								<< pars.groupMeta << std::endl;
			for(const auto & samp : analysisSetup.groupMetaData_->missingMetaForSamples_){
				errorOutput << "\tSample: " << samp << "\n";
			}
		}
	}

	if("" != pars.refSeqsDir){
		if(!analysisSetup.forRefSeqs_.missing_.empty()){
			foundErrors = true;
			errorOutput
								<< "Warning, reference sequences were supplied but reference files for the following targets couldn't be found" << "\n";
			for(const auto & tar : analysisSetup.forRefSeqs_.missing_){
				errorOutput << "\tTarget: " << tar << "\n";
			}
		}
		if(!analysisSetup.forRefSeqs_.notMatching_.empty()){
			foundErrors = true;
			errorOutput
								<< "Warning, reference sequences were supplied but there are no matching reference targets in project for the following" << "\n";
			for(const auto & tar : analysisSetup.forRefSeqs_.notMatching_){
				errorOutput << "\tTarget: " << tar << "\n";
			}
			errorOutput << "\tOptions are: " << bib::conToStr(analysisSetup.getTargets(), ", ") << "\n";
		}
	}

	if("" != pars.lenCutOffsFnp){
		if(!analysisSetup.forLenCutOffs_.missing_.empty()){
			foundErrors = true;
			errorOutput
								<< "Warning, length cut offs were supplied but cut offs for the following targets couldn't be found" << "\n";
			for(const auto & tar : analysisSetup.forLenCutOffs_.missing_){
				errorOutput << "\tTarget: " << tar << "\n";
			}
		}
		if(!analysisSetup.forLenCutOffs_.notMatching_.empty()){
			foundErrors = true;
			errorOutput
								<< "Warning, length cut offs were supplied but cut offs for targets didn't match targets in project for the following targets" << "\n";
			for(const auto & tar : analysisSetup.forRefSeqs_.notMatching_){
				errorOutput << "\tTarget: " << tar << "\n";
			}
			errorOutput << "\tOptions are: " << bib::conToStr(analysisSetup.getTargets(), ", ") << "\n";
		}
	}
	//now write id files
	analysisSetup.writeOutIdFiles();



	auto files = bib::files::listAllFiles(pars.inputDir.string(), false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;

	auto reps = analysisSetup.getExpectantInputNames();

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
	const std::string zcatTempCmdR1 = analysisSetup.pars_.zcatCmd + " " + "{FILES} > "
			+  "{OUTPUT}_R1.fastq";
	const std::string zcatTempCmdR2 = analysisSetup.pars_.zcatCmd + " " + "{FILES} > "
			+  "{OUTPUT}_R2.fastq";
	const std::string stitchCmd = analysisSetup.pars_.stitcherCmd + " " +
			"{OUTPUT}_R1.fastq " +
			"{OUTPUT}_R2.fastq " +
			"-o {OUTPUT} " +
			"--max-overlap " + estd::to_string(analysisSetup.pars_.maxOverlap);

	auto extractStitchFiles =
			[&analysisSetup ,&samplesExtracted, &samplesEmpty,
			 &zcatTempCmdR1, &zcatTempCmdR2, &stitchCmd,
			 &filesKeys, &readsByPairs,
			 &logs, &logsMut,
			 &pars]() {
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
					if(std::numeric_limits<uint32_t>::max() != pars.r1Trim){
						auto finalFnp = bib::files::make_path(analysisSetup.dir_, key).string() + "_R1.fastq";
						auto tempFnp = bib::files::make_path(analysisSetup.dir_, "tempTrim_" + key).string() + "_R1.fastq";
						auto opts = SeqIOOptions::genFastqInOut(finalFnp,tempFnp,false);
						SeqIO reader(opts);
						reader.openIn();
						reader.openOut();
						seqInfo r1seq;
						while(reader.readNextRead(r1seq)){
							readVecTrimmer::trimToMaxLength(r1seq, pars.r1Trim);
							reader.write(r1seq);
						}
						reader.closeIn();
						reader.closeOut();
						bfs::remove(finalFnp);
						bfs::rename(tempFnp, finalFnp);
					}

					if(std::numeric_limits<uint32_t>::max() != pars.r2Trim){
						auto finalFnp = bib::files::make_path(analysisSetup.dir_, key).string() + "_R2.fastq";
						auto tempFnp =  bib::files::make_path(analysisSetup.dir_, "tempTrim_" + key).string() + "_R2.fastq";
						auto opts = SeqIOOptions::genFastqInOut(finalFnp,tempFnp,false);
						SeqIO reader(opts);
						reader.openIn();
						reader.openOut();
						seqInfo r2seq;
						while(reader.readNextRead(r2seq)){
							readVecTrimmer::trimToMaxLength(r2seq, pars.r2Trim);
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
	for(uint32_t tNum = 0; tNum < analysisSetup.pars_.numThreads; ++tNum){
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
		errorOutput << "The following files were found but didn't match any input sample names in " << analysisSetup.pars_.samplesNamesFnp << std::endl;
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
	analysisSetup.setUpPopClusteringDirs(setUp.pars_.verbose_);



	VecStr extractorCmds;
	VecStr qlusterCmds;
	if (analysisSetup.pars_.byIndex) {
		//extractor cmds
		std::string extractorCmdTemplate = setUp.commands_.masterProgram_ + " extractor --fastq {INDEX}.extendedFrags.fastq --dout {INDEX}_extraction --overWriteDir --illumina";
		std::string qualTrimTemplate = "--trimAtQual 2";
		std::string lenCutOffsTemplate = "--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt";
		std::string idTemplate = "--id info/ids/{TARS}.id.txt --multiplex ";
		//qluster cmds;
		auto extractionDirs = bib::files::make_path(bfs::absolute(analysisSetup.dir_), "{INDEX}_extraction");
		std::string qlusterCmdTemplate = "cd \"" + extractionDirs.string() + "\" && "
				+ "if [ -f {TARGET}{MIDREP}.fastq  ]; then "
				+ setUp.commands_.masterProgram_ + " qluster "
				"--fastq {TARGET}{MIDREP}.fastq  --illumina --qualThres 25,20 "
				"--alnInfoDir {TARGET}{MIDREP}_alnCache --overWriteDir "
				"--additionalOut \"../popClustering/{TARGET}/locationByIndex/{INDEX}.tab.txt\" "
				"--overWrite --dout {TARGET}{MIDREP}_qlusterOut ";
		auto indexes = analysisSetup.getIndexes();
		for(const auto & index : indexes){
			if(bib::in(index, samplesExtracted) &&
					!bib::in(index, noReadsStitched)){
				auto tarsNames = bib::conToStr(analysisSetup.indexToTars_[index], "_");
				auto cmds = VecStr { extractorCmdTemplate, idTemplate };
				if (!analysisSetup.pars_.noQualTrim) {
					cmds.emplace_back(qualTrimTemplate);
				}
				if (bfs::exists(
						bib::files::make_path(analysisSetup.idsDir_,
								tarsNames + "_lenCutOffs.tab.txt"))) {
					cmds.emplace_back(lenCutOffsTemplate);
				}
				if (analysisSetup.indexToTars_[index].size() > 1) {
					cmds.emplace_back("--multipleTargets");
				}
				if ("" != analysisSetup.pars_.extraExtractorCmds) {
					cmds.emplace_back(analysisSetup.pars_.extraExtractorCmds);
				}
				auto currentExtractCmd = bib::conToStr(cmds, " ");
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{INDEX}",
						index);
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{TARS}",
						tarsNames);
				extractorCmds.emplace_back(currentExtractCmd);
				for(const auto & mid : analysisSetup.idsMids_->getMids()){
					for (const auto & tar : analysisSetup.indexToTars_[index]) {
						auto currentQlusterCmdTemplate = qlusterCmdTemplate;
						if ("" != analysisSetup.pars_.extraQlusterCmds) {
							currentQlusterCmdTemplate  += " "
									+ analysisSetup.pars_.extraQlusterCmds;
						}
						currentQlusterCmdTemplate += "; fi";
						currentQlusterCmdTemplate = bib::replaceString(currentQlusterCmdTemplate,
								"{INDEX}", index);
						currentQlusterCmdTemplate = bib::replaceString(currentQlusterCmdTemplate,
								"{TARGET}", tar);
						currentQlusterCmdTemplate = bib::replaceString(currentQlusterCmdTemplate,
								"{MIDREP}", mid);
						qlusterCmds.emplace_back(currentQlusterCmdTemplate);
					}
				}
			}
		}
	} else {
		//extractor cmds
		std::string extractorCmdTemplate = setUp.commands_.masterProgram_ + " extractor --fastq {REP}.extendedFrags.fastq --dout {REP}_extraction --overWriteDir --illumina";
		std::string qualTrimTemplate = "--trimAtQual 2";
		std::string lenCutOffsTemplate = "--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt";
		std::string idTemplate = "--id info/ids/{TARS}.id.txt";
		std::string sampleNameTemplate = "--sampleName {REP}";
		//qluster cmds;
		auto extractionDirs = bib::files::make_path(bfs::absolute(analysisSetup.dir_), "{REP}_extraction");

		std::string qlusterCmdTemplate = "cd " + extractionDirs.string() + " && " + setUp.commands_.masterProgram_ + " qluster "
				"--fastq {TARGET}{MIDREP}.fastq  --illumina --qualThres 25,20 "
				"--alnInfoDir {TARGET}{MIDREP}_alnCache --overWriteDir "
				"--additionalOut ../popClustering/locationByIndex/{TARGET}.tab.txt "
				"--overWrite --dout {TARGET}{MIDREP}_qlusterOut";
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
			if ((bib::in(rep.first, samplesExtracted)
					|| bib::in(bib::replaceString(rep.first, "MID", ""), samplesExtracted))
					&& !(bib::in(rep.first, noReadsStitched)
							|| bib::in(bib::replaceString(rep.first, "MID", ""),
									noReadsStitched))) {
				std::string fName = rep.first;
				if (!bib::in(rep.first, samplesExtracted)) {
					fName = bib::replaceString(rep.first, "MID", "");
				}

				std::string sampName = rep.first;
				if (!bib::beginsWith(rep.first, "MID")) {
					sampName = "MID" + rep.first;
				}

				auto tarsNames = bib::conToStr(rep.second, "_");
				auto currentSampTemp = bib::replaceString(sampleNameTemplate, "{REP}",
						sampName);
				auto cmds = VecStr { extractorCmdTemplate, idTemplate, currentSampTemp };
				if (!analysisSetup.pars_.noQualTrim) {
					cmds.emplace_back(qualTrimTemplate);
				}
				if (bfs::exists(
						bib::files::make_path(analysisSetup.idsDir_,
								tarsNames + "_lenCutOffs.tab.txt"))) {
					cmds.emplace_back(lenCutOffsTemplate);
				}
				if (rep.second.size() > 1) {
					cmds.emplace_back("--multipleTargets");
				}
				if ("" != analysisSetup.pars_.extraExtractorCmds) {
					cmds.emplace_back(analysisSetup.pars_.extraExtractorCmds);
				}
				auto currentExtractCmd = bib::conToStr(cmds, " ");
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{REP}",
						fName);
				currentExtractCmd = bib::replaceString(currentExtractCmd, "{TARS}",
						tarsNames);
				extractorCmds.emplace_back(currentExtractCmd);
				for (const auto & tar : rep.second) {
					std::string currentQlusterCmdTemplate = qlusterCmdTemplate;
					if ("" != analysisSetup.pars_.extraQlusterCmds) {
						currentQlusterCmdTemplate += " "
								+ analysisSetup.pars_.extraQlusterCmds;
					}
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
	if("" != analysisSetup.pars_.extraProcessClusterCmds){
		processClusterTemplate += " " + analysisSetup.pars_.extraProcessClusterCmds;
	}
	VecStr processClusterCmds;
	if(nullptr != analysisSetup.groupMetaData_){
		processClusterTemplate += " --groupingsFile " + bib::files::make_path(bfs::absolute(analysisSetup.infoDir_), "groupMeta.tab.txt").string();
	}
	auto popDir = bib::files::make_path(bfs::absolute(analysisSetup.dir_),
			"popClustering");

	for (const auto & tar : targets) {
		processClusterCmds.emplace_back(
				"cd " + bib::files::make_path(popDir, tar).string() + " && "
						+ processClusterTemplate);
	}
	OutOptions processClusterCmdsOpts(bib::files::make_path(analysisSetup.dir_, "processClusterCmds.txt").string());
	std::ofstream processClusterCmdsFile;
	openTextFile(processClusterCmdsFile, processClusterCmdsOpts);
	printVector(processClusterCmds, "\n", processClusterCmdsFile);

	//gen analysis configs
	std::string genConfigTemplate = setUp.commands_.masterProgram_ + " genProjectConfig --projectName {TARGET} --out serverConfigs/{TARGET}.config";


	VecStr genConfigCmds;
	for (const auto & tar : targets) {
		genConfigCmds.emplace_back(
				bib::replaceString(
						genConfigTemplate + " --mainDir popClustering/{TARGET}/analysis",
						"{TARGET}", tar));
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
			<< "	echo \"Illegal number of parameters, needs either 0 or 2 arguments, if 2 args 1) port number to server on 2) the name to serve on\""
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
	runAnalysisFile << "##except start the server which is done by running ./startServerCmd.sh" << std::endl;
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
