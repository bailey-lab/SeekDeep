/*
 * SeekDeepUtilsRunner_setUpTarAmpAnalysis.cpp
 *
 *  Created on: Jul 1, 2018
 *      Author: nick
 */

#include "SeekDeepUtilsRunner.hpp"

#include <bibseq/ProgramRunners.h>

namespace bibseq {



template<typename CON>
std::string conToStrEndSpecial(const CON & container,
		const std::string & delim, const std::string & lastDelim) {
	std::string ret = "";
	if (container.size() > 1) {
		uint32_t count = 0;
		for (const auto & element : container) {
			if (count + 1 == container.size()) {
				ret += lastDelim;
			} else if (count > 0) {
				ret += delim;
			}
			ret += estd::to_string(element);
			++count;
		}
	} else {
		ret = bib::conToStr(container, delim);
	}
	return ret;
}

int SeekDeepUtilsRunner::setupTarAmpAnalysis(
		const bib::progutils::CmdArgs & inputCommands) {
	VecStr acceptableTechs{"454", "IonTorrent", "Illumina", "Illumina-SingleEnd"};

	TarAmpAnalysisSetup::TarAmpPars pars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.debug = setUp.pars_.debug_;
	setUp.setOption(pars.technology, "--technology",
			"Sequencing Technology (should be " + conToStrEndSpecial(acceptableTechs, ", ", " or ") + ")",
			false, "Technology");
	bib::for_each(acceptableTechs, [](std::string & tech){
		stringToLower(tech);
	});
	stringToLower(pars.technology);


	if (!bib::in(pars.technology, acceptableTechs)) {
		setUp.failed_ = true;
		std::stringstream ss;
		ss
				<< "Error in setting technology, should be "
				<< conToStrEndSpecial(acceptableTechs, ", ", " or ")
				<< " not "
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
		if (analysisSetup.pars_.techIsIllumina()) {
			extractorCmdTemplate = setUp.commands_.masterProgram_
					+ " extractorPairedEnd --dout {INDEX}_extraction --overWriteDir  ";
		} else {
			extractorCmdTemplate = setUp.commands_.masterProgram_
					+ " extractor --dout {INDEX}_extraction --overWriteDir  ";
			if (analysisSetup.pars_.techIsIlluminaSingleEnd()) {
				extractorCmdTemplate += " --illumina ";
			}
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
		if (analysisSetup.pars_.techIsIllumina() || analysisSetup.pars_.techIsIlluminaSingleEnd()) {
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
									+ " extractor --dout {INDEX}_extraction --overWriteDir ";
			 if(analysisSetup.pars_.techIsIlluminaSingleEnd()){
				 extractorCmdTemplate += " --illumina ";
			 }
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
		if (analysisSetup.pars_.techIsIllumina() || analysisSetup.pars_.techIsIlluminaSingleEnd()) {
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


}  // namespace bibseq
