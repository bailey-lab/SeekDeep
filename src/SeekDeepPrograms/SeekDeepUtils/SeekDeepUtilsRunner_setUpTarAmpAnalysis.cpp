/*
 * SeekDeepUtilsRunner_setUpTarAmpAnalysis.cpp
 *
 *  Created on: Jul 1, 2018
 *      Author: nick
 */

#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/ProgramRunners.h>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {



int SeekDeepUtilsRunner::setupTarAmpAnalysis(
		const njh::progutils::CmdArgs & inputCommands) {
	VecStr acceptableTechs{"454", "IonTorrent", "Illumina", "Illumina-SingleEnd", "Nanopore", "Pacbio"};
	TarAmpAnalysisSetup::TarAmpPars pars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.debug = setUp.pars_.debug_;
	setUp.setOption(pars.technology, "--technology",
			"Sequencing Technology (should be " + njh::conToStrEndSpecial(acceptableTechs, ", ", " or ") + ")",
			false, "Technology");
	njh::for_each(acceptableTechs, [](std::string & tech){
		stringToLower(tech);
	});
	stringToLower(pars.technology);



	if (!njh::in(pars.technology, acceptableTechs)) {
		setUp.failed_ = true;
		std::stringstream ss;
		ss
				<< "Error in setting technology, should be "
				<< njh::conToStrEndSpecial(acceptableTechs, ", ", " or ")
				<< " not "
				<< pars.technology << "\n";
		setUp.addWarning(ss.str());
	}
	setUp.setOption(pars.development,   "--development",   "development mode, will generate a lot more output", false, "Output");
	setUp.setOption(pars.outDir,   "--outDir",   "Directory to setup for analysis", true, "Output");
	setUp.setOption(pars.inputDir, "--inputDir", "Input Directory of raw data files", true, "Input");
	setUp.setOption(pars.idFile,   "--idFile",   "ID file containing primer and possible additional MIDs", true, "ID Files");

	setUp.setOption(pars.samplesNamesFnp, "--samples", "A file containing the samples names, columns should go target,sample,pcr1,pcr2(optional)",
			false, "ID Files");

  setUp.setOption(pars.uniqueKmersPerTarget, "--uniqueKmersPerTarget", "unique Kmers Per Target", pars.techIsNanoporeOrPacbio());
	std::string overlapStatus{"auto"};
	std::set<std::string> allowableOverlapStatuses{"AUTO", "R1BEGINSINR2", "R1ENDSINR2", "NOOVERLAP", "ALL"};

	std::function<njh::progutils::ProgramSetUp::FlagCheckResult(const std::string&)> overlapStatusCheck = [allowableOverlapStatuses](const std::string & flagSet){
		bool success = true;
		std::string mess = "";
		std::string upper = njh::strToUpperRet(flagSet);
		if(!njh::in(upper, allowableOverlapStatuses)){
			success = false;
			mess = njh::pasteAsStr("--defaultOverlapStatus needs to be one of the following (case insensitive) ", njh::conToStr(allowableOverlapStatuses, ","), " not ", upper);
		}
		return njh::progutils::ProgramSetUp::FlagCheckResult{success, mess};
	};

	bool setDefaultOverlapStatus = setUp.setOption(overlapStatus,
			"--defaultOverlapStatus",
			"Set a overlap status for all targets, can be 1 of 5 values(case insensitive), AUTO, R1BEGINSINR2, R1ENDSINR2, NOOVERLAP, ALL. ALL=(R1BEGINSINR2 and R1ENDSINR2). Setting to auto will go with the status was most commonly found for a target, this can be dangerous as with unspecific amplification this can end up being set as the incorrect status",
			false, overlapStatusCheck);
	if(setDefaultOverlapStatus){
		std::string upper = njh::strToUpperRet(overlapStatus);
		if("ALL" == upper){
			pars.defaultStatuses_.emplace_back(
					PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2);
			pars.defaultStatuses_.emplace_back(
					PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2);
			pars.defaultStatuses_.emplace_back(
					PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP);
		} else if("AUTO" == upper){
			pars.defaultStatuses_.emplace_back(PairedReadProcessor::ReadPairOverLapStatus::AUTO);
		} else if("R1BEGINSINR2" == upper){
			pars.defaultStatuses_.emplace_back(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2);
		}else if("R1ENDSINR2" == upper){
			pars.defaultStatuses_.emplace_back(PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2);
		}else if("NOOVERLAP" == upper){
			pars.defaultStatuses_.emplace_back(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP);
		}
	}

	setUp.setOption(pars.overlapStatusFnp, "--overlapStatusFnp",
			"A file with two columns, target,status; status column should contain 1 of 3 values (capitalization doesn't matter): r1BegOverR2End,r1EndOverR2Beg,NoOverlap. r1BegOverR2End=target size < read length (causes read through),r1EndOverR2Beg= target size > read length less than 2 x read length, NoOverlap=target size > 2 x read length",
			pars.techIsIllumina() && !setDefaultOverlapStatus, "Illumina Stitching");

	setUp.setOption(pars.targetsToIndexFnp, "--targetsToIndex",
			"A tsv file with two columns named targets and index where targets in a comma separated value with the targets for the index in index",
			false, "ID Files");

	setUp.setOption(pars.numThreads, "--numThreads", "Number of CPUs to use");

	// Extra arguments to give the sub programs
	setUp.setOption(pars.extraExtractorCmds, "--extraExtractorCmds",
			"Extra extractor cmds to add to the defaults", false, "Extra Commands");
	setUp.setOption(pars.extraQlusterCmds, "--extraQlusterCmds,--extraKlusterCmds",
			"Extra qluster/kluster cmds to add to the defaults", false, "Extra Commands");
	if(setUp.setOption(pars.extraQlusterCmdsPerTargetFnp, "--extraQlusterCmdsPerTarget,--extraKlusterCmdsPerTarget",
			"Extra qluster/kluster cmds to add to the defaults per target, table needs cols, 1)target,2)extraQlusterCmds", false, "Extra Commands")){
		if(!bfs::exists(pars.extraQlusterCmdsPerTargetFnp)){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr(pars.extraQlusterCmdsPerTargetFnp, " needs to exist for --extraQlusterCmdsPerTarget,--extraKlusterCmdsPerTarget"));
		}
		table extraQlusterCmdsPerTargetTab(pars.extraQlusterCmdsPerTargetFnp, "\t", true);
		extraQlusterCmdsPerTargetTab.changeHeaderToLowerCase();
		VecStr requiredColumns{"target", "extraqlustercmds"};
		VecStr missing = extraQlusterCmdsPerTargetTab.getMissingHeaders(requiredColumns);
		if(!missing.empty()){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr("table ", pars.extraQlusterCmdsPerTargetFnp, " is missing columns ", njh::conToStr(requiredColumns, ",")));
		}
		for(const auto & row : extraQlusterCmdsPerTargetTab){
			pars.extraQlusterCmdsPerTarget[row[extraQlusterCmdsPerTargetTab.getColPos("target")]] = row[extraQlusterCmdsPerTargetTab.getColPos("extraqlustercmds")];
		}
	}

	setUp.setOption(pars.extraProcessClusterCmds, "--extraProcessClusterCmds",
			"Extra process clusters cmds to add to the defaults", false, "Extra Commands");
	if(setUp.setOption(pars.extraProcessClusterCmdsPerTargetFnp, "--extraProcessClusterCmdsPerTarget",
			"Extra process clusters cmds to add to the defaults, table needs cols, 1)target,2)extraProcessClusterCmds", false, "Extra Commands")){
		if(!bfs::exists(pars.extraProcessClusterCmdsPerTargetFnp)){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr(pars.extraProcessClusterCmdsPerTargetFnp, " needs to exist for --extraProcessClusterCmdsPerTarget"));
		}
		table extraProcessClusterCmdsPerTargetTab(pars.extraProcessClusterCmdsPerTargetFnp, "\t", true);
		extraProcessClusterCmdsPerTargetTab.changeHeaderToLowerCase();
		VecStr requiredColumns{"target", "extraprocessclustercmds"};
		VecStr missing = extraProcessClusterCmdsPerTargetTab.getMissingHeaders(requiredColumns);
		if(!missing.empty()){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr("table ", pars.extraProcessClusterCmdsPerTargetFnp, " is missing columns ", njh::conToStr(requiredColumns, ",")));
		}
		for(const auto & row : extraProcessClusterCmdsPerTargetTab){
			pars.extraProcessClusterCmdsPerTarget[row[extraProcessClusterCmdsPerTargetTab.getColPos("target")]] = row[extraProcessClusterCmdsPerTargetTab.getColPos("extraprocessclustercmds")];
		}
	}




	setUp.setOption(pars.useKlusterClustering_, "--useKCrushClustering",
                  "Use kmer clustering for high error rate sequences like nanopore and pacbio");


	setUp.setOption(pars.conservative, "--conservativePopClus",
			"Do conservative population clustering which skips possible artifact cleanup step", false, "processClusters");
	setUp.setOption(pars.rescueFilteredHaplotypes, "--rescueFilteredHaplotypes", "Add on resuce of haplotypes that filtered due to low frequency or chimera filtering if it appears as a major haplotype in another sample", false, "processClusters");
	setUp.setOption(pars.groupMeta, "--groupMeta", "Group Metadata", false, "Meta");
	setUp.setOption(pars.lenCutOffsFnp, "--lenCutOffs",
			"A file with 3 columns, 1)target, 2)minlen, 3)maxlen to supply length cut off specifically for each target", pars.techIsNanoporeOrPacbio(), "Extractor");
	setUp.setOption(pars.refSeqsDir, "--refSeqsDir",
			"A directory of fasta files where each file is named with the input target names", false, "Extractor");
  setUp.setOption(pars.previousPopSeqsDir, "--previousPopSeqsDir",
                  "A directory of fasta files where each file is named with the input target names to rename population names", false, "ProcessClusters");


	if (pars.techIs454() || pars.techIsIonTorrent()) {
		pars.inputFilePat = ".*.fastq";
	}

	setUp.setOption(pars.inputFilePat, "--inputFilePat",
			"The input file pattern in the input directory to work on", false, "Input");

	setUp.setOption(pars.noGuessSampNames, "--noGuessSampNames",
			"Don't guess the sample names from the raw fastq directory input", false, "Input");
	if("" == pars.samplesNamesFnp){
		if(!setUp.setOption(pars.samplesNamesWithBarcodeInfoFnp, "--samplesNamesWithBarcodeInfo",
				"Sample file 3 or 4 required columns 1)library,2)sample,3)fbarcode,4(optional))rbarcode."
			"\n\t\t\t1) name of input file without extension/illumina info e.g. Sample1 for Sample1_S2_R1_001.fastq.gz"
			"\n\t\t\t2) sample name for this barcode in this library file, sample names can be duplicate to indicate replicates"
			"\n\t\t\t3) barcode sequence associated with forward primer"
			"\n\t\t\t4) if sample is dual barcoded, barcode associated with reverse primer", false, "Input")){
			setUp.setOption(pars.ignoreSamples, "--ignoreSamplesWhenGuessing",
					"Ignore Samples with these names", false, "Input");
			setUp.setOption(pars.replicatePattern, "--replicatePatternWhenGuessing",
							"Replicate name regex pattern when guessing samples to order samples into replicates, should be two regex group, the first being sample, the 2nd being the replicate, e.g. --replicatePatternWhenGuessing\"(.*)(-rep.*)\"", false, "Input");
		}
		setUp.setOption(pars.samplesNamesByLibraryNameFnp, "--samplesNamesBySeqFileName",
				"Sample file 2 required columns 1)library,2)sample."
			"\n\t\t\t1) name of input file without extension/illumina info e.g. Sample1 for Sample1_S2_R1_001.fastq.gz, this cannot be duplicated as only one sample can be assoicated with an input file when there are no barcodes"
			"\n\t\t\t2) sample name for library file, sample names can be duplicate to indicate replicates ", false, "Input");
	}

	setUp.setOption(pars.doNotGuessRecFlags, "--doNotGuessRecFlags", "Don't guess at additional SeekDeep extractor/extratorPairedEnd flags by investigating input sequence files", false, "Extra Commands");
	setUp.setOption(pars.numberOfFilesToInvestigate, "--numberOfFilesToInvestigate", "Number of files to investigate when adding additional recommended flags", false, "Extra Commands");
	setUp.setOption(pars.testNumberOfReadsToInvestigate, "--testNumberOfReadsToInvestigate", "Number of reads per file to investigate when adding additional recommended flags", false, "Extra Commands");

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
	setUp.startARunLog(njh::appendAsNeededRet(analysisSetup.dir_.string(), "/"));
	auto targets = analysisSetup.idsMids_->getTargets();
	auto sampNamesTargets = analysisSetup.getTargets();
	njh::sort(targets);
	njh::sort(sampNamesTargets);
	VecStr idMissingTargets;
	VecStr sampNamesTargetsTarMissing;
	for (const auto & tar : sampNamesTargets) {
		if (!njh::in(tar, targets)) {
			idMissingTargets.emplace_back(tar);
		}
	}
	for (const auto & tar : targets) {
		if (!njh::in(tar, sampNamesTargets)) {
			sampNamesTargetsTarMissing.emplace_back(tar);
		}
	}

	if (!idMissingTargets.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": error, missing the following targets from the id file "
				<< pars.idFile << "\n";
		ss << "Targets: " << njh::conToStr(idMissingTargets, ", ") << "\n";
		ss << "AvailableTargets: " << njh::conToStr(targets, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}

	if (!sampNamesTargetsTarMissing.empty()) {
		foundErrors = true;
		errorOutput << __PRETTY_FUNCTION__
				<< ": warning, missing the following targets from the sample name file"
				<< pars.samplesNamesFnp << "\n";
		errorOutput << "Targets: "
				<< njh::conToStr(sampNamesTargetsTarMissing, ", ") << "\n";
		errorOutput << "AvailableTargets: " << njh::conToStr(sampNamesTargets, ", ")
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

	if ("" != pars.refSeqsDir || !pars.previousPopSeqsDir.empty()) {
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
					<< njh::conToStr(analysisSetup.getTargets(), ", ") << "\n";
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
					<< njh::conToStr(analysisSetup.getTargets(), ", ") << "\n";
		}
	}

	//now write id files
	analysisSetup.writeOutIdFiles();

	auto files = njh::files::listAllFiles(pars.inputDir.string(), false, {
			std::regex ( analysisSetup.pars_.inputFilePat ) });

	if (setUp.pars_.debug_) {
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}
	auto expectedSamples = analysisSetup.getExpectantInputNames();
	if(setUp.pars_.debug_){
		std::cout << "Expected samples: " << std::endl;
		std::cout << njh::conToStr(expectedSamples, "\n") << std::endl;
	}
	VecStr unrecognizedInput;
	VecStr sampleFilesFound;
	VecStr sampleFilesNotFound;
	//VecStr samplesExtracted;
	//VecStr samplesEmpty;
	Json::Value logs;
	std::mutex logsMut;
	std::map<std::string, std::pair<VecStr, VecStr>> readsByPairs ;
	std::map<std::string, bfs::path> filesByPossibleName;

	ReadPairsOrganizer rpOrganizer(expectedSamples);
	rpOrganizer.doNotGuessSampleNames_ = pars.noGuessSampNames;

	if (analysisSetup.pars_.techIsIllumina()) {
		rpOrganizer.processFiles(files);
		readsByPairs = rpOrganizer.processReadPairs();
		auto keys = getVectorOfMapKeys(readsByPairs);
		njh::sort(keys);
		sampleFilesFound = keys;
		if(setUp.pars_.debug_){
			std::cout << "Samples found: " << std::endl;
			std::cout << njh::conToStr(keys, "\n") << std::endl;
		}
		for(const auto & expected : expectedSamples){
			if(!njh::in(expected, sampleFilesFound)){
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
		// need to find samples that are empty
		// add to inputPassed
		// just checking for possible compression with file ending, might consider changing to libmagic or something to make sure
		bool compressed = false;
		if (njh::endsWith(analysisSetup.pars_.inputFilePat, ".gz")) {
			compressed = true;
		}
		for (const auto & file : files) {
			auto fNameNoExt = file.first.filename().replace_extension("");
			if (compressed) {
				fNameNoExt.replace_extension("");
			}
			if (njh::in(fNameNoExt.string(), expectedSamples)) {
				filesByPossibleName[fNameNoExt.string()] = file.first;
			} else {
				unrecognizedInput.emplace_back(file.first.string());
			}
		}
		auto keys = njh::getVecOfMapKeys(filesByPossibleName);
		njh::sort(keys);
		sampleFilesFound = keys;
		VecStr samplesNotFound;
		//std::cout <<"expectedSamples.size(): " << expectedSamples.size() << std::endl;
		for (const auto & expected : expectedSamples) {
			if (!njh::in(expected, keys)
					&& !njh::in(njh::replaceString(expected, "MID", ""),
							keys)) {
				samplesNotFound.emplace_back(expected);
			}
		}

		if (!samplesNotFound.empty()) {
			foundErrors = true;
			errorOutput << "The following input files were not found" << std::endl;
			for (const auto & samp : samplesNotFound) {
				errorOutput << "\tSample: " << samp << std::endl;
			}
		}
	}



	OutOptions wrningsOpts(njh::files::make_path(analysisSetup.dir_, "WARNINGS_PLEASE_READ.txt"));
	if (foundErrors) {
		std::ofstream outWarnings;
		openTextFile(outWarnings, wrningsOpts);
		outWarnings << errorOutput.str() << std::endl;
	}

	//make population clustering directory
	analysisSetup.setUpPopClusteringDirs(setUp.pars_.verbose_);

	VecStr extractorCmds;
	VecStr qlusterCmds;


	//investigate input seq files to recommend
	if(!pars.doNotGuessRecFlags){
		std::vector<SeqIOOptions> filesToInvestigate;
		/**@todo consider doing random selection here*/
		if(pars.techIsIllumina()){
			for(const auto & pair : readsByPairs){
				if(njh::endsWith(pair.second.first.front(), ".gz")){
					filesToInvestigate.emplace_back(SeqIOOptions::genPairedInGz(bfs::path(pair.second.first.front()), bfs::path(pair.second.second.front())));
				}else{
					filesToInvestigate.emplace_back(SeqIOOptions::genPairedIn(bfs::path(pair.second.first.front()), bfs::path(pair.second.second.front())));
				}
				if(filesToInvestigate.size() > pars.numberOfFilesToInvestigate){
					break;
				}
			}
		}else{
			for(const auto & file : filesByPossibleName){
				filesToInvestigate.emplace_back(SeqIOOptions(file.second, SeqIOOptions::getInFormatFromFnp(file.second), false));
				if(filesToInvestigate.size() > pars.numberOfFilesToInvestigate){
					break;
				}
			}
		}
		TarAmpSeqInvestigator::TarAmpSeqInvestigatorPars investPars;
		investPars.idFnp = analysisSetup.pars_.idFile;
		investPars.testNumber = analysisSetup.pars_.testNumberOfReadsToInvestigate;
		investPars.pars.corePars_.pDetPars.primerWithin_ = 40;

		TarAmpSeqInvestigator masterInvestigator(investPars);
		std::mutex masterInvesMut;
		if(setUp.pars_.verbose_){
			std::cout << "Investigating files for additional recommended extractor flags which will be automatically added" << std::endl;
		}

		njh::concurrent::LockableQueue<SeqIOOptions> optsQueue(filesToInvestigate);
		bool verbose = setUp.pars_.verbose_;
		std::function<void()> investigateFile = [&optsQueue,&investPars,&masterInvesMut,&masterInvestigator,&verbose](){
			SeqIOOptions seqOpts;
			TarAmpSeqInvestigator investigator(investPars);

			while(optsQueue.getVal(seqOpts)){
				if(verbose){
					std::cout << "Investigating " << seqOpts.firstName_ << " " << ("" ==seqOpts.secondName_ ? std::string("") : seqOpts.secondName_.string()) << std::endl;
				}
				investigator.investigateFile(seqOpts, verbose);
			}
			{
				std::lock_guard<std::mutex> lock(masterInvesMut);
				masterInvestigator.addOtherCounts(investigator);
			}
		};
		njh::concurrent::runVoidFunctionThreaded(investigateFile, analysisSetup.pars_.numThreads);

		masterInvestigator.processCounts();


		auto recFlags = masterInvestigator.recommendSeekDeepExtractorFlags();

		VecStr flagsToAdd;
		auto currentExtraExtractorCmds = njh::strToLowerRet(analysisSetup.pars_.extraExtractorCmds);

		for(const auto & recFlag : recFlags){
			auto rflag = njh::strToLowerRet(recFlag);
			trimAtFirstWhitespace(rflag);
			njh::lstrip(rflag, '-');
			if(std::string::npos == currentExtraExtractorCmds.find(rflag)){
				flagsToAdd.emplace_back(recFlag);
			} else {
				if(setUp.pars_.verbose_){
					std::cout << "Already have " << recFlag << " no need to add" << std::endl;
				}
			}
		}
		if(!flagsToAdd.empty()){
			auto addingStr = njh::conToStr(flagsToAdd, " ");
			if(setUp.pars_.verbose_){
				std::cout << "Adding " << addingStr << std::endl;
			}
			analysisSetup.pars_.extraExtractorCmds.append(" ");
			analysisSetup.pars_.extraExtractorCmds.append(addingStr);
		}
		auto investDir = njh::files::make_path(analysisSetup.reportsDir_, "fileInvestigationCounts");
		njh::files::makeDir(njh::files::MkdirPar{investDir});
		masterInvestigator.writeOutTables(investDir, true);
	}

	if (analysisSetup.pars_.byIndex) {
		if(setUp.pars_.debug_){
			std::cout << "Samples:" << std::endl;
			std::cout << njh::conToStr(getVectorOfMapKeys(analysisSetup.samples_), "\n") << std::endl;
		}


		//extractor cmds
		std::string extractorCmdTemplate;
		if (analysisSetup.pars_.techIsIllumina()) {
			extractorCmdTemplate = setUp.commands_.masterProgram_
					+ " extractorPairedEnd --dout {INDEX}_extraction --overWriteDir  "  ;
		} else if (analysisSetup.pars_.techIsNanoporeOrPacbio()) {
      extractorCmdTemplate = setUp.commands_.masterProgram_
                             + " extractorByKmerMatching --dout {INDEX}_extraction --overWriteDir  "
                               " --uniqueKmersPerTarget info/ids/{TARS}_uniqueKmers.tab.txt.gz "
                               " --rename "
                               ;
    } else {
			extractorCmdTemplate = setUp.commands_.masterProgram_
					+ " extractor --dout {INDEX}_extraction --overWriteDir  ";
			if (analysisSetup.pars_.techIsIlluminaSingleEnd()) {
				extractorCmdTemplate += " --illumina ";
			}
		}
		if (analysisSetup.pars_.techIsIllumina()) {
			extractorCmdTemplate += "--fastq1gz \"{INDEX_R1}\" --fastq2gz \"{INDEX_R2}\" ";
		} else {
			extractorCmdTemplate += "--fastqgz \"{INDEX_InFile}\" ";
		}
		std::string lenCutOffsTemplate = "--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt ";
		std::string idTemplate = "--id info/ids/{TARS}.id.txt ";
		std::string overLapStatusTemplate = "--overlapStatusFnp info/ids/{TARS}_overlapStatus.tab.txt ";
		std::string refSeqsDir = "--compareSeq info/refs/ ";
		//qluster cmds;
		auto extractionDirs = njh::files::make_path(
				bfs::absolute(analysisSetup.dir_), "{INDEX}_extraction");
		std::string qlusterCmdTemplate =
				"cd \"" + extractionDirs.string() + "\" && "
						+ "if [ -f {TARGET}{MIDREP}.fastq.gz  ]; then "
						+ setUp.commands_.masterProgram_
						+ " qluster "
								"--fastqgz \"{TARGET}{MIDREP}.fastq.gz\" ";
		if(pars.development){
			qlusterCmdTemplate = qlusterCmdTemplate +  "--alnInfoDir {TARGET}{MIDREP}_alnCache --development";
		}
		qlusterCmdTemplate = qlusterCmdTemplate +
								" --overWriteDir  --additionalOut \"../popClustering/{TARGET}/locationByIndex/{INDEX}.tab.txt\" "
								"--overWrite --dout {TARGET}{MIDREP}_qlusterOut ";

		std::string qlusterCmdSepMatesTemplate =
				"cd \"" + extractionDirs.string() + "\" && "
						+ "if [ -f {TARGET}{MIDREP}_{MATEFILE}.fastq.gz  ]; then "
						+ setUp.commands_.masterProgram_
						+ " qluster "
								"--fastqgz \"{TARGET}{MIDREP}_{MATEFILE}.fastq.gz\" ";
		if(pars.development){
			qlusterCmdSepMatesTemplate = qlusterCmdSepMatesTemplate +  " --alnInfoDir {TARGET}{MIDREP}_{MATEFILE}_alnCache --development";
		}
		qlusterCmdSepMatesTemplate = qlusterCmdSepMatesTemplate + " --overWriteDir --additionalOut \"../popClustering/{TARGET}-{MATEFILE}/locationByIndex/{INDEX}.tab.txt\" "
								"--overWrite --dout {TARGET}{MIDREP}_{MATEFILE}_qlusterOut "
						    "--illumina --qualThres 25,20"
						;

		if (pars.useKlusterClustering_ || pars.techIsNanoporeOrPacbio()) {
			qlusterCmdTemplate =
							"cd \"" + extractionDirs.string() + "\" && "
							+ "if [ -f {TARGET}{MIDREP}.fastq.gz  ]; then "
							+ " SeekDeep "
							+ " kluster "
								"--fastqgz \"{TARGET}{MIDREP}.fastq.gz\" ";
			if (pars.development) {
				qlusterCmdTemplate = qlusterCmdTemplate + "--alnInfoDir {TARGET}{MIDREP}_alnCache --development ";
			}
			qlusterCmdTemplate = qlusterCmdTemplate +
													 "--overWriteDir  --additionalOut \"../popClustering/{TARGET}/locationByIndex/{INDEX}.tab.txt\" "
													 "--overWrite --dout {TARGET}{MIDREP}_klusterOut ";
      //add in current defaults commonly used
      //qlusterCmdTemplate += " --cutOff 0.05  --sizeCutOff 1%,3  --map --recalcConsensus --writeInitialClusters --qualThres 15,10 -checkIndelsWhenMapping --checkChimeras  ";
      qlusterCmdTemplate += " --sizeCutOff 1%,3  --map --recalcConsensus --writeInitialClusters --qualThres 15,10 -checkIndelsWhenMapping --checkChimeras  ";

//			qlusterCmdTemplate = "cd \"" + extractionDirs.string() + "\" && "
//					+ " if [ -f {TARGET}{MIDREP}.fastq.gz  ]; then "
//											+ " SeekDeep "
//											+ " kluster "
//							"--fastqgz \"{TARGET}{MIDREP}.fastq.gz\" "
//							"--alnInfoDir {TARGET}{MIDREP}_alnCache --overWriteDir "
//							"--additionalOut ../popClustering/locationByIndex/{TARGET}.tab.txt "
//							"--overWrite --dout {TARGET}{MIDREP}_klusterOut ";
		}

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
			//if (njh::in(index, inputPassed)) {
			if (true) {

				auto tarsNames = njh::conToStr(analysisSetup.indexToTars_[index], "_");
				auto cmds = VecStr { extractorCmdTemplate, idTemplate };
				if (bfs::exists(
						njh::files::make_path(analysisSetup.idsDir_,
								analysisSetup.tarsToTargetSubSets_[tarsNames] + "_lenCutOffs.tab.txt"))) {
					cmds.emplace_back(lenCutOffsTemplate);
				}
				bool anyRefs = false;
				for (const auto & tar : analysisSetup.indexToTars_[index]) {
					if (bfs::exists(
							njh::files::make_path(analysisSetup.refsDir_, tar + ".fasta"))) {
						anyRefs = true;
						break;
					}
				}
				if(anyRefs && !pars.refSeqsDir.empty()){
					cmds.emplace_back(refSeqsDir);
				}

				if ("" != analysisSetup.pars_.extraExtractorCmds) {
					cmds.emplace_back(analysisSetup.pars_.extraExtractorCmds);
				}
				if(pars.techIsIllumina()){
					cmds.emplace_back(overLapStatusTemplate);
				}
				auto currentExtractCmd = njh::conToStr(cmds, " ");
				currentExtractCmd = njh::replaceString(currentExtractCmd, "{INDEX}",
						index);
				if(pars.techIsIllumina()){
					//INDEX_R1, INDEX_R2
					auto searchForPairs = readsByPairs.find(index);
					if(searchForPairs != readsByPairs.end()){
						currentExtractCmd = njh::replaceString(currentExtractCmd, "{INDEX_R1}",
								bfs::absolute(searchForPairs->second.first.front()).string() );
						currentExtractCmd = njh::replaceString(currentExtractCmd, "{INDEX_R2}",
								bfs::absolute(searchForPairs->second.second.front()).string() );
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error index: " << index << " was not found in readsByPairs" << "\n";
						ss << "Options: " << njh::conToStr(njh::getVecOfMapKeys(readsByPairs), ", ");
					}
				}else{
					//INDEX_InFile
					auto searchForFile = filesByPossibleName.find(index);
					if(searchForFile != filesByPossibleName.end()){
						currentExtractCmd = njh::replaceString(currentExtractCmd, "{INDEX_InFile}",
								bfs::absolute(searchForFile->second).string());
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error index: " << index << " was not found in filesByPossibleName" << "\n";
						ss << "Options: " << njh::conToStr(njh::getVecOfMapKeys(filesByPossibleName), ", ");
					}
				}
				currentExtractCmd = njh::replaceString(currentExtractCmd, "{TARS}",
						analysisSetup.tarsToTargetSubSets_[tarsNames]);
				extractorCmds.emplace_back(currentExtractCmd);

				std::set<std::string> midsForIndex;
				for(const auto & sampleMid : analysisSetup.samples_.at(index).samples_ ){
					midsForIndex.insert(sampleMid.second.reps_.begin(),sampleMid.second.reps_.end() );
				}
				for (const auto & mid : midsForIndex) {
				//for (const auto & mid : analysisSetup.idsMids_->getMids()) {

					for (const auto & tar : analysisSetup.indexToTars_[index]) {
						auto currentQlusterCmdTemplate = qlusterCmdTemplate;
						if (1
								== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
								&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
										== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front() &&
										analysisSetup.pars_.techIsIllumina()) {
							//by default will analyze R1 and R2 as two separate targets
							currentQlusterCmdTemplate = qlusterCmdSepMatesTemplate;
						}

						if ("" != analysisSetup.pars_.extraQlusterCmds) {
							currentQlusterCmdTemplate += " "
									+ analysisSetup.pars_.extraQlusterCmds;
						}

						if (1 == analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size() &&
								PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP == analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front() &&
								analysisSetup.pars_.techIsIllumina()
										) {
							if(njh::in(tar + "-R1", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
								currentQlusterCmdTemplate += " {EXTRAR1ARGS}";
							}
							if(njh::in(tar + "-R2", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
								currentQlusterCmdTemplate += " {EXTRAR2ARGS}";
							}
						}else{
							if(njh::in(tar, analysisSetup.pars_.extraQlusterCmdsPerTarget)){
								currentQlusterCmdTemplate += " "
										+ analysisSetup.pars_.extraQlusterCmdsPerTarget[tar];
							}
						}

						currentQlusterCmdTemplate += "; fi";

						currentQlusterCmdTemplate = njh::replaceString(
								currentQlusterCmdTemplate, "{INDEX}", index);
						currentQlusterCmdTemplate = njh::replaceString(
								currentQlusterCmdTemplate, "{TARGET}", tar);
						currentQlusterCmdTemplate = njh::replaceString(
								currentQlusterCmdTemplate, "{MIDREP}", mid);
						if (1
								== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
								&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
										== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front() &&
										analysisSetup.pars_.techIsIllumina()) {
							auto r1TemplateCmd = njh::replaceString(
									currentQlusterCmdTemplate, "{MATEFILE}", "R1");
							auto r2TemplateCmd = njh::replaceString(
									currentQlusterCmdTemplate, "{MATEFILE}", "R2");
							if(njh::in(tar + "-R1", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
								r1TemplateCmd= njh::replaceString(
										r1TemplateCmd, "{EXTRAR1ARGS}", analysisSetup.pars_.extraQlusterCmdsPerTarget[tar + "-R1"]);
								r2TemplateCmd = njh::replaceString(r2TemplateCmd, "{EXTRAR1ARGS}","");
							}
							if(njh::in(tar + "-R2", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
								r2TemplateCmd= njh::replaceString(
										r2TemplateCmd, "{EXTRAR2ARGS}", analysisSetup.pars_.extraQlusterCmdsPerTarget[tar + "-R2"]);
								r1TemplateCmd = njh::replaceString(r1TemplateCmd, "{EXTRAR2ARGS}","");
							}
							qlusterCmds.emplace_back(r1TemplateCmd);
							qlusterCmds.emplace_back(r2TemplateCmd);
						} else {
							qlusterCmds.emplace_back(currentQlusterCmdTemplate);
						}
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
    } else if (analysisSetup.pars_.techIsNanoporeOrPacbio()) {
      extractorCmdTemplate = setUp.commands_.masterProgram_
                             + " extractorByKmerMatching --dout {REP}_extraction --overWriteDir  "
                               " --uniqueKmersPerTarget info/ids/{TARS}_uniqueKmers.tab.txt.gz "
                               " --rename ";
    } else {
			 extractorCmdTemplate =
							setUp.commands_.masterProgram_
									+ " extractor --dout {REP}_extraction --overWriteDir ";
			 if(analysisSetup.pars_.techIsIlluminaSingleEnd()){
				 extractorCmdTemplate += " --illumina ";
			 }
		}
		if (analysisSetup.pars_.techIsIllumina()) {
			extractorCmdTemplate += "--fastq1gz \"{REP_R1}\" --fastq2gz \"{REP_R2}\" ";
		} else {
			extractorCmdTemplate += "--fastqgz \"{REP_InFile}\" ";
		}

		std::string lenCutOffsTemplate ="--lenCutOffs info/ids/{TARS}_lenCutOffs.tab.txt ";
		std::string idTemplate = "--id info/ids/{TARS}.id.txt ";
		std::string overLapStatusTemplate = "--overlapStatusFnp info/ids/{TARS}_overlapStatus.tab.txt ";
		std::string refSeqsDir = "--compareSeq info/refs/ ";

		std::string sampleNameTemplate = "--sampleName {REP}";
		//qluster cmds;
		auto extractionDirs = njh::files::make_path(
				bfs::absolute(analysisSetup.dir_), "{REP}_extraction");

		std::string qlusterCmdTemplate = "cd \"" + extractionDirs.string() + "\" && "
				+ " if [ -f {TARGET}{MIDREP}.fastq.gz  ]; then "
										+ setUp.commands_.masterProgram_
										+ " qluster "
						"--fastqgz \"{TARGET}{MIDREP}.fastq.gz\" ";
		if (pars.development) {
			qlusterCmdTemplate = qlusterCmdTemplate +	"--alnInfoDir {TARGET}{MIDREP}_alnCache --development";
		}
		qlusterCmdTemplate = qlusterCmdTemplate +	" --overWriteDir  --additionalOut ../popClustering/locationByIndex/{TARGET}.tab.txt "
						"--overWrite --dout {TARGET}{MIDREP}_qlusterOut ";

		std::string qlusterCmdSepMatesTemplate = "cd \"" + extractionDirs.string() + "\" && "
				+ " if [ -f {TARGET}{MIDREP}_{MATEFILE}.fastq.gz  ]; then "
										+ setUp.commands_.masterProgram_
										+ " qluster "
						"--fastqgz \"{TARGET}{MIDREP}_{MATEFILE}.fastq.gz\" ";
		if (pars.development) {
			qlusterCmdSepMatesTemplate = qlusterCmdSepMatesTemplate + "--alnInfoDir {TARGET}{MIDREP}_{MATEFILE}_alnCache --development ";
		}
		qlusterCmdSepMatesTemplate = qlusterCmdSepMatesTemplate + "--overWriteDir  --additionalOut ../popClustering/locationByIndex/{TARGET}-{MATEFILE}.tab.txt "
						"--overWrite --dout {TARGET}{MIDREP}_{MATEFILE}_qlusterOut "
						"--illumina --qualThres 25,20";

		if(pars.useKlusterClustering_ || pars.techIsNanoporeOrPacbio()){
      qlusterCmdTemplate =
          "cd \"" + extractionDirs.string() + "\" && "
          + "if [ -f {TARGET}{MIDREP}.fastq.gz  ]; then "
          + " SeekDeep "
          + " kluster "
            "--fastqgz \"{TARGET}{MIDREP}.fastq.gz\" ";
			if (pars.development) {
				qlusterCmdTemplate = qlusterCmdTemplate + "--alnInfoDir {TARGET}{MIDREP}_alnCache --development ";
			}
			qlusterCmdTemplate = qlusterCmdTemplate + " --overWriteDir  --additionalOut \"../popClustering/locationByIndex/{TARGET}.tab.txt\" "
            "--overWrite --dout {TARGET}{MIDREP}_klusterOut ";
      //add in current defaults commonly used
      qlusterCmdTemplate += " --sizeCutOff 1%,3  --map --recalcConsensus --writeInitialClusters --qualThres 15,10 -checkIndelsWhenMapping --checkChimeras  ";
      //qlusterCmdTemplate += " --cutOff 0.05  --sizeCutOff 1%,3  --map --recalcConsensus --writeInitialClusters --qualThres 15,10 -checkIndelsWhenMapping --checkChimeras  ";
		}
		if (analysisSetup.pars_.techIsIllumina() || analysisSetup.pars_.techIsIlluminaSingleEnd()) {
			qlusterCmdTemplate += "--illumina --qualThres 25,20";
		} else if (analysisSetup.pars_.techIs454()) {
			qlusterCmdTemplate += "--454";
		} else if (analysisSetup.pars_.techIsIonTorrent()) {
			qlusterCmdTemplate += "--ionTorrent";
		}
		if(setUp.pars_.debug_){
			std::cout << "Sample Files Found: " << std::endl;
			std::cout << njh::conToStr(sampleFilesFound, "\n") << std::endl;

			std::cout << "Targets:" << std::endl;
			std::cout << njh::conToStr(getVectorOfMapKeys(analysisSetup.samples_), "\n") << std::endl;
		}

		std::unordered_map<std::string, VecStr> targetsForReps;
		for (const auto & tars : analysisSetup.samples_) {
			for (const auto & rep : tars.second.getReps()) {
				targetsForReps[rep].emplace_back(tars.first);
			}
		}
		for (auto & rep : targetsForReps) {
			njh::sort(rep.second);
			if(setUp.pars_.debug_){
				std::cout << "Sample: " << rep.first << std::endl;
				std::cout << "\t" << njh::conToStr(rep.second, ", ") << std::endl;
			}
		}

		for (const auto & rep : targetsForReps) {
			if (njh::in(rep.first, sampleFilesFound)
					|| njh::in(njh::replaceString(rep.first, "MID", ""), sampleFilesFound)) {
				std::string fName = rep.first;
//				if (!njh::in(rep.first, samplesExtracted)) {
//					fName = njh::replaceString(rep.first, "MID", "");
//				}

				std::string sampName = rep.first;
				if (!njh::beginsWith(rep.first, "MID")) {
					sampName = "MID" + rep.first;
				}

				auto tarsNames = njh::conToStr(rep.second, "_");
				auto currentSampTemp = njh::replaceString(sampleNameTemplate, "{REP}",
						sampName);
				auto cmds = VecStr { extractorCmdTemplate, idTemplate, currentSampTemp };
				if (bfs::exists(
						njh::files::make_path(analysisSetup.idsDir_,
								analysisSetup.tarsToTargetSubSets_[tarsNames] + "_lenCutOffs.tab.txt"))) {
					cmds.emplace_back(lenCutOffsTemplate);
				}
				if ("" != analysisSetup.pars_.extraExtractorCmds) {
					cmds.emplace_back(analysisSetup.pars_.extraExtractorCmds);
				}
				//add ref files
				bool anyRefs = false;
				for (const auto & tar : rep.second) {
					if (bfs::exists(
							njh::files::make_path(analysisSetup.refsDir_, tar + ".fasta"))) {
						anyRefs = true;
						break;
					}
				}
				if(anyRefs && !pars.refSeqsDir.empty()){
					cmds.emplace_back(refSeqsDir);
				}
				// add overlap status
				if(pars.techIsIllumina()){
					cmds.emplace_back(overLapStatusTemplate);
				}


				auto currentExtractCmd = njh::conToStr(cmds, " ");
				currentExtractCmd = njh::replaceString(currentExtractCmd, "{REP}",
						fName);
				if(pars.techIsIllumina()){
					//INDEX_R1, INDEX_R2
					auto searchForPairs = readsByPairs.find(fName);
					if(searchForPairs != readsByPairs.end()){
						currentExtractCmd = njh::replaceString(currentExtractCmd, "{REP_R1}",
								bfs::absolute(searchForPairs->second.first.front()).string() );
						currentExtractCmd = njh::replaceString(currentExtractCmd, "{REP_R2}",
								bfs::absolute(searchForPairs->second.second.front()).string() );
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error REP: " << fName << " was not found in readsByPairs" << "\n";
						ss << "Options: " << njh::conToStr(njh::getVecOfMapKeys(readsByPairs), ", ");
					}
				}else{
					//INDEX_InFile
					auto searchForFile = filesByPossibleName.find(fName);
					if(searchForFile != filesByPossibleName.end()){
						currentExtractCmd = njh::replaceString(currentExtractCmd, "{REP_InFile}",
								bfs::absolute(searchForFile->second).string());
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error REP: " << fName << " was not found in filesByPossibleName" << "\n";
						ss << "Options: " << njh::conToStr(njh::getVecOfMapKeys(filesByPossibleName), ", ");
					}
				}
				currentExtractCmd = njh::replaceString(currentExtractCmd, "{TARS}",
						analysisSetup.tarsToTargetSubSets_[tarsNames]);
				extractorCmds.emplace_back(currentExtractCmd);

				for (const auto & tar : rep.second) {
					std::string currentQlusterCmdTemplate = qlusterCmdTemplate;
					if(1 == analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
							&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
									== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front()&&
									analysisSetup.pars_.techIsIllumina()){
						currentQlusterCmdTemplate = qlusterCmdSepMatesTemplate;
					}

					if ("" != analysisSetup.pars_.extraQlusterCmds) {
						currentQlusterCmdTemplate += " "
								+ analysisSetup.pars_.extraQlusterCmds;
					}

					if (1
							== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
							&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
									== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front()&&
									analysisSetup.pars_.techIsIllumina()) {
						if(njh::in(tar + "-R1", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
							currentQlusterCmdTemplate += " {EXTRAR1ARGS}";
						}
						if(njh::in(tar + "-R2", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
							currentQlusterCmdTemplate += " {EXTRAR2ARGS}";
						}
					}else{
						if(njh::in(tar, analysisSetup.pars_.extraQlusterCmdsPerTarget)){
							currentQlusterCmdTemplate += " "
									+ analysisSetup.pars_.extraQlusterCmdsPerTarget[tar];
						}
					}

					currentQlusterCmdTemplate += "; fi";
					currentQlusterCmdTemplate = njh::replaceString(
							currentQlusterCmdTemplate, "{REP}", fName);
					currentQlusterCmdTemplate = njh::replaceString(
							currentQlusterCmdTemplate, "{MIDREP}", sampName);
					currentQlusterCmdTemplate = njh::replaceString(
							currentQlusterCmdTemplate, "{TARGET}", tar);
					if(1 == analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
												&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
														== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front()&&
														analysisSetup.pars_.techIsIllumina()){
						auto r1TemplateCmd = njh::replaceString(
								currentQlusterCmdTemplate, "{MATEFILE}", "R1");
						auto r2TemplateCmd = njh::replaceString(
								currentQlusterCmdTemplate, "{MATEFILE}", "R2");
						if(njh::in(tar + "-R1", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
							r1TemplateCmd= njh::replaceString(
									r1TemplateCmd, "{EXTRAR1ARGS}", analysisSetup.pars_.extraQlusterCmdsPerTarget[tar + "-R1"]);
							r2TemplateCmd = njh::replaceString(r2TemplateCmd, "{EXTRAR1ARGS}","");
						}
						if(njh::in(tar + "-R2", analysisSetup.pars_.extraQlusterCmdsPerTarget)){
							r2TemplateCmd= njh::replaceString(
									r2TemplateCmd, "{EXTRAR2ARGS}", analysisSetup.pars_.extraQlusterCmdsPerTarget[tar + "-R2"]);
							r1TemplateCmd = njh::replaceString(r1TemplateCmd, "{EXTRAR2ARGS}","");
						}
						qlusterCmds.emplace_back(r1TemplateCmd);
						qlusterCmds.emplace_back(r2TemplateCmd);
					}else{
						qlusterCmds.emplace_back(currentQlusterCmdTemplate);
					}
				}
			}
		}
	}

	OutputStream extractorCmdsFile(
			njh::files::make_path(analysisSetup.dir_, "extractorCmds.txt"));
	njh::sort(extractorCmds);
	printVector(extractorCmds, "\n", extractorCmdsFile);

	OutputStream qlusterCmdsFile(
			njh::files::make_path(analysisSetup.dir_, "qlusterCmds.txt"));
	njh::sort(qlusterCmds);
	printVector(qlusterCmds, "\n", qlusterCmdsFile);

	//process cluster cmds
	std::string processClusterTemplate =
			setUp.commands_.masterProgram_
					+ " processClusters "
							"--alnInfoDir alnCache --strictErrors --dout analysis --fastqgz output.fastq.gz --overWriteDir ";

  if(pars.techIsNanoporeOrPacbio()){
    processClusterTemplate += " --allowHomopolymerCollapse --qualThres 15,10 ";
  }

	if(pars.techIsIllumina() || pars.techIsIlluminaSingleEnd()){
		auto lowerCaseExtracProcessArgs = stringToLowerReturn(analysisSetup.pars_.extraProcessClusterCmds);
		if(std::string::npos == lowerCaseExtracProcessArgs.find("--illumina")){
			processClusterTemplate  += " --illumina";
		}
	}
	if (!analysisSetup.pars_.conservative) {

		processClusterTemplate += " --removeOneSampOnlyOneOffHaps --excludeCommonlyLowFreqHaplotypes --excludeLowFreqOneOffs --rescueExcludedOneOffLowFreqHaplotypes";
		//processClusterTemplate += " --excludeCommonlyLowFreqHaplotypes --excludeLowFreqOneOffs --fracCutOff 0 ";
	}

	if (!analysisSetup.pars_.conservative && analysisSetup.pars_.rescueFilteredHaplotypes) {
		//processClusterTemplate += " --rescueExcludedOneOffLowFreqHaplotypes --rescueMatchingExpected --rescueExcludedChimericHaplotypes";
		processClusterTemplate += " --rescueMatchingExpected --rescueExcludedChimericHaplotypes";
	}

	if ("" != analysisSetup.pars_.extraProcessClusterCmds) {
		processClusterTemplate += " " + analysisSetup.pars_.extraProcessClusterCmds;
	}



	VecStr processClusterCmds;
	if (nullptr != analysisSetup.groupMetaData_) {
		processClusterTemplate += " --groupingsFile \""
				+ njh::files::make_path(bfs::absolute(analysisSetup.infoDir_),
						"groupMeta.tab.txt").string() + "\"";
	}

	auto popDir = njh::files::make_path(bfs::absolute(analysisSetup.dir_),
			"popClustering");

	for (const auto & tar : targets) {
		if (1 == analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
				&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
						== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front()&&
						analysisSetup.pars_.techIsIllumina()) {
			//by default will analyze R1 and R2 as two separate targets
			for(const auto & mate : VecStr{"-R1", "-R2"}){
				if(njh::in(tar + mate, analysisSetup.pars_.extraProcessClusterCmdsPerTarget)){
					processClusterTemplate += " " + analysisSetup.pars_.extraProcessClusterCmdsPerTarget[tar + mate];
				}
				std::stringstream processClustersCmdsStream;
				processClustersCmdsStream
						<< "cd \"" + njh::files::make_path(popDir, tar + mate ).string() + "\" && "
								+ processClusterTemplate + " --experimentName " + tar+ mate;
				auto refSeqFnp = njh::files::make_path(pars.outDir, "info/refs/" + tar+ mate + ".fasta");
				if(bfs::exists(refSeqFnp) && !pars.refSeqsDir.empty()){
					processClustersCmdsStream << " --ref " << bfs::absolute(refSeqFnp);
				}
        if(bfs::exists(refSeqFnp) && !pars.previousPopSeqsDir.empty()){
          processClustersCmdsStream << " --previousPop " << bfs::absolute(refSeqFnp);
        }
				processClusterCmds.emplace_back(processClustersCmdsStream.str());
			}
		}else{
			if(njh::in(tar, analysisSetup.pars_.extraProcessClusterCmdsPerTarget)){
				processClusterTemplate += " " + analysisSetup.pars_.extraProcessClusterCmdsPerTarget[tar];
			}
			std::stringstream processClustersCmdsStream;
			processClustersCmdsStream
					<< "cd \"" + njh::files::make_path(popDir, tar).string() + "\" && "
							+ processClusterTemplate + " --experimentName " + tar;
			auto refSeqFnp = njh::files::make_path(pars.outDir, "info/refs/" + tar + ".fasta");
      if(bfs::exists(refSeqFnp) && !pars.refSeqsDir.empty()){
        processClustersCmdsStream << " --ref " << bfs::absolute(refSeqFnp);
      }
      if(bfs::exists(refSeqFnp) && !pars.previousPopSeqsDir.empty()){
        processClustersCmdsStream << " --previousPop " << bfs::absolute(refSeqFnp);
      }
			processClusterCmds.emplace_back(processClustersCmdsStream.str());
		}
	}
	OutOptions processClusterCmdsOpts(
			njh::files::make_path(analysisSetup.dir_, "processClusterCmds.txt"));
	std::ofstream processClusterCmdsFile;
	openTextFile(processClusterCmdsFile, processClusterCmdsOpts);
	printVector(processClusterCmds, "\n", processClusterCmdsFile);

	//gen analysis configs
	std::string genConfigTemplate =
			setUp.commands_.masterProgram_
					+ " genProjectConfig --projectName {TARGET} --out serverConfigs/{TARGET}.config";

	VecStr genConfigCmds;
	for (const auto & tar : targets) {
		if (1 == analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.size()
				&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
						== analysisSetup.idsMids_->targets_.at(tar).overlapStatuses_.front()&&
						analysisSetup.pars_.techIsIllumina()) {
			//by default will analyze R1 and R2 as two separate targets
			genConfigCmds.emplace_back(
					njh::replaceString(
							genConfigTemplate + " --overWrite --mainDir popClustering/{TARGET}/analysis",
							"{TARGET}", tar + "-R1"));
			genConfigCmds.emplace_back(
					njh::replaceString(
							genConfigTemplate + " --overWrite --mainDir popClustering/{TARGET}/analysis",
							"{TARGET}", tar + "-R2"));
		} else {
			genConfigCmds.emplace_back(
					njh::replaceString(
							genConfigTemplate + " --overWrite --mainDir popClustering/{TARGET}/analysis",
							"{TARGET}", tar));
		}
	}
	OutOptions genConfigCmdsOpts(
			njh::files::make_path(analysisSetup.dir_, "genConfigCmds.txt"));
	OutputStream genConfigCmdsFile(genConfigCmdsOpts);
	printVector(genConfigCmds, "\n", genConfigCmdsFile);

	//start server config
	OutOptions startServerCmdOpts(
			njh::files::make_path(analysisSetup.dir_, "startServerCmd.sh"));
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
			<< " popClusteringViewer --verbose --configDir \"$(pwd)/serverConfigs\" --port $1 --name $2 &"
			<< std::endl;
	startServerCmdFile << "else" << std::endl;
	startServerCmdFile << "	nohup " << setUp.commands_.masterProgram_
			<< " popClusteringViewer --verbose --configDir \"$(pwd)/serverConfigs\" & "
			<< std::endl;
	startServerCmdFile << "fi" << std::endl;
	//make file executable
	chmod(startServerCmdOpts.outFilename_.c_str(),
			S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	//combining extraction counts
	OutOptions combineExtractionCmdOpts(
			njh::files::make_path(analysisSetup.dir_, "combineExtractionCountsCmd.sh"));
	OutputStream combineExtractionCmdOut(combineExtractionCmdOpts);
	combineExtractionCmdOut << "#!/usr/bin/env bash" << std::endl;
	//combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains allPrimerCounts.tab.txt --delim tab --header --out reports/combinedAllPrimerCounts.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains extractionProfile.tab.txt --delim tab --header --out reports/allExtractionProfile.tab.txt  --overWrite" << std::endl;
	combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains extractionStats.tab.txt --delim tab --header --out reports/allExtractionStats.tab.txt  --overWrite" << std::endl;
	if(analysisSetup.pars_.techIsIllumina()){
		combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains processPairsCounts.tab.txt --delim tab --header --out reports/allProcessPairsCounts.tab.txt  --overWrite" << std::endl;
		combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains allFailedPrimerCounts.tab.txt --delim tab --header --out reports/combinedAllFailedPrimerCounts.tab.txt  --overWrite" << std::endl;
	}
	if(analysisSetup.idsMids_->containsMids()){
		combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains midCounts.tab.txt --delim tab --header --out reports/allMidCounts.tab.txt  --overWrite" << std::endl;
		combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains top_mostCommonR1Starts_for_unrecognizedBarcodes.tab.txt --delim tab --header --out reports/allTop_mostCommonR1Starts_for_unrecognizedBarcodes.tab.txt  --overWrite" << std::endl;
		combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains top_mostCommonR2Starts_for_unrecognizedBarcodes.tab.txt --delim tab --header --out reports/allTop_mostCommonR2Starts_for_unrecognizedBarcodes.tab.txt  --overWrite" << std::endl;
		combineExtractionCmdOut << setUp.commands_.masterProgram_ << " rBind --recursive --depth 1 --contains top_mostCommonR1AndR2Starts_for_unrecognizedBarcodes.tab.txt --delim tab --header --out reports/allTop_mostCommonR1AndR2Starts_for_unrecognizedBarcodes.tab.txt  --overWrite" << std::endl;
	}
	chmod(combineExtractionCmdOpts.outFilename_.c_str(),
			S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	//run file analysis file
	OutOptions runAnalysisOpts(
			njh::files::make_path(analysisSetup.dir_, "runAnalysis.sh"));
	std::ofstream runAnalysisFile;
	openTextFile(runAnalysisFile, runAnalysisOpts);
	runAnalysisFile << "#!/usr/bin/env bash" << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "##run all parts of the pipeline" << std::endl;
	runAnalysisFile
			<< "##except start the server which is done by running ./startServerCmd.sh"
			<< std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "numThreads=" << pars.numThreads << std::endl;
	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "if [[ $# -eq 1 ]]; then" << std::endl;
	runAnalysisFile << "	numThreads=$1" << std::endl;
	runAnalysisFile << "fi\n" << std::endl;
	runAnalysisFile << R"(if [[ $# -gt 1 ]]; then
  echo "Illegal number of parameters, should either be no arguments or a number to indicate the number of threads"
  echo "Examples:"
  echo "Example 1"
  echo "Run with the default number of threads, which was set when running setupTarAmpAnalysis"
  echo "./runAnalysis.sh"
  echo ""
  echo "Example 2"
  echo "Run with 7 threads"
  echo "./runAnalysis.sh 7"
  exit
fi)" << std::endl;

	runAnalysisFile << "" << std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile extractorCmds.txt      --numThreads $numThreads --raw --logDir logs "
			<< std::endl;
	runAnalysisFile << "./combineExtractionCountsCmd.sh"<< std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile qlusterCmds.txt        --numThreads $numThreads --raw --logDir logs "
			<< std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile processClusterCmds.txt --numThreads $numThreads --raw --logDir logs "
			<< std::endl;
	runAnalysisFile << "" << setUp.commands_.masterProgram_
			<< " runMultipleCommands --cmdFile genConfigCmds.txt      --numThreads $numThreads --raw --logDir logs "
			<< std::endl;
  runAnalysisFile << "/bin/ls popClustering/*/analysis/selectedClustersInfo.tab.txt.gz | elucidator rBind --files STDIN --delim tab --header --overWrite --out popClustering/allSelectedClustersInfo.tab.txt.gz" << std::endl;
	runAnalysisFile << "" << std::endl;
	//make file executable
	chmod(runAnalysisOpts.outFilename_.c_str(),
			S_IWUSR | S_IRUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IEXEC | S_IXGRP);

	if (foundErrors) {
		std::cerr << njh::bashCT::flashing << njh::bashCT::red << njh::bashCT::bold
				<< "ERRORS FOUND!!" << njh::bashCT::reset << std::endl;
		std::cerr << "Read Warnings in "
				<< njh::bashCT::boldRed(wrningsOpts.outFilename_.string()) << std::endl;
	}

	return 0;
}


}  // namespace njhseq
