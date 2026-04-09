/*
 * gatherInfoOnTargetedAmpliconSeqFile.cpp
 *
 *  Created on: Oct 16, 2019
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"

#include "SeekDeep/objects.h"
#include "SeekDeep/parameters.h"


namespace njhseq {



int SeekDeepUtilsRunner::gatherInfoOnTargetedAmpliconSeqFile(
		const njh::progutils::CmdArgs & inputCommands) {
	TarAmpAnalysisSetup::TarAmpPars tar_amp_pars;
	TarAmpSeqInvestigator::TarAmpSeqInvestigatorPars investPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	tar_amp_pars.numberOfFilesToInvestigate = 20;
	investPars.testNumber = 5000;
  setUp.setOption(investPars.max_len_to_investigate, "--max_len_to_investigate", "don't investigate sequences above this length");

  setUp.setOption(investPars.testNumber, "--testNumber", "Just use this number of reads of the top of the file");
	setUp.setOption(tar_amp_pars.numberOfFilesToInvestigate, "--numberOfFilesToInvestigate", "Number of files to investigate when adding additional recommended flags", false, "Extra Commands");

	setUp.setOption(investPars.dontCollapsePossibleMIDs, "--dontCollapsePossibleMIDs",
			"Don't Collapse Possible MIDs", false);
	setUp.setOption(investPars.unrecogBaseSampling, "--unrecogBaseSampling",
			"Number of bases to sample from file for unrecognized sequences", false);
	setUp.setOption(tar_amp_pars.extraExtractorCmds, "--extraExtractorCmds",
			"Extra extractor cmds to add to the defaults", false, "Extra Commands");
	setUp.setOption(investPars.precdingBaseFreqCutOff, "--precdingBaseFreqCutOff", "Preceding Base Freq Cut Off", false);
	setUp.setOption(tar_amp_pars.numThreads, "--numThreads", "number of Threads to utilize", false);

	investPars.pars.corePars_.pDetPars.primerWithin_ = 30;
	setUp.setOption(investPars.pars.corePars_.pDetPars.primerWithin_, "--primerWithin", "Primer Within bases search", false, "Primer");
	setUp.setOption(investPars.fracUndeterminedToTriggerRecount_, "--fracUndeterminedToTriggerRecount", "fraction of undetermined primers to trigger recount", false, "Primer");

	bool primerToUpperCase = false;
	setUp.setOption(primerToUpperCase, "--primerUpper",
			"Leave primers in upper case", false, "Primer");
	investPars.pars.corePars_.pDetPars.primerToLowerCase_ = !primerToUpperCase;
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.distances_.query_.coverage_, "--primerCoverage",
			"Amount of primers found", false, "Primer");
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.hqMismatches_, "--primerNumOfMismatches",
			"Number of Mismatches to allow in primers", false, "Primer");
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.oneBaseIndel_, "--primerOneBaseIndels",
			"Number Of One base indels to allow in primers", false, "Primer");
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.twoBaseIndel_, "--primerTwoBaseIndels",
			"Number Of Two base indels to allow in primers", false, "Primer");

	setUp.pars_.gapInfo_.gapOpen_ = 5;
	setUp.pars_.gapInfo_.gapExtend_ = 1;
	setUp.pars_.gap_ = "5,1";
	setUp.pars_.gapInfo_.gapRightQueryOpen_ = 0;
	setUp.pars_.gapInfo_.gapRightQueryExtend_ = 0;
	setUp.pars_.gapInfo_.gapRightRefOpen_ = 0;
	setUp.pars_.gapInfo_.gapRightRefExtend_ = 0;
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapInfo_.gapLeftQueryOpen_ = 0;
	setUp.pars_.gapInfo_.gapLeftQueryExtend_ = 0;
	setUp.pars_.gapInfo_.gapLeftRefOpen_ = 0;
	setUp.pars_.gapInfo_.gapLeftRefExtend_ = 0;
	setUp.pars_.gapLeft_ = "0,0";
	setUp.processGap();
	investPars.verbose_ = setUp.pars_.verbose_;
	investPars.gapInfo_ = setUp.pars_.gapInfo_;
	setUp.setOption(investPars.idFnp, "--id", "SeekDeep primers file", true);
	setUp.setOption(tar_amp_pars.inputDir, "--reads_dir", "a directory of sequence files to investigate instead of a single input read, files to investigate will depend on the technology set", false);
	VecStr acceptableTechs{"454", "IonTorrent", "Illumina", "Illumina-SingleEnd", "Nanopore", "Pacbio"};

	setUp.setOption(tar_amp_pars.technology, "--technology",
			"Sequencing Technology (should be " + njh::conToStrEndSpecial(acceptableTechs, ", ", " or ") + ")",
			!tar_amp_pars.inputDir.empty(), "Technology");
	njh::for_each(acceptableTechs, [](std::string & tech){
		stringToLower(tech);
	});
	stringToLower(tar_amp_pars.technology);
	if (!tar_amp_pars.inputDir.empty() && !njh::in(tar_amp_pars.technology, acceptableTechs)) {
		setUp.failed_ = true;
		std::stringstream ss;
		ss
				<< "Error in setting technology, should be "
				<< njh::conToStrEndSpecial(acceptableTechs, ", ", " or ")
				<< " not "
				<< tar_amp_pars.technology << "\n";
		setUp.addWarning(ss.str());
	}
	setUp.processReadInNames(VecStr{"--fastq1", "--fastq", "--fasta", "--fastq1gz", "--fastqgz", "--fastagz"}, tar_amp_pars.inputDir.empty());

	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);

	std::function investigate = [](const SeqIOOptions& seqOpts, const TarAmpSeqInvestigator::TarAmpSeqInvestigatorPars& investPars) {
		std::shared_ptr<TarAmpSeqInvestigator> investigator = std::make_shared<TarAmpSeqInvestigator>(investPars);
		auto prepCounts = investigator->prepareForInvestiagteFile(seqOpts);
		investigator->investigateFile(seqOpts, prepCounts);
		investigator->processCounts();
		if (investPars.verbose_) {
			std::cout << "investigator->getFractionOfUnrecognizedPrimers(): " << investigator->getFractionOfUnrecognizedPrimers() << std::endl;
		}
		if (investigator->getFractionOfUnrecognizedPrimers() >= investPars.fracUndeterminedToTriggerRecount_) {
			investigator->pars_.pars.corePars_.pDetPars.primerWithin_ = static_cast<uint32_t>(std::round(prepCounts.readMedian));
			if (investigator->ids_.containsMids()) {
				investigator->pars_.pars.corePars_.primIdsPars.mPars_.searchStop_ = static_cast<uint32_t>(std::round(prepCounts.readMedian));
				investigator->ids_.initMidDeterminator(investigator->pars_.midPars);
			}
			if (investPars.verbose_) {
				std::cout << "Fraction of reads with undetermined primers " << investigator->getFractionOfUnrecognizedPrimers() <<
						" is more than fracUndeterminedToTriggerRecount " << investPars.fracUndeterminedToTriggerRecount_ << " so will set primer within to " <<
						investigator->pars_.pars.corePars_.pDetPars.primerWithin_ << " and recount" << std::endl;
			}

			investigator->resetCounts();
			investigator->investigateFile(seqOpts, prepCounts);
			investigator->processCounts();
			if (investPars.verbose_) {
				std::cout << "Fraction of reads with undetermined primers after recount is " << investigator->getFractionOfUnrecognizedPrimers() << std::endl;
			}
		}
		return investigator;
	};
	auto masterInvestigator = std::make_shared<TarAmpSeqInvestigator>(investPars);
	std::mutex masterInvesMut;
	if (tar_amp_pars.inputDir.empty()) {
		masterInvestigator = investigate(setUp.pars_.ioOptions_, investPars);
	} else {
		std::map<std::string, std::pair<VecStr, VecStr>> readsByPairs ;
		std::map<std::string, bfs::path> filesByPossibleName;
		auto guessedSamples = GuessPossibleSamps(tar_amp_pars);
		auto expectedSamples = guessedSamples.getColumn("sample");
		std::regex inputFilePat( tar_amp_pars.inputFilePat);
		auto files = njh::files::listAllFilesThrowOnDupSymlink(tar_amp_pars.inputDir.string(), false, {inputFilePat});
		ReadPairsOrganizer rpOrganizer{expectedSamples};
		rpOrganizer.illuminaPat_ = tar_amp_pars.illuminaInputFilePat;
		if (tar_amp_pars.techIsIllumina()) {
			rpOrganizer.processFiles(files);
			readsByPairs = rpOrganizer.processReadPairs();
		} else {
			std::regex filePatReg{tar_amp_pars.inputFilePat};
			for (const auto & file : files) {
				auto fNameNoExt = njh::files::removeExtension(file.first.filename());
				if (njh::in(fNameNoExt, expectedSamples)) {
					filesByPossibleName[fNameNoExt] = file.first;
				}
			}
		}
		//investigate input seq files to recommend
		std::vector<SeqIOOptions> filesToInvestigate;
		njh::randomGenerator rgen;

		if (tar_amp_pars.techIsIllumina()) {
			double fractionToBeat = tar_amp_pars.numberOfFilesToInvestigate/static_cast<double>(readsByPairs.size());
			for (const auto& pair: readsByPairs) {
				if (rgen.unifRand() <= fractionToBeat) {
					if (njh::endsWith(pair.second.first.front(), ".gz")) {
						filesToInvestigate.emplace_back(
							SeqIOOptions::genPairedInGz(bfs::path(pair.second.first.front()), bfs::path(pair.second.second.front())));
					} else {
						filesToInvestigate.emplace_back(
							SeqIOOptions::genPairedIn(bfs::path(pair.second.first.front()), bfs::path(pair.second.second.front())));
					}
					if (filesToInvestigate.size() > tar_amp_pars.numberOfFilesToInvestigate) {
						break;
					}
				}
			}
		} else {
			double fractionToBeat = tar_amp_pars.numberOfFilesToInvestigate/static_cast<double>(filesByPossibleName.size());
			for (const auto& file: filesByPossibleName) {
				if (rgen.unifRand() <= fractionToBeat) {
					filesToInvestigate.emplace_back(file.second, SeqIOOptions::getInFormatFromFnpExcludePaired(file.second), false);
					if (filesToInvestigate.size() > tar_amp_pars.numberOfFilesToInvestigate) {
						break;
					}
				}
			}
		}
		njh::concurrent::LockableQueue<SeqIOOptions> optsQueue(filesToInvestigate);
		std::function investigateFile = [&optsQueue,&investPars,&masterInvesMut,&masterInvestigator,&investigate](){
			SeqIOOptions seqOpts;
			TarAmpSeqInvestigator current_masterInvestigator(investPars);
			while(optsQueue.getVal(seqOpts)){
				if(investPars.verbose_){
					std::cout << "Investigating " << seqOpts.firstName_ << " " << (seqOpts.secondName_.empty() ? std::string("") : seqOpts.secondName_.string()) << std::endl;
				}
				auto investigator = investigate(seqOpts, investPars);
				current_masterInvestigator.addOtherCounts(*investigator);
			}
			{
				std::lock_guard<std::mutex> lock(masterInvesMut);
				masterInvestigator->addOtherCounts(current_masterInvestigator);
			}
		};
		njh::concurrent::runVoidFunctionThreaded(investigateFile, tar_amp_pars.numThreads);
		masterInvestigator->processCounts();
	}

	masterInvestigator->writeOutTables(setUp.pars_.directoryName_, true);

	std::stringstream ss;
	auto possibleRevComp = masterInvestigator->reverseComplementLikely();
	auto possiblePrecedingRandomeBases = masterInvestigator->hasPossibleRandomPrecedingBases(masterInvestigator->ids_.getMaxMIDSize());
	ss << "Has Possible Reverse Complement directed reads: " << njh::boolToStr(possibleRevComp) << std::endl;
	ss << "Has Possible Random Preceding bases: " << njh::boolToStr(possiblePrecedingRandomeBases) << std::endl;
	auto recFlags = masterInvestigator->recommendSeekDeepExtractorFlags();
	VecStr flagsToAdd;
	auto currentExtraExtractorCmds = njh::strToLowerRet(tar_amp_pars.extraExtractorCmds);

	for(const auto & recFlag : recFlags){
		auto rflag = njh::strToLowerRet(recFlag);
		trimAtFirstWhitespace(rflag);
		njh::lstrip(rflag, '-');
		if(std::string::npos == currentExtraExtractorCmds.find(rflag)){
			// std::cout << "recFlag: " << recFlag << std::endl;
			// std::cout << "recFlag: " << recFlag << std::endl;

			if(!(tar_amp_pars.techIsNanoporeOrPacbio() && "checkrevcomplementforprimers" == rflag) &&
			!(tar_amp_pars.techIsNanoporeOrPacbio() && "checkrevcomplementformids" == rflag) &&
			!(tar_amp_pars.techIsNanoporeOrPacbio() && njh::beginsWith(rflag, "midwithinstart") )) {
				// std::cout << "rflag: " << rflag << std::endl;
				flagsToAdd.emplace_back(recFlag);
			}
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
		tar_amp_pars.extraExtractorCmds.append(" ");
		tar_amp_pars.extraExtractorCmds.append(addingStr);
	}
	if(!recFlags.empty()){
		ss << "Recommended SeekDeep extractor additional flags: " << std::endl;
		ss << njh::conToStr(recFlags, " ")<< std::endl;
	} else {
		ss << "No additional recommended SeekDeep extractor flags" << std::endl;
	}
	OutputStream outSeekDeepExtractorFlags(njh::files::make_path(setUp.pars_.directoryName_, "outSeekDeepExtractorFlags.txt"));
	outSeekDeepExtractorFlags << tar_amp_pars.extraExtractorCmds << std::endl;

	OutputStream raw_outSeekDeepExtractorFlags(njh::files::make_path(setUp.pars_.directoryName_, "raw_outSeekDeepExtractorFlags.txt"));
	raw_outSeekDeepExtractorFlags << njh::conToStr(recFlags, "\n") << std::endl;
	//
	OutputStream outSeekDeepMessage(njh::files::make_path(setUp.pars_.directoryName_, "message.txt"));
	outSeekDeepMessage << ss.str();
	if(setUp.pars_.verbose_){
		std::cout << ss.str();
	}
	return 0;

}




}  // namespace njhseq

