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

	TarAmpSeqInvestigator::TarAmpSeqInvestigatorPars investPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(investPars.testNumber, "--testNumber", "Just use this number of reads of the top of the file");
	setUp.setOption(investPars.dontCollapsePossibleMIDs, "--dontCollapsePossibleMIDs",
			"Don't Collapse Possible MIDs", false);
	setUp.setOption(investPars.unrecogBaseSampling, "--unrecogBaseSampling",
			"Number of bases to sample from file for unrecognized sequences", false);

	setUp.setOption(investPars.precdingBaseFreqCutOff, "--precdingBaseFreqCutOff", "Preceding Base Freq Cut Off", false);

	investPars.pars.corePars_.pDetPars.primerWithin_ = 30;
	setUp.setOption(investPars.pars.corePars_.pDetPars.primerWithin_, "--primerWithin", "Primer Within bases search", false, "Primer");


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
	investPars.gapInfo_ = setUp.pars_.gapInfo_;
	setUp.setOption(investPars.idFnp, "--id", "SeekDeep primers file", true);
	setUp.processReadInNames(VecStr{"--fastq1", "--fastq", "--fasta", "--fastq1gz", "--fastqgz", "--fastagz"});
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);



	TarAmpSeqInvestigator investigator(investPars);
	investigator.investigateFile(setUp.pars_.ioOptions_, setUp.pars_.verbose_);
	investigator.processCounts();
	investigator.writeOutTables(setUp.pars_.directoryName_, true);


	std::stringstream ss;
	auto possibleRevComp = investigator.reverseComplementLikely();
	auto possiblePrecedingRandomeBases = investigator.hasPossibleRandomPrecedingBases(investigator.ids_.getMaxMIDSize());
	ss << "Has Possible Reverse Complement directed reads: " << njh::boolToStr(possibleRevComp) << std::endl;
	ss << "Has Possible Random Preceding bases: " << njh::boolToStr(possiblePrecedingRandomeBases) << std::endl;
	auto recFlags = investigator.recommendSeekDeepExtractorFlags();
	if(!recFlags.empty()){
		ss << "Recommended SeekDeep extractor additional flags: " << std::endl;
		ss << njh::conToStr(recFlags, " ")<< std::endl;
	} else {
		ss << "No additional recommended SeekDeep extractor flags" << std::endl;
	}

	OutputStream seekdeepFlagRecs(njh::files::make_path(setUp.pars_.directoryName_, "message.txt"));
	seekdeepFlagRecs << ss.str();
	if(setUp.pars_.verbose_){
		std::cout << ss.str();
	}
	return 0;

}




}  // namespace njhseq

