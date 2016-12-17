/*
 * extractorSetUp.cpp
 *
 *  Created on: Dec 17, 2016
 *      Author: nick
 */

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
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {

void SeekDeepSetUp::setUpExtractor(extractorPars & pars) {
	//ioOptions_.lowerCaseBases = "remove";
	if (needsHelp()) {
		std::stringstream tempOut;
		tempOut << "extractor" << std::endl;
		tempOut << "Extracts sequences from fasta/fastq  files and some "
				"pre-filtering" << std::endl;
		tempOut << "Commands, order not necessary and case insensitive"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Required commands" << bib::bashCT::reset
				<< std::endl;
		tempOut << "Input can be fasta or fastq" << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		printInputUsage(std::cout);
		tempOut
				<< "2) -id [option]: Full name of the id file with primers and Ids info"
				<< std::endl;

		tempOut << "\tThe id file should be tab delimited and contains the "
				"primer and reverse primer and MIDs if the data is multiplex, "
				"example below" << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		std::cout << "\tgene	forwardPrimer	reversePrimer" << std::endl;
		std::cout << "\tPFMSP1	AACTAGAAGCTTTAGAAGATGCA	ACATATGATTGGTTAAATCAAAG"
				<< std::endl;
		std::cout << "\tPFMSP2	AGATGCAGCTTTAACTAGAAGAA	TTAAAACATATGATTGGTCAAAG"
				<< std::endl;
		std::cout << "\tid	barcode" << std::endl;
		std::cout << "\tMID01	ACGAGTGCGT" << std::endl;
		std::cout << "\tMID02	ACGCTCGACA" << std::endl;
		tempOut << bib::bashCT::bold << "Optional commands" << bib::bashCT::reset
				<< std::endl;
		tempOut << "-dout [option]: Name of an output directory, will default "
				"to the name of the input file plus the date if none given"
				<< std::endl;
		tempOut << "-multiplex : Whether the id file is multiplex with "
				"multiple MIDs, the program will assume only primers are "
				"given if this is not switched on" << std::endl;
		tempOut << "-checkComplement : Check the complement of the sequences "
				"for the barcodes as well" << std::endl;
		tempOut
				<< "-rename : Rename the sequence to the associated primers and mids "
				<< std::endl;
		tempOut << bib::bashCT::boldBlack("Quality filtering options") << std::endl;
		tempOut << " -minLen [option]: The minimum length for the read to be "
				"extracted including primers, defaults to " << pars.minLen << std::endl;
		tempOut << " -maxLen [option]: The maximum length for the read to be "
				"extracted including primers, defaults to " << pars.maxLength
				<< std::endl;
		tempOut << " -qualWindow [option]: Quality window checking options, "
				"given separated by commas and in the order of windowSize, "
				"windowStep, and windowMinimumThres, ex. 50,5,20 (default)"
				<< std::endl;
		tempOut << "-qualCheck	Qual Check Level; default" << pars.qualCheck
				<< std::endl;
		tempOut
				<< "-qualCheckCutOff	Cut Off for fraction of bases above qual check of "
				<< pars.qualCheck << "; default=" << pars.qualCheckCutOff << std::endl;
		tempOut << " -numberOfNs [option]: the number of N's(the symbol for "
				"unknown base) to be allowed in sequence, defaults to 0" << std::endl;
		tempOut << " -samllFragmentCutOff [option]: the size cutoff for a "
				"sequence to be considered simply a fragment, defaults to 50"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Options for looking for primers"
				<< bib::bashCT::reset << std::endl;
		tempOut << "-noReverse: if you don't want to look for the reverse "
				"primer as well as the forward primer, will default to "
				"looking for both" << std::endl;
		tempOut << "-rPrimerCoverage [option] : fraction of the reverse "
				"primer that needs to be present, defaults to "
				<< pars.rPrimerErrors.distances_.query_.coverage_ << "  so at least "
						"half of the primer needs to be present" << std::endl;
		tempOut
				<< "-rPrimerMismatches [option] : Number of mismatches to allow in reverse primer "
						", defaults to " << pars.rPrimerErrors.hqMismatches_ << std::endl;
		tempOut
				<< "-rUpper: set this to make the reverse primer into "
						"upper case, the default is to make it and everything that follows it lower case"
				<< std::endl;
		tempOut
				<< "-noForwardPrimer: Don't look for the forward primer, this means only one primer pair can be put in the id file"
				<< std::endl;
		tempOut << "-fPrimerCoverage [option] : fraction of the forward "
				"primer that needs to be present, defaults to "
				<< pars.fPrimerErrors.distances_.query_.coverage_ << std::endl;
		tempOut
				<< "-fPrimerMismatches [option] : Number of mismatches to allow in forward primer "
						", defaults to " << pars.fPrimerErrors.hqMismatches_ << std::endl;
		tempOut
				<< "-fUpper: set this to make the forward primer into "
						"upper case, the default is to make everything up to and including to it lower case"
				<< std::endl;
		tempOut << bib::bashCT::bold
				<< "Options for screening for possible contamination"
				<< bib::bashCT::reset << std::endl;
		tempOut
				<< "-contamination : Where to screen for possible "
						"contamination, a comparison sequence needs to be supplied if "
						"this is turned on, this option is supported for when using multiple primer pairs yet"
				<< std::endl;
		tempOut
				<< "-compareSeq [option]: A comparison sequence to which the "
						"reads are compared and and if they are dissimilar enough they will be removed as contamination"
				<< std::endl;
		tempOut
				<< "-contaminationKLen [option]:	Contamination Kmer Length, defaults to "
				<< pars.contaminationKLen << std::endl;
		tempOut
				<< "-kmerCutOff [option]:	Kmer cut off for contamination check, defaults to "
				<< pars.kmerCutOff << std::endl;
		tempOut << bib::bashCT::bold
				<< "Options for when extracting with multiple targets"
				<< bib::bashCT::reset << std::endl;
		tempOut
				<< "--multipleTargets : A flag to indicate you are extracting on multiple targets"
				<< std::endl;
		tempOut
				<< "--lenCutOffs [option]: a filename containing specific cut offs for each target need three columns, target,minlen,and maxlen."
						"  If target is not found here length cut offs will default to options assoiciated with -minLen and -maxLen"
				<< std::endl;
		tempOut
				<< "-compareSeq [option]: When the --multipleTargets flag is on this can be a fasta or fastq file containing a different seq for each target"
				<< std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		std::cout << "examples: " << std::endl
				<< "\tSeekDeep extractor -fasta reads.fasta -id "
						"idFile.tab.txt -dout outPutDir " << std::endl
				<< "\tSeekDeep extractor -stub reads -id "
						"idFile.tab.txt -dout outPutDir " << std::endl
				<< "\tSeekDeep extractor "
						"-fastq reads.fastq -id idFile.tab.txt" << std::endl;
		tempOut.str(std::string());
		tempOut << bib::bashCT::bold << "Output Files:" << bib::bashCT::reset
				<< std::endl;
		tempOut
				<< "extractionProfile.tab.txt: This breaks down the filtering per final extraction sequence file"
				<< std::endl;
		tempOut
				<< "extractionStats.tab.txt: This has info how the whole extraction went"
				<< std::endl;
		tempOut
				<< "seqFiles: Extracted files will be named [PrimerName][MIDNAME].fastq, if not multiplexed then just [PrimerName][MIDNAME].fastq"
				<< std::endl;
		tempOut
				<< "filteredOff: Reads that failed filtering will be found in this directory along with reads that don't match any supplied barcodes and/or primers"
				<< std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		exit(0);
	}
	processVerbose();
	processDebug();
	setOption(pars.noForwardPrimer, "-noForwardPrimer",
			"No Forward Primer Required");
	setOption(pars.forwardPrimerToUpperCase, "-fUpper",
			"Leave reverse primer upper case");
	setOption(pars.fPrimerErrors.distances_.query_.coverage_, "-fCoverage",
			"Amount of Forward Primer to find");
	setOption(pars.fPrimerErrors.hqMismatches_, "-fNumOfMismatches",
			"Number of Mismatches to allow in Forward Primer");
	setOption(pars.fPrimerErrors.oneBaseIndel_, "-fOneBaseIndels",
			"Number Of One base indels to allow in forward primer");
	setOption(pars.fPrimerErrors.twoBaseIndel_, "-fTwoBaseIndels",
			"Number Of Two base indels to allow in forward primer");
	//fPrimerErrors.largeBaseIndel_ = .99;
	setOption(pars.noReversePrimer, "-noReverse",
			"Don't look for reverse Primer");
	setOption(pars.reversePrimerToUpperCase, "-rUpper",
			"Leave reverse primer upper case");
	setOption(pars.rPrimerErrors.distances_.query_.coverage_, "-rPrimerCoverage",
			"Amount Of Reverse Primer Required");
	setOption(pars.rPrimerErrors.hqMismatches_, "-rNumOfMismatches",
			"Number of Mismatches to allow in Reverse Primer");
	setOption(pars.rPrimerErrors.oneBaseIndel_, "-rOneBaseIndels",
			"Number Of One base indels to allow in reverse primer");
	setOption(pars.rPrimerErrors.twoBaseIndel_, "-rTwoBaseIndels",
			"Number Of Two base indels to allow in reverse primer");
	//rPrimerErrors.largeBaseIndel_ = .99;

	setOption(pars.filterOffSmallReadCounts, "--filterOffSmallReadCounts",
			"Whether to Filter Off Extraction of Small size (size set by --smallExtractReadCount) ");
	setOption(pars.smallExtractReadCount, "--smallExtractReadCount",
			"Filter Off Extraction of This Size or Smaller");

	setOption(pars.sampleName, "--sampleName",
			"A name to append to the output files");
	processAlnInfoInput();

	setOption(pars.rename, "-rename", "Rename With Barcode Names");
	setOption(pars.barcodeErrors, "-barcodeErrors", "Errors Allowed in Barcode");
	setOption(pars.trimTcag, "-trimTcag", "trim TCAG(454 tag) if it is present ");
	if (setOption(pars.HMP, "-HMP", "HMP")) {
		setOption(pars.primerLen, "-len,-primerLen", "PrimerLen");
	}

	setOption(pars.mDetPars.barcodesBothEnds_, "-barcodeBothEnds",
			"Look for Barcodes in Both Primers");
	setOption(pars.midEndsRevComp, "--midEndsRevComp",
			"Barcodes on both ends are in the reverse complement of each other");
	setOption(pars.mDetPars.checkForShorten_, "-checkShortenBars",
			"Check for shorten Barcodes if the first base may have been trimmed off");

	setOption(pars.mDetPars.variableStop_, "-variableStart",
			"variableStart");
	if (setOption(pars.qualCheck, "-qualCheck", "Qual Check Level")) {
		pars.checkingQCheck = true;
	}
	if (setOption(pars.qualCheckCutOff, "-qualCheckCutOff",
			"Cut Off for fraction of bases above qual check of "
					+ estd::to_string(pars.qualCheck))) {
		pars.checkingQCheck = true;
	}
	setOption(pars.illumina, "--illumina", "If input reads are from Illumina");
	if(pars.illumina){
		pars.checkingQCheck = true;
	}
	processDefaultReader(true);
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);
	if (setOption(pars.idFilename, "-id", "The name of the ID file")) {
		if (!fexists(pars.idFilename)) {
			failed_ = true;
			warnings_.emplace_back("error in opening " + pars.idFilename);
		}
	}

	setOption(pars.idFileDelim, "-idFileDelim", "Id File Delim");
	// unknown primers and ids
	setOption(pars.multiplex, "-multiplex",
			"Indicates that the Reads are multiplex barcoded");
	setOption(pars.mDetPars.checkComplement_, "-checkComplement",
			"Check the Complement of the Seqs As Well");
	// setOption(within, "-within");

	setOption(pars.minLen, "-minlen", "Minimum read length");
	setOption(pars.maxLength, "-maxlen", "Maximum read length");

	setOption(pars.numberOfNs, "-numberofns", "Number Of Ns Cut Off");
	std::string qualWindow = "";
	if (setOption(qualWindow, "-qualWindow",
			"Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold")) {
		seqUtil::processQualityWindowString(qualWindow, pars.qualityWindowLength,
				pars.qualityWindowStep, pars.qualityWindowThres);
	} else {
		pars.qualityWindowLength = 50;
		pars.qualityWindowStep = 5;
		pars.qualityWindowThres = 25;
	}

	bool compareSeqSuccess = processSeq(pars.compareSeq, "-compareSeq",
			"Comparison Sequence For Contamination ");
	if(compareSeqSuccess){
		commands_.lookForOptionCaseInsen(pars.compareSeqFilename, "-compareSeq");
	}
	setOption(pars.contaminationKLen, "-contaminationKLen",
			"Contamination Kmer Length comparison");
	setOption(pars.kmerCutOff, "-kmerCutOff",
			"Kmer cut off for contamination check");
	if (setOption(pars.screenForPossibleContamination, "-contamination",
			"Screening For Contamination")) {
		if (!compareSeqSuccess) {
			warnings_.emplace_back(
					"Need to have -compareSeq if checking for contamination");
		}
	}
	setOption(pars.contaminationMutlipleCompare, "-contaminationMutlipleCompare",
			"Compare to all the sequences within the compare seq file if a file is given");
	setOption(pars.multipleTargets, "--multipleTargets",
			"Id file contains multiple targets");
	if (pars.multipleTargets) {
		if (pars.screenForPossibleContamination) {
			commands_.lookForOptionCaseInsen(pars.compareSeqFilename, "-compareSeq");
		}
	}
	setOption(pars.multipleLenCutOffFilename, "--lenCutOffs",
			"Length cut offs for when extracting multiple targets");
	setOption(pars.qualWindowTrim, "-qualWindowTrim", "Trim To Qual Window");
	setOption(pars.smallFragmentCutoff, "-smallFragmentCutOff",
			"Small Fragment Cut Off Size");

	setOption(pars.mothurExtract, "-mothurExtract",
			"Extract Data Files for mothur analysis");
	setOption(pars.pyroExtract, "-pyro",
				"Do Pyro extract for files necessary for AmpliconNoise");
	if (pars.mothurExtract || pars.pyroExtract) {
		if (!setOption(pars.maxFlowCutoff, "-maxFlows", "Max Flows", true)) {
			addWarning("Need to supply maxFlows if doing -pyro or -mothurExtract");
			failed_ = true;
		}
		if (SeqIOOptions::inFormats::SFFBIN != pars_.ioOptions_.inFormat_
				&& SeqIOOptions::inFormats::SFFTXT != pars_.ioOptions_.inFormat_) {
			addWarning(
					"Seq input format needs to be sff or sffbin if using -pyro or -mothurExtract");
			failed_ = true;
		}
	}

	pars.trimAtQual = setOption(pars.trimAtQualCutOff, "--trimAtQual",
			"Trim Reads at first occurrence of quality score");

	pars_.gapInfo_.gapOpen_ = 5;
	pars_.gapInfo_.gapExtend_ = 1;
	pars_.gap_ = "5,1";
	pars_.gapInfo_.gapRightOpen_ = 0;
	pars_.gapInfo_.gapRightExtend_ = 0;
	pars_.gapRight_ = "0,0";
	pars_.gapInfo_.gapLeftOpen_ = 0;
	pars_.gapInfo_.gapLeftExtend_ = 0;
	pars_.gapLeft_ = "0,0";
	processGap();
	finishSetUp(std::cout);
}

}  // namespace bibseq


