/*
 * extractorSetUp.cpp
 *
 *  Created on: Dec 17, 2016
 *      Author: nick
 */

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
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {


void SeekDeepSetUp::setUpExtractorPairedEnd(ExtractorPairedEndPars & pars) {
	if (needsHelp()) {
		commands_.arguments_["-h"] = "";
	}
	//ioOptions_.lowerCaseBases = "remove";
	description_ = "Extract sequences from various sequences input types (fastq,fasta,sff,etc.) with primers and barcodes plus some filtering";
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fasta reads.fasta --id idFile.tab.txt --dout outPutDir");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq reads.fastq --id idFile.tab.txt --dout outPutDir --illumina");


	processVerbose();
	processDebug();
	//core
	pars.corePars_.setCorePars(*this);
	processAlnInfoInput();
	processReadInNames(VecStr{"--fastq1", "--fastq1gz", "--fastq2", "--fastq2gz"},true);
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);

	//paired end specific stuff
	pars.pairProcessorParams_.verbose_ = pars_.verbose_;
	setOption(pars.corePars_.primIdsPars.noOverlapProcessForNoOverlapStatusTargets_, "--noOverlapProcessForNoOverlapStatusTargets_", "noOverlapProcessForNoOverlapStatusTargets",
			false, "Post-Processing-PairProcessing");
	setOption(pars.pairProcessorParams_.errorAllowed_, "--overLapErrorAllowed",
			"The amount of error to allow in the overlap processing", false, "Post-Processing-PairProcessing");
	setOption(pars.pairProcessorParams_.hardMismatchCutOff_, "--hardMismatchCutOff",
			 "A hard cut off for number of mismatches between overlapping sequences in pair processing", false, "Post-Processing-PairProcessing");
	setOption(pars.pairProcessorParams_.minOverlap_, "--minOverlap",
			"The minimal amount of over lap in pair processing", false, "Post-Processing-PairProcessing");
	setOption(pars.pairProcessorParams_.writeOverHangs_, "--writeOverHangs",
			"Write out the overhang for sequences that have read through", false, "Post-Processing-PairProcessing");
	setOption(pars.r1Trim_, "--r1Trim",
			"Remove this many sequences off of the end of r1 reads", false, "Post Processing");
	setOption(pars.r2Trim_, "--r2Trim",
			"Remove this many sequences off of the end of r2 reads", false, "Post Processing");
	setOption(pars.corePars_.primIdsPars.overlapStatusFnp_, "--overlapStatusFnp",
			"A file with two columns, target,status; status column should contain 1 of 3 values (capitalization doesn't matter): r1BegOverR2End,r1EndOverR2Beg,NoOverlap. r1BegOverR2End=target size < read length (causes read through),r1EndOverR2Beg= target size > read length less than 2 x read length, NoOverlap=target size > 2 x read length", true, "Post Processing");



	pars_.gapInfo_.gapOpen_ = 5;
	pars_.gapInfo_.gapExtend_ = 1;
	pars_.gap_ = "5,1";
	pars_.gapInfo_.gapRightQueryOpen_ = 0;
	pars_.gapInfo_.gapRightQueryExtend_ = 0;
	pars_.gapInfo_.gapRightRefOpen_ = 0;
	pars_.gapInfo_.gapRightRefExtend_ = 0;
	pars_.gapRight_ = "0,0";
	pars_.gapInfo_.gapLeftQueryOpen_ = 0;
	pars_.gapInfo_.gapLeftQueryExtend_ = 0;
	pars_.gapInfo_.gapLeftRefOpen_ = 0;
	pars_.gapInfo_.gapLeftRefExtend_ = 0;
	pars_.gapLeft_ = "0,0";
	processGap();
	if (needsHelp()) {
		printFlags(std::cout);
		std::cout << "The id file should be tab delimited and contains the "
				"primer and reverse primer and MIDs if the data is multiplex, "
				"example below" << std::endl;
		std::cout << "\ttarget\tforwardPrimer\treversePrimer" << std::endl;
		std::cout << "\tPFMSP1\tAACTAGAAGCTTTAGAAGATGCA\tACATATGATTGGTTAAATCAAAG"
				<< std::endl;
		std::cout << "\tPFMSP2\tAGATGCAGCTTTAACTAGAAGAA\tTTAAAACATATGATTGGTCAAAG"
				<< std::endl;
		std::cout << "\tid	barcode" << std::endl;
		std::cout << "\tMID01\t	ACGAGTGCGT" << std::endl;
		std::cout << "\tMID02	\tACGCTCGACA" << std::endl;
		std::cout << bib::bashCT::bold << "Output Files:" << bib::bashCT::reset
				<< std::endl;
		std::cout
				<< "extractionProfile.tab.txt: This breaks down the filtering per final extraction sequence file"
				<< std::endl;
		std::cout
				<< "extractionStats.tab.txt: This has info how the whole extraction went"
				<< std::endl;
		std::cout
				<< "seqFiles: Extracted files will be named [PrimerName][MIDNAME].fast(a/q), if not multiplexed then just [PrimerName].fast(a/q)"
				<< std::endl;
		std::cout
				<< "filteredOff: Reads that failed filtering will be found in this directory along with reads that don't match any supplied barcodes and/or primers"
				<< std::endl;
		exit(1);
	}
	finishSetUp(std::cout);
}

void SeekDeepSetUp::setUpExtractor(extractorPars & pars) {
	if (needsHelp()) {
		commands_.arguments_["-h"] = "";
	}
	//ioOptions_.lowerCaseBases = "remove";
	description_ = "Extract sequences from various sequences input types (fastq,fasta,etc.) with primers and barcodes plus some filtering";
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fasta reads.fasta --id idFile.tab.txt --dout outPutDir");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq reads.fastq --id idFile.tab.txt --dout outPutDir --illumina");


	processVerbose();
	processDebug();
	setOption(pars.illumina, "--illumina", "If input reads are from Illumina",false, "Technology");
	if(pars.illumina){
		pars.corePars_.qPars_.checkingQFrac_ = true;

		pars.corePars_.pDetPars.allowable_.hqMismatches_ = 2;
		pars.corePars_.pDetPars.allowable_.lqMismatches_ = 5;
		pars.corePars_.pDetPars.allowable_.distances_.query_.coverage_ = 1;
		pars.corePars_.pDetPars.allowable_.largeBaseIndel_ = 0.99;
		pars.corePars_.pDetPars.allowable_.oneBaseIndel_ = 0.5;
		pars.corePars_.pDetPars.allowable_.twoBaseIndel_ = 0.5;

	}

	//core
	pars.corePars_.setCorePars(*this);
	setOption(pars.corePars_.noReversePrimer_, "--noReversePrimer", "Don't check for reverse primer");

	//copy over the regular determined primers to the back end primers
	pars.corePars_.backEndpDetPars.primerWithin_ = pars.corePars_.pDetPars.primerWithin_;
	pars.corePars_.backEndpDetPars.primerToLowerCase_ = pars.corePars_.pDetPars.primerToLowerCase_ ;

	pars.corePars_.backEndpDetPars.allowable_.hqMismatches_ = pars.corePars_.pDetPars.allowable_.hqMismatches_;
	pars.corePars_.backEndpDetPars.allowable_.lqMismatches_ = pars.corePars_.pDetPars.allowable_.lqMismatches_;
	pars.corePars_.backEndpDetPars.allowable_.distances_.query_.coverage_ = pars.corePars_.pDetPars.allowable_.distances_.query_.coverage_;
	pars.corePars_.backEndpDetPars.allowable_.largeBaseIndel_ = pars.corePars_.pDetPars.allowable_.largeBaseIndel_;
	pars.corePars_.backEndpDetPars.allowable_.oneBaseIndel_ = pars.corePars_.pDetPars.allowable_.oneBaseIndel_;
	pars.corePars_.backEndpDetPars.allowable_.twoBaseIndel_ = pars.corePars_.pDetPars.allowable_.twoBaseIndel_;
	setOption(pars.corePars_.pDetPars.primerWithin_, "--frontEndPrimerWithinStart",
			"By default the primer or barcodes are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives, this is for the front primer",
			false, "Primers");
	setOption(pars.corePars_.backEndpDetPars.primerWithin_, "--backEndPrimerWithinStart",
			"By default the primer or barcodes are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives, this is for the back primer",
			false, "Primers");



	//fPrimerErrors.largeBaseIndel_ = .99;
//	setOption(pars.noReversePrimer, "--noReverse",
//			"Don't look for reverse Primer", false, "Primer");
//	setOption(pars.reversePrimerToUpperCase, "--rUpper",
//			"Leave reverse primer upper case", false, "Primer");
//	setOption(pars.rPrimerErrors.distances_.query_.coverage_, "--rPrimerCoverage",
//			"Amount Of Reverse Primer Required", false, "Primer");
//	setOption(pars.rPrimerErrors.hqMismatches_, "--rNumOfMismatches",
//			"Number of Mismatches to allow in Reverse Primer", false, "Primer");
//	setOption(pars.rPrimerErrors.oneBaseIndel_, "--rOneBaseIndels",
//			"Number Of One base indels to allow in reverse primer", false, "Primer");
//	setOption(pars.rPrimerErrors.twoBaseIndel_, "--rTwoBaseIndels",
//			"Number Of Two base indels to allow in reverse primer", false, "Primer");
	//rPrimerErrors.largeBaseIndel_ = .99;

	setOption(pars.filterOffSmallReadCounts, "--filterOffSmallReadCounts", "Whether to Filter Off Extraction of Small size (size set by --smallExtractReadCount)", false, "Post Extraction Filtering");
	setOption(pars.smallExtractReadCount, "--smallExtractReadCount", "Filter Off Extraction of This Size or Smaller", false, "Post Extraction Filtering");

	processAlnInfoInput();


	processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz"},true);
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);

	// setOption(within, "-within");
	setOption(pars.minLen, "--minlen", "Minimum read length", false, "Filtering");
	setOption(pars.maxLength, "--maxlen", "Maximum read length", false, "Filtering");
	std::string qualWindow = "";
	if (setOption(qualWindow, "--qualWindow",
			"Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold", false, "Filtering")) {
		seqUtil::processQualityWindowString(qualWindow, pars.corePars_.qPars_.qualityWindowLength_,
				pars.corePars_.qPars_.qualityWindowStep_, pars.corePars_.qPars_.qualityWindowThres_);
	} else {
		pars.corePars_.qPars_.qualityWindowLength_ = 50;
		pars.corePars_.qPars_.qualityWindowStep_ = 5;
		pars.corePars_.qPars_.qualityWindowThres_ = 25;
	}
	setOption(pars.qualWindowTrim, "--qualWindowTrim", "Trim To Qual Window", false, "Post Processing");
	pars.trimAtQual = setOption(pars.trimAtQualCutOff, "--trimAtQual", "Trim Reads at first occurrence of quality score", false, "Post Processing");

	setOption(pars.trimToMaxLength, "--trimToMaxLength", "Trim sequences to max expected length to improve primer determination for mixed target datasets");
	pars_.gapInfo_.gapOpen_ = 5;
	pars_.gapInfo_.gapExtend_ = 1;
	pars_.gap_ = "5,1";
	pars_.gapInfo_.gapRightQueryOpen_ = 0;
	pars_.gapInfo_.gapRightQueryExtend_ = 0;
	pars_.gapInfo_.gapRightRefOpen_ = 0;
	pars_.gapInfo_.gapRightRefExtend_ = 0;
	pars_.gapRight_ = "0,0";
	pars_.gapInfo_.gapLeftQueryOpen_ = 0;
	pars_.gapInfo_.gapLeftQueryExtend_ = 0;
	pars_.gapInfo_.gapLeftRefOpen_ = 0;
	pars_.gapInfo_.gapLeftRefExtend_ = 0;
	pars_.gapLeft_ = "0,0";
	processGap();

	if (pars_.verbose_) {
		if (pars.corePars_.qPars_.checkingQFrac_) {
			std::cout << "Quality Check: " << pars.corePars_.qPars_.qualCheck_ << std::endl;
			std::cout << "Quality Check Cut Off: " << pars.corePars_.qPars_.qualCheckCutOff_<< std::endl;
		} else {
			std::cout << "Quality Window Length: " << pars.corePars_.qPars_.qualityWindowLength_
					<< std::endl;
			std::cout << "Quality Window Step: " << pars.corePars_.qPars_.qualityWindowStep_
					<< std::endl;
			std::cout << "Quality Window Threshold: " << pars.corePars_.qPars_.qualityWindowThres_
					<< std::endl;
		}
	}
	if (needsHelp()) {
		printFlags(std::cout);
		std::cout << "The id file should be tab delimited and contains the "
				"primer and reverse primer and MIDs if the data is multiplex, "
				"example below" << std::endl;
		std::cout << "\ttarget\tforwardPrimer\treversePrimer" << std::endl;
		std::cout << "\tPFMSP1\tAACTAGAAGCTTTAGAAGATGCA\tACATATGATTGGTTAAATCAAAG"
				<< std::endl;
		std::cout << "\tPFMSP2\tAGATGCAGCTTTAACTAGAAGAA\tTTAAAACATATGATTGGTCAAAG"
				<< std::endl;
		std::cout << "\tid	barcode" << std::endl;
		std::cout << "\tMID01\t	ACGAGTGCGT" << std::endl;
		std::cout << "\tMID02	\tACGCTCGACA" << std::endl;
		std::cout << bib::bashCT::bold << "Output Files:" << bib::bashCT::reset
				<< std::endl;
		std::cout
				<< "extractionProfile.tab.txt: This breaks down the filtering per final extraction sequence file"
				<< std::endl;
		std::cout
				<< "extractionStats.tab.txt: This has info how the whole extraction went"
				<< std::endl;
		std::cout
				<< "seqFiles: Extracted files will be named [PrimerName][MIDNAME].fast(a/q), if not multiplexed then just [PrimerName].fast(a/q)"
				<< std::endl;
		std::cout
				<< "filteredOff: Reads that failed filtering will be found in this directory along with reads that don't match any supplied barcodes and/or primers"
				<< std::endl;
		exit(1);
	}
	finishSetUp(std::cout);
}

}  // namespace bibseq


