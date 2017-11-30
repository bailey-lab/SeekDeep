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
	setOption(pars.noForwardPrimer, "--noForwardPrimer",
			"No Forward Primer Required", false, "Primer");
	setOption(pars.forwardPrimerToUpperCase, "--fUpper",
			"Leave forward primer upper case", false, "Primer");
	setOption(pars.fPrimerErrors.distances_.query_.coverage_, "--fCoverage",
			"Amount of Forward Primer to find", false, "Primer");
	setOption(pars.fPrimerErrors.hqMismatches_, "--fNumOfMismatches",
			"Number of Mismatches to allow in Forward Primer", false, "Primer");
	setOption(pars.fPrimerErrors.oneBaseIndel_, "--fOneBaseIndels",
			"Number Of One base indels to allow in forward primer", false, "Primer");
	setOption(pars.fPrimerErrors.twoBaseIndel_, "--fTwoBaseIndels",
			"Number Of Two base indels to allow in forward primer", false, "Primer");
	//fPrimerErrors.largeBaseIndel_ = .99;
	setOption(pars.noReversePrimer, "--noReverse",
			"Don't look for reverse Primer", false, "Primer");
	setOption(pars.reversePrimerToUpperCase, "--rUpper",
			"Leave reverse primer upper case", false, "Primer");
	setOption(pars.rPrimerErrors.distances_.query_.coverage_, "--rPrimerCoverage",
			"Amount Of Reverse Primer Required", false, "Primer");
	setOption(pars.rPrimerErrors.hqMismatches_, "--rNumOfMismatches",
			"Number of Mismatches to allow in Reverse Primer", false, "Primer");
	setOption(pars.rPrimerErrors.oneBaseIndel_, "--rOneBaseIndels",
			"Number Of One base indels to allow in reverse primer", false, "Primer");
	setOption(pars.rPrimerErrors.twoBaseIndel_, "--rTwoBaseIndels",
			"Number Of Two base indels to allow in reverse primer", false, "Primer");


	setOption(pars.sampleName, "--sampleName",
			"A name to append to the output files",
			false, "Output Naming");
	processAlnInfoInput();

	setOption(pars.rename, "--rename", "Rename Sequences With Barcode Names",
			false, "Output Naming");
	setOption(pars.barcodeErrors, "--barcodeErrors", "Errors Allowed in Barcode", false, "Barcodes");


	setOption(pars.mDetPars.barcodesBothEnds_, "--barcodeBothEnds",
			"Look for Barcodes in Both Primers", false, "Barcodes");
	setOption(pars.midEndsRevComp, "--midEndsRevComp",
			"Barcodes on both ends are in the reverse complement of each other", false, "Barcodes");
	setOption(pars.mDetPars.checkForShorten_, "--checkShortenBars",
			"Check for shorten Barcodes if the first base may have been trimmed off", false, "Barcodes");

	setOption(pars.mDetPars.variableStop_, "--variableStart",
			"By default the primer or barcodes are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives",
			false, "Barcodes");

	processReadInNames(VecStr{"--fastq1", "--fastq1gz", "--fastq2", "--fastq2gz"},true);
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);
	if (setOption(pars.idFilename, "--id", "The name of the ID file", true, "ID File")) {
		if (!bfs::exists(pars.idFilename)) {
			failed_ = true;
			warnings_.emplace_back("Error the id file doesn't exist: " + pars.idFilename);
		}
	}

	setOption(pars.idFileDelim, "--idFileDelim", "Id File Delim", false, "ID File");
	// unknown primers and ids
	setOption(pars.multiplex, "--multiplex",
			"Indicates that the Reads are multiplex barcoded", false, "Barcodes");
	setOption(pars.mDetPars.checkComplement_, "--checkComplement",
			"Check the Complement of the Seqs As Well", false, "Complement");

	setOption(pars.smallFragmentCutoff, "--smallFragmentCutOff",
			"Remove sequences smaller than this length", false, "Pre Processing");


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
	description_ = "Extract sequences from various sequences input types (fastq,fasta,sff,etc.) with primers and barcodes plus some filtering";
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fasta reads.fasta --id idFile.tab.txt --dout outPutDir");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq reads.fastq --id idFile.tab.txt --dout outPutDir --illumina");


	processVerbose();
	processDebug();
	setOption(pars.noForwardPrimer, "--noForwardPrimer",
			"No Forward Primer Required", false, "Primer");
	setOption(pars.forwardPrimerToUpperCase, "--fUpper",
			"Leave forward primer upper case", false, "Primer");
	setOption(pars.fPrimerErrors.distances_.query_.coverage_, "--fCoverage",
			"Amount of Forward Primer to find", false, "Primer");
	setOption(pars.fPrimerErrors.hqMismatches_, "--fNumOfMismatches",
			"Number of Mismatches to allow in Forward Primer", false, "Primer");
	setOption(pars.fPrimerErrors.oneBaseIndel_, "--fOneBaseIndels",
			"Number Of One base indels to allow in forward primer", false, "Primer");
	setOption(pars.fPrimerErrors.twoBaseIndel_, "--fTwoBaseIndels",
			"Number Of Two base indels to allow in forward primer", false, "Primer");
	//fPrimerErrors.largeBaseIndel_ = .99;
	setOption(pars.noReversePrimer, "--noReverse",
			"Don't look for reverse Primer", false, "Primer");
	setOption(pars.reversePrimerToUpperCase, "--rUpper",
			"Leave reverse primer upper case", false, "Primer");
	setOption(pars.rPrimerErrors.distances_.query_.coverage_, "--rPrimerCoverage",
			"Amount Of Reverse Primer Required", false, "Primer");
	setOption(pars.rPrimerErrors.hqMismatches_, "--rNumOfMismatches",
			"Number of Mismatches to allow in Reverse Primer", false, "Primer");
	setOption(pars.rPrimerErrors.oneBaseIndel_, "--rOneBaseIndels",
			"Number Of One base indels to allow in reverse primer", false, "Primer");
	setOption(pars.rPrimerErrors.twoBaseIndel_, "--rTwoBaseIndels",
			"Number Of Two base indels to allow in reverse primer", false, "Primer");
	//rPrimerErrors.largeBaseIndel_ = .99;

	setOption(pars.filterOffSmallReadCounts, "--filterOffSmallReadCounts",
			"Whether to Filter Off Extraction of Small size (size set by --smallExtractReadCount)",
			false, "Post Extraction Filtering");
	setOption(pars.smallExtractReadCount, "--smallExtractReadCount",
			"Filter Off Extraction of This Size or Smaller",
			false, "Post Extraction Filtering");

	setOption(pars.sampleName, "--sampleName",
			"A name to append to the output files",
			false, "Output Naming");
	processAlnInfoInput();

	setOption(pars.rename, "--rename", "Rename Sequences With Barcode Names",
			false, "Output Naming");
	setOption(pars.barcodeErrors, "--barcodeErrors", "Errors Allowed in Barcode", false, "Barcodes");
	setOption(pars.trimTcag, "--trimTcag", "trim TCAG(454 tag) if it is present",
			false, "Preprocessing");
	if (setOption(pars.HMP, "--HMP", "HMP")) {
		setOption(pars.primerLen, "--len,--primerLen", "PrimerLen");
	}

	setOption(pars.mDetPars.barcodesBothEnds_, "--barcodeBothEnds",
			"Look for Barcodes in Both Primers", false, "Barcodes");
	setOption(pars.midEndsRevComp, "--midEndsRevComp",
			"Barcodes on both ends are in the reverse complement of each other", false, "Barcodes");
	setOption(pars.mDetPars.checkForShorten_, "--checkShortenBars",
			"Check for shorten Barcodes if the first base may have been trimmed off", false, "Barcodes");

	setOption(pars.mDetPars.variableStop_, "--variableStart",
			"By default the primer or barcodes are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives",
			false, "Barcodes");
	if (setOption(pars.qualCheck, "--qualCheck", "Check the fraction of qualities above a threshold set by --qualCheckCutOff",
			false, "Filtering")) {
		pars.checkingQCheck = true;
	}
	if (setOption(pars.qualCheckCutOff, "--qualCheckCutOff",
			"Cut Off for fraction of bases above qual check of "
					+ estd::to_string(pars.qualCheck), false, "Filtering")) {
		pars.checkingQCheck = true;
	}
	setOption(pars.illumina, "--illumina", "If input reads are from Illumina",
			false, "Technology");
	if(pars.illumina){
		pars.checkingQCheck = true;
	}
	processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz"},true);
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);
	if (setOption(pars.idFilename, "--id", "The name of the ID file", true, "ID File")) {
		if (!bfs::exists(pars.idFilename)) {
			failed_ = true;
			warnings_.emplace_back("Error the id file doesn't exist: " + pars.idFilename);
		}
	}

	setOption(pars.idFileDelim, "--idFileDelim", "Id File Delim", false, "ID File");
	// unknown primers and ids
	setOption(pars.multiplex, "--multiplex",
			"Indicates that the Reads are multiplex barcoded", false, "Barcodes");
	setOption(pars.mDetPars.checkComplement_, "--checkComplement",
			"Check the Complement of the Seqs As Well", false, "Complement");
	// setOption(within, "-within");
	setOption(pars.minLen, "--minlen", "Minimum read length", false, "Filtering");
	setOption(pars.maxLength, "--maxlen", "Maximum read length", false, "Filtering");
	setOption(pars.numberOfNs, "--numberofns", "Number Of Ns Cut Off", false, "Filtering");
	std::string qualWindow = "";
	if (setOption(qualWindow, "--qualWindow",
			"Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold", false, "Filtering")) {
		seqUtil::processQualityWindowString(qualWindow, pars.qualityWindowLength,
				pars.qualityWindowStep, pars.qualityWindowThres);
	} else {
		pars.qualityWindowLength = 50;
		pars.qualityWindowStep = 5;
		pars.qualityWindowThres = 25;
	}

	bool compareSeqSuccess = processSeq(pars.compareSeq, "--compareSeq",
			"Comparison Sequence For Contamination", false, "Contamination Filtering");
	if(compareSeqSuccess){
		commands_.lookForOptionDashCaseInsen(pars.compareSeqFilename, "--compareSeq");
	}
	setOption(pars.contaminationKLen, "--contaminationKLen",
			"Contamination Kmer Length comparison", false, "Contamination Filtering");
	setOption(pars.kmerCutOff, "--kmerCutOff",
			"Kmer cut off for contamination check", false, "Contamination Filtering");
	if (setOption(pars.screenForPossibleContamination, "--contamination",
			"Screening For Contamination", false, "Contamination Filtering")) {
		if (!compareSeqSuccess) {
			warnings_.emplace_back(
					"Need to have --compareSeq if checking for contamination");
		}
	}
	setOption(pars.contaminationMutlipleCompare, "--contaminationMutlipleCompare",
			"Compare to all the sequences within the compare seq file if a file is given", false, "Contamination Filtering");
	setOption(pars.multipleTargets, "--multipleTargets",
			"Id file contains multiple targets", false, "Primers");
	if (pars.multipleTargets) {
		if (pars.screenForPossibleContamination) {
			commands_.lookForOptionDashCaseInsen(pars.compareSeqFilename, "--compareSeq");
		}
	}
	setOption(pars.multipleLenCutOffFilename, "--lenCutOffs",
			"Length cut offs for when extracting multiple targets", false, "Filtering");
	setOption(pars.qualWindowTrim, "--qualWindowTrim", "Trim To Qual Window", false, "Post Processing");
	setOption(pars.smallFragmentCutoff, "--smallFragmentCutOff",
			"Remove sequences smaller than this length", false, "Pre Processing");

	setOption(pars.mothurExtract, "--mothurExtract",
			"Extract Data Files for mothur analysis");
	setOption(pars.pyroExtract, "--pyro",
				"Do Pyro extract for files necessary for AmpliconNoise");
	if (pars.mothurExtract || pars.pyroExtract) {
		if (!setOption(pars.maxFlowCutoff, "--maxFlows", "Max Flows", true)) {
			addWarning("Need to supply maxFlows if doing --pyro or --mothurExtract");
			failed_ = true;
		}
		if (SeqIOOptions::inFormats::SFFBIN != pars_.ioOptions_.inFormat_
				&& SeqIOOptions::inFormats::SFFTXT != pars_.ioOptions_.inFormat_) {
			addWarning(
					"Seq input format needs to be sff or sffbin if using --pyro or --mothurExtract");
			failed_ = true;
		}
	}

	pars.trimAtQual = setOption(pars.trimAtQualCutOff, "--trimAtQual",
			"Trim Reads at first occurrence of quality score", false, "Post Processing");

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

}  // namespace bibseq


