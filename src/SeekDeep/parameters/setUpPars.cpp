/*
 * setUpPars.cpp
 *
 *  Created on: Nov 25, 2016
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


#include "setUpPars.hpp"

namespace bibseq {

void CoreExtractorPars::setCorePars(seqSetUp & setUp){
	bool primerToUpperCase = false;
	setUp.setOption(primerToUpperCase, "--primerUpper",
			"Leave primers in upper case", false, "Primer");
	pDetPars.primerToLowerCase_ = !primerToUpperCase;
	setUp.setOption(pDetPars.allowable_.distances_.query_.coverage_, "--primerCoverage",
			"Amount of primers found", false, "Primer");
	setUp.setOption(pDetPars.allowable_.hqMismatches_, "--primerNumOfMismatches",
			"Number of Mismatches to allow in primers", false, "Primer");
	setUp.setOption(pDetPars.allowable_.oneBaseIndel_, "--primerOneBaseIndels",
			"Number Of One base indels to allow in primers", false, "Primer");
	setUp.setOption(pDetPars.allowable_.twoBaseIndel_, "--primerTwoBaseIndels",
			"Number Of Two base indels to allow in primers", false, "Primer");
	setUp.setOption(sampleName, "--sampleName",
			"A name to append to the output files",
			false, "Output Naming");
	setUp.setOption(rename, "--rename", "Rename Sequences With Barcode/Primer Names",
			false, "Output Naming");
	setUp.setOption(primIdsPars.barcodeErrors_ , "--barcodeErrors", "Errors Allowed in Barcode", false, "Barcodes");
	setUp.setOption(mDetPars.barcodesBothEnds_, "--barcodeBothEnds",
			"Look for Barcodes in Both Primers", false, "Barcodes");
	setUp.setOption(primIdsPars.midEndsRevComp_, "--midEndsRevComp",
			"Barcodes on both ends are in the reverse complement of each other", false, "Barcodes");
	setUp.setOption(mDetPars.checkForShorten_, "--checkShortenBars",
			"Check for shorten Barcodes if the first base may have been trimmed off", false, "Barcodes");

	setUp.setOption(mDetPars.variableStop_, "--midWithinStart",
			"By default the primer or barcodes are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives",
			false, "Barcodes");

	setUp.setOption(pDetPars.primerWithin_, "--primerWithinStart",
			"By default the primer or barcodes are searched at the very beginning of seq, use this flag to extended the search, should be kept low to cut down on false positives",
			false, "Primers");
	if (setUp.setOption(primIdsPars.idFile_, "--id", "The name of the ID file", true, "ID File")) {
		if (!bfs::exists(primIdsPars.idFile_)) {
			setUp.failed_ = true;

			setUp.addWarning("Error the id file doesn't exist: " + primIdsPars.idFile_.string());
		}
	}
	setUp.setOption(idFileDelim, "--idFileDelim", "Id File Delim", false, "ID File");
	setUp.setOption(mDetPars.checkComplement_, "--checkRevComplementForMids", "Check the Reverse Complement of the Seqs As Well For MIDs", false, "Complement");
	setUp.setOption(pDetPars.checkComplement_, "--checkRevComplementForPrimers", "Check the Reverse Complement of the Seqs As Well For Primers", false, "Complement");

	setUp.setOption(smallFragmentCutoff, "--smallFragmentCutOff",
			"Remove sequences smaller than this length", false, "Pre Processing");
	setUp.setOption(primIdsPars.lenCutOffFilename_, "--lenCutOffs",
			"A file with at least three columns, target,minlen,maxlen the target column should match up with the first column in the id file", false, "Post Processing");
	setUp.setOption(primIdsPars.comparisonSeqFnp_, "--compareSeq",
			"A fasta file or a directory to fasta files, with references to check against, if file record name need to match target name, if directory file name should be TARGET.fasta", false, "Post Processing");
	if(setUp.setOption(qPars_.qualCheck_, "--qualCheckLevel",
			"Bin qualities at this quality to do filtering on fraction above this", false, "Post Processing")){
		qPars_.checkingQFrac_ = true;
	}
	if(setUp.setOption(qPars_.qualCheckCutOff_, "--qualCheckCutOff",
			"The fractions of bases that have to be above the qualCheckLevel to be kept", false, "Post Processing")){
		qPars_.checkingQFrac_ = true;
	}
	setUp.setOption(numberOfNs, "--numberOfNs", "Number Of Ns Cut Off", false, "Filtering");

	setUp.setOption(noPrimers_, "--noPrimers", "If no primers is set, only one line can be found under the targets/gene headers, the sequences in the forward/reverse primers will be ignored", false, "Primer");

	setUp.setOption(keepUnfilteredReads, "--keepUnfilteredReads", "Keep the unfiltered reads for debugging purposes", false);

}


extractorPars::extractorPars(){
//  rPrimerErrors.hqMismatches_ = 4;
//  rPrimerErrors.distances_.query_.coverage_ = .50;
//  rPrimerErrors.largeBaseIndel_ = .99;
//  rPrimerErrors.oneBaseIndel_ = 2;
//  rPrimerErrors.twoBaseIndel_ = 1;

  corePars_.pDetPars.allowable_.hqMismatches_ = 2;
  corePars_.pDetPars.allowable_.distances_.query_.coverage_ = .90;
  corePars_.pDetPars.allowable_.largeBaseIndel_ = .99;
  corePars_.pDetPars.allowable_.oneBaseIndel_ = 2;
  corePars_.pDetPars.allowable_.twoBaseIndel_ = 1;
}


ExtractorPairedEndPars::ExtractorPairedEndPars(){

	corePars_.pDetPars.allowable_.hqMismatches_ = 2;
	corePars_.pDetPars.allowable_.lqMismatches_ = 5; /**@todo incorporate this*/
	corePars_.pDetPars.allowable_.distances_.query_.coverage_ = 1;
	corePars_.pDetPars.allowable_.largeBaseIndel_ = 0.99;
	corePars_.pDetPars.allowable_.oneBaseIndel_ = 0.5;
	corePars_.pDetPars.allowable_.twoBaseIndel_ = 0.5;

	pairProcessorParams_.r1Trim_ = 1;
	pairProcessorParams_.r2Trim_ = 1;

}

}  // namespace bibseq
