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
#include "SeekDeepSetUp.hpp"
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

void SeekDeepSetUp::setUpClusterDown(clusterDownPars & pars) {
	// input file info
	pars_.ioOptions_.out_.outFilename_ = "output";
	pars_.ioOptions_.lowerCaseBases_ = "remove";
	if (needsHelp()) {
		std::stringstream tempOut;
		tempOut << commands_.subProgram_ << std::endl;
		tempOut << "Iteratively clusters reads by using the allowable errors given "
				"in a parameters file" << std::endl;
		tempOut << "Commands, order not necessary and case insensitive"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Required commands" << bib::bashCT::reset
				<< std::endl;
		// std::cout << cleanOut(tempOut.str(), width_, indent_);
		// tempOut.str(std::string());
		printInputUsage(tempOut);
		tempOut << "-par [option]: parameter file" << std::endl << std::endl;

		tempOut
				<< "Parameter file is set up as where each line is one iteration and that line indicates "
						"which errors to allow on that iteration " << std::endl;
		tempOut << "The errors allowed are separated by a colon : " << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		std::cout
				<< "Header should be as follows, followed by the iteration information "
				<< std::endl;
		std::cout
				<< "stopCheck:smallCutoff:1baseIndel:2baseIndel:>2baseIndel:HQMismatches:LQMismatches:LKMismatches"
				<< std::endl;
		std::cout << "100:3:1:0:0:0:0:0" << std::endl;
		std::cout << "100:3:2:0:0:0:0:0" << std::endl;
		std::cout << "100:3:3:0:0:0:1:0" << std::endl;
		std::cout << "100:3:4:0:0:0:2:0" << std::endl;
		std::cout << "100:0:1:0:0:0:0:0" << std::endl;
		std::cout << "100:0:2:0:0:0:0:0" << std::endl;
		std::cout << "100:0:3:0:0:0:1:0" << std::endl;
		std::cout << "100:0:4:0:0:0:2:0" << std::endl;
		tempOut
				<< "So here there would be eight iteration, in the first iteration reads"
						" will only cluster if they differ only by 1 1base indel"
				<< std::endl;
		tempOut
				<< "The first two columns control how many reads to compare."
						"  The first column will control the number of reads check, so in this"
						" instance only the first 100 clusters will be checked.  The second number"
						" controls size of clusters to compare against, in the first iteration a 3 "
						"means clusters smaller than 3 will be allowed to have clusters collapsed "
						"into them though they can still collapse into larger clusters"
				<< std::endl;
		tempOut
				<< "Alternatively if the -onPerId is used, this parameter flag can contain percent idendities instead"
				<< std::endl;
		tempOut << "example: " << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		tempOut << "100:3:.99" << std::endl;
		tempOut << "100:3:.99" << std::endl;
		tempOut << "100:3:.98" << std::endl;
		tempOut << "100:3:.98" << std::endl;
		tempOut
				<< "Here, four iteration where clusters were clustered at first .99 percent identity and then .98"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Optional commands" << bib::bashCT::reset
				<< std::endl;
		tempOut << "-dout [option]: Name of an output directory, will default "
				"to the stub name plus the date" << std::endl;
		tempOut << "-additionaOut [option]: Name of a tab delimited file with "
				"two columns, column one is the name of the ID and column two "
				"is the file path location of additional output directory"
				", such a file can be made by makeSampleDirectories" << std::endl;
		printAdditionalClusteringUsage(tempOut);
		tempOut << bib::bashCT::bold << "Alignment comparison options"
				<< bib::bashCT::reset << std::endl;
		printQualThresUsage(tempOut);
		printAlignmentUsage(tempOut);
		printKmerProfilingUsage(tempOut);
		printAlnInfoDirUsage(tempOut);
		printAdditionaInputUsage(tempOut, pars_.ioOptions_.lowerCaseBases_);
		tempOut
				<< "-qualRep [option] : Sets the quality for the identical clusters, "
						"Options are median (default), average, bestQual, or worst"
				<< std::endl;
		tempOut << "-qualTrim [option]: Will trim off any bases below this "
				"quality, mainly used to trim bases with a quality of 1 (80% "
				"chance of error) by setting this to 2" << std::endl;
		tempOut << "-adjustHomopolyerRuns : Will take average quality across "
				"a homopolymer run and set all the quality to this average, "
				"mainly used with Ion Torrent reads due to the very low "
				"quality of the last base in a long homopolymer run" << std::endl;
		tempOut << bib::bashCT::bold << "Additional Processing options"
				<< bib::bashCT::reset << std::endl;
		tempOut << "-markChimeras : Have the program mark possible chimeras "
				"before outputting" << std::endl;
		tempOut << "-parFreqs [option] : The minimum frequency multiplier "
				"reads must have in order to be considered to be a possible "
				"parent of a chimera or for sequences to be collapsed on gaps "
				"in tandem repeats, defaults to 2" << std::endl;
		tempOut << bib::bashCT::bold << "Technology Specific flags"
				<< bib::bashCT::reset << std::endl;
		tempOut
				<< "-useCompPerCutOff: This turns on filtering off clustering that have reads that only come from one direction, only should be used if reads were extracted in both direction like in ion torrent"
				<< std::endl;
		tempOut
				<< "-ionTorrent : This turns on several of the previously mentioned flags as they are beneficial for ion torrent data, turns on -qualTrim,-adjustHomopolyerRuns, and -useCompPerCutOff"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Percent Identity Falgs"
				<< bib::bashCT::reset << std::endl;
		tempOut
				<< "-onPerId: cluster on percent identity rather than specific errors"
				<< std::endl;

		printReferenceComparisonUsage(tempOut);
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		std::cout << "examples" << std::endl;
		std::cout << "\tSeekDeep " << commands_.subProgram_
				<< " -fasta MID1.fasta -qual "
						"MID1.fasta.qual -par par" << std::endl;
		std::cout << "\tSeekDeep " << commands_.subProgram_
				<< " -fastq MID1.fastq -par par -local" << std::endl;
		std::cout << "\tSeekDeep " << commands_.subProgram_
				<< " -stub MID1 -par par -local "
						"-qualThres 25,20" << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		tempOut << bib::bashCT::bold << "Output Files:" << bib::bashCT::reset
				<< std::endl;
		tempOut
				<< "output.fastq: This is the final clusters with their consensus sequence, the sequences are named so that the suffif _t[NUM] where NUM is the number of reads that fell into that cluster"
				<< std::endl;
		tempOut
				<< "outputInfo.tab.txt: This contains cluster number info for each cluster"
				<< std::endl;
		tempOut
				<< "clusters: A seq file is available for each final cluster that contains all the sequences that clustered into this cluster"
				<< std::endl;
		tempOut
				<< "internalSnpInfo: A table of internal snps frequencies is available for each cluster in case overcollapsing is suspected"
				<< std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		exit(0);
	}
	processDefaultReader(true);

	setOption(pars.binParameters, "-binPar",
			"Parameters Filename for when reads are binned");
	setOption(pars.ionTorrent, "-ionTorrent",
			"Flag to indicate reads are ion torrent and therefore turns on -useCompPerCutOff,-adjustHomopolyerRuns, and -qualTrim");
	setOption(pars.illumina, "-illumina",
				"Flag to indicate reads are illumina");
	bool needsParFlag = !pars.illumina && !pars.ionTorrent;
	if (pars.ionTorrent) {
		pars_.colOpts_.iTOpts_.removeLowQualityBases_ = true;
		pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_ = true;
		//pars.useCompPerCutOff = true;
		auto temp = CollapseIterations::gen454ItDefaultPars(100);
		pars.iteratorMap = temp;
	}

	if(pars.illumina){
		pars.iteratorMap = CollapseIterations::genIlluminaDefaultPars(100);
		pars_.colOpts_.iTOpts_.weighHomopolyer_ = false;
	}

	setOption(pars.parameters, "-par", "ParametersFilename", needsParFlag);
	setOption(pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_, "-adjustHomopolyerRuns",
				"Adjust homopolymer runs to be same quality scores");

	setOption(pars.compPerCutOff, "-compPerCutOff",
			"Percentage of reads in one direction Cut Off");
	setOption(pars.useCompPerCutOff, "-useCompPerCutOff",
			"Throw out clusters that are made up reads of > "
					+ estd::to_string(pars.useCompPerCutOff * 100)
					+ "% of reads in only one direction, used with Ion Torrent Reads");

	setOption(pars.diffCutOffStr, "-diffCutOffs",
			"Difference Cutoff to form Nuc Comp Clusters");

	pars_.colOpts_.nucCompBinOpts_.diffCutOffVec_ = vecStrToVecNum<double>(
			tokenizeString(pars.diffCutOffStr, ","));

	setOption(pars_.colOpts_.nucCompBinOpts_.findBestNuc_, "-findBestNuc", "Find Best Nucleotide Cluster");
	setOption(pars_.colOpts_.nucCompBinOpts_.useNucComp_, "-useNucComp",
			"Cluster on Nucleotide Composition First");

	setOption(pars_.colOpts_.nucCompBinOpts_.useMinLenNucComp_, "-useMinLenNucComp",
			"Use Nucleotide Composition of Only Front of seqs");

	setOption(pars.startWithSingles, "-startWithSingles",
			"Start The Clustering With Singletons, rather then adding them afterwards");
	setOption(pars.mapBackSinglets, "-mapBackSinglets",
				"If Singlets are left out, map them back in for frequency estimates");
	setOption(pars.singletCutOff, "-singletCutOff",
				"Naturally the cut off for being a singlet is by default 1 but can use -singletCutOff to raise the number");

	setOption(pars.createMinTree, "--createMinTree",
			"Create Psudo minimum Spanning Trees For Mismatches for Final Clusters");
	setOption(pars_.colOpts_.kmerBinOpts_.useKmerBinning_, "-useKmerBinning", "useKmerBinning");
	setOption(pars_.colOpts_.kmerBinOpts_.kmerCutOff_, "-kmerCutOff", "kmerCutOff");
	setOption(pars_.colOpts_.kmerBinOpts_.kCompareLen_, "-kCompareLen", "kCompareLen");

	setOption(pars.leaveOutSinglets, "-leaveOutSinglets",
			"Leave out singlet clusters out of all analysis");
	setOption(pars.onPerId, "-onPerId", "Cluster on Percent Identity Instead");

	pars_.colOpts_.iTOpts_.removeLowQualityBases_= setOption(pars_.colOpts_.iTOpts_.lowQualityBaseTrim_, "-qualTrim",
			"LowQualityCutOff");
	setOption(pars.qualRep, "-qualRep",
			"QualityRepresentative_for_unique_clusters");
	setOption(pars.extra, "-extra", "Extra");
	setOption(pars_.chiOpts_.checkChimeras_, "-markChimeras", "MarkChimeras");
	setOption(pars_.chiOpts_.parentFreqs_, "-parfreqs", "Parent_freq_multiplier_cutoff");

	setOption(pars.snapShotsOpts_.snapShots_, "-snapshots", "OutputSnapShots");
	setOption(pars.sortBy, "-sortBy", "SortClustersBy");
	pars.additionalOut = setOption(pars.additionalOutLocationFile,
			"-additionalOut", "AdditionalOutFilename");
	setOption(pars.collapsingTandems, "-collapseTandems", "CollapsingTandems");

	setOption(pars_.colOpts_.alignOpts_.noAlign_, "--noAlignCompare",
			"Do comparisons without aligning");
	processSkipOnNucComp();
	processVerbose();
	processDebug();
	pars_.colOpts_.verboseOpts_.verbose_ = pars_.verbose_;
	pars_.colOpts_.verboseOpts_.debug_ = pars_.debug_;
	processRefFilename();
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "5,1";
	processAlignerDefualts();
	pars.smallReadSize = pars_.colOpts_.kmerOpts_.kLength_ * 2;
	setOption(pars.smallReadSize, "--smallReadSize",
			"A cut off to remove small reads");
	if (!failed_) {
		if ("" != pars.parameters) {
			if (pars.onPerId) {
				pars.iteratorMap = processIteratorMapOnPerId(pars.parameters);
			} else {
				//std::cout << __PRETTY_FUNCTION__ << 2 << std::endl;
				pars.iteratorMap = processIteratorMap(pars.parameters);
				//std::cout << __PRETTY_FUNCTION__ << 3 << std::endl;
			}
		}
		if ("" != pars.binParameters) {
			if (pars.onPerId) {
				pars.binIteratorMap = processIteratorMapOnPerId(pars.binParameters);
			} else {
				pars.binIteratorMap = processIteratorMap(pars.binParameters);
			}
		} else {
			pars.binIteratorMap = pars.iteratorMap;
		}
		if (pars_.verbose_) {
			std::cout << "p: " << pars_.qScorePars_.primaryQual_ << std::endl;
			std::cout << "s: " << pars_.qScorePars_.secondaryQual_ << std::endl;
			std::cout << "go: " << pars_.gapInfo_.gapOpen_ << std::endl;
			std::cout << "ge: " << pars_.gapInfo_.gapExtend_ << std::endl;
		}
	}
	finishSetUp(std::cout);
}

void SeekDeepSetUp::setUpMultipleSampleCluster(processClustersPars & pars) {
	// parse the command line options
	if (needsHelp()) {
		std::stringstream tempOut;
		tempOut << "processClusters" << std::endl;
		tempOut << "Ran from inside directory tree set up such that "
				"currentDir/SampleDir/RunDir/seqFiles" << std::endl;
		tempOut << "Example dir tree: " << std::endl;
		tempOut << "Samp01/Run1/output.fastq" << std::endl;
		tempOut << "Samp01/Run2/output.fastq" << std::endl;
		tempOut << "Samp02/Run1/output.fastq" << std::endl;
		tempOut << "Samp02/Run2/output.fastq" << std::endl;
		tempOut << "currentDir/" << std::endl;
		tempOut << "----Samp01/" << std::endl;
		tempOut << "--------Run1/" << std::endl;
		tempOut << "------------output.fastq" << std::endl;
		tempOut << "--------Run2/" << std::endl;
		tempOut << "------------output.fastq" << std::endl;
		tempOut << "----Samp02/" << std::endl;
		tempOut << "--------Run1/" << std::endl;
		tempOut << "------------output.fastq" << std::endl;
		tempOut << "--------Run2/" << std::endl;
		tempOut << "------------output.fastq" << std::endl;

		tempOut << "program would be run from the directory that"
				" contains Samp01 and Samp02 directories (here currentDir)"
				<< std::endl;
		tempOut
				<< "sequences in the seq files should have a suffix of _[NUM] or _t[NUM] "
						"where [NUM] is the number of reads associated with that sequence"
				<< std::endl;
		tempOut << "Final results are named with the name of the Sample Directories"
				<< std::endl;
		tempOut << "Commands, order not necessary and case insensitive"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Required commands" << bib::bashCT::reset
				<< std::endl;
		printInputUsage(tempOut);
		tempOut
				<< "in the previous example tree the flag would be -fastq output.fastq and all seq "
						"files should be named output.fastq, others will be ignored"
				<< std::endl;
		tempOut
				<< "-par : parameters file, see SeekDeep qluster --help for details on the format of the file"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Optional commands" << bib::bashCT::reset
				<< std::endl;
		tempOut << "-dout [option]: Name of an output directory, will default "
				"to cluster_CURRENT_DATE" << std::endl;
		tempOut << "-onPerId : Cluster reads on percent identity rather"
				" than specific errors, format of parameters file would have to change,"
				" see SeekDeep qluster --help for details" << std::endl;
		tempOut
				<< "-experimentName [option] : Name prefix to give to the final population haplotypes, defaults to PopUID"
				<< std::endl;

		tempOut << bib::bashCT::boldBlack("Final analysis inclusion criteria")
				<< std::endl;
		tempOut << "-fracCutOff [option] : The fraction threshold to be "
				"included in final analysis, threshold is compared to the averaged "
				"fraction across runs, defaults to 0.005 " << std::endl;
		tempOut << "-runsRequired [option] : Number of runs a cluster has to "
				"appear in to be"
				" included in final analysis, defaults to all the runs included "
				"in the sample" << std::endl;
		tempOut << bib::bashCT::bold << "Population clustering options"
				<< bib::bashCT::reset << std::endl;
		tempOut << "-popPar [option] : Separate population paramters for "
				"clustering will default to the parameters given for the "
				"between sample clustering" << std::endl;

		tempOut << bib::bashCT::bold << "Additional Pre-processing options"
				<< bib::bashCT::reset << std::endl;
		tempOut << "-cutoff [option]: Size of cluster not to include in "
				"clustering, defaults to 1" << std::endl;
		tempOut << "-markChimeras : Have the program mark possible chimeras "
				"before starting clustering" << std::endl;
		tempOut << "-parFreqs [option] : The minimum frequency multiplier "
				"reads must have in order to be considered to be a possible "
				"parrent of a chimera, defaults to 2" << std::endl;
		tempOut << "-keepChimeras : Whether to include chimeras in the final "
				"analysis,"
				" defaults to excluded any cluster made of at least half of "
				"reads with CHI_ flag in"
				" in their name " << std::endl;
		printReferenceComparisonUsage(tempOut);
		printAdditionalClusteringUsage(tempOut);
		printAlignmentUsage(tempOut);
		printQualThresUsage(tempOut);
		printAlnInfoDirUsage(tempOut);
		tempOut << "examples" << std::endl;
		tempOut << "SeekDeep processClusters -fasta output.fasta "
				"-par pars.txt" << std::endl;
		tempOut << "SeekDeep processClusters -fasta output.fasta "
				"-par pars.txt" << std::endl;
		tempOut << "SeekDeep processClusters -fasta output.fasta -par "
				"pars.txt -popPar otherPars.txt" << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		tempOut << bib::bashCT::bold << "Output Files:" << bib::bashCT::reset
				<< std::endl;
		tempOut
				<< "selectedClustersInfo.tab.txt: This contains the final haplotype information and replicate comparison results, it is a very large table, consult the SeekDeep (http://bib2.umassmed.edu/~hathawan/SeekDeep.html) website for details on what each column means"
				<< std::endl;
		tempOut << "allClustersInfo.tab.txt: " << std::endl;
		tempOut << "dotFiles: " << std::endl;
		tempOut << "final: " << std::endl;
		tempOut << "population/populationCluster.tab.txt: " << std::endl;
		tempOut << "population/PopUID.fastq: " << std::endl;
		tempOut << "" << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		exit(0);
	}

	setOption(pars.ionTorrent, "-ionTorrent",
			"Flag to indicate reads are ion torrent and therefore turns on -useCompPerCutOff,-adjustHomopolyerRuns, and -qualTrim");
	if (pars.ionTorrent) {
		pars.removeLowQualBases = true;
		pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_ = true;
	}

	setOption(pars.masterDir, "--masterDir", "Master directory containing sample sequence files");
	pars.masterDir = bfs::absolute(pars.masterDir).string();
	pars.removeLowQualBases = setOption(pars.lowQualityCutOff, "-qualTrim",
			"LowQualityCutOff");
	setOption(pars.writeExcludedOriginals, "--writeExcludedOriginals",
			"Write out the excluded Originals");
	setOption(pars.chiCutOff, "--chiCutOff",
			"The Fraction of a cluster to determine if it chimeric");
	setOption(pars.recheckChimeras, "--recheckChimeras",
			"Re Check chimeras after replicate comparison");
	setOption(pars.eventBasedRef, "-eventBasedRef", "Do Event Based Ref Count");
	setOption(pars.customCutOffs, "--custumCutOffs",
			"Two Column Table, first column is sample name, second is a custom frac cut off, if sample not found will default to -fracCutOff");
	setOption(pars.previousPopFilename, "-previousPop", "previousPopFilename");
	processComparison(pars.previousPopErrors, "previousPop");
	setOption(pars.groupingsFile, "--groupingsFile",
			"A file to sort samples into different groups");
	setOption(pars.investigateChimeras, "--investigateChimeras",
			"Check to see if a chimera appears as a high variant in another sample");
	processDebug();
	processVerbose();
	pars_.colOpts_.verboseOpts_.verbose_ = pars_.verbose_;
	pars_.colOpts_.verboseOpts_.debug_ = pars_.debug_;
	setOption(pars.onPerId, "-onPerId", "Cluster on Percent Identity Instead");
	processSkipOnNucComp();
	// input file info
	pars_.ioOptions_.lowerCaseBases_ = "upper";
	pars_.ioOptions_.processed_ = true;

	processDefaultReader(true);
	bool noErrorsSet  = false;
	setOption(noErrorsSet, "--noErrors", "Parameters with no errors");
	bool strictErrorsSet = false;
	setOption(strictErrorsSet, "--strictErrors", "Parameters with a several low quality mismatches");
	bool strictErrorsSetHq1 = false;
	setOption(strictErrorsSetHq1, "--strictErrors-hq1", "Parameters with a several low quality mismatches and 1 high quality mismatch");

	uint32_t hqMismatches = 0;
	setOption(hqMismatches, "--hq", "Number of high quality mismatches to allow");

	setOption(pars.parameters, "-par", "ParametersFileName", !noErrorsSet && !strictErrorsSet && !strictErrorsSetHq1);

	setOption(pars.binParameters, "-binPar", "bin Parameters Filename");

	setOption(pars.runsRequired, "-runsRequired", "runsRequired");
	setOption(pars.experimentName, "-experimentName", "ExperimentName");
	if (bib::containsSubString(pars.experimentName, ".")) {
		addWarning("Error in populationCollapse::populationCollapse, populationName can't contain '.', "
						+ pars.experimentName);
		failed_ = true;
	}
	setOption(pars.clusterCutOff, "-clusterCutOff", "Cluster Size Cut Off");
	setOption(pars.noTrees, "--noTrees", "Don't generate html difference trees");
	processDirectoryOutputName("clusters_" + getCurrentDate(), true);

	setOption(pars.extra, "-extra", "ExtraOutput");
	processRefFilename();
	setOption(pars.noPopulation, "-noPopulation",
			"Don't do Population Clustering");

	setOption(pars.fracCutoff, "-fracCutOff",
			"PopulationClusteringFractionCutoff");
	pars.differentPar = setOption(pars.parametersPopulation, "-popPar",
			"ParametersForPopulationCollapse");
	setOption(pars_.chiOpts_.checkChimeras_, "-markChimeras", "MarkChimeras");
	setOption(pars.keepChimeras, "-keepChimeras", "KeepChimeras");
	setOption(pars_.chiOpts_.parentFreqs_, "-parFreqs", "ParentFrequence_multiplier_cutoff");

	setOption(pars.numThreads, "--numThreads", "Number of threads to use");

	setOption(pars.plotRepAgreement, "-plotRepAgreement", "Plot Rep Agreement");
	processAlignerDefualts();
	if (pars_.debug_ && !failed_) {
		std::cout << "p: " << pars_.qScorePars_.primaryQual_ << std::endl;
		std::cout << "s: " << pars_.qScorePars_.secondaryQual_ << std::endl;
		std::cout << "go: " << pars_.gapInfo_.gapOpen_ << std::endl;
		std::cout << "ge: " << pars_.gapInfo_.gapExtend_ << std::endl;
	}
	if (!failed_ ) {
		if("" != pars.parameters){
			// read in the parameters from the parameters file
			if (pars.onPerId) {
				pars.iteratorMap = processIteratorMapOnPerId(pars.parameters);
			} else {
				pars.iteratorMap = processIteratorMap(pars.parameters);
			}
		}else if(noErrorsSet){
			pars.iteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(100);
		}else if(strictErrorsSet){
			if (hqMismatches > 0) {
				pars.iteratorMap = CollapseIterations::genStrictDefaultParsWithHqs(100,
						hqMismatches);
			} else {
				pars.iteratorMap = CollapseIterations::genStrictDefaultPars(100);
			}
		} else if (strictErrorsSetHq1) {
			pars.iteratorMap = CollapseIterations::genStrictDefaultParsHq1(100);
		}

		if (pars.binParameters != "") {
			if (pars.onPerId) {
				pars.binIteratorMap = processIteratorMapOnPerId(pars.binParameters);
			} else {
				pars.binIteratorMap = processIteratorMap(pars.binParameters);
			}
		}
		if (pars.differentPar) {
			if (pars.onPerId) {
				pars.popIteratorMap = processIteratorMapOnPerId(pars.parametersPopulation);
			} else {
				pars.popIteratorMap = processIteratorMap(pars.parametersPopulation);
			}
		} else {
			pars.popIteratorMap = pars.iteratorMap;
		}
	}
	finishSetUp(std::cout);

}

void SeekDeepSetUp::setUpMakeSampleDirectories(
		makeSampleDirectoriesPars & pars) {
	if (needsHelp()) {
		std::stringstream tempStream;
		tempStream << "makeSampleDirectoires" << std::endl;
		tempStream << "Set up a directory tree for processClusters" << std::endl;
		tempStream << "Commands, order not necessary, flags are case insensitive"
				<< std::endl;
		tempStream << "Required commands" << std::endl;
		tempStream << "-file [option], name of the file of sample names to read in"
				<< std::endl;
		tempStream << "-dout [option], name of the main directory to create"
				<< std::endl;
		tempStream << "File should be tab delimited and a few examples are below"
				<< std::endl;
		tempStream << "File should have at least three columns" << std::endl;
		tempStream << "Where first column is the name of the index or sff file "
				"used, second column is the sample names, and all following "
				"columns are the MIDs for that samples " << std::endl;
		std::cout << cleanOut(tempStream.str(), width_, indent_);
		tempStream.str(std::string());
		std::cout << "Example with two replicates and two separate master indexes"
				<< std::endl;
		std::cout << "1\t090-00\tMID01\tMID02" << std::endl;
		std::cout << "1\t090-24\tMID03\tMID04" << std::endl;
		std::cout << "1\t090-48\tMID05\tMID06" << std::endl;
		std::cout << "1\t090-72\tMID07\tMID08" << std::endl;
		std::cout << "2\t095-00\tMID01\tMID02" << std::endl;
		std::cout << "2\t095-24\tMID03\tMID04" << std::endl;
		std::cout << "2\t095-48\tMID05\tMID06" << std::endl;
		std::cout << "2\t095-72\tMID07\tMID08" << std::endl;
		std::cout << std::endl;
		std::cout << "Example with one replicate" << std::endl;
		std::cout << "1\t090-00\tMID01" << std::endl;
		std::cout << "1\t090-24\tMID02" << std::endl;
		std::cout << "2\t095-00\tMID01" << std::endl;
		std::cout << "2\t095-24\tMID02" << std::endl;
		std::cout << std::endl;
		std::cout << "Example with with mix amount of replicates" << std::endl;
		std::cout << "1\t090-00\tMID01" << std::endl;
		std::cout << "1\t090-24\tMID02" << std::endl;
		std::cout << "2\t095-00\tMID01\tMID02" << std::endl;
		std::cout << "2\t095-24\tMID03\tMID04" << std::endl;

		std::cout << "examples, SeekDeep makesampleDirectories -file "
				"names.tab.txt -dout clustering" << std::endl;
		exit(0);
	}
	processVerbose();
	setOption(pars.separatedDirs, "--separatedDirs",
				"Create a separate directory for each index");
	setOption(pars.sampleNameFilename, "--file", "Sample Names Filename", true);
	setOption(pars_.directoryName_, "-dout", "Main Out Directory Name", true);
	setOption(pars_.overWriteDir_, "--overWriteDir",
			"If the directory already exists over write it");
	if (!failed_) {
		std::string newDirectoryName = "./"
				+ bib::replaceString(bib::replaceString(pars_.directoryName_, "./", ""), "TODAY",
						getCurrentDate()) + "/";
		pars_.directoryName_ = bib::files::makeDir("./",
				bib::files::MkdirPar(bib::replaceString(pars_.directoryName_, "./", ""),
						pars_.overWriteDir_));
	}
	finishSetUp(std::cout);
}

}  // namespace bibseq
