#include "SeekDeepSetUp.hpp"
#include <bibcpp/bashUtils.h>
namespace bibseq {

void SeekDeepSetUp::setUpExtractor(
    std::string& idFileName, bool& multiplex, bool& condensed, int& minLen,
    int& maxLength, int& within, int& qualityWindowLength,
    int& qualityWindowStep, int& qualityWindowThres, bool& unknown,
    int& unknownPrimerSize, int& unknownMinSize, bool& findReversePrimer,
    double& queryCoverageCutoff, double& percentIdentityCutoff, int& numberOfNs,
    bool& pyroExtract, bool& reversePrimerToLowerCase, bool& flowFiltering,
    int& maxFlows, bool& checkComplement, bool& screenForPossibleContaimination,
    std::string& compareSeq, std::string& idFileDelim, int& smallFragmentCutOff,
    bool& qualWindowTrim) {
  //ioOptions_.lowerCaseBases = "remove";
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "sffExtractor" << std::endl;
    tempOut << "Extracts sequences from  sff.txt files and some "
               "pre-filtering" << std::endl;
    tempOut << "Commands, order not necessary and case insensitive"
            << std::endl;
    tempOut << bib::bashCT::bold << "Required commands"<< bib::bashCT::reset << std::endl;
    tempOut << "Input can be fasta, fastq, or sff.txt" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    printInputUsage(std::cout);
    tempOut << "1e) -sff [option]: Full name of sff.txt file converted from "
               ".sff by sffinfo" << std::endl;
    tempOut
        << "2) -id [option]: Full name of the id file with primers and Ids info"
        << std::endl;

    tempOut << "\tThe id file should be whitespace delimited (tabs or just simply spaces) and contains the "
               "primer and reverse primer and MIDs if the data is multiplex, "
               "example bellow" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    std::cout << "\tgene	forwardPrimer	reversePrimer" << std::endl;
    std::cout << "\tPFMSP1	AACTAGAAGCTTTAGAAGATGCA	ACATATGATTGGTTAAATCAAAG"
              << std::endl;
    std::cout << "\tid	barcode" << std::endl;
    std::cout << "\tMID01	ACGAGTGCGT" << std::endl;
    std::cout << "\tMID02	ACGCTCGACA" << std::endl;
    tempOut << bib::bashCT::bold << "Optional commands"<< bib::bashCT::reset << std::endl;
    tempOut << "1) -dout [option]: Name of an output directory, will default "
               "to the name of the input file plus the date if none given"
            << std::endl;
    tempOut << "2) -multiplex : Whether the id file is multiplex with "
               "multiple MIDs, the program will assume only primers are "
               "given if this is not switched on" << std::endl;
    tempOut << "3) -checkComplement : Check the complement of the sequences "
               "for the barcodes as well" << std::endl;
    tempOut << "Quality filtering options" << std::endl;
    tempOut << "1) -minLen [option]: The minimum length for the read to be "
               "extracted including primers, defaults to 200" << std::endl;
    tempOut << "2) -maxLen [option]: The maximum length for the read to be "
               "extracted including primers, defaults to 300" << std::endl;
    tempOut << "3) -qualWindow [option]: Quality window checking options, "
               "given separated by commas and in the order of windowSize, "
               "windowStep, and windowMinimumThres, ex. 50,5,20 (default)"
            << std::endl;
    tempOut << "4) -numberOfNs [option]: the number of N's(the symbol for "
               "unknown base) to be allowed in sequence, defaults to 0"
            << std::endl;
    tempOut << "5) -samllFragmentCutOff [option]: the size cutoff for a "
               "sequence to be considered simply a fragment, defaults to 50"
            << std::endl;
    tempOut << bib::bashCT::bold <<
                   "Options if flow values are supplied as well from 454 data"<< bib::bashCT::reset            << std::endl;
    tempOut << "1) -flowFiltering : Whether to filter on the first noisy "
               "flow, defined as a flow between 0.5 and 0.7 " << std::endl;
    tempOut << "2) -pyro : Where to output flow gram files in a format that "
               "pyronoise/ampiconnoise uses" << std::endl;
    tempOut << "3) -maxFlows [option]: The number of maxFlows to be allowed, "
               "titantium should be 720 and pre-titantium should be 260"
            << std::endl;
    tempOut << bib::bashCT::bold << "Options for looking for primers"<< bib::bashCT::reset << std::endl;
    tempOut << "1) -noReverse: if you don't want to look for the reverse "
               "primer as well as the forward primer, will default to "
               "looking for both" << std::endl;
    // tempOut<<"2) -reversepar: parameters for allowing errors in the reverse
    // primer, will default to letting one 1 base indel and some homopolymer
    // difference"<<std::endl;
    tempOut << "2) -rPrimerCoverage [option] : fraction of the reverse "
               "primer that needs to be present, defaults to 0.5 so at least "
               "half of the primer needs to be present" << std::endl;
    tempOut << "3) -rPrimerIdCutoff [option] : of the coverage how much "
               "needs to match, defaults to 0.8 so of the primer present, at "
               "least 80% of it has to match the reverse primer" << std::endl;
    tempOut << "4) -reverseLower: set this to make the reverse primer into "
               "lower case, the default is to keep as upper case" << std::endl;
    tempOut << bib::bashCT::bold << "Options for screening for possible contamination"<< bib::bashCT::reset            << std::endl;
    tempOut << "1) -contamination : Where to screen for possible "
               "contamination, a comparison sequence needs to be supplied if "
               "this is turned on" << std::endl;
    tempOut << "2) -compareSeq [option]: A comparison sequence to which the "
               "reads are compared and if the alignment score falls bellow "
               "zero the read is considered contamination, therefore this "
               "flag only works if you are extracting sequences that are at "
               "least similar to each other and won't be confused as "
               "contamination" << std::endl;
    tempOut << bib::bashCT::bold << "Options for when primers are unknown"<< bib::bashCT::reset            << std::endl;
    tempOut << "1) -unknown : switch for when the primer info is unknown"
            << std::endl;
    tempOut << "2) -unknownPrimerSize [option]: Size of unknown primers, "
               "defaults to 10 which means the reads will be extracted using "
               "the first 10 bases of the reads" << std::endl;
    tempOut << "3) -unknownMinSize [option]: Minimum number of reads to be "
               "extracted for each unknown primer, defaults to 100 which "
               "means after the reads are organized by the unknownPrimerSize "
               "(see above), this is the minimum number of reads that have "
               "that primer in order to be extracted" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    std::cout << "examples: " << std::endl
              << "\tSeekDeep sffExtractor -sff reads.sff.txt -id "
                 "idFile.tab.txt -dout outPutDir " << std::endl
              << "\tSeekDeep sffExtractor "
                 "-fastq reads.fastq -id idFile.tab.txt" << std::endl;

    exit(0);
  }
  processDefaultReader(true);
  bool mustMakeDirectory = true;
  processDirectoryOutputName(mustMakeDirectory);
  if (!setOption(idFileName, "-id", "The name of the ID file")) {
    if (!setOption(unknown, "-unknown", "If the id's are unknown and the first part of seq is used instead")) {
      warnings_.emplace_back(
          "need to have -id or -unknown, see SeekDeep extractor "
          "-help for details");
      failed_ = true;
    } else {
      setOption(unknownPrimerSize, "-unknownPrimerSize", "The size of the first part of the sequence that should be used");
      setOption(unknownMinSize, "-unknownMinSize", "The size cutoff to for the unknown id cluster size to be kept");
    }
  }else{
  	if(!fexists(idFileName)){
  		failed_ = true;
  		 warnings_.emplace_back(
  		          "error in opening " + idFileName);
  	}
  }
  setOption(queryCoverageCutoff, "-rPrimerCoverage",
            "Amount Of Reverse Primer Required");
  setOption(percentIdentityCutoff, "-rPrimerIdCutoff",
            "Percent Identity Of Reverser Primer Required");
  setOption(idFileDelim, "-idFileDelim", "Id File Delim");
  // unknown primers and ids
  setOption(multiplex, "-multiplex", "Indicates that the Reads are multiplex barcoded");
  setOption(condensed, "-condensed", "Check_with_condensed_seq");
  setOption(checkComplement, "-checkComplement", "Check the Complement of the Seqs As Well");
  // setOption(within, "-within");

  setOption(minLen, "-minlen", "Minimum read length");
  setOption(maxLength, "-maxlen", "Maximum read length");
  setBoolOptionFalse(findReversePrimer, "-noReverse", "Don't look for reverse Primer");
  setBoolOptionFalse(reversePrimerToLowerCase, "-reverseUpper","Leave reverse primer upper case");
  setOption(numberOfNs, "-numberofns", "Number Of Ns Cut Off");
  std::string qualWindow = "";
  if (setOption(qualWindow, "-qualWindow", "Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold")) {
    seqUtil::processQualityWindowString(qualWindow, qualityWindowLength,
                                        qualityWindowStep, qualityWindowThres);
  } else {
    qualityWindowLength = 50;
    qualityWindowStep = 5;
    qualityWindowThres = 25;
  }
  if (setOption(flowFiltering, "-flowFiltering", "Filter on flowgram (454)")) {
    if (!setOption(maxFlows, "-maxFlows", "Max Flows", true)) {
      warnings_.emplace_back("Need to supply maxFlows if doing flowFiltering");
    }
  }
  // printOutMapContents(commands.arguments, std::cout);
  // std::cout<<std::endl;
  if (setOption(pyroExtract, "-pyro", "Do Pyro extract for files neccessary for Pyro Noise")) {
    if (!setOption(maxFlows, "-maxFlows", "Max Flows", true)) {
      warnings_.emplace_back("Need to supply maxFlows if doing pyroExtract");
      // printOutMapContents(commands.arguments, std::cout);
    }
  }
  if (setOption(screenForPossibleContaimination, "-contamination",
                "Screening For Contamination")) {
  	if (! processSeq(compareSeq, "-compareSeq,-seq", "Comparison Sequence For Contamination ", true) ) {
  	      warnings_.emplace_back(
  	          "Need to have -compareSeq if checking for contamination");
  	 }
  }
  setOption(qualWindowTrim, "-qualWindowTrim", "Trim To Qual Window");
  setOption(smallFragmentCutOff, "-smallFragmentCutOff", "Small Fragment Cut Off Size");
  gapInfo_.gapOpen_ = 7.0;
  gapInfo_.gapExtend_ = 1;
  gap_ = "7.0,0.5";
  processGap();
  finishSetUp(std::cout);
}

void SeekDeepSetUp::setUpExtractor(extractorPars & pars) {
  //ioOptions_.lowerCaseBases = "remove";
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "extractor" << std::endl;
    tempOut << "Extracts sequences from fasta/fastq  files and some "
               "pre-filtering" << std::endl;
    tempOut << "Commands, order not necessary and case insensitive"
            << std::endl;
    tempOut << bib::bashCT::bold << "Required commands"<< bib::bashCT::reset << std::endl;
    tempOut << "Input can be fasta or fastq" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    printInputUsage(std::cout);
    tempOut
        << "2) -id [option]: Full name of the id file with primers and Ids info"
        << std::endl;

    tempOut << "\tThe id file should be tab delimited and contains the "
               "primer and reverse primer and MIDs if the data is multiplex, "
               "example bellow" << std::endl;
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
    tempOut << bib::bashCT::bold << "Optional commands"<< bib::bashCT::reset << std::endl;
    tempOut << "-dout [option]: Name of an output directory, will default "
               "to the name of the input file plus the date if none given"
            << std::endl;
    tempOut << "-multiplex : Whether the id file is multiplex with "
               "multiple MIDs, the program will assume only primers are "
               "given if this is not switched on" << std::endl;
    tempOut << "-variableStart[option]: If the primers and barcodes aren't directly at the beginning of"
    		" the sequence this can be given to extend the search up this position" << std::endl;
    tempOut << "-checkComplement : Check the complement of the sequences "
               "for the barcodes and primers as well" << std::endl;
    tempOut << "-rename : Rename the sequence to the associated primers and mids "<< std::endl;
    tempOut << bib::bashCT::boldBlack("Quality filtering options") << std::endl;
    tempOut << " -minLen [option]: The minimum length for the read to be "
               "extracted including primers, defaults to " << pars.minLen << std::endl;
    tempOut << " -maxLen [option]: The maximum length for the read to be "
               "extracted including primers, defaults to " << pars.maxLength << std::endl;
    tempOut << " -qualWindow [option]: Quality window checking options, "
               "given separated by commas and in the order of windowSize, "
               "windowStep, and windowMinimumThres, ex. 50,5,20 (default)"
            << std::endl;
    tempOut << "-qualCheck [option]:	Alternative to using using a sliding quality window threshold and is good "
    		"with illumina, this methods filters reads where the fraction of their qualities fall bellow this value, "
    		"use in conjunction with -qualcheckCutOff, eg. setting -qualCheck 25 -qualCheckCutOff 0.75 would filter"
    		" off reads if less than 75% of their qualities were above 25; default " << pars.qualCheck << std::endl;
    tempOut << "-qualCheckCutOff [option]: Cut Off for fraction of bases above qual "
    		"check of "<< pars.qualCheck << "; default=" << pars.qualCheckCutOff << std::endl;
    tempOut << " -numberOfNs [option]: the number of N's(the symbol for "
               "unknown base) to be allowed in sequence, defaults to 0"
            << std::endl;
    tempOut << " -samllFragmentCutOff [option]: the size cutoff for a "
               "sequence to be considered simply a fragment, defaults to 50"
            << std::endl;

    tempOut << bib::bashCT::bold << "Options for looking for primers"<< bib::bashCT::reset << std::endl;
    tempOut << "-noReverse: if you don't want to look for the reverse "
               "primer as well as the forward primer, will default to "
               "looking for both" << std::endl;
    tempOut << "-rPrimerCoverage [option] : fraction of the reverse "
               "primer that needs to be present, defaults to " << pars.rPrimerErrors.distances_.queryCoverage_
							 <<"  so at least "
               "half of the primer needs to be present" << std::endl;
    tempOut << "-rPrimerMismatches [option] : Number of mismatches to allow in reverse primer "
               ", defaults to " << pars.rPrimerErrors.hqMismatches_ << std::endl;
    tempOut << "-rUpper: set this to make the reverse primer into "
               "upper case, the default is to make it and everything that follows it lower case" << std::endl;
    tempOut << "-noForwardPrimer: Don't look for the forward primer, this means only one primer pair can be put in the id file" << std::endl;
    tempOut << "-fPrimerCoverage [option] : fraction of the forward "
               "primer that needs to be present, defaults to " << pars.fPrimerErrors.distances_.queryCoverage_ << std::endl;
    tempOut << "-fPrimerMismatches [option] : Number of mismatches to allow in forward primer "
               ", defaults to " << pars.fPrimerErrors.hqMismatches_ << std::endl;
    tempOut << "-fUpper: set this to make the forward primer into "
               "upper case, the default is to make everything up to and including to it lower case" << std::endl;
    tempOut << bib::bashCT::bold << "Options for screening for possible contamination"<< bib::bashCT::reset << std::endl;
    tempOut << "-contamination : Where to screen for possible "
               "contamination, a comparison sequence needs to be supplied if "
               "this is turned on, this option is supported for when using multiple primer pairs yet" << std::endl;
    tempOut << "-compareSeq [option]: A comparison sequence to which the "
               "reads are compared and and if they are dissimilar enough they will be removed as contamination" << std::endl;
    tempOut << "-contaminationKLen [option]:	Contamination Kmer Length, defaults to " <<  pars.contaminationKLen << std::endl;
    tempOut << "-kmerCutOff [option]:	Kmer cut off for contamination check, defaults to " <<  pars.kmerCutOff  << std::endl;

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
    tempOut << bib::bashCT::bold << "Output Files:"<< bib::bashCT::reset << std::endl;
    tempOut << bib::bashCT::boldBlack("extractionProfile.tab.txt") << ": This breaks down the filtering per final extraction sequence file" << std::endl;
    tempOut << bib::bashCT::boldBlack("extractionStats.tab.txt") << ": This has info how the whole extraction went" << std::endl;
    tempOut << bib::bashCT::boldBlack("seqFiles") << ": Extracted files will be named [PrimerName][MIDNAME].fastq, if not multiplexed then just [PrimerName].fastq" << std::endl;
    tempOut << bib::bashCT::boldBlack("filteredOff") << ": Reads that failed filtering will be found in this directory along with reads that don't match any supplied barcodes and/or primers" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    exit(0);
  }

  setOption(pars.noForwardPrimer, "-noForwardPrimer", "No Forward Primer Required");
  setOption(pars.forwardPrimerToUpperCase, "-fUpper","Leave reverse primer upper case");
  setOption(pars.fPrimerErrors.distances_.queryCoverage_, "-fCoverage", "Amount of Forward Primer to find");
  setOption(pars.fPrimerErrors.hqMismatches_, "-fNumOfMismatches", "Number of Mismatches to allow in Forward Primer");
  setOption(pars.fPrimerErrors.oneBaseIndel_, "-fOneBaseIndels", "Number Of One base indels to allow in forward primer");
  setOption(pars.fPrimerErrors.twoBaseIndel_, "-fTwoBaseIndels", "Number Of Two base indels to allow in forward primer");
  //fPrimerErrors.largeBaseIndel_ = .99;
  setOption(pars.noReversePrimer, "-noReverse", "Don't look for reverse Primer");
  setOption(pars.reversePrimerToUpperCase, "-rUpper","Leave reverse primer upper case");
  setOption(pars.rPrimerErrors.distances_.queryCoverage_, "-rPrimerCoverage", "Amount Of Reverse Primer Required");
  setOption(pars.rPrimerErrors.hqMismatches_, "-rNumOfMismatches", "Number of Mismatches to allow in Reverse Primer");
  setOption(pars.rPrimerErrors.oneBaseIndel_, "-rOneBaseIndels", "Number Of One base indels to allow in reverse primer");
  setOption(pars.rPrimerErrors.twoBaseIndel_, "-rTwoBaseIndels", "Number Of Two base indels to allow in reverse primer");
  //rPrimerErrors.largeBaseIndel_ = .99;



  setOption(pars.rename, "-rename", "Rename With Barcode Names");
  setOption(pars.barcodeErrors, "-barcodeErrors", "Errors Allowed in Barcode");
  setOption(pars.trimTcag, "-trimTcag", "trim TCAG(454 tag) if it is present ");
  if(setOption(pars.HMP, "-HMP", "HMP")){
  	setOption(pars.primerLen, "-len,-primerLen", "PrimerLen");
  }



  setOption(pars.barcodesBothEnds, "-barcodeBothEnds", "Look for Barcodes in Both Primers");
  processDebug();

  pars.variableStart = setOption(pars.variableStop, "-variableStart", "variableStart");
  if(setOption(pars.qualCheck, "-qualCheck", "Qual Check Level")){
  	pars.checkingQCheck = true;
  }
  if(setOption(pars.qualCheckCutOff, "-qualCheckCutOff",
  		"Cut Off for fraction of bases above qual check of " + estd::to_string(pars.qualCheck))){
  	pars.checkingQCheck = true;
  }
  processDefaultReader(true);
  bool mustMakeDirectory = true;
  processDirectoryOutputName(mustMakeDirectory);
  if(setOption(pars.idFilename, "-id", "The name of the ID file")){
  	if(!fexists(pars.idFilename)){
  		failed_ = true;
  		 warnings_.emplace_back(
  		          "error in opening " + pars.idFilename);
  	}
  }


  setOption(pars.idFileDelim, "-idFileDelim", "Id File Delim");
  // unknown primers and ids
  setOption(pars.multiplex, "-multiplex", "Indicates that the Reads are multiplex barcoded");
  setOption(pars.checkComplement, "-checkComplement", "Check the Complement of the Seqs As Well");
  // setOption(within, "-within");

  setOption(pars.minLen, "-minlen", "Minimum read length");
  setOption(pars.maxLength, "-maxlen", "Maximum read length");


  setOption(pars.numberOfNs, "-numberofns", "Number Of Ns Cut Off");
  std::string qualWindow = "";
  if (setOption(qualWindow, "-qualWindow", "Sliding Quality Window, goes WindowSize,WindowStep,Thresdhold")) {
    seqUtil::processQualityWindowString(qualWindow, pars.qualityWindowLength,
    		pars.qualityWindowStep, pars.qualityWindowThres);
  } else {
  	pars.qualityWindowLength = 50;
  	pars.qualityWindowStep = 5;
  	pars.qualityWindowThres = 25;
  }
  bool compareSeqSuccess = processSeq(pars.compareSeq, "-compareSeq", "Comparison Sequence For Contamination ");
  setOption(pars.contaminationKLen, "-contaminationKLen", "Contamination Kmer Length comparison");
  setOption(pars.kmerCutOff, "-kmerCutOff", "Kmer cut off for contamination check");
	if (setOption(pars.screenForPossibleContamination, "-contamination",
			"Screening For Contamination")) {
		if (!compareSeqSuccess) {
			warnings_.emplace_back(
					"Need to have -compareSeq if checking for contamination");
		}
	}
  setOption(pars.qualWindowTrim, "-qualWindowTrim", "Trim To Qual Window");
  setOption(pars.smallFragmentCutoff, "-smallFragmentCutOff", "Small Fragment Cut Off Size");
  gapInfo_.gapOpen_ = 7;
  gapInfo_.gapExtend_ = 1;
  gap_ = "7,1";
  processGap();
  finishSetUp(std::cout);
}

void SeekDeepSetUp::setUpClusterDown(
    std::string& qualRep, std::string& parameters, bool& extra,
    std::map<int, std::vector<double>>& iteratorMap, bool& smallestFirst,
    bool& markChimeras, double& parFreqs, bool& bestMatch, int& bestMatchCheck,
    bool& snapShots, std::string& sortBy, bool& additionalOut,
    std::string& additonalOutLocationFile, bool& collapsingTandems,
    bool& kmerCheckingOnAlignProfile, bool& condensedCollapse,
    bool& removeLowQualBases, int& lowQualityCutOff,
    bool& adjustHomopolyerRuns, bool & onPerId) {
  // input file info
  ioOptions_.outFilename_ = "output";
  ioOptions_.lowerCaseBases_ = "remove";
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << commands_["-program"] << std::endl;
    tempOut << "Iteratively clusters reads by using the allowable errors given "
               "in a parameters file" << std::endl;
    tempOut << "Commands, order not necessary and case insensitive"
            << std::endl;
    tempOut << bib::bashCT::bold << "Required commands"<< bib::bashCT::reset << std::endl;
    // std::cout << cleanOut(tempOut.str(), width_, indent_);
    // tempOut.str(std::string());
    printInputUsage(tempOut);
    tempOut << "-par [option]: parameter file" << std::endl << std::endl;

    tempOut << "Parameter file is set up as where each line is one iteration and that line indicates "
    		"which errors to allow on that iteration " << std::endl;
    tempOut << "The errors allowed are separated by a colon : " << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    std::cout << "Header should be as follows, followed by the iteration information " << std::endl;
    std::cout << "stopCheck:smallCutoff:1baseIndel:2baseIndel:>2baseIndel:HQMismatches:LQMismatches:LKMismatches" << std::endl;
    std::cout << "100:3:1:0:0:0:0:0" << std::endl;
    std::cout << "100:3:2:0:0:0:0:0" << std::endl;
    std::cout << "100:3:3:0:0:0:1:0" << std::endl;
    std::cout << "100:3:4:0:0:0:2:0" << std::endl;
    std::cout << "100:0:1:0:0:0:0:0" << std::endl;
    std::cout << "100:0:2:0:0:0:0:0" << std::endl;
    std::cout << "100:0:3:0:0:0:1:0" << std::endl;
    std::cout << "100:0:4:0:0:0:2:0" << std::endl;
    tempOut << "So here there would be eight iteration, in the first iteration reads"
    		" will only cluster if they differ only by 1 1base indel" << std::endl;
    tempOut << "The first two columns control how many reads to compare."
    		"  The first column will control the number of reads check, so in this"
    		" instance only the first 100 clusters will be checked.  The second number"
    		" controls size of clusters to compare against, in the first iteration a 3 "
    		"means clusters smaller than 3 will be allowed to have clusters collapsed "
    		"into them though they can still collapse into larger clusters" << std::endl;
    tempOut << "Alternatively if the -onPerId is used, this parameter flag can contain percent idendities instead" << std::endl;
    tempOut << "example: " << std::endl;
    std::cout << "100:3:.99" << std::endl;
    std::cout << "100:3:.99" << std::endl;
    std::cout << "100:3:.98" << std::endl;
    std::cout << "100:3:.98" << std::endl;
    tempOut << "Here, four iteration where clusters were clustered at first .99 percent identity and then .98" << std::endl;
    tempOut << bib::bashCT::bold << "Optional commands"<< bib::bashCT::reset << std::endl;
    tempOut << "-dout [option]: Name of an output directory, will default "
               "to the stub name plus the date" << std::endl;
    tempOut << "-additionaOut [option]: Name of a tab delimited file with "
               "two columns, column one is the name of the ID and column two "
               "is the file path location of additional output directory"
               ", such a file can be made by makeSampleDirectories"
            << std::endl;
    printAdditionalClusteringUsage(tempOut);
    tempOut << bib::bashCT::bold << "Alignment comparison options"<< bib::bashCT::reset << std::endl;
    printQualThresUsage(tempOut);
    printAlignmentUsage(tempOut);
    printKmerUsage(tempOut);
    printAlnInfoDirUsage(tempOut);
    printAdditionaInputUsage(tempOut, ioOptions_.lowerCaseBases_);
    tempOut << "-qualRep [option] : Sets the quality for the identical clusters, "
               "Options are median (default), average, bestQual, or worst"
            << std::endl;
    tempOut << "-qualTrim [option]: Will trim off any bases bellow this "
               "quality, mainly used to trim bases with a quality of 1 (80% "
               "chance of error) by setting this to 2" << std::endl;
    tempOut << "-adjustHomopolyerRuns : Will take average quality across "
               "a homopolymer run and set all the quality to this average, "
               "mainly used with Ion Torrent reads due to the very low "
               "quality of the last base in a long homopolymer run"
            << std::endl;
    tempOut << bib::bashCT::bold << "Additional Processing options"<< bib::bashCT::reset << std::endl;
    tempOut << "-markChimeras : Have the program mark possible chimeras "
               "before outputting" << std::endl;
    tempOut << "-parFreqs [option] : The minimum frequency multiplier "
               "reads must have in order to be considered to be a possible "
               "parent of a chimera or for sequences to be collapsed on gaps "
               "in tandem repeats, defaults to 2" << std::endl;
    tempOut << bib::bashCT::bold << "Technology Specific flags"<< bib::bashCT::reset << std::endl;
    tempOut << "-useCompPerCutOff: This turns on filtering off clustering that have reads that only come from one direction, only should be used if reads were extracted in both direction like in ion torrent" << std::endl;
    tempOut << "-ionTorrent : This turns on several of the previously mentioned flags as they are beneficial for ion torrent data, turns on -qualTrim,-adjustHomopolyerRuns, and -useCompPerCutOff" << std::endl;
    tempOut << bib::bashCT::bold << "Percent Identity Flags"<< bib::bashCT::reset << std::endl;
    tempOut << "-onPerId: cluster on percent identity rather than specific errors" << std::endl;

    printReferenceComparisonUsage(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    std::cout << "examples" << std::endl;
    std::cout << "\tSeekDeep " << commands_["-program"] <<" -fasta MID1.fasta -qual "
                 "MID1.fasta.qual -par par" << std::endl;
    std::cout << "\tSeekDeep " << commands_["-program"] <<" -fastq MID1.fastq -par par"
              << std::endl;
    std::cout << "\tSeekDeep " << commands_["-program"] <<" -fastq MID1.fastq -par par -runCutOff 0.2%"
              << std::endl;
    std::cout << "\tSeekDeep " << commands_["-program"] <<" -fastq MID1.fastq -par par -runCutOff 5"
              << std::endl;
    std::cout << "\tSeekDeep " << commands_["-program"] <<" -fastq MID1.fastq -par par -runCutOff 0.1%,2"
              << std::endl;
    std::cout << "\tSeekDeep " << commands_["-program"] <<" -stub MID1 -par par "
                 "-qualThres 25,20" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    tempOut << bib::bashCT::bold << "Output Files:"<< bib::bashCT::reset << std::endl;
    tempOut << bib::bashCT::boldBlack("output.fastq") << ": This is the final clusters with their consensus sequence, the sequences are named so that the suffif _t[NUM] where NUM is the number of reads that fell into that cluster" << std::endl;
    tempOut << bib::bashCT::boldBlack("outputInfo.tab.txt") << ": This contains cluster number info for each cluster" << std::endl;
    tempOut << bib::bashCT::boldBlack("clusters") << ": A seq file is available for each final cluster that contains all the sequences that clustered into this cluster" << std::endl;
    tempOut << bib::bashCT::boldBlack("internalSnpInfo") << ": A table of internal snps frequencies is available for each cluster in case overcollapsing is suspected" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    exit(0);
  }
  setOption(onPerId, "-onPerId", "Cluster on Percent Identity Instead");
  processDefaultReader(true);
  if(setOption(parameters, "-par", "ParametersFilename", true)){
  	if(onPerId){
  		processIteratorMapOnPerId(parameters, iteratorMap);
  	}else{
  		processIteratorMap(parameters, iteratorMap);
  	}
  }

  //setOption(condensedCollapse, "-condensedCollapse", "condensedCollapse");
  setBoolOptionFalse(condensedCollapse, "-noCondensedCollapse", "condensedCollapse");
  if (setOption(removeLowQualBases, "-qualTrim", "Remove Low Quality Bases")) {
    setOption(lowQualityCutOff, "-qualTrim", "LowQualityCutOff");
  }
  setOption(adjustHomopolyerRuns,
            "-adjustHomopolyerRuns,-adjustHRuns,-adjustHPRuns",
            "AdjustHomopolyerRunsToBeSameQual");
  processVerbose();
  processScoringPars();
  setOption(qualRep, "-qualRep", "QualityRepresentative_for_unique_clusters");
  setOption(extra, "-extra", "Extra");
  setBoolOptionFalse(smallestFirst, "-largestfirst",
                     "CompareSmallestClustersFirst");
  setOption(bestMatch, "-bestMatch", "FindBestMatch");
  setBoolOptionFalse(bestMatch, "-matchFirst,-firstMatch", "FindBestMatch");
  setOption(bestMatchCheck, "-bestmatchcheck", "BestMatchCheckNumber");

  setOption(markChimeras, "-markChimeras", "MarkChimeras");
  setOption(parFreqs, "-parfreqs", "Parent_freq_multiplier_cutoff");
  processRefFilename();
  setOption(snapShots, "-snapshots", "OutputSnapShots");
  setOption(sortBy, "-sortBy", "SortClustersBy");
  if (setOption(additionalOut, "-additionalOut", "AdditionalOutputing")) {
    setOption(additonalOutLocationFile, "-additionalOut",
              "AdditionalOutFilename", true);
  }
  setOption(collapsingTandems, "-collapseTandems", "CollapsingTandems");
  setOption(kmerCheckingOnAlignProfile, "-checkKmerAlnProfile",
            "Leaving_out_low_kmersMismatches");

  bool mustMakeDirectory = true;
  processDirectoryOutputName(mustMakeDirectory);
  processQualThres();
  processGap();
  if (!failed_ && !commands_.containsFlagCaseInsensitive("-flags") &&
      !commands_.containsFlagCaseInsensitive("-getFlags")) {
    std::cout << "p: " << primaryQual_ << std::endl;
    std::cout << "s: " << secondaryQual_ << std::endl;
    std::cout << "go: " << gapInfo_.gapOpen_ << std::endl;
    std::cout << "ge: " << gapInfo_.gapExtend_ << std::endl;
  }
  finishSetUp(std::cout);

}

void SeekDeepSetUp::setUpMultipleSampleCluster(
    std::string& parameters, bool& extra, int& cutOff,
    std::map<int, std::vector<double>>& iteratorMap, bool& population,
    double& fracCutoff, bool& smallestFirst, bool& bestMatch,
    int& bestMatchCheck, bool& checkChimeras, double& parFreqs,
    std::string& parametersPopulation, bool& differentPar,
    std::map<int, std::vector<double>>& popIteratorMap, bool& popBoth,
    bool& keepChimeras, std::string& experimentName, uint32_t& runsRequired, bool & onPerId) {
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
               " contains Samp01 and Samp02 directories (here currentDir)" << std::endl;
    tempOut << "sequences in the seq files should have a suffix of _[NUM] or _t[NUM] "
    		"where [NUM] is the number of reads associated with that sequence" << std::endl;
    tempOut << "Final results are named with the name of the Sample Directories" << std::endl;
    tempOut << "Commands, order not necessary and case insensitive"
            << std::endl;
    tempOut << bib::bashCT::bold << "Required commands"<< bib::bashCT::reset << std::endl;
    printInputUsage(tempOut);
    tempOut << "in the previous example tree the flag would be -fastq output.fastq and all seq "
    		"files should be named output.fastq, others will be ignored" << std::endl;
    tempOut << "-par : parameters file, see SeekDeep qluster -help for details on the format of the file" << std::endl;
    tempOut << bib::bashCT::bold << "Optional commands"<< bib::bashCT::reset << std::endl;
    tempOut << "-dout [option]: Name of an output directory, will default "
               "to cluster_CURRENT_DATE" << std::endl;
    tempOut << "-onPerId : Cluster reads on percent identity rather"
    		" than specific errors, format of parameters file would have to change,"
    		" see SeekDeep qluster -help for details" << std::endl;
    tempOut << "-experimentName [option] : Name prefix to give to the final population haplotypes, "
    		"defaults to PopUID" << std::endl;
    tempOut << "--groupingsFile [option] : This is an optional file that give groupings information "
    		"for the samples in the analysis to output additional information.  The first column contains"
    		" the sample names and each additional column is another grouping category" << std::endl;
    tempOut << "For example if the analysis contained Control1,Control2,Child1,Child2,Adult1,Adult2 "
    		"a groupings file might look like bellow" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    std::cout << "samples  type    age" << std::endl;
    std::cout << "control1 control control" << std::endl;
    std::cout << "control2 control control" << std::endl;
    std::cout << "Child1   patient child" << std::endl;
    std::cout << "Child2   patient child" << std::endl;
    std::cout << "Adult1   patient adult" << std::endl;
    std::cout << "Adult2   patient adult" << std::endl;
    tempOut << "This will create a directory called groups and two additional group directories "
    		"type and age and then further subgroup directories in there.  Each final directory will "
    		"contain info on the specifc samples for that group and the population haplotypes found "
    		"in that group." << std::endl;
    tempOut << "This is still a fairly new feature and will be expanded upon in future releases, "
    		"things to come would be outputting statically analysis between the different groups" << std::endl;
    tempOut << bib::bashCT::boldBlack("Final analysis inclusion criteria") << std::endl;
    tempOut
        << "-fracCutOff [option] : The fraction threshold to be "
           "included in final analysis, threshold is compared to the averaged "
           "fraction across runs, defaults to 0.005 " << std::endl;
    tempOut << "-runsRequired [option] : Number of runs a cluster has to "
               "appear in to be"
               " included in final analysis, defaults to all the runs included "
               "in the sample" << std::endl;
    tempOut << bib::bashCT::bold << "Population clustering options"<< bib::bashCT::reset << std::endl;
    tempOut << "-population : do a population clustering across the samples "
               "as well" << std::endl;
    tempOut << "-popPar [option] : Separate population paramters for "
               "clustering will default to the parameters given for the "
               "between sample clustering" << std::endl;

    tempOut << bib::bashCT::bold << "Additional Pre-processing options"<< bib::bashCT::reset << std::endl;
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
               "-par pars.txt -population" << std::endl;
    tempOut << "SeekDeep processClusters -fasta output.fasta -par "
               "pars.txt -popPar otherPars.txt -population" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    tempOut << bib::bashCT::bold << "Output Files:"<< bib::bashCT::reset << std::endl;
    tempOut << bib::bashCT::boldBlack("selectedClustersInfo.tab.txt") << ": This contains the final haplotype information and replicate comparison results, it is a very large table, consult the SeekDeep (http://bib2.umassmed.edu/~hathawan/SeekDeep.html) website for details on what each column means" << std::endl;
    tempOut << bib::bashCT::boldBlack("allClustersInfo.tab.txt") << ": This table contains information on the final clustering process before final filtering is done (removal of chimrias, frac cut off, etc" << std::endl;
    tempOut << bib::bashCT::boldBlack("dotFiles") << ": This contains .dot files for the final haplotypes that can be converted by graphviz to a pdf of a psuedo-minimum spanning tree of the differences in the final haplotypes (cmd: neato -T pdf -o [FILE.dot]" << std::endl;
    tempOut << bib::bashCT::boldBlack("final") << ": This contains the final consensus of the final haplotypes after final filtering has been done" << std::endl;
    tempOut << "If population clustering was done there will be a directly called population" << std::endl;
    tempOut << bib::bashCT::boldBlack("population/populationCluster.tab.txt") << ": This contains information on final haplotypes across the population (number of samples it appears in, etc.) see SeekDeep website (http://bib2.umassmed.edu/~hathawan/SeekDeep.html) for detail on all columns " << std::endl;
    tempOut << bib::bashCT::boldBlack("population/PopUID.fastq") << ": The consensus sequences of the final haplotypes in the population, the prefix name is determine by the -experimentName flag (default PopUID), haplotypes are sorted by the number of samples they appear in, suffic is _f[FRAC] where [FRAC] is the average fraction it appears at in the samples" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    tempOut.str(std::string());
    exit(0);
  }
  setOption(onPerId, "-onPerId", "Cluster on Percent Identity Instead");
  // input file info
  ioOptions_.lowerCaseBases_ = "upper";
  processDefaultReader(true);
  setOption(parameters, "-par", "ParametersFileName", true);
  processAlnInfoInput();
  setOption(runsRequired, "-runsRequired", "runsRequired");
  setOption(experimentName, "-experimentName", "ExperimentName");
  setOption(cutOff, "-cutoff", "Cluster Size Cut Off");
  processDirectoryOutputName("clusters_" + getCurrentDate(), true);
  processVerbose();
  processScoringPars();
  setOption(extra, "-extra", "ExtraOutput");
  processRefFilename();
  setOption(population, "-population", "PopulationClustering");
  setOption(fracCutoff, "-fraccutoff", "PopulationClusteringFractionCutoff");
  if (setOption(parametersPopulation, "-poppar",
                "ParametersForPopulationCollapse")) {
    differentPar = true;
  }

  if (!setOption(popBoth, "-popBoth", "PopBoth")) {
    setOption(popBoth, "-popAll", "PopAll");
  }

  setBoolOptionFalse(smallestFirst, "-largestFirst",
                     "Check_smallestClustersFirst");
  setOption(bestMatch, "-bestMatch", "FindingBestMatch");
  setOption(bestMatchCheck, "-bestMatchCheck", "BestMatchCheck");
  setOption(checkChimeras, "-markChimeras", "MarkChimeras");
  setOption(keepChimeras, "-keepChimeras", "KeepChimeras");
  setOption(parFreqs, "-parfreqs", "ParentFrequence_multiplier_cutoff");

  processKmerOptions();
  // get the qualities
  processQualThres();
  processGap();
  if (!failed_ && !commands_.containsFlagCaseInsensitive("-flags") &&
      !commands_.containsFlagCaseInsensitive("-getFlags")) {
    std::cout << "p: " << primaryQual_ << std::endl;
    std::cout << "s: " << secondaryQual_ << std::endl;
    std::cout << "go: " << gapInfo_.gapOpen_ << std::endl;
    std::cout << "ge: " << gapInfo_.gapExtend_ << std::endl;
  }
  // read in the paramteres from the parameters file
  finishSetUp(std::cout);
  if(onPerId){
  	 processIteratorMapOnPerId(parameters, iteratorMap);
  }else{
  	 processIteratorMap(parameters, iteratorMap);
  }
  if(differentPar){
    if(onPerId){
    	processIteratorMapOnPerId(parametersPopulation, popIteratorMap);
    }else{
    	processIteratorMap(parametersPopulation, popIteratorMap);
    }
  }else{
  	popIteratorMap = iteratorMap;
  }
}



void SeekDeepSetUp::setUpMakeSampleDirectories(
    std::string& sampleNameFilename) {
  if (needsHelp()) {
  	std::stringstream tempStream;
    tempStream << "makeSampleDirectoires" << std::endl;
    tempStream << "Set up a directory tree for processClusters" << std::endl;
    tempStream << "Commands, order not necessary, flags are case insensitive" << std::endl;
    tempStream << "Required commands" << std::endl;
    tempStream << "-file [option], name of the file of sample names to read in"
              << std::endl;
    tempStream << "-dout [option], name of the main directory to create"
              << std::endl;
    tempStream << "File should be tab delimited and a few examples are bellow" << std::endl;
    tempStream << "File should have at least three columns" << std::endl;
    tempStream << "Where first column is the name of the index or sff file "
                 "used, second column is the sample names, and all following "
                 "columns are the MIDs for that samples " << std::endl;
    std::cout << cleanOut(tempStream.str(), width_, indent_	);
    tempStream.str(std::string());
    std::cout << "Example with two replicates and two separate master indexes" << std::endl;
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
  setOption(sampleNameFilename, "-file", "SampleNamesFileName", true);
  setOption(directoryName_, "-dout", "MainOutDirectoryName", true);
  processDirectoryOutputName(true);
  finishSetUp(std::cout);
}

}  // namespace bib
