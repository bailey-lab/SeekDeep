/*
 * clusterDownSetUp.cpp
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
//
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {

void SeekDeepSetUp::setUpClusterDown(clusterDownPars & pars) {
	// input file info
	pars_.ioOptions_.out_.outFilename_ = "output";
	pars_.ioOptions_.lowerCaseBases_ = "remove";

	if (false) {
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
		//printInputUsage(tempOut);
		tempOut << "--par [option]: parameter file" << std::endl << std::endl;
		//" | "
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
				<< "Alternatively if the --onPerId is used, this parameter flag can contain percent identities instead"
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
		tempOut << "--dout [option]: Name of an output directory, will default "
				"to the stub name plus the date" << std::endl;
		tempOut << "--additionaOut [option]: Name of a tab delimited file with "
				"two columns, column one is the name of the ID and column two "
				"is the file path location of additional output directory"
				", such a file can be made by makeSampleDirectories" << std::endl;
		tempOut
				<< "--qualRep [option] : Sets the quality for the identical clusters, "
						"Options are median (default), average, bestQual, or worst"
				<< std::endl;
		tempOut << "--qualTrim [option]: Will trim off any bases below this "
				"quality, mainly used to trim bases with a quality of 1 (80% "
				"chance of error) by setting this to 2" << std::endl;
		tempOut << "--adjustHomopolyerRuns : Will take average quality across "
				"a homopolymer run and set all the quality to this average, "
				"mainly used with Ion Torrent reads due to the very low "
				"quality of the last base in a long homopolymer run" << std::endl;
		tempOut << bib::bashCT::bold << "Additional Processing options"
				<< bib::bashCT::reset << std::endl;
		tempOut << "--markChimeras : Have the program mark possible chimeras "
				"before outputting" << std::endl;
		tempOut << "--parFreqs [option] : The minimum frequency multiplier "
				"reads must have in order to be considered to be a possible "
				"parent of a chimera or for sequences to be collapsed on gaps "
				"in tandem repeats, defaults to 2" << std::endl;
		tempOut << bib::bashCT::bold << "Technology Specific flags"
				<< bib::bashCT::reset << std::endl;
		tempOut
				<< "--useCompPerCutOff: This turns on filtering off clustering that have reads that only come from one direction, only should be used if reads were extracted in both direction like in ion torrent"
				<< std::endl;
		tempOut
				<< "--ionTorrent : This turns on several of the previously mentioned flags as they are beneficial for ion torrent data, turns on --qualTrim,--adjustHomopolyerRuns, and --useCompPerCutOff"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Percent Identity Falgs"
				<< bib::bashCT::reset << std::endl;
		tempOut
				<< "--onPerId: cluster on percent identity rather than specific errors"
				<< std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		std::cout << "examples" << std::endl;
		std::cout << "\tSeekDeep " << commands_.subProgram_
				<< " --fasta MID1.fasta --qual "
						"MID1.fasta.qual --par par" << std::endl;
		std::cout << "\tSeekDeep " << commands_.subProgram_
				<< " --fastq MID1.fastq --par par --local" << std::endl;
		std::cout << "\tSeekDeep " << commands_.subProgram_
				<< " --stub MID1 --par par --local "
						"--qualThres 25,20" << std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		tempOut << bib::bashCT::bold << "Output Files:" << bib::bashCT::reset
				<< std::endl;
		tempOut
				<< "outputDirectory/output.fastq: This is the final clusters with their consensus sequence, the sequences are named so that the suffix _t[NUM] where NUM is the number of reads that fell into that cluster"
				<< std::endl;
		tempOut
				<< "outputDirectory/outputInfo.tab.txt: This contains cluster number info for each cluster"
				<< std::endl;
		std::cout << cleanOut(tempOut.str(), width_, indent_);
		tempOut.str(std::string());
		exit(0);
	}
	description_ = "Cluster input sequences by collapsing on set differences between the sequences";
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq inputSeqs.fastq --par parFile.txt");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq inputSeqs.fastq --illumina   #use default Illumina parameters for collapsing instead of supplying a parameters file");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq inputSeqs.fastq --ionTorrent #use default IonTorrent parameters for collapsing instead of supplying a parameters file");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq inputSeqs.fastq --454        #use default 454 parameters for collapsing instead of supplying a parameters file");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq inputSeqs.fastq --otu .99    #use default parameters to OTU clustering at an OTU threshold of .99 instead of supplying a parameters file");
	processVerbose();
	processDebug();
	//processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz", "--fastq1", "--fastq2", "--fastq1gz", "fastq2gz"}, true);
	processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz"}, true);

	setOption(pars.binParameters, "--binPar",
			"Parameters Filename for when reads are binned");
	setOption(pars.binParameters, "--initialPar",
				"Parameters Filename for when doing non singlet analysis first, otherwise defaults to regular parameters");
	setOption(pars.ionTorrent, "--ionTorrent", "Flag to indicate reads are IonTorrent and therefore turns on --adjustHomopolyerRuns and --qualTrim", false, "Technology");
	setOption(pars.tech454, "--454", "Flag to indicate reads are 454", false, "Technology");
	setOption(pars.illumina, "--illumina","Flag to indicate reads are Illumina", false, "Technology");
	bool needsParFlag = true;

	setOption(pars.hq, "--hq",
			"When the --illumina,--454,or --ionTorrent flag is on also allow this many high quality mismatch to the defaults for those technology", false, "Technology");
	if (pars.ionTorrent) {
		pars_.colOpts_.iTOpts_.removeLowQualityBases_ = true;
		pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_ = true;
		//pars.useCompPerCutOff = true;
		needsParFlag = false;
	}

	if((pars.tech454 || pars.ionTorrent) && pars.illumina){
		failed_ = true;
		std::stringstream ss;
		if(pars.tech454){
			ss << "Error, both --454 and --illumina where used" << "\n";
		}else{
			ss << "Error, both --ionTorrent and --illumina where used" << "\n";
		}
		addWarning(ss.str());
	}

	if(pars.tech454 || pars.ionTorrent){
		if(pars.hq  > 0){
			pars.binIteratorMap = CollapseIterations::gen454ItDefaultParsWithHqs(100, pars.hq);
			pars.intialParameters = CollapseIterations::gen454ItDefaultParsWithHqs(100, 0);
			pars.iteratorMap = CollapseIterations::gen454ItDefaultParsWithHqs(100, pars.hq);
		}else{
			pars.binIteratorMap = CollapseIterations::gen454ItDefaultPars(100);
			pars.intialParameters = CollapseIterations::gen454ItDefaultPars(100);
			pars.iteratorMap = CollapseIterations::gen454ItDefaultPars(100);
		}
		needsParFlag = false;
	}

	if(pars.illumina){
		if(pars.hq  > 0){
			pars.binIteratorMap = CollapseIterations::genIlluminaDefaultParsWithHqs(100, pars.hq);
			pars.intialParameters = CollapseIterations::genIlluminaDefaultParsWithHqs(100, 0);
			pars.iteratorMap = CollapseIterations::genIlluminaDefaultParsWithHqs(100, pars.hq);
		}else{
			pars.binIteratorMap =  CollapseIterations::genIlluminaDefaultPars(100);
			pars.intialParameters = CollapseIterations::genIlluminaDefaultPars(100);
			pars.iteratorMap = CollapseIterations::genIlluminaDefaultPars(100);
		}
		pars_.colOpts_.iTOpts_.weighHomopolyer_ = false;
		needsParFlag = false;
	}
	bool otuSet = setOption(pars.otuPerc, "--otu",
			"Collapse on this OTU percentage, should be between (0,1)", false, "OTU Clustering");
	bool regOtu = false;
	setOption(regOtu, "--regOtu",
				"Collapse on percent identity calculated using regular otu percent identity rather than using only high quality differences",
				false, "OTU Clustering");
	pars.onHqPerId = !regOtu;
	if (otuSet) {
		if ((pars.otuPerc <= 0) || (pars.otuPerc >= 1)) {
			failed_ = true;
			std::stringstream ss;
			ss
					<< "Error, otu percentage should be between 0 and 1 (non inclusive), not"
					<< pars.otuPerc << "\n";
			addWarning(ss.str());
		}
		pars.onPerId = true;
		pars.binIteratorMap = CollapseIterations::genOtuPars(100, pars.otuPerc, pars.onHqPerId);
		pars.iteratorMap = CollapseIterations::genOtuPars(100, pars.otuPerc, pars.onHqPerId);
		pars.intialParameters = CollapseIterations::genOtuPars(100, pars.otuPerc, pars.onHqPerId);
		needsParFlag = false;
		pars_.colOpts_.kmerOpts_.runCutOff_ = 1;
		pars_.colOpts_.kmerOpts_.runCutOffString_ = "1";
	}


	setOption(pars.parameters, "--par", "Parameters filename to supply differences to allow to cluster, required unless technology flag used", needsParFlag, "Clustering");
	setOption(pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_, "--adjustHomopolyerRuns",
				"Adjust homopolymer runs to be same quality scores");

	setOption(pars.compPerCutOff, "--compPerCutOff",
			"Percentage of reads in one direction Cut Off");
	setOption(pars.useCompPerCutOff, "--useCompPerCutOff",
			"Throw out clusters that are made up reads of > "
					+ estd::to_string(pars.useCompPerCutOff * 100)
					+ "% of reads in only one direction, used with Ion Torrent Reads");

	setOption(pars.diffCutOffStr, "--diffCutOffs",
			"Difference Cutoff to form Nuc Comp Clusters when --useNucComp is used", false, "Clustering");

	pars_.colOpts_.nucCompBinOpts_.diffCutOffVec_ = vecStrToVecNum<double>(
			tokenizeString(pars.diffCutOffStr, ","));

	setOption(pars_.colOpts_.nucCompBinOpts_.findBestNuc_, "--findBestNuc", "Find Best Nucleotide Cluster when --useNucComp is used", false, "Clustering");
	setOption(pars_.colOpts_.nucCompBinOpts_.useNucComp_, "--useNucComp",
			"Cluster on Nucleotide Composition First", false, "Clustering");

	setOption(pars_.colOpts_.nucCompBinOpts_.useMinLenNucComp_, "--useMinLenNucComp",
			"Use Nucleotide Composition of Only Front of seqs");

	setOption(pars.startWithSingles, "--startWithSingles",
			"Start The Clustering With Singletons, rather then adding them afterwards", false, "Clustering");
//	setOption(pars.mapBackSinglets, "--mapBackSinglets",
//				"If Singlets are left out, map them back in for frequency estimates");
	setOption(pars.singletCutOff, "--singletCutOff",
				"Naturally the cut off for being a singlet is by default 1 but can use --singletCutOff to raise the number", false, "Clustering");
	setOption(pars.createMinTree, "--createMinTree",
			"Create Pseudo minimum Spanning Trees For Mismatches for Final Clusters", false, "Additional Output");
	setOption(pars_.colOpts_.kmerBinOpts_.useKmerBinning_, "--useKmerBinning", "Use Kmer Binning for initial clustering to speed up clustering", false, "Clustering");
	setOption(pars_.colOpts_.kmerBinOpts_.kmerCutOff_, "--kmerCutOff", "kmer Cut Off for when --useKmerBinning is used", false, "Clustering");
	setOption(pars_.colOpts_.kmerBinOpts_.kCompareLen_, "--kCompareLen", "kmer Compare Length for when bining by kmers first for when --useKmerBinning is used", false, "Clustering");
	setOption(pars.leaveOutSinglets, "--leaveOutSinglets",
			"Leave out singlet clusters out of all analysis", false, "Clustering");
	setOption(pars.onPerId, "--onPerId", "Cluster on Percent Identity Instead", false, "OTU Clustering");
	pars_.colOpts_.iTOpts_.removeLowQualityBases_= setOption(pars_.colOpts_.iTOpts_.lowQualityBaseTrim_, "--qualTrim",
			"Low Quality Cut Off", false, "Preprocessing");
	setOption(pars.qualRep, "--qualRep",
			"Per base quality score calculation for initial unique clusters collapse", false, "Preprocessing");
	//setOption(pars.extra, "--extra", "Extra");
	setOption(pars.writeOutFinalInternalSnps, "--writeOutFinalInternalSnps", "Write out Internal (within the clusters) SNP class, useful for debugging if over collapsing is happening", false, "Additional Output");

	pars_.chiOpts_.checkChimeras_ = true;
	pars_.chiOpts_.parentFreqs_ = 2;
	bool dontMarkChimeras = false;
	setOption(dontMarkChimeras, "--noMarkChimeras", "Don't do chimera marking", false, "Chimeras");
	pars_.chiOpts_.checkChimeras_ = !dontMarkChimeras;
	setOption(pars_.chiOpts_.parentFreqs_, "--parFreqs",
			"Parent freq multiplier cutoff", false, "Chimeras");
	setOption(pars_.chiOpts_.keepLowQaulityMismatches_, "--chiKeepLowQaulityMismatches", "When marking chimeras also consider low quality mismatches, these are removed by default to prevent errors from incorrectly mising chimeras but could also hurt chimera marking", false, "Chimeras");

	setOption(pars.snapShotsOpts_.snapShots_, "--snapShots", "Output Snap Shots of clustering results after each iteration", false, "Additional Output");
	setOption(pars.sortBy, "--sortBy", "Sort Clusters By");
	pars.additionalOut = setOption(pars.additionalOutLocationFile,
			"--additionalOut", "Additional out filename for sorting final results", false, "Additional Output");
	//setOption(pars.collapsingTandems, "--collapseTandems", "Collapsing Tandems");

	setOption(pars_.colOpts_.alignOpts_.noAlign_, "--noAlignCompare",
			"Do comparisons without globally aligning", false, "Alignment");
	processSkipOnNucComp();
	setOption(pars_.colOpts_.clusOpts_.converge_, "--converge", "Keep clustering at each iteration until there is no more collapsing, could increase run time significantly", false, "Clustering");
	setOption(pars.writeOutInitalSeqs, "--writeOutInitalSeqs", "Write out the sequences that make up each cluster", false, "Additional Output");
	pars_.colOpts_.verboseOpts_.verbose_ = pars_.verbose_;
	pars_.colOpts_.verboseOpts_.debug_ = pars_.debug_;
	processRefFilename();
	bool mustMakeDirectory = true;
	processDirectoryOutputName(mustMakeDirectory);
	processWritingOptions();
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "5,1";
	processAlignerDefualts();
	pars.smallReadSize = pars_.colOpts_.kmerOpts_.kLength_ * 2;
	setOption(pars.smallReadSize, "--smallReadSize",
			"A cut off to remove small reads", false, "Preprocessing");
	if(pars.smallReadSize <= pars_.colOpts_.kmerOpts_.kLength_){
		pars.smallReadSize = 	pars_.colOpts_.kmerOpts_.kLength_ + 1;
	}
	if (!failed_) {
		if ("" != pars.parameters) {
			if (pars.onPerId) {
				pars.iteratorMap = processIteratorMapOnPerId(pars.parameters);
			} else {
				pars.iteratorMap = processIteratorMap(pars.parameters);
			}
		}
		if ("" != pars.binParameters) {
			if (pars.onPerId) {
				pars.binIteratorMap = processIteratorMapOnPerId(pars.binParameters);
			} else {
				pars.binIteratorMap = processIteratorMap(pars.binParameters);
			}
		} else if("" != pars.parameters){
			//set to reg iter map if it was set by --par flag
			pars.binIteratorMap = pars.iteratorMap;
		}
		if ("" != pars.initalParsFnp) {
			if (pars.onPerId) {
				pars.intialParameters = processIteratorMapOnPerId(pars.initalParsFnp.string());
			} else {
				pars.binIteratorMap = processIteratorMap(pars.initalParsFnp.string());
			}
		} else if("" != pars.parameters){
			//set to reg iter map if it was set by --par flag
			pars.intialParameters = pars.iteratorMap;
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

}   // namespace bibseq

