/*
 * processClustersSetUp.cpp
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
				<< "in the previous example tree the flag would be --fastq output.fastq and all seq "
						"files should be named output.fastq, others will be ignored"
				<< std::endl;
		tempOut
				<< "--par : parameters file, see SeekDeep qluster --help for details on the format of the file"
				<< std::endl;
		tempOut << bib::bashCT::bold << "Optional commands" << bib::bashCT::reset
				<< std::endl;
		tempOut << "--dout [option]: Name of an output directory, will default "
				"to cluster_CURRENT_DATE" << std::endl;
		tempOut << "--onPerId : Cluster reads on percent identity rather"
				" than specific errors, format of parameters file would have to change,"
				" see SeekDeep qluster --help for details" << std::endl;
		tempOut
				<< "--experimentName [option] : Name prefix to give to the final population haplotypes, defaults to PopUID"
				<< std::endl;

		tempOut << bib::bashCT::boldBlack("Final analysis inclusion criteria")
				<< std::endl;
		tempOut << "--fracCutOff [option] : The fraction threshold to be "
				"included in final analysis, threshold is compared to the averaged "
				"fraction across runs, defaults to 0.005 " << std::endl;
		tempOut << "--runsRequired [option] : Number of runs a cluster has to "
				"appear in to be"
				" included in final analysis, defaults to all the runs included "
				"in the sample" << std::endl;
		tempOut << bib::bashCT::bold << "Population clustering options"
				<< bib::bashCT::reset << std::endl;
		tempOut << "--popPar [option] : Separate population paramters for "
				"clustering will default to the parameters given for the "
				"between sample clustering" << std::endl;

		tempOut << bib::bashCT::bold << "Additional Pre-processing options"
				<< bib::bashCT::reset << std::endl;
		tempOut << "--clusterCutOff [option]: Size of cluster not to include in "
				"clustering, defaults to 1" << std::endl;
		tempOut << "--markChimeras : Have the program mark possible chimeras "
				"before starting clustering" << std::endl;
		tempOut << "--parFreqs [option] : The minimum frequency multiplier "
				"reads must have in order to be considered to be a possible "
				"parrent of a chimera, defaults to 2" << std::endl;
		tempOut << "--keepChimeras : Whether to include chimeras in the final "
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
		tempOut << "SeekDeep processClusters --fasta output.fasta "
				"--par pars.txt" << std::endl;
		tempOut << "SeekDeep processClusters --fasta output.fasta "
				"--par pars.txt" << std::endl;
		tempOut << "SeekDeep processClusters --fasta output.fasta --par "
				"pars.txt --popPar otherPars.txt" << std::endl;
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

	setOption(pars.illumina, "--illumina",
			"Flag to indicate reads are Illumina");

	setOption(pars.ionTorrent, "--ionTorrent",
			"Flag to indicate reads are ion torrent and therefore turns on --useCompPerCutOff,--adjustHomopolyerRuns, and --qualTrim");

	if (pars.ionTorrent) {
		pars.removeLowQualBases = true;
		pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_ = true;
	}

	setOption(pars.masterDir, "--masterDir", "Master directory containing sample sequence files");
	pars.masterDir = bfs::absolute(pars.masterDir);
	pars.removeLowQualBases = setOption(pars.lowQualityCutOff, "--qualTrim",
			"LowQualityCutOff");
	setOption(pars.writeExcludedOriginals, "--writeExcludedOriginals",
			"Write out the excluded Originals");
	setOption(pars.chiCutOff, "--chiCutOff",
			"The Fraction of a cluster to determine if it chimeric");
	setOption(pars.recheckChimeras, "--recheckChimeras",
			"Re Check chimeras after replicate comparison");
	setOption(pars.eventBasedRef, "--eventBasedRef", "Do Event Based Ref Count");
	setOption(pars.customCutOffs, "--custumCutOffs",
			"Two Column Table, first column is sample name, second is a custom frac cut off, if sample not found will default to --fracCutOff");
	setOption(pars.previousPopFilename, "--previousPop", "previousPopFilename");
	processComparison(pars.previousPopErrors, "previousPop");
	setOption(pars.groupingsFile, "--groupingsFile",
			"A file to sort samples into different groups");
	setOption(pars.investigateChimeras, "--investigateChimeras",
			"Check to see if a chimera appears as a high variant in another sample");
	processDebug();
	processVerbose();
	pars_.colOpts_.verboseOpts_.verbose_ = pars_.verbose_;
	pars_.colOpts_.verboseOpts_.debug_ = pars_.debug_;
	setOption(pars.onPerId, "--onPerId", "Cluster on Percent Identity Instead");
	processSkipOnNucComp();
	// input file info
	pars_.ioOptions_.lowerCaseBases_ = "upper";
	pars_.ioOptions_.processed_ = true;

	processDefaultReader(true);

	setOption(pars.noErrorsSet, "--noErrors", "Collapse parameters with no errors");


	setOption(pars.strictErrorsSetHq1, "--strictErrors-hq1", "Collapse parameters with a several low quality mismatches and 1 high quality mismatch");
	if(pars.strictErrorsSetHq1){
		pars.strictErrorsSet = true;
		pars.hqMismatches = 1;
	}
	setOption(pars.strictErrorsSet, "--strictErrors", "Collapse parameters with a several low quality mismatches");
	setOption(pars.hqMismatches, "--hq", "Number of high quality mismatches to allow");
	setOption(pars.stopAfter, "--stopAfter", "Number of top haplotypes to check");

	setOption(pars.parameters, "--par", "ParametersFileName", !pars.noErrorsSet && !pars.strictErrorsSet && !pars.strictErrorsSetHq1);

	setOption(pars.binParameters, "--binPar", "bin Parameters Filename");
	setOption(pars.sampleMinTotalReadCutOff, "--sampleMinTotalReadCutOff", "Sample Minimum Total Read Cut Off, if the total read count for the sample is below this it will be thrown out");

	setOption(pars.runsRequired, "--runsRequired", "Number of PCR runs Required for a haplotype to be kept");
	setOption(pars.experimentName, "--experimentName", "ExperimentName");
	if (bib::containsSubString(pars.experimentName, ".")) {
		addWarning("Error in populationCollapse::populationCollapse, populationName can't contain '.', "
						+ pars.experimentName);
		failed_ = true;
	}
	setOption(pars.clusterCutOff, "--clusterCutOff", "Cluster Size Cut Off");
	setOption(pars.noTrees, "--noTrees", "Don't generate html difference trees");
	processDirectoryOutputName("clusters_" + getCurrentDate(), true);

	setOption(pars.extra, "--extra", "Extra Output");
	processRefFilename();
	setOption(pars.noPopulation, "--noPopulation",
			"Don't do Population Clustering");

	setOption(pars.fracCutoff, "--fracCutOff",
			"PopulationClusteringFractionCutoff");
	pars.differentPar = setOption(pars.parametersPopulation, "--popPar",
			"Parameters For Population Collapse");
	struct ClusteringParametersPars {
		bool noErrorsSet = false;
		bool strictErrorsSet = false;
		uint32_t stopAfter = 100;
		uint32_t hqMismatches = 0;
	};
	ClusteringParametersPars popClusParsPars;
	setOption(popClusParsPars.noErrorsSet, "--pop-noErrors", "Collapse parameters with no errors in population clustering");
	setOption(popClusParsPars.strictErrorsSet, "--pop-strictErrors", "Collapse parameters with a several low quality mismatches in population clustering");
	setOption(popClusParsPars.hqMismatches, "--pop-hq", "Number of high quality mismatches to allow in population clustering");
	setOption(popClusParsPars.stopAfter, "--pop-stopAfter", "Number of top haplotypes to check in population clustering");

	setOption(pars_.chiOpts_.checkChimeras_, "--markChimeras", "Mark Chimeras");
	setOption(pars.keepChimeras, "--keepChimeras", "KeepChimeras");
	setOption(pars_.chiOpts_.parentFreqs_, "--parFreqs", "Chimeric Parent Frequency multiplier cutoff");

	setOption(pars.numThreads, "--numThreads", "Number of threads to use");
	setOption(pars_.colOpts_.clusOpts_.converge_, "--converge", "Keep clustering at each iteration until there is no more collapsing, could increase run time significantly");
	//setOption(pars.plotRepAgreement, "--plotRepAgreement", "Plot Rep Agreement");
	processAlignerDefualts();
	if (pars_.debug_ && !failed_) {
		std::cout << "p: " << pars_.qScorePars_.primaryQual_ << std::endl;
		std::cout << "s: " << pars_.qScorePars_.secondaryQual_ << std::endl;
		std::cout << "go: " << pars_.gapInfo_.gapOpen_ << std::endl;
		std::cout << "ge: " << pars_.gapInfo_.gapExtend_ << std::endl;
	}
	if (!failed_) {
		if ("" != pars.parameters) {
			// read in the parameters from the parameters file
			if (pars.onPerId) {
				pars.iteratorMap = processIteratorMapOnPerId(pars.parameters);
			} else {
				pars.iteratorMap = processIteratorMap(pars.parameters);
			}
		} else if (pars.noErrorsSet) {
			if (pars.hqMismatches > 0) {
				pars.iteratorMap =
						CollapseIterations::genStrictNoErrorsDefaultParsWithHqs(
								pars.stopAfter, pars.hqMismatches);
			} else {
				pars.iteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(
						pars.stopAfter);
			}
		} else if (pars.strictErrorsSet) {
			pars.iteratorMap = CollapseIterations::genStrictDefaultParsWithHqs(
					pars.stopAfter, pars.hqMismatches, pars.illumina);
		}

		if (pars.binParameters != "") {
			if (pars.onPerId) {
				pars.binIteratorMap = processIteratorMapOnPerId(pars.binParameters);
			} else {
				pars.binIteratorMap = processIteratorMap(pars.binParameters);
			}
		}
		if (pars.differentPar || popClusParsPars.noErrorsSet
				|| popClusParsPars.strictErrorsSet) {
			if (popClusParsPars.noErrorsSet) {
				if (popClusParsPars.hqMismatches > 0) {
					pars.popIteratorMap =
							CollapseIterations::genStrictNoErrorsDefaultParsWithHqs(
									popClusParsPars.stopAfter, popClusParsPars.hqMismatches);
				} else {
					pars.popIteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(
							popClusParsPars.stopAfter);
				}
			} else if (popClusParsPars.strictErrorsSet) {
				pars.popIteratorMap = CollapseIterations::genStrictDefaultParsWithHqs(
						popClusParsPars.stopAfter, popClusParsPars.hqMismatches, pars.illumina);
			} else {
				if (pars.onPerId) {
					pars.popIteratorMap = processIteratorMapOnPerId(
							pars.parametersPopulation);
				} else {
					pars.popIteratorMap = processIteratorMap(pars.parametersPopulation);
				}
			}

		} else {
			pars.popIteratorMap = pars.iteratorMap;
		}
	}
	finishSetUp(std::cout);

}


}  // namespace bibseq
