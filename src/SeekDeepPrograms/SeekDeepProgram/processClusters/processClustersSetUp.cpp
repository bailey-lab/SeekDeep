/*
 * processClustersSetUp.cpp
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
#include <njhcpp/bashUtils.h>

namespace njhseq {

void SeekDeepSetUp::setUpMultipleSampleCluster(processClustersPars & pars) {
	// parse the command line options
	description_ = "Ran from inside directory tree set up such that currentDir/SampleDir/RunDir/seqFiles";
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fasta output.fasta --par pars.txt");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fastq output.fastq --strictErrors #just collapse on a few low quality errors ");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fasta output.fasta --strictErrors --pop-noErrors #collapse across replicates within sample on low quality errors but allow no errors in population compare");
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --fasta output.fasta --strictErrors -hq 1 --pop-noErrors #collapse across replicates within sample on 1 high quality error but allow no errors in population compare");
	if (needsHelp()) {
		commands_.arguments_["-h"] = "";
	}

	setOption(pars.illumina, "--illumina",
			"Flag to indicate reads are Illumina", false, "Technology");

	setOption(pars.ionTorrent, "--ionTorrent",
			"Flag to indicate reads are ion torrent and therefore turns on --adjustHomopolyerRuns and --qualTrim",
			false, "Technology");

	if (pars.ionTorrent) {
		pars.removeLowQualBases = true;
		pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_ = true;
	}

	setOption(pars.masterDir, "--masterDir", "Master directory containing sample sequence files", false, "Input");
	pars.masterDir = bfs::absolute(pars.masterDir);
	pars.removeLowQualBases = setOption(pars.lowQualityCutOff, "--qualTrim",
			"Low Quality Cut Off", false, "Pre Processing");
	setOption(pars.writeExcludedOriginals, "--writeExcludedOriginals",
			"Write out the excluded Originals", false, "Additional Output");
	setOption(pars.chiCutOff, "--chiCutOff",
			"The Fraction of a cluster to determine if it chimeric", false, "Chimeras");
	setOption(pars.recheckChimeras, "--recheckChimeras",
			"Re Check chimeras after replicate comparison", false, "Chimeras");
	setOption(pars.eventBasedRef, "--eventBasedRef", "Do Event Based Ref Count");
	setOption(pars.customCutOffs, "--custumCutOffs",
			"Two Column Table, first column is sample name, second is a custom frac cut off, if sample not found will default to --fracCutOff", false, "Filtering");
	setOption(pars.previousPopFilename, "--previousPop", "previousPopFilename", false, "Population");
	processComparison(pars.previousPopErrors, "previousPop");
	setOption(pars.groupingsFile, "--groupingsFile",
			"A file to sort samples into different groups", false, "Meta");
	setOption(pars.investigateChimeras, "--investigateChimeras",
			"Check to see if a chimera appears as a high variant in another sample", false, "Chimeras");
	processDebug();
	processVerbose();
	pars_.colOpts_.verboseOpts_.verbose_ = pars_.verbose_;
	pars_.colOpts_.verboseOpts_.debug_ = pars_.debug_;
	setOption(pars.onPerId, "--onPerId", "Cluster on Percent Identity Instead", false, "Clustering");
	processSkipOnNucComp();
	// input file info
	pars_.ioOptions_.lowerCaseBases_ = "upper";
	pars_.ioOptions_.processed_ = true;

	processReadInNames(VecStr{"--fasta", "--fastagz", "--fastq", "--fastqgz"}, true);

	setOption(pars.noErrorsSet, "--noErrors", "Collapse parameters with no errors", false, "Clustering");


	setOption(pars.strictErrorsSetHq1, "--strictErrors-hq1", "Collapse parameters with a several low quality mismatches and 1 high quality mismatch", false, "Clustering");
	if(pars.strictErrorsSetHq1){
		pars.strictErrorsSet = true;
		pars.hqMismatches = 1;
	}
	setOption(pars.strictErrorsSet, "--strictErrors", "Collapse parameters with a several low quality mismatches", false, "Clustering");
	setOption(pars.hqMismatches, "--hq", "Number of high quality mismatches to allow", false, "Clustering");
	setOption(pars.stopAfter, "--stopAfter", "Number of top haplotypes to check", false, "Clustering");

	setOption(pars.parameters, "--par", "ParametersFileName", !pars.noErrorsSet && !pars.strictErrorsSet && !pars.strictErrorsSetHq1, "Clustering");

	setOption(pars.binParameters, "--binPar", "bin Parameters Filename", false, "Clustering");
	setOption(pars.preFiltCutOffs.sampleMinReadCount, "--sampleMinTotalReadCutOff",
			"Sample Minimum Total Read Cut Off, if the total read count for the sample is below this it will be thrown out", false, "Filtering");
	setOption(pars.preFiltCutOffs.replicateMinReadCount, "--replicateMinTotalReadCutOff",
				"Replicate Minimum Total Read Cut Off, if the total read count for the replicate is below this it will be thrown out", false, "Filtering");

	setOption(pars.runsRequired, "--runsRequired", "Number of PCR runs Required for a haplotype to be kept", false, "Filtering");
	setOption(pars.experimentName, "--experimentName", "Name given to the final population haplotypes", false, "Population");
	if (njh::containsSubString(pars.experimentName, ".")) {
		addWarning("Error in populationCollapse::populationCollapse, populationName can't contain '.', "
						+ pars.experimentName);
		failed_ = true;
	}
	setOption(pars.preFiltCutOffs.clusterSizeCutOff, "--clusterCutOff", "Input Cluster Size Cut Off", false, "Filtering");
	setOption(pars.fracExcludeOnlyInFinalAverageFrac, "--fracExcludeOnlyInFinalAverageFrac", "By default fraction exclusion is done per rep, use fracExcludeOnlyInFinalAverageFrac to exclude only on the final averaged frac", false, "Filtering");


	setOption(pars.removeCommonlyLowFreqHaplotypes_, "--excludeCommonlyLowFreqHaplotypes", "Remove Commonly Low Freq Haplotypes", false, "Filtering");
	setOption(pars.lowFreqHaplotypeFracCutOff_, "--lowFreqHaplotypeFracCutOff", "Low Freq Haplotype Frac Cut Off", false, "Filtering");



	setOption(pars.noTrees, "--noTrees", "Don't generate html difference trees");
	processDirectoryOutputName("clusters_" + getCurrentDate(), true);

	setOption(pars.extra, "--extra", "Extra Output", false, "Additional Output");
	processRefFilename();
	setOption(pars.noPopulation, "--noPopulation",
			"Don't do Population Clustering", false, "Population");

	setOption(pars.controlSamples, "--controlSamples", "Samples that shouldn't be included in frequency filtering calcs", false, "Filtering");


	setOption(pars.collapseLowFreqOneOffs, "--excludeLowFreqOneOffs",
			"Collapse any haplotypes that are low frequency compared to another haplotype (determined by lowFreqMultiplier) and only differs by 1 base", false, "Filtering");
	setOption(pars.lowFreqMultiplier, "--lowFreqMultiplier",
			"Low Freq Multiplier used for --excludeLowFreqOneOffs, considered low frequency if haplotype frac is less than its fraction times this number than the other haplotype", false, "Filtering");

	setOption(pars.rescueExcludedChimericHaplotypes, "--rescueExcludedChimericHaplotypes", "Rescue Excluded chimeric Haplotypes if they appear as a major haplotype in another sample");
	setOption(pars.rescueExcludedOneOffLowFreqHaplotypes, "--rescueExcludedOneOffLowFreqHaplotypes", "Rescue Excluded chimeric Haplotypes if they appear as a major haplotype in another sample");

	setOption(pars.fracCutoff, "--fracCutOff",
			"Final cluster Fraction Cut off", false, "Filtering");
	pars.differentPar = setOption(pars.parametersPopulation, "--popPar",
			"Parameters For Population Collapse", false, "Population");
	struct ClusteringParametersPars {
		bool noErrorsSet = false;
		bool strictErrorsSet = false;
		uint32_t stopAfter = 100;
		uint32_t hqMismatches = 0;
	};
	ClusteringParametersPars popClusParsPars;
	setOption(popClusParsPars.noErrorsSet, "--pop-noErrors", "Collapse parameters with no errors in population clustering", false, "Population");
	setOption(popClusParsPars.strictErrorsSet, "--pop-strictErrors", "Collapse parameters with a several low quality mismatches in population clustering", false, "Population");
	setOption(popClusParsPars.hqMismatches, "--pop-hq", "Number of high quality mismatches to allow in population clustering", false, "Population");
	setOption(popClusParsPars.stopAfter, "--pop-stopAfter", "Number of top haplotypes to check in population clustering", false, "Population");

	setOption(pars_.chiOpts_.checkChimeras_, "--markChimeras", "Check Input sequences for possible Chimeras", false, "Chimeras");
	setOption(pars.keepChimeras, "--keepChimeras", "KeepChimeras", false, "Chimeras");
	setOption(pars_.chiOpts_.parentFreqs_, "--parFreqs", "Chimeric Parent Frequency multiplier cutoff", false, "Chimeras");

	setOption(pars.numThreads, "--numThreads", "Number of threads to use");
	setOption(pars.writeOutAllInfoFile, "--writeOutAllInfoFile", "Write Out All Info File that contains information on all clusters including excluded ones");

	setOption(pars_.colOpts_.clusOpts_.converge_, "--converge", "Keep clustering at each iteration until there is no more collapsing, could increase run time significantly", false, "Clustering");
	//setOption(pars.plotRepAgreement, "--plotRepAgreement", "Plot Rep Agreement");
	processAlignerDefualts();
	if (needsHelp()) {
		printFlags(std::cout);
		std::cout << "Example input directory tree set up: " << std::endl;
		std::cout << "Samp01/Run1/output.fastq" << std::endl;
		std::cout << "Samp01/Run2/output.fastq" << std::endl;
		std::cout << "Samp02/Run1/output.fastq" << std::endl;
		std::cout << "Samp02/Run2/output.fastq" << std::endl;
		std::cout << "currentDir/" << std::endl;
		std::cout << "----Samp01/" << std::endl;
		std::cout << "--------Run1/" << std::endl;
		std::cout << "------------output.fastq" << std::endl;
		std::cout << "--------Run2/" << std::endl;
		std::cout << "------------output.fastq" << std::endl;
		std::cout << "----Samp02/" << std::endl;
		std::cout << "--------Run1/" << std::endl;
		std::cout << "------------output.fastq" << std::endl;
		std::cout << "--------Run2/" << std::endl;
		std::cout << "------------output.fastq" << std::endl;
		std::cout << "program would be run from the directory that"
				" contains Samp01 and Samp02 directories (here currentDir)"
				<< std::endl;
		std::cout
				<< "sequences in the seq files should have a suffix of _[NUM] or _t[NUM] "
						"where [NUM] is the number of reads associated with that sequence"
				<< std::endl;
		std::cout << "Final results are named with the name of the Sample Directories"
				<< std::endl;
		std::cout << njh::bashCT::bold << "Output Files:" << njh::bashCT::reset
				<< std::endl;
		std::cout
				<< "selectedClustersInfo.tab.txt: This contains the final haplotype information and replicate comparison results, it is a very large table, consult the SeekDeep (http://njh2.umassmed.edu/~hathawan/SeekDeep.html) website for details on what each column means"
				<< std::endl;
		std::cout << "allClustersInfo.tab.txt: " << std::endl;
		std::cout << "dotFiles: " << std::endl;
		std::cout << "final: " << std::endl;
		std::cout << "population/populationCluster.tab.txt: " << std::endl;
		std::cout << "population/PopUID.fastq: " << std::endl;
		exit(0);
	}
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


}  // namespace njhseq
