/*
 * SeekDeepUtilsRunner_benchmarkControlMixtures.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/ControlBenchmarking.h"
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {

std::string plotMultiTargetResultsQmd = R"(---
title: "Process Results"
format:
  html:
    theme: cosmo
    fig-height: 10
    fig-width: 10
    fig-align: center
    toc: true
    toc-depth: 4 # default is 3
    toc-title: Contents
    number-sections: true
    anchor-sections: true
    smooth-scroll: true
    highlight-style: dracula
    page-layout: full
    code-fold: true
    code-tools: true
    code-link: true
---


```{r setup, echo=FALSE, message=FALSE}
require(knitr)
require(DT)
require(tidyverse)
require(plotly)
#turn off messages and warnings and make it so output isn't prefixed by anything,
#default is to put "##" in front of all output for some reason
#also set tidy to true so code is wrapped properly
opts_chunk$set(message=FALSE, warning=FALSE, comment = "", cache = F)
options(width = 200)

myFormula= x~y
require(ggpmisc)
`%!in%` <- Negate(`%in%`)

sofonias_theme = theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
   theme(legend.position = "bottom") +
   theme(plot.title = element_text(hjust = 0.5))

sofonias_theme_xRotate = sofonias_theme +
  theme(axis.text.x = element_text(size=12, angle = -90, vjust = 0.5, hjust = 0))
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))),
                filter = "top")
}
colorPalette_08 = c("#2271B2","#F748A5","#359B73","#F0E442","#D55E00","#3DB7E9","#E69F00","#000000")
colorPalette_12 = c("#E20134","#FF6E3A","#008DF9","#8400CD","#FFC33B","#9F0162","#009F81","#FF5AAF","#00FCCF","#00C2F9","#FFB2FD","#A40122")
colorPalette_15 = c("#F60239","#003C86","#EF0096","#9400E6","#009FFA","#008169","#68023F","#00DCB5","#FFCFE2","#FF71FD","#7CFFFA","#6A0213","#008607","#00E307","#FFDC3D")


```


```{r}
samplesToMixs = readr::read_tsv("samplesToMix.tsv")
mixSetUps = readr::read_tsv("mixSetUps.tsv") %>%
  group_by(MixName) %>%
  mutate(COI = n())
create_dt(mixSetUps)

samplesToMixs = samplesToMixs %>%
  left_join(mixSetUps)

samplesToMixs_mixtures = samplesToMixs %>%
  filter(COI > 1)
```


```{r}

samplesToMixs_mod = samplesToMixs %>%
  mutate(strain = factor(strain)) %>%
  mutate(label = paste0(sample, "\n", MixName)) %>%
  arrange(MixName) %>%
  ungroup() %>%
  mutate(label = factor(label, levels = unique(.$label)))


mixturesSampleNames_mix_mod_plot = ggplot(samplesToMixs_mod) +
  geom_tile(aes(x = label, y = strain, fill = relative_abundance),
            color = "black", width = 0.8) +
  geom_text(aes(x = label, y = strain, label = paste0(relative_abundance *100, "%"), color = ifelse(relative_abundance > 50, "black", "white")), angle = "-90", fontsize = 12) +
  sofonias_theme_xRotate + labs(x = "Sample/Mixture", fill = "Within Sample\nAbunance") +
  scale_fill_gradient(high = "#fed976", low = "#bd0026",
                      limits=c(0,1)) +
  scale_color_identity()
```

```{r}
#| column: screen-inset-shaded
#| fig-align: center
#| fig-width: 20
print(mixturesSampleNames_mix_mod_plot)
```


```{r}
performancePerTarget = readr::read_tsv("performancePerTarget.tsv")
create_dt(performancePerTarget)
```


```{r}
classifiedHaplotypes = readr::read_tsv("classifiedHaplotypes.tsv")
create_dt(classifiedHaplotypes)
```

```{r}
#| fig-column: screen-inset-shaded
#| fig-align: center
ggplotly(ggplot(performancePerTarget) +
  geom_boxplot(aes(y = expectedHapRecoveryRate, x = AnalysisName, group = AnalysisName, AnalysisName = AnalysisName)) +
  sofonias_theme_xRotate)

```


```{r}
#| fig-column: screen-inset-shaded
#| fig-align: center
ggplotly(ggplot(performancePerTarget) +
  geom_boxplot(aes(y = 1 - falseHapsRate, x = AnalysisName, group = AnalysisName, AnalysisName = AnalysisName)) +
  sofonias_theme_xRotate)
```

```{r}
#| fig-column: screen-inset-shaded
#| fig-align: center

performancePerTarget_sum = performancePerTarget %>%
  group_by(AnalysisName) %>%
  summarise(expectedHapRecoveryRate = mean(expectedHapRecoveryRate),
            falseHapsRate = mean(falseHapsRate))

ggplotly(ggplot(performancePerTarget_sum) +
  geom_point(aes(y = expectedHapRecoveryRate, x = 1 - falseHapsRate,
                 AnalysisName = AnalysisName)) +
  sofonias_theme)
```

```{r}
#| fig-column: screen-inset-shaded
#| fig-align: center


classifiedHaplotypes_mix_expected = classifiedHaplotypes %>%
  filter(sample %in% samplesToMixs_mixtures$sample) %>%
  filter(matchExpected | is.na(matchExpected)) %>%
  mutate(frac = ifelse(is.na(frac), 0, frac))

ggplot(classifiedHaplotypes_mix_expected %>%
         filter(expectedFrac < 1, frac > 0), aes(x = expectedFrac, y = frac)) +
  # geom_boxplot(aes(group = expectedFrac)) +
  # geom_point(aes(group = expectedFrac)) +
  geom_bin2d(aes(fill=..ndensity.., group = expectedFrac), binwidth = c(.05, .05), color = "black") +
  stat_poly_eq(formula = myFormula, aes(x = expectedFrac, y =frac,  label = paste(..rr.label.., sep = "~~~")),size=3, parse = TRUE,
               label.y = "bottom", label.x = "right") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  sofonias_theme_xRotate +
  # scale_x_continuous(limits = c(0, 1)) +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_fill_gradient(low = "#e5f5f9", high = "#2ca25f")

```
)";

int SeekDeepUtilsRunner::benchmarkMultiTarAmpControlMixtures(
				const njh::progutils::CmdArgs &inputCommands) {

	ControlBencher::ControlBencherPars conBenchPars;

	bfs::path resultsFnp;
	std::string sampleColName = "s_Sample";
	std::string withinSampleFreqColName = "c_AveragedFrac";
	std::string withinSampleReadCntColName = "c_ReadCnt";
	std::string popHapIdColName = "h_popUID";
	std::string popHapSeqColName = "h_Consensus";

	std::string targetNameColName = "p_name";

	bfs::path expectedSeqsDirFnp = "";
	bfs::path popSeqsDirFnp = "";

	bfs::path metaFnp = "";
	bool skipMissingSamples = false;
	bool fillInMissingSamples = false;
	std::string popSeqsRegexPatRemoval = R"(_([tf])?\d+(\.\d+)?$)";
	std::string elementStr;
	std::string columnName;


	comparison allowableError;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(resultsFnp, "--resultsFnp",
									"results tab delimited file, each row is a haplotype, should have at least 5 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), 5)target name column (--targetNameColName), optionally 4th col with hap sequence (--popHapSeqColName) or read in from --popSeqsDirFnp",
									true);

	setUp.setOption(sampleColName, "--sampleColName", "sample Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleFreqColName, "--withinSampleFreqColName",
									"within Sample haplotype frequency Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleReadCntColName, "--withinSampleReadCntColName", "within Sample Read Cnt Col Column Name", false, "Results Column Names");
	setUp.setOption(popHapIdColName, "--popHapIdColName", "popHapIdColName", false, "Results Column Names");
	setUp.setOption(popHapSeqColName, "--popHapSeqColName",
									"population Haplotype Sequence Column Name, the seq to compare to expected", false, "Results Column Names");
	setUp.setOption(targetNameColName, "--targetNameColName",
									"target Name Column Name, the column name in the table which indicates the different targets", false, "Results Column Names");

	setUp.setOption(popSeqsDirFnp, "--popSeqsDirFnp",
									"Population Sequences, in this directory should be a fasta file with the name of each target");
	setUp.setOption(expectedSeqsDirFnp, "--expectedSeqsDirFnp",
									"Expected Seqs directory of fasta files, in this directory should be a fasta file with the name of each target",
									true);
	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");
	setUp.setOption(popSeqsRegexPatRemoval, "--popSeqsRegexPatRemoval",
									"optional regex pattern to process the input pop sequences from --popSeqsFnp to make it match up with the input results folder");

	setUp.setOption(conBenchPars.samplesToMixFnp_, "--sampleToMixture",
									"Sample To Mixture, 2 columns 1)sample, 2)MixName", true);
	setUp.setOption(conBenchPars.mixSetUpFnp_, "--mixtureSetUp",
									"Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);
	setUp.setOption(skipMissingSamples, "--skipMissingSamples", "Skip Samples if they are missing");
	setUp.setOption(fillInMissingSamples, "--fillInMissingSamples", "Fill in missing Samples with placeholders");

	setUp.setOption(columnName, "--newColumnName",
									"Name of a new Column to add to table, can be comma delimited to add multiple");
	setUp.setOption(elementStr, "--newColumnElement",
									"What to Add to the Table Under new Column, can be comma delimited when adding multiple new columns");

	setUp.processComparison(allowableError);
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName("benchmarkMultiTarAmpControlMixtures_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow({resultsFnp, expectedSeqsDirFnp,
																	 conBenchPars.samplesToMixFnp_, conBenchPars.mixSetUpFnp_},
																	__PRETTY_FUNCTION__);


	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double>>> readCountsPerHapPerSample;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> cNameToPopUID;
	std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> hPopUIDPopSamps;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hPopUID_to_hConsensus;

	std::set<std::string> targetNamesSet;

	VecStr requiredColumns{sampleColName, withinSampleReadCntColName,
												 withinSampleFreqColName, withinSampleReadCntColName,
												 popHapIdColName,
												 targetNameColName};
	if (!exists(popSeqsDirFnp)) {
		requiredColumns.emplace_back(popHapSeqColName);
	}
	uint64_t maxLen = 0;
	//population seqs;
	std::unordered_map<std::string, std::vector<seqInfo>> popSeqs;
	auto sampInfoFnp = resultsFnp;
	TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));
	sampInfoReader.header_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);
	VecStr row;
	std::unordered_map<std::string, std::unordered_set<std::string>> allSamplesInOutput;
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> allResultSeqs;
	while (sampInfoReader.getNextRow(row)) {
		auto target = row[sampInfoReader.header_.getColPos(targetNameColName)];
		targetNamesSet.emplace(target);
		auto sample = row[sampInfoReader.header_.getColPos(sampleColName)];

		auto c_ReadCnt = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);
		auto c_AveragedFrac = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleFreqColName)]);
		auto readCnt = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);

		auto h_popUID = row[sampInfoReader.header_.getColPos(popHapIdColName)];
		auto hapName = njh::pasteAsStr(sample, "__", h_popUID);
		std::string hapSeq;
		if (popSeqsDirFnp.empty()) {
			hPopUID_to_hConsensus[target][h_popUID] = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
			hapSeq = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
		} else {
			hapSeq = njh::mapAt(njh::mapAt(hPopUID_to_hConsensus, target), h_popUID);
		}
		seqInfo clus(hapName, hapSeq);
		clus.cnt_ = c_ReadCnt;
		clus.frac_ = c_AveragedFrac;
		bool add = true;
		for (auto &seq: allResultSeqs[target][sample]) {
			if (seq.seq_ == clus.seq_) {
				add = false;
				seq.cnt_ += c_ReadCnt;
				seq.frac_ += c_AveragedFrac;
				readCountsPerHapPerSample[target][sample][hapName] += readCnt;
				break;
			}
		}
		if (add) {
			readCountsPerHapPerSample[target][sample][hapName] = readCnt;
			allResultSeqs[target][sample].emplace_back(clus);
		}

		readVec::getMaxLength(clus, maxLen);
		cNameToPopUID[target][hapName] = h_popUID;
		hPopUIDPopSamps[target][h_popUID].emplace(sample);
		allSamplesInOutput[target].emplace(sample);
	}
	if (!exists(popSeqsDirFnp)) {
		for (const auto &tarPopHaps: hPopUID_to_hConsensus) {
			for (const auto &popHap: tarPopHaps.second) {
				popSeqs[tarPopHaps.first].emplace_back(popHap.first, popHap.second);
			}
		}
	}

	auto targetNames = std::vector<std::string>(targetNamesSet.begin(), targetNamesSet.end());
	njh::naturalSortNameSet(targetNames);

	//benchers
	std::map<std::string, std::unique_ptr<ControlBencher>> benchers;

	table sampleToMixTab(conBenchPars.samplesToMixFnp_, "\t", true);
	table mixSetUpTab(conBenchPars.mixSetUpFnp_, "\t", true);

	std::map<std::string, table> sampleToMixSplitTab;
	std::map<std::string, table> mixSetUpSplitTab;

	if (njh::in(targetNameColName, sampleToMixTab.columnNames_) &&
			njh::in(targetNameColName, mixSetUpTab.columnNames_)) {
		sampleToMixSplitTab = sampleToMixTab.splitTableOnColumn(targetNameColName);
		mixSetUpSplitTab = mixSetUpTab.splitTableOnColumn(targetNameColName);
		auto sampleToMixSplitTab_targetNames = njh::getVecOfMapKeys(sampleToMixSplitTab);
		auto mixSetUpSplitTab_targetNames = njh::getVecOfMapKeys(mixSetUpSplitTab);

		VecStr inSampToMix;
		VecStr inMixSetup;
		VecStr shared;
		njh::decompose_sets(sampleToMixSplitTab_targetNames.begin(), sampleToMixSplitTab_targetNames.end(),
												mixSetUpSplitTab_targetNames.begin(), mixSetUpSplitTab_targetNames.end(),
												std::back_inserter(inSampToMix),
												std::back_inserter(inMixSetup),
												std::back_inserter(shared));
		VecStr warnings;
		if (!inSampToMix.empty()) {
			warnings.emplace_back(
							njh::pasteAsStr("error, found the following targets in ", conBenchPars.samplesToMixFnp_, " but not in ",
															conBenchPars.mixSetUpFnp_,
															"\n", njh::conToStr(inSampToMix, ",")));
		}
		if (!inMixSetup.empty()) {
			warnings.emplace_back(
							njh::pasteAsStr("error, found the following targets in ", conBenchPars.mixSetUpFnp_, " but not in ",
															conBenchPars.samplesToMixFnp_,
															"\n", njh::conToStr(inMixSetup, ",")));
		}
		if (!warnings.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " found the following errors" << "\n";
			throw std::runtime_error{ss.str()};
		}
		targetNames = shared;
		njh::naturalSortNameSet(targetNames);
	}
	if (mixSetUpSplitTab.empty()) {
		for (const auto &target: targetNames) {
			benchers[target] = std::make_unique<ControlBencher>(conBenchPars);
		}
	} else {
		for (const auto &target: targetNames) {
			benchers[target] = std::make_unique<ControlBencher>(
							njh::mapAt(mixSetUpSplitTab, target),
							njh::mapAt(sampleToMixSplitTab, target)
			);
		}
	}



	if (exists(popSeqsDirFnp)) {
		seqInfo seq;
		for (const auto &target: targetNames) {
			if (!bfs::exists(njh::files::make_path(popSeqsDirFnp, target + ".fasta")) &&
					!bfs::exists(njh::files::make_path(popSeqsDirFnp, target + ".fasta.gz"))
							) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "error, could not find population seqs file for target: " << target
					 << ", should be " << njh::files::make_path(popSeqsDirFnp, target + ".fasta") << " or "
					 << njh::files::make_path(popSeqsDirFnp, target + ".fasta.gz") << "\n";
				throw std::runtime_error{ss.str()};
			}
			auto popSeqsFnp = njh::files::make_path(popSeqsDirFnp, target + ".fasta");
			if (!exists(popSeqsFnp)) {
				popSeqsFnp = njh::files::make_path(popSeqsDirFnp, target + ".fasta.gz");
			}
			SeqInput reader{SeqIOOptions(popSeqsFnp, SeqIOOptions::getInFormatFromFnp(popSeqsFnp))};
			reader.openIn();
			while (reader.readNextRead(seq)) {
				readVec::getMaxLength(seq, maxLen);
				seq.name_ = std::regex_replace(seq.name_, std::regex{popSeqsRegexPatRemoval}, "");
				hPopUID_to_hConsensus[target][seq.name_] = seq.seq_;
			}
		}
	}

	//check for samples in analysis
	std::unordered_map<std::string, std::set<std::string>> missingSamples;
	std::unordered_map<std::string, std::set<std::string>> controlSamples;
	for (const auto &target: targetNames) {
		auto currentControlSamples = benchers[target]->getSamples();
		controlSamples[target].insert(currentControlSamples.begin(), currentControlSamples.end());
		for (const auto &sname: currentControlSamples) {
			//if (!analysisMaster.hasSample(sname)) {
			if (!njh::in(sname, allSamplesInOutput[target])) {
				missingSamples[target].emplace(sname);
			}
		}
	}

	if (!skipMissingSamples && !missingSamples.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
			 << "missing the following samples in " << resultsFnp << " for targets" << "\n";
		for (const auto &missing: missingSamples) {
			ss << "tar: " << missing.first << ", samples: " << njh::conToStrEndSpecial(missing.second, ", ", " and ") << "\n";
		}

		throw std::runtime_error{ss.str()};
	} else if (skipMissingSamples && !missingSamples.empty()) {
		OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples.txt"));
		OutputStream outMissing(outOptsMissing);
		outMissing << "target\tsample" << std::endl;
		for (const auto &missing: missingSamples) {
			for (const auto &samp: missing.second) {
				outMissing << missing.first << "\t" << samp << std::endl;
			}
		}
	}



	//read in expected seqs
	std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> initialExpSeqs;
	for (const auto &target: targetNames) {
		if (!bfs::exists(njh::files::make_path(expectedSeqsDirFnp, target + ".fasta")) &&
				!bfs::exists(njh::files::make_path(expectedSeqsDirFnp, target + ".fasta.gz"))
						) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "error, could not find population seqs file for target: " << target
				 << ", should be " << njh::files::make_path(expectedSeqsDirFnp, target + ".fasta") << " or "
				 << njh::files::make_path(expectedSeqsDirFnp, target + ".fasta.gz") << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto expectedSeqsFnp = njh::files::make_path(expectedSeqsDirFnp, target + ".fasta");
		if (!exists(expectedSeqsFnp)) {
			expectedSeqsFnp = njh::files::make_path(expectedSeqsDirFnp, target + ".fasta.gz");
		}
		SeqInput expSeqsSeqIn{SeqIOOptions(expectedSeqsFnp, SeqIOOptions::getInFormatFromFnp(expectedSeqsFnp))};
		initialExpSeqs[target] = expSeqsSeqIn.readAllReadsPtrs<seqInfo>();
	}



	//check the needed expected names
	std::map<std::string, std::map<std::string, uint32_t>> initialExpSeqsPositions;
	std::map<std::string, std::map<std::string, std::string>> expSeqsKey;
	std::map<std::string, std::map<std::string, VecStr>> matchingExpectedSeqs;
	std::map<std::string, std::vector<std::shared_ptr<seqInfo>>> expSeqs;
	std::map<std::string, std::map<std::string, uint32_t> > finalExpSeqsPositions;

	for (const auto &target: targetNames) {
		VecStr removeSeqs;
		std::set<std::string> expNames;

		for (const auto &expSeq: initialExpSeqs[target]) {

			if (std::all_of(expSeq->seq_.begin(), expSeq->seq_.end(),
											[](const char base) {
												return 'N' == base || 'n' == base;
											})) {
				removeSeqs.emplace_back(expSeq->name_);
				continue;
			}

			if (njh::in(expSeq->name_, expNames)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " already have expected sequence named " << expSeq->name_ << "\n";
				throw std::runtime_error{ss.str()};
			}

			expNames.emplace(expSeq->name_);
			bool found = false;
			for (const auto otherSeqPos: iter::range(expSeqs[target].size())) {
				const auto &otherSeq = expSeqs[target][otherSeqPos];
				if (otherSeq->seq_ == expSeq->seq_) {
					found = true;
					otherSeq->name_ += "," + expSeq->name_;
					initialExpSeqsPositions[target][expSeq->name_] = otherSeqPos;
					break;
				}
			}

			if (!found) {
				//if not found add
				//add position of the expected seq
				initialExpSeqsPositions[target][expSeq->name_] = expSeqs[target].size();
				expSeqs[target].emplace_back(std::make_shared<seqInfo>(*expSeq));
			}
		}

		//remove ref seqs that were zero;
		benchers[target]->removeStrains(removeSeqs);



		//

		for (const auto expSeqPos: iter::range(expSeqs[target].size())) {
			readVec::getMaxLength(expSeqs[target][expSeqPos], maxLen);
			expSeqsKey[target][expSeqs[target][expSeqPos]->name_] = expSeqs[target][expSeqPos]->name_;
			finalExpSeqsPositions[target][expSeqs[target][expSeqPos]->name_] = expSeqPos;
		}

		benchers[target]->checkForStrainsThrow(expNames, __PRETTY_FUNCTION__);

	}
	bfs::copy(njh::files::normalize(conBenchPars.samplesToMixFnp_), njh::files::make_path(setUp.pars_.directoryName_, "samplesToMix.tsv"));

	if (njh::in(targetNameColName, sampleToMixTab.columnNames_) &&
			njh::in(targetNameColName, mixSetUpTab.columnNames_)) {
		//all benches are basically the same so just use the beginning
		OutputStream mixOut(njh::files::make_path(setUp.pars_.directoryName_, "mixSetUps.tsv"));
		table allMixTab;
		auto tarBenchNames = njh::getVecOfMapKeys(benchers);
		njh::naturalSortNameSet(tarBenchNames);
		for (const auto &tarBenchName: tarBenchNames) {

			auto mixTab = benchers.at(tarBenchName)->genMixSetUpsInSamplesTab();
			mixTab.addColumn(VecStr{tarBenchName}, targetNameColName);
			if (allMixTab.empty()) {
				allMixTab = mixTab;
			} else {
				allMixTab.rbind(mixTab, false);
			}
		}
		allMixTab.outPutContents(mixOut, "\t");
	} else {
		benchers.begin()->second->writeMixSetUpsInSamples(
						njh::files::make_path(setUp.pars_.directoryName_, "mixSetUps.tsv"));
	}


	OutputStream falseHaplotypesToExpClassified(
					njh::files::make_path(setUp.pars_.directoryName_, "falseHaplotypesComparedToExpected.tsv"));
	falseHaplotypesToExpClassified
					<< "AnalysisName\tsample\tmix\tInputReadName\treadCnt\tfrac\tRefName\tExpectedRefFreq\tExpectedMajorOrMinor\tIdentityScore\tbestMatchScore\tmismatches\toneBaseIndels\ttwoBaseIndels\tlargeIndels\ttotalErrors";
	OutputStream falseHaplotypesToOtherResultsClassified(
					njh::files::make_path(setUp.pars_.directoryName_, "falseHaplotypesComparedToOthers.tsv"));
	falseHaplotypesToOtherResultsClassified
					<< "AnalysisName\tsample\tmix\tInputReadName\treadCnt\tfrac\tOtherName\tOtherReadCnt\tOtherFrac\tratio\tOtherMajor\tOtherMatchesExpected\tOtherExpectedMatchName\tIdentityScore\tbestMatchScore\tmismatches\toneBaseIndels\ttwoBaseIndels\tlargeIndels\ttotalErrors";


	OutputStream haplotypesClassified(njh::files::make_path(setUp.pars_.directoryName_, "classifiedHaplotypes.tsv"));
	//haplotypesClassified << "AnalysisName\tsample\tmix\tInputReadName\tHapPopUID\tHapSampleCount\treadCnt\tfrac\tmatchExpected\texpectedRef\texpectedFrac\tMajorOrMinor\tmatchingPopulation\tPopName";
	haplotypesClassified
					<< "AnalysisName\tsample\tmix\tInputReadName\tHapPopUID\tHapSampleCount\treadCnt\tfrac\tmatchExpected\texpectedRef\texpectedFrac\tMajorOrMinor";
	OutputStream performanceOut(njh::files::make_path(setUp.pars_.directoryName_, "performancePerTarget.tsv"));
	performanceOut
					<< "AnalysisName\tsample\tmix\ttotalReads\trecoveredExpectedHaps\tfalseHaps\ttotalHaps\ttotalExpectedHaps\texpectedHapRecoveryRate\tfalseHapsRate\texpectedMissing\tRMSE";
	VecStr metalevels;

	VecStr newColumnEleToks;
	VecStr newColumnToks;
	if (!columnName.empty()) {
		newColumnEleToks = tokenizeString(elementStr, ",");
		newColumnToks = tokenizeString(columnName, ",");
		if (newColumnEleToks.size() != newColumnToks.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
				 << ", error, when adding multiple columns the number of elements need to match number of new columns" << "\n"
				 << "New Columns #: " << newColumnToks.size() << ", New Elements #: " << newColumnEleToks.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	if (nullptr != metaGroupData) {
		metalevels = getVectorOfMapKeys(metaGroupData->groupData_);
		njh::sort(metalevels);
		for (const auto &meta: metalevels) {
			haplotypesClassified << "\t" << meta;
			performanceOut << "\t" << meta;
			falseHaplotypesToExpClassified << "\t" << meta;
			falseHaplotypesToOtherResultsClassified << "\t" << meta;
		}
	}
	if (!newColumnToks.empty()) {
		for (const auto &col: newColumnToks) {
			haplotypesClassified << "\t" << col;
			performanceOut << "\t" << col;
			falseHaplotypesToExpClassified << "\t" << col;
			falseHaplotypesToOtherResultsClassified << "\t" << col;
		}
	}

	haplotypesClassified << std::endl;
	performanceOut << std::endl;
	falseHaplotypesToExpClassified << std::endl;
	falseHaplotypesToOtherResultsClassified << std::endl;
	aligner alignerObj(maxLen, gapScoringParameters(5, 1, 0, 0, 0, 0));
	alignerObj.weighHomopolymers_ = setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_;
	alignerObj.countEndGaps_ = true;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, false);
	for (const auto &target: targetNames) {

		std::cout << "target: " << target << std::endl;

		for (const auto &sname: njh::mapAt(controlSamples,target )) {
			//skip completely missing
			if (skipMissingSamples && njh::in(sname, missingSamples)) {
				continue;
			}

			const std::vector<seqInfo> &resultSeqs = allResultSeqs[target][sname];


			std::map<std::string, uint32_t> resSeqToPos;
			double maxResFrac = 0;
			for (const auto pos: iter::range(resultSeqs.size())) {
				resSeqToPos[resultSeqs[pos].name_] = pos;
				if (resultSeqs[pos].frac_ > maxResFrac) {
					maxResFrac = resultSeqs[pos].frac_;
				}
			}


			//get current expected seqs
			std::map<std::string, double> currentExpectedSeqsFrac;
			double maxExpFrac = 0;
			std::map<std::string, VecStr> expectedSeqNameToCurrentSeqsKey;


			for (const auto &expSeqFrac: benchers[target]->mixSetups_.at(
							benchers[target]->samplesToMix_.at(sname)).relativeAbundances_) {
				currentExpectedSeqsFrac[expSeqs[target][initialExpSeqsPositions[target][expSeqFrac.first]]->name_] += expSeqFrac.second;
				expectedSeqNameToCurrentSeqsKey[expSeqs[target][initialExpSeqsPositions[target][expSeqFrac.first]]->name_].emplace_back(
								expSeqFrac.first);
			}
			for (auto &key: expectedSeqNameToCurrentSeqsKey) {
//			std::cout << "key: " << key.first << std::endl;
//			std::cout << "\t" << njh::conToStr(key.second, ",") << std::endl;
				njh::sort(key.second);
			}


			for (const auto &expFrac: currentExpectedSeqsFrac) {
				if (expFrac.second > maxExpFrac) {
					maxExpFrac = expFrac.second;
				}
			}


			std::map<std::string, std::string> expectedToMajorClass;
			for (const auto &expFrac: currentExpectedSeqsFrac) {
				if (expFrac.second == maxExpFrac) {
					expectedToMajorClass[expFrac.first] = "major";
				} else {
					expectedToMajorClass[expFrac.first] = "minor";
				}
			}
			std::map<std::string, std::string> resultsToMajorClass;
			for (const auto &res: resultSeqs) {
				if (res.frac_ == maxResFrac) {
					resultsToMajorClass[res.name_] = "major";
				} else {
					resultsToMajorClass[res.name_] = "minor";
				}
			}


			std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
			currentExpectedSeqs.reserve(currentExpectedSeqsFrac.size());
			for (const auto &finalExp: currentExpectedSeqsFrac) {
				currentExpectedSeqs.emplace_back(expSeqs[target][finalExpSeqsPositions[target][finalExp.first]]);
			}

			auto res = njhseq::ControlBencher::benchmark(resultSeqs, currentExpectedSeqs,
																									 currentExpectedSeqsFrac,
																									 expSeqsKey[target], alignerObj, allowableError);
			double total = 0;



			//haplotype classification
			for (const auto &seq: resultSeqs) {
				total += readCountsPerHapPerSample[target][sname][seq.name_];
				haplotypesClassified << target
														 << "\t" << sname
														 << "\t" << benchers[target]->samplesToMix_[sname]
														 << "\t" << seq.name_
														 << "\t" << cNameToPopUID[target][seq.name_]
														 << "\t" << hPopUIDPopSamps[target][cNameToPopUID[target][seq.name_]].size()
														 << "\t" << readCountsPerHapPerSample[target][sname][seq.name_]
														 << "\t" << seq.frac_
														 << "\t" << (res.resSeqToExpSeq_[seq.name_].empty() ? "FALSE" : "TRUE")
														 << "\t"
														 << njh::conToStr(expectedSeqNameToCurrentSeqsKey[res.resSeqToExpSeq_[seq.name_]], ",")
														 << "\t" << currentExpectedSeqsFrac[res.resSeqToExpSeq_[seq.name_]]
														 << "\t" << expectedToMajorClass[res.resSeqToExpSeq_[seq.name_]];

				std::string matchingPop;
				for (const auto &popSeq: popSeqs[target]) {
					if (popSeq.seq_ == seq.seq_) {
						matchingPop = popSeq.name_;
						break;
					}
				}
//			haplotypesClassified
//							<< "\t" << (matchingPop.empty() ? "TRUE" : "FALSE")
//							<< "\t" << matchingPop;
				if (nullptr != metaGroupData) {
					for (const auto &meta: metalevels) {
						haplotypesClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
					}
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						haplotypesClassified << "\t" << ele;
					}
				}
				haplotypesClassified << std::endl;
			}
			for (const auto &missing: res.missingExpecteds_) {
				haplotypesClassified << target
														 << "\t" << sname
														 << "\t" << benchers[target]->samplesToMix_[sname]
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ",")
														 << "\t" << currentExpectedSeqsFrac[missing]
														 << "\t" << expectedToMajorClass[missing];

//			haplotypesClassified
//							<< "\t" << "NA"
//							<< "\t" << "NA";
				if (nullptr != metaGroupData) {
					for (const auto &meta: metalevels) {
						haplotypesClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
					}
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						haplotypesClassified << "\t" << ele;
					}
				}
				haplotypesClassified << std::endl;
			}
			//performance
			VecStr missingExpectedDecoded;
			for (const auto &missing: res.missingExpecteds_) {
				missingExpectedDecoded.emplace_back(njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ","));
			}
			performanceOut << target
										 << "\t" << sname
										 << "\t" << benchers[target]->samplesToMix_[sname]
										 << "\t" << total
										 << "\t" << res.recoveredHaps_
										 << "\t" << res.falseHaps_
										 << "\t" << res.totalHaps()
										 << "\t" << res.expectedHapCnt_
										 << "\t" << res.hapRecoveryRate()
										 << "\t" << res.falseHapRate()
										 << "\t" << njh::conToStr(missingExpectedDecoded, ";")
										 << "\t" << res.RMSE();
			if (nullptr != metaGroupData) {
				for (const auto &meta: metalevels) {
					performanceOut << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
				}
			}
			if (!newColumnToks.empty()) {
				for (const auto &ele: newColumnEleToks) {
					performanceOut << "\t" << ele;
				}
			}
			performanceOut << std::endl;



			//further classification of false haplotypes to expected
			for (const auto &falseHap: res.falseHapsCompsToExpected) {
				double bestScore = 0;
				for (const auto &exp: falseHap.second) {
					if (exp.second.distances_.eventBasedIdentityHq_ > bestScore) {
						bestScore = exp.second.distances_.eventBasedIdentityHq_;
					}
				}
				for (const auto &exp: falseHap.second) {
					falseHaplotypesToExpClassified << target
																				 << "\t" << sname
																				 << "\t" << benchers[target]->samplesToMix_[sname]
																				 << "\t" << falseHap.first
																				 << "\t" << readCountsPerHapPerSample[target][sname][falseHap.first]
																				 << "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
																				 << "\t" << exp.first
																				 << "\t" << currentExpectedSeqsFrac[exp.first]
																				 << "\t" << expectedToMajorClass[exp.first]
																				 << "\t" << exp.second.distances_.eventBasedIdentityHq_
																				 << "\t"
																				 << (exp.second.distances_.eventBasedIdentityHq_ == bestScore ? "TRUE"
																																																			: "FALSE")
																				 << "\t" << exp.second.hqMismatches_ + exp.second.lqMismatches_ +
																										exp.second.lowKmerMismatches_
																				 << "\t" << exp.second.oneBaseIndel_
																				 << "\t" << exp.second.twoBaseIndel_
																				 << "\t" << exp.second.largeBaseIndel_
																				 << "\t" << exp.second.distances_.getNumOfEvents(true);
					if (nullptr != metaGroupData) {
						for (const auto &meta: metalevels) {
							falseHaplotypesToExpClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
						}
					}
					if (!newColumnToks.empty()) {
						for (const auto &ele: newColumnEleToks) {
							falseHaplotypesToExpClassified << "\t" << ele;
						}
					}
					falseHaplotypesToExpClassified << std::endl;
				}
			}

			//further classification of false haplotypes to other sequences
			for (const auto &falseHap: res.falseHapsCompsToOthers) {
				double bestScore = 0;
				for (const auto &other: falseHap.second) {
					if (other.second.distances_.eventBasedIdentityHq_ > bestScore) {
						bestScore = other.second.distances_.eventBasedIdentityHq_;
					}
				}
				for (const auto &other: falseHap.second) {
					falseHaplotypesToOtherResultsClassified << target
																									<< "\t" << sname
																									<< "\t" << benchers[target]->samplesToMix_[sname]
																									<< "\t" << falseHap.first
																									<< "\t" << readCountsPerHapPerSample[target][sname][falseHap.first]
																									<< "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
																									<< "\t" << other.first
																									<< "\t" << readCountsPerHapPerSample[target][sname][other.first]
																									<< "\t" << resultSeqs[resSeqToPos[other.first]].frac_
																									<< "\t" << resultSeqs[resSeqToPos[other.first]].frac_ /
																														 resultSeqs[resSeqToPos[falseHap.first]].frac_
																									<< "\t" << resultsToMajorClass[other.first]
																									<< "\t"
																									<< (res.resSeqToExpSeq_[other.first].empty() ? "FALSE" : "TRUE")
																									<< "\t" << res.resSeqToExpSeq_[other.first]
																									<< "\t" << other.second.distances_.eventBasedIdentityHq_
																									<< "\t"
																									<< (other.second.distances_.eventBasedIdentityHq_ == bestScore
																											? "TRUE"
																											: "FALSE")
																									<< "\t" << other.second.hqMismatches_ + other.second.lqMismatches_ +
																														 other.second.lowKmerMismatches_
																									<< "\t" << other.second.oneBaseIndel_
																									<< "\t" << other.second.twoBaseIndel_
																									<< "\t" << other.second.largeBaseIndel_
																									<< "\t" << other.second.distances_.getNumOfEvents(true);
					if (nullptr != metaGroupData) {
						for (const auto &meta: metalevels) {
							falseHaplotypesToOtherResultsClassified << "\t"
																											<< metaGroupData->groupData_[meta]->getGroupForSample(sname);
						}
					}
					if (!newColumnToks.empty()) {
						for (const auto &ele: newColumnEleToks) {
							falseHaplotypesToOtherResultsClassified << "\t" << ele;
						}
					}
					falseHaplotypesToOtherResultsClassified << std::endl;
				}
			}
		}


		if (fillInMissingSamples) {
			for (const auto &sname: missingSamples[target]) {

				const std::vector<seqInfo> &resultSeqs = allResultSeqs[target][sname];

				std::map<std::string, uint32_t> resSeqToPos;
				double maxResFrac = 0;
				for (const auto pos: iter::range(resultSeqs.size())) {
					resSeqToPos[resultSeqs[pos].name_] = pos;
					if (resultSeqs[pos].frac_ > maxResFrac) {
						maxResFrac = resultSeqs[pos].frac_;
					}
				}
				//get current expected seqs
				std::map<std::string, double> currentExpectedSeqsFrac;
				double maxExpFrac = 0;
				std::map<std::string, VecStr> expectedSeqNameToCurrentSeqsKey;


				for (const auto &expSeqFrac: benchers[target]->mixSetups_.at(
								benchers[target]->samplesToMix_.at(sname)).relativeAbundances_) {
					currentExpectedSeqsFrac[expSeqs[target][initialExpSeqsPositions[target][expSeqFrac.first]]->name_] += expSeqFrac.second;
					expectedSeqNameToCurrentSeqsKey[expSeqs[target][initialExpSeqsPositions[target][expSeqFrac.first]]->name_].emplace_back(
									expSeqFrac.first);
				}
				for (auto &key: expectedSeqNameToCurrentSeqsKey) {
					njh::sort(key.second);
				}
				for (const auto &expFrac: currentExpectedSeqsFrac) {
					if (expFrac.second > maxExpFrac) {
						maxExpFrac = expFrac.second;
					}
				}

				std::map<std::string, std::string> expectedToMajorClass;
				for (const auto &expFrac: currentExpectedSeqsFrac) {
					if (expFrac.second == maxExpFrac) {
						expectedToMajorClass[expFrac.first] = "major";
					} else {
						expectedToMajorClass[expFrac.first] = "minor";
					}
				}
				std::map<std::string, std::string> resultsToMajorClass;
				for (const auto &res: resultSeqs) {
					if (res.frac_ == maxResFrac) {
						resultsToMajorClass[res.name_] = "major";
					} else {
						resultsToMajorClass[res.name_] = "minor";
					}
				}



				std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
				currentExpectedSeqs.reserve(currentExpectedSeqsFrac.size());
				for (const auto &finalExp: currentExpectedSeqsFrac) {
					currentExpectedSeqs.push_back(expSeqs[target][finalExpSeqsPositions[target][finalExp.first]]);
				}
				VecStr missingExpectedDecoded;
				for (const auto &missing: currentExpectedSeqs) {
					missingExpectedDecoded.emplace_back(njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing->name_], ","));
				}
				//performance
				performanceOut << target
											 << "\t" << sname
											 << "\t" << benchers[target]->samplesToMix_[sname]
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << currentExpectedSeqs.size()
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << njh::conToStr(missingExpectedDecoded, ";")
											 << "\t" << 0;
				if (nullptr != metaGroupData) {
					for (const auto &meta: metalevels) {
						performanceOut << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
					}
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						performanceOut << "\t" << ele;
					}
				}
				performanceOut << std::endl;
				for (const auto &missing: readVec::getNames(currentExpectedSeqs)) {
					haplotypesClassified << target
															 << "\t" << sname
															 << "\t" << benchers[target]->samplesToMix_[sname]
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ",")
															 << "\t" << currentExpectedSeqsFrac[missing]
															 << "\t" << expectedToMajorClass[missing];

//				haplotypesClassified
//								<< "\t" << "NA"
//								<< "\t" << "NA";
					if (nullptr != metaGroupData) {
						for (const auto &meta: metalevels) {
							haplotypesClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
						}
					}
					if (!newColumnToks.empty()) {
						for (const auto &ele: newColumnEleToks) {
							haplotypesClassified << "\t" << ele;
						}
					}
					haplotypesClassified << std::endl;
				}
			}
		}
	}


	OutputStream qmdOut(njh::files::make_path(setUp.pars_.directoryName_, "process.qmd"));
	qmdOut << plotMultiTargetResultsQmd << std::endl;
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, false);

	return 0;

}


std::string plotSingleTargetResultsQmd = R"(---
title: "Process Results"
format:
  html:
    theme: cosmo
    fig-height: 10
    fig-width: 10
    fig-align: center
    toc: true
    toc-depth: 4 # default is 3
    toc-title: Contents
    number-sections: true
    anchor-sections: true
    smooth-scroll: true
    highlight-style: dracula
    page-layout: full
    code-fold: true
    code-tools: true
    code-link: true
---


```{r setup, echo=FALSE, message=FALSE}
require(knitr)
require(DT)
require(tidyverse)
#turn off messages and warnings and make it so output isn't prefixed by anything,
#default is to put "##" in front of all output for some reason
#also set tidy to true so code is wrapped properly
opts_chunk$set(message=FALSE, warning=FALSE, comment = "", cache = F)
options(width = 200)

myFormula= x~y
require(ggpmisc)
`%!in%` <- Negate(`%in%`)

sofonias_theme = theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
   theme(legend.position = "bottom") +
   theme(plot.title = element_text(hjust = 0.5))

sofonias_theme_xRotate = sofonias_theme +
  theme(axis.text.x = element_text(size=12, angle = -90, vjust = 0.5, hjust = 0))
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))),
                filter = "top")
}
colorPalette_08 = c("#2271B2","#F748A5","#359B73","#F0E442","#D55E00","#3DB7E9","#E69F00","#000000")
colorPalette_12 = c("#E20134","#FF6E3A","#008DF9","#8400CD","#FFC33B","#9F0162","#009F81","#FF5AAF","#00FCCF","#00C2F9","#FFB2FD","#A40122")
colorPalette_15 = c("#F60239","#003C86","#EF0096","#9400E6","#009FFA","#008169","#68023F","#00DCB5","#FFCFE2","#FF71FD","#7CFFFA","#6A0213","#008607","#00E307","#FFDC3D")


```


```{r}
samplesToMixs = readr::read_tsv("samplesToMix.tsv")
mixSetUps = readr::read_tsv("mixSetUps.tsv") %>%
  group_by(MixName) %>%
  mutate(COI = n())
create_dt(mixSetUps)

samplesToMixs = samplesToMixs %>%
  left_join(mixSetUps)

samplesToMixs_mixtures = samplesToMixs %>%
  filter(COI > 1)
```


```{r}

samplesToMixs_mod = samplesToMixs %>%
  mutate(strain = factor(strain)) %>%
  mutate(label = paste0(sample, "\n", MixName)) %>%
  arrange(MixName) %>%
  ungroup() %>%
  mutate(label = factor(label, levels = unique(.$label)))


mixturesSampleNames_mix_mod_plot = ggplot(samplesToMixs_mod) +
  geom_tile(aes(x = label, y = strain, fill = relative_abundance),
            color = "black", width = 0.8) +
  geom_text(aes(x = label, y = strain, label = paste0(relative_abundance *100, "%"), color = ifelse(relative_abundance > 50, "black", "white")), angle = "-90", fontsize = 12) +
  sofonias_theme_xRotate + labs(x = "Sample/Mixture", fill = "Within Sample\nAbunance") +
  scale_fill_gradient(high = "#fed976", low = "#bd0026",
                      limits=c(0,1)) +
  scale_color_identity()
```

```{r}
#| column: screen-inset-shaded
#| fig-align: center
#| fig-width: 20
print(mixturesSampleNames_mix_mod_plot)
```


```{r}
performancePerTarget = readr::read_tsv("performancePerTarget.tsv")
create_dt(performancePerTarget)
```


```{r}
classifiedHaplotypes = readr::read_tsv("classifiedHaplotypes.tsv")
create_dt(classifiedHaplotypes)
```

```{r}
ggplot(performancePerTarget) +
  geom_boxplot(aes(y = expectedHapRecoveryRate)) +
  sofonias_theme

ggplot(performancePerTarget) +
  geom_boxplot(aes(y = 1 - falseHapsRate)) +
  sofonias_theme

ggplot(performancePerTarget) +
  geom_bin2d(aes(y = expectedHapRecoveryRate, x = 1 - falseHapsRate)) +
  sofonias_theme
```

```{r}
#| fig-column: screen-inset-shaded
#| fig-align: center


classifiedHaplotypes_mix_expected = classifiedHaplotypes %>%
  filter(sample %in% samplesToMixs_mixtures$sample) %>%
  filter(matchExpected | is.na(matchExpected)) %>%
  mutate(frac = ifelse(is.na(frac), 0, frac))

ggplot(classifiedHaplotypes_mix_expected %>%
         filter(expectedFrac < 1, frac > 0), aes(x = expectedFrac, y = frac)) +
  # geom_boxplot(aes(group = expectedFrac)) +
  # geom_point(aes(group = expectedFrac)) +
  geom_bin2d(aes(fill=..ndensity.., group = expectedFrac), binwidth = c(.05, .05), color = "black") +
  stat_poly_eq(formula = myFormula, aes(x = expectedFrac, y =frac,  label = paste(..rr.label.., sep = "~~~")),size=3, parse = TRUE,
               label.y = "bottom", label.x = "right") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  sofonias_theme_xRotate +
  # scale_x_continuous(limits = c(0, 1)) +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_fill_gradient(low = "#e5f5f9", high = "#2ca25f")

```
)";

int SeekDeepUtilsRunner::benchmarkTarAmpControlMixtures(
				const njh::progutils::CmdArgs &inputCommands) {
	ControlBencher::ControlBencherPars conBenchPars;

	bfs::path resultsFnp;
	std::string sampleColName = "s_Sample";
	std::string withinSampleFreqColName = "c_AveragedFrac";
	std::string withinSampleReadCntColName = "c_ReadCnt";
	std::string popHapIdColName = "h_popUID";
	std::string popHapSeqColName = "h_Consensus";

	bfs::path expectedSeqsFnp = "";
	std::string name;
	bfs::path metaFnp = "";
	bool skipMissingSamples = false;
	bool fillInMissingSamples = false;
	bfs::path popSeqsFnp = "";
	std::string popSeqsRegexPatRemoval = R"(_([tf])?\d+(\.\d+)?$)";
	std::string elementStr;
	std::string columnName;


	comparison allowableError;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(resultsFnp, "--resultsFnp",
									"results tab delimited file, each row is a haplotype, should have at least 4 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), optionally 4th col with hap sequence (--popHapSeqColName) or read in from --popSeqsFnp",
									true);

	setUp.setOption(sampleColName, "--sampleColName", "sample Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleFreqColName, "--withinSampleFreqColName",
									"within Sample haplotype frequency Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleReadCntColName, "--withinSampleReadCntColName", "within Sample Read Cnt Col Column Name", false, "Results Column Names");
	setUp.setOption(popHapIdColName, "--popHapIdColName", "popHapIdColName", false, "Results Column Names");
	setUp.setOption(popHapSeqColName, "--popHapSeqColName",
									"population Haplotype Sequence Column Name, the seq to compare to expected", false, "Results Column Names");


	setUp.setOption(expectedSeqsFnp, "--expectedSeqsFnp", "Expected Seqs fasta file", true);
	setUp.setOption(name, "--name", "Name to give the current analysis", true);
	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");
	setUp.setOption(popSeqsFnp, "--popSeqsFnp", "Population Sequences");
	setUp.setOption(popSeqsRegexPatRemoval, "--popSeqsRegexPatRemoval",
									"optional regex pattern to process the input pop sequences from --popSeqsFnp to make it match up with the input results folder");

	setUp.setOption(conBenchPars.samplesToMixFnp_, "--sampleToMixture",
									"Sample To Mixture, 2 columns 1)sample, 2)MixName", true);
	setUp.setOption(conBenchPars.mixSetUpFnp_, "--mixtureSetUp",
									"Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);
	setUp.setOption(skipMissingSamples, "--skipMissingSamples", "Skip Samples if they are missing");
	setUp.setOption(fillInMissingSamples, "--fillInMissingSamples", "Fill in missing Samples with placeholders");

	setUp.setOption(columnName, "--newColumnName",
									"Name of a new Column to add to table, can be comma delimited to add multiple");
	setUp.setOption(elementStr, "--newColumnElement",
									"What to Add to the Table Under new Column, can be comma delimited when adding multiple new columns");

	setUp.processComparison(allowableError);
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName(name + "_benchmarkControlMixtures_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow({resultsFnp, expectedSeqsFnp,
																	 conBenchPars.samplesToMixFnp_, conBenchPars.mixSetUpFnp_},
																	__PRETTY_FUNCTION__);

	ControlBencher bencher(conBenchPars);
	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
	std::map<std::string, std::map<std::string, double>> readCountsPerHapPerSample;
	std::map<std::string, std::string> cNameToPopUID;
	std::map<std::string, std::set<std::string>> hPopUIDPopSamps;
	std::map<std::string, std::string> hPopUID_to_hConsensus;
	//population seqs;
	std::vector<seqInfo> popSeqs;

	VecStr requiredColumns{sampleColName, withinSampleReadCntColName,
												 withinSampleFreqColName, withinSampleReadCntColName,
												 popHapIdColName};

	uint64_t maxLen = 0;
	if (exists(popSeqsFnp)) {
		seqInfo seq;
		SeqInput reader{SeqIOOptions(popSeqsFnp, SeqIOOptions::getInFormatFromFnp(popSeqsFnp))};
		reader.openIn();
		while (reader.readNextRead(seq)) {
			readVec::getMaxLength(seq, maxLen);
			seq.name_ = std::regex_replace(seq.name_, std::regex{popSeqsRegexPatRemoval}, "");
			hPopUID_to_hConsensus[seq.name_] = seq.seq_;
			popSeqs.emplace_back(seq);
		}
	} else {
		requiredColumns.emplace_back(popHapSeqColName);
	}

	auto sampInfoFnp = resultsFnp;
	TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));

	sampInfoReader.header_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);

	VecStr row;
	std::unordered_set<std::string> allSamplesInOutput;

	std::map<std::string, std::vector<seqInfo>> allResultSeqs;
	while (sampInfoReader.getNextRow(row)) {
		auto sample = row[sampInfoReader.header_.getColPos(sampleColName)];

		auto c_ReadCnt = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);
		auto c_AveragedFrac = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleFreqColName)]);
		auto readCnt = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);

		auto h_popUID = row[sampInfoReader.header_.getColPos(popHapIdColName)];
		auto hapName = njh::pasteAsStr(sample, "__", h_popUID);
		std::string hapSeq;
		if (popSeqsFnp.empty()) {
			hPopUID_to_hConsensus[h_popUID] = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
			hapSeq = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
		} else {
			hapSeq = njh::mapAt(hPopUID_to_hConsensus, h_popUID);
		}
		seqInfo clus(hapName, hapSeq);
		clus.cnt_ = c_ReadCnt;
		clus.frac_ = c_AveragedFrac;
		bool add = true;
		for (auto &seq: allResultSeqs[sample]) {
			if (seq.seq_ == clus.seq_) {
				add = false;
				seq.cnt_ += c_ReadCnt;
				seq.frac_ += c_AveragedFrac;
				readCountsPerHapPerSample[sample][hapName] += readCnt;
				break;
			}
		}
		if (add) {
			readCountsPerHapPerSample[sample][hapName] = readCnt;
			allResultSeqs[sample].emplace_back(clus);
		}

		readVec::getMaxLength(clus, maxLen);
		cNameToPopUID[hapName] = h_popUID;
		hPopUIDPopSamps[h_popUID].emplace(sample);
		allSamplesInOutput.emplace(sample);
	}
	if (!exists(popSeqsFnp)) {
		for (const auto &popHap: hPopUID_to_hConsensus) {
			popSeqs.emplace_back(popHap.first, popHap.second);
		}
	}


	//check for samples in analysis
	auto controlSamples = bencher.getSamples();
	VecStr missingSamples;
	for (const auto &sname: controlSamples) {
		//if (!analysisMaster.hasSample(sname)) {
		if (!njh::in(sname, allSamplesInOutput)) {
			missingSamples.emplace_back(sname);
		}
	}

	if (!skipMissingSamples && !missingSamples.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
			 << "missing the following samples in " << resultsFnp << ": "
			 << njh::conToStrEndSpecial(missingSamples, ", ", " and ") << "\n";
		throw std::runtime_error{ss.str()};
	} else if (skipMissingSamples && !missingSamples.empty()) {
		OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples.txt"));
		OutputStream outMissing(outOptsMissing);
		outMissing << njh::conToStr(missingSamples, "\n") << std::endl;
	}
	bfs::copy(njh::files::normalize(conBenchPars.samplesToMixFnp_), njh::files::make_path(setUp.pars_.directoryName_, "samplesToMix.tsv"));
	bencher.writeMixSetUpsInSamples(njh::files::make_path(setUp.pars_.directoryName_, "mixSetUps.tsv"));

	//read in expected seqs
	SeqInput expSeqsSeqIn(SeqIOOptions::genFastaIn(expectedSeqsFnp));
	auto initialExpSeqs = expSeqsSeqIn.readAllReadsPtrs<seqInfo>();


	//check the needed expected names
	std::set<std::string> expNames;
	std::map<std::string, uint32_t> initialExpSeqsPositions;
	std::map<std::string, std::string> expSeqsKey;
	std::map<std::string, VecStr> matchingExpectedSeqs;

	std::vector<std::shared_ptr<seqInfo>> expSeqs;

	VecStr removeSeqs;
	for (const auto &expSeq: initialExpSeqs) {
		if (std::all_of(expSeq->seq_.begin(), expSeq->seq_.end(),
										[](const char base) {
											return 'N' == base || 'n' == base;
										})) {
			removeSeqs.emplace_back(expSeq->name_);
			continue;
		}
		if (njh::in(expSeq->name_, expNames)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have expected sequence named " << expSeq->name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		expNames.emplace(expSeq->name_);
		bool found = false;
		for (const auto otherSeqPos: iter::range(expSeqs.size())) {
			const auto &otherSeq = expSeqs[otherSeqPos];
			if (otherSeq->seq_ == expSeq->seq_) {
				found = true;
				otherSeq->name_ += "," + expSeq->name_;
				initialExpSeqsPositions[expSeq->name_] = otherSeqPos;
				break;
			}
		}
		if (!found) {
			//if not found add
			//add position of the expected seq
			initialExpSeqsPositions[expSeq->name_] = expSeqs.size();
			expSeqs.emplace_back(std::make_shared<seqInfo>(*expSeq));
		}
	}
	//remove ref seqs that were zero;
	bencher.removeStrains(removeSeqs);


	std::map<std::string, uint32_t> finalExpSeqsPositions;

	//

	for (const auto expSeqPos: iter::range(expSeqs.size())) {
		readVec::getMaxLength(expSeqs[expSeqPos], maxLen);
		expSeqsKey[expSeqs[expSeqPos]->name_] = expSeqs[expSeqPos]->name_;
		finalExpSeqsPositions[expSeqs[expSeqPos]->name_] = expSeqPos;
	}

	bencher.checkForStrainsThrow(expNames, __PRETTY_FUNCTION__);


	OutputStream falseHaplotypesToExpClassified(
					njh::files::make_path(setUp.pars_.directoryName_, "falseHaplotypesComparedToExpected.tsv"));
	falseHaplotypesToExpClassified
					<< "AnalysisName\tsample\tmix\tInputReadName\treadCnt\tfrac\tRefName\tExpectedRefFreq\tExpectedMajorOrMinor\tIdentityScore\tbestMatchScore\tmismatches\toneBaseIndels\ttwoBaseIndels\tlargeIndels\ttotalErrors";
	OutputStream falseHaplotypesToOtherResultsClassified(
					njh::files::make_path(setUp.pars_.directoryName_, "falseHaplotypesComparedToOthers.tsv"));
	falseHaplotypesToOtherResultsClassified
					<< "AnalysisName\tsample\tmix\tInputReadName\treadCnt\tfrac\tOtherName\tOtherReadCnt\tOtherFrac\tratio\tOtherMajor\tOtherMatchesExpected\tOtherExpectedMatchName\tIdentityScore\tbestMatchScore\tmismatches\toneBaseIndels\ttwoBaseIndels\tlargeIndels\ttotalErrors";


	OutputStream haplotypesClassified(njh::files::make_path(setUp.pars_.directoryName_, "classifiedHaplotypes.tsv"));
	//haplotypesClassified << "AnalysisName\tsample\tmix\tInputReadName\tHapPopUID\tHapSampleCount\treadCnt\tfrac\tmatchExpected\texpectedRef\texpectedFrac\tMajorOrMinor\tmatchingPopulation\tPopName";
	haplotypesClassified
					<< "AnalysisName\tsample\tmix\tInputReadName\tHapPopUID\tHapSampleCount\treadCnt\tfrac\tmatchExpected\texpectedRef\texpectedFrac\tMajorOrMinor";

	OutputStream performanceOut(njh::files::make_path(setUp.pars_.directoryName_, "performancePerTarget.tsv"));
	performanceOut
					<< "AnalysisName\tsample\tmix\ttotalReads\trecoveredExpectedHaps\tfalseHaps\ttotalHaps\ttotalExpectedHaps\texpectedHapRecoveryRate\tfalseHapsRate\texpectedMissing\tRMSE";
	VecStr metalevels;

	VecStr newColumnEleToks;
	VecStr newColumnToks;
	if (!columnName.empty()) {
		newColumnEleToks = tokenizeString(elementStr, ",");
		newColumnToks = tokenizeString(columnName, ",");
		if (newColumnEleToks.size() != newColumnToks.size()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
				 << ", error, when adding multiple columns the number of elements need to match number of new columns" << "\n"
				 << "New Columns #: " << newColumnToks.size() << ", New Elements #: " << newColumnEleToks.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	if (nullptr != metaGroupData) {
		metalevels = getVectorOfMapKeys(metaGroupData->groupData_);
		njh::sort(metalevels);
		for (const auto &meta: metalevels) {
			haplotypesClassified << "\t" << meta;
			performanceOut << "\t" << meta;
			falseHaplotypesToExpClassified << "\t" << meta;
			falseHaplotypesToOtherResultsClassified << "\t" << meta;
		}
	}
	if (!newColumnToks.empty()) {
		for (const auto &col: newColumnToks) {
			haplotypesClassified << "\t" << col;
			performanceOut << "\t" << col;
			falseHaplotypesToExpClassified << "\t" << col;
			falseHaplotypesToOtherResultsClassified << "\t" << col;
		}
	}

	haplotypesClassified << std::endl;
	performanceOut << std::endl;
	falseHaplotypesToExpClassified << std::endl;
	falseHaplotypesToOtherResultsClassified << std::endl;
	aligner alignerObj(maxLen, gapScoringParameters(5, 1, 0, 0, 0, 0));
	alignerObj.weighHomopolymers_ = setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_;
	alignerObj.countEndGaps_ = true;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, false);
	for (const auto &sname: controlSamples) {
		//skip completely missing
		if (skipMissingSamples && njh::in(sname, missingSamples)) {
			continue;
		}
		//read in result sequences
		//auto resultsSeqsFnp = analysisMaster.getSampleFinalHapsPath(sname);
		if (!njh::in(sname, allSamplesInOutput) && skipMissingSamples) {
			OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples.txt"));
			outOptsMissing.append_ = true;
			OutputStream outMissing(outOptsMissing);
			outMissing << sname << std::endl;
			continue;
		} else if (!njh::in(sname, allSamplesInOutput)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "missing results for the following sample: " << sname << "\n";
			throw std::runtime_error{ss.str()};
		}


		std::vector<seqInfo> resultSeqs = allResultSeqs[sname];

		std::map<std::string, uint32_t> resSeqToPos;
		double maxResFrac = 0;
		for (const auto pos: iter::range(resultSeqs.size())) {
			resSeqToPos[resultSeqs[pos].name_] = pos;
			if (resultSeqs[pos].frac_ > maxResFrac) {
				maxResFrac = resultSeqs[pos].frac_;
			}
		}
		//get current expected seqs
		std::map<std::string, double> currentExpectedSeqsFrac;
		double maxExpFrac = 0;
		std::map<std::string, VecStr> expectedSeqNameToCurrentSeqsKey;

		for (const auto &expSeqFrac: bencher.mixSetups_.at(bencher.samplesToMix_.at(sname)).relativeAbundances_) {
			currentExpectedSeqsFrac[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_] += expSeqFrac.second;
			expectedSeqNameToCurrentSeqsKey[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_].emplace_back(
							expSeqFrac.first);
		}
		for (auto &key: expectedSeqNameToCurrentSeqsKey) {
//			std::cout << "key: " << key.first << std::endl;
//			std::cout << "\t" << njh::conToStr(key.second, ",") << std::endl;
			njh::sort(key.second);
		}
		for (const auto &expFrac: currentExpectedSeqsFrac) {
			if (expFrac.second > maxExpFrac) {
				maxExpFrac = expFrac.second;
			}
		}

		std::map<std::string, std::string> expectedToMajorClass;
		for (const auto &expFrac: currentExpectedSeqsFrac) {
			if (expFrac.second == maxExpFrac) {
				expectedToMajorClass[expFrac.first] = "major";
			} else {
				expectedToMajorClass[expFrac.first] = "minor";
			}
		}
		std::map<std::string, std::string> resultsToMajorClass;
		for (const auto &res: resultSeqs) {
			if (res.frac_ == maxResFrac) {
				resultsToMajorClass[res.name_] = "major";
			} else {
				resultsToMajorClass[res.name_] = "minor";
			}
		}


		std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
		currentExpectedSeqs.reserve(currentExpectedSeqsFrac.size());
		for (const auto &finalExp: currentExpectedSeqsFrac) {
			currentExpectedSeqs.emplace_back(expSeqs[finalExpSeqsPositions[finalExp.first]]);
		}

		auto res = njhseq::ControlBencher::benchmark(resultSeqs, currentExpectedSeqs,
																								 currentExpectedSeqsFrac,
																								 expSeqsKey, alignerObj, allowableError);
		double total = 0;


		//haplotype classification
		for (const auto &seq: resultSeqs) {
			total += readCountsPerHapPerSample[sname][seq.name_];
			haplotypesClassified << name
													 << "\t" << sname
													 << "\t" << bencher.samplesToMix_[sname]
													 << "\t" << seq.name_
													 << "\t" << cNameToPopUID[seq.name_]
													 << "\t" << hPopUIDPopSamps[cNameToPopUID[seq.name_]].size()
													 << "\t" << readCountsPerHapPerSample[sname][seq.name_]
													 << "\t" << seq.frac_
													 << "\t" << (res.resSeqToExpSeq_[seq.name_].empty() ? "FALSE" : "TRUE")
													 << "\t"
													 << njh::conToStr(expectedSeqNameToCurrentSeqsKey[res.resSeqToExpSeq_[seq.name_]], ",")
													 << "\t" << currentExpectedSeqsFrac[res.resSeqToExpSeq_[seq.name_]]
													 << "\t" << expectedToMajorClass[res.resSeqToExpSeq_[seq.name_]];

			std::string matchingPop;
			for (const auto &popSeq: popSeqs) {
				if (popSeq.seq_ == seq.seq_) {
					matchingPop = popSeq.name_;
					break;
				}
			}
//			haplotypesClassified
//							<< "\t" << (matchingPop.empty() ? "TRUE" : "FALSE")
//							<< "\t" << matchingPop;
			if (nullptr != metaGroupData) {
				for (const auto &meta: metalevels) {
					haplotypesClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
				}
			}
			if (!newColumnToks.empty()) {
				for (const auto &ele: newColumnEleToks) {
					haplotypesClassified << "\t" << ele;
				}
			}
			haplotypesClassified << std::endl;
		}
		for (const auto &missing: res.missingExpecteds_) {
			haplotypesClassified << name
													 << "\t" << sname
													 << "\t" << bencher.samplesToMix_[sname]
													 << "\t" << "NA"
													 << "\t" << "NA"
													 << "\t" << "NA"
													 << "\t" << "NA"
													 << "\t" << "NA"
													 << "\t" << "NA"
													 << "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ",")
													 << "\t" << currentExpectedSeqsFrac[missing]
													 << "\t" << expectedToMajorClass[missing];

//			haplotypesClassified
//							<< "\t" << "NA"
//							<< "\t" << "NA";
			if (nullptr != metaGroupData) {
				for (const auto &meta: metalevels) {
					haplotypesClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
				}
			}
			if (!newColumnToks.empty()) {
				for (const auto &ele: newColumnEleToks) {
					haplotypesClassified << "\t" << ele;
				}
			}
			haplotypesClassified << std::endl;
		}
		//performance
		VecStr missingExpectedDecoded;
		for (const auto &missing: res.missingExpecteds_) {
			missingExpectedDecoded.emplace_back(njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ","));
		}
		performanceOut << name
									 << "\t" << sname
									 << "\t" << bencher.samplesToMix_[sname]
									 << "\t" << total
									 << "\t" << res.recoveredHaps_
									 << "\t" << res.falseHaps_
									 << "\t" << res.totalHaps()
									 << "\t" << res.expectedHapCnt_
									 << "\t" << res.hapRecoveryRate()
									 << "\t" << res.falseHapRate()
									 << "\t" << njh::conToStr(missingExpectedDecoded, ";")
									 << "\t" << res.RMSE();
		if (nullptr != metaGroupData) {
			for (const auto &meta: metalevels) {
				performanceOut << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
			}
		}
		if (!newColumnToks.empty()) {
			for (const auto &ele: newColumnEleToks) {
				performanceOut << "\t" << ele;
			}
		}
		performanceOut << std::endl;


		//further classification of false haplotypes to expected
		for (const auto &falseHap: res.falseHapsCompsToExpected) {
			double bestScore = 0;
			for (const auto &exp: falseHap.second) {
				if (exp.second.distances_.eventBasedIdentityHq_ > bestScore) {
					bestScore = exp.second.distances_.eventBasedIdentityHq_;
				}
			}
			for (const auto &exp: falseHap.second) {
				falseHaplotypesToExpClassified << name
																			 << "\t" << sname
																			 << "\t" << bencher.samplesToMix_[sname]
																			 << "\t" << falseHap.first
																			 << "\t" << readCountsPerHapPerSample[sname][falseHap.first]
																			 << "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
																			 << "\t" << exp.first
																			 << "\t" << currentExpectedSeqsFrac[exp.first]
																			 << "\t" << expectedToMajorClass[exp.first]
																			 << "\t" << exp.second.distances_.eventBasedIdentityHq_
																			 << "\t"
																			 << (exp.second.distances_.eventBasedIdentityHq_ == bestScore ? "TRUE" : "FALSE")
																			 << "\t" << exp.second.hqMismatches_ + exp.second.lqMismatches_ +
																									exp.second.lowKmerMismatches_
																			 << "\t" << exp.second.oneBaseIndel_
																			 << "\t" << exp.second.twoBaseIndel_
																			 << "\t" << exp.second.largeBaseIndel_
																			 << "\t" << exp.second.distances_.getNumOfEvents(true);
				if (nullptr != metaGroupData) {
					for (const auto &meta: metalevels) {
						falseHaplotypesToExpClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
					}
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						falseHaplotypesToExpClassified << "\t" << ele;
					}
				}
				falseHaplotypesToExpClassified << std::endl;
			}
		}

		//further classification of false haplotypes to other sequences
		for (const auto &falseHap: res.falseHapsCompsToOthers) {
			double bestScore = 0;
			for (const auto &other: falseHap.second) {
				if (other.second.distances_.eventBasedIdentityHq_ > bestScore) {
					bestScore = other.second.distances_.eventBasedIdentityHq_;
				}
			}
			for (const auto &other: falseHap.second) {
				falseHaplotypesToOtherResultsClassified << name
																								<< "\t" << sname
																								<< "\t" << bencher.samplesToMix_[sname]
																								<< "\t" << falseHap.first
																								<< "\t" << readCountsPerHapPerSample[sname][falseHap.first]
																								<< "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
																								<< "\t" << other.first
																								<< "\t" << readCountsPerHapPerSample[sname][other.first]
																								<< "\t" << resultSeqs[resSeqToPos[other.first]].frac_
																								<< "\t" << resultSeqs[resSeqToPos[other.first]].frac_ /
																													 resultSeqs[resSeqToPos[falseHap.first]].frac_
																								<< "\t" << resultsToMajorClass[other.first]
																								<< "\t" << (res.resSeqToExpSeq_[other.first].empty() ? "FALSE" : "TRUE")
																								<< "\t" << res.resSeqToExpSeq_[other.first]
																								<< "\t" << other.second.distances_.eventBasedIdentityHq_
																								<< "\t"
																								<< (other.second.distances_.eventBasedIdentityHq_ == bestScore ? "TRUE"
																																																							 : "FALSE")
																								<< "\t" << other.second.hqMismatches_ + other.second.lqMismatches_ +
																													 other.second.lowKmerMismatches_
																								<< "\t" << other.second.oneBaseIndel_
																								<< "\t" << other.second.twoBaseIndel_
																								<< "\t" << other.second.largeBaseIndel_
																								<< "\t" << other.second.distances_.getNumOfEvents(true);
				if (nullptr != metaGroupData) {
					for (const auto &meta: metalevels) {
						falseHaplotypesToOtherResultsClassified << "\t"
																										<< metaGroupData->groupData_[meta]->getGroupForSample(sname);
					}
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						falseHaplotypesToOtherResultsClassified << "\t" << ele;
					}
				}
				falseHaplotypesToOtherResultsClassified << std::endl;
			}
		}
	}

	if (fillInMissingSamples) {
		for (const auto &sname: missingSamples) {

			std::vector<seqInfo> resultSeqs = allResultSeqs[sname];

			std::map<std::string, uint32_t> resSeqToPos;
			double maxResFrac = 0;
			for (const auto pos: iter::range(resultSeqs.size())) {
				resSeqToPos[resultSeqs[pos].name_] = pos;
				if (resultSeqs[pos].frac_ > maxResFrac) {
					maxResFrac = resultSeqs[pos].frac_;
				}
			}
			//get current expected seqs
			std::map<std::string, double> currentExpectedSeqsFrac;
			double maxExpFrac = 0;
			std::map<std::string, VecStr> expectedSeqNameToCurrentSeqsKey;


			for (const auto &expSeqFrac: bencher.mixSetups_.at(bencher.samplesToMix_.at(sname)).relativeAbundances_) {
				currentExpectedSeqsFrac[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_] += expSeqFrac.second;
				expectedSeqNameToCurrentSeqsKey[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_].emplace_back(
								expSeqFrac.first);
			}
			for (auto &key: expectedSeqNameToCurrentSeqsKey) {
				njh::sort(key.second);
			}
			for (const auto &expFrac: currentExpectedSeqsFrac) {
				if (expFrac.second > maxExpFrac) {
					maxExpFrac = expFrac.second;
				}
			}

			std::map<std::string, std::string> expectedToMajorClass;
			for (const auto &expFrac: currentExpectedSeqsFrac) {
				if (expFrac.second == maxExpFrac) {
					expectedToMajorClass[expFrac.first] = "major";
				} else {
					expectedToMajorClass[expFrac.first] = "minor";
				}
			}
			std::map<std::string, std::string> resultsToMajorClass;
			for (const auto &res: resultSeqs) {
				if (res.frac_ == maxResFrac) {
					resultsToMajorClass[res.name_] = "major";
				} else {
					resultsToMajorClass[res.name_] = "minor";
				}
			}


			std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
			currentExpectedSeqs.reserve(currentExpectedSeqsFrac.size());
			for (const auto &finalExp: currentExpectedSeqsFrac) {
				currentExpectedSeqs.push_back(expSeqs[finalExpSeqsPositions[finalExp.first]]);
			}
			VecStr missingExpectedDecoded;
			for (const auto &missing: currentExpectedSeqs) {
				missingExpectedDecoded.emplace_back(njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing->name_], ","));
			}
			//performance
			performanceOut << name
										 << "\t" << sname
										 << "\t" << bencher.samplesToMix_[sname]
										 << "\t" << 0
										 << "\t" << 0
										 << "\t" << 0
										 << "\t" << 0
										 << "\t" << currentExpectedSeqs.size()
										 << "\t" << 0
										 << "\t" << 0
										 << "\t" << njh::conToStr(missingExpectedDecoded, ";")
										 << "\t" << 0;
			if (nullptr != metaGroupData) {
				for (const auto &meta: metalevels) {
					performanceOut << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
				}
			}
			if (!newColumnToks.empty()) {
				for (const auto &ele: newColumnEleToks) {
					performanceOut << "\t" << ele;
				}
			}
			performanceOut << std::endl;
			for (const auto &missing: readVec::getNames(currentExpectedSeqs)) {
				haplotypesClassified << name
														 << "\t" << sname
														 << "\t" << bencher.samplesToMix_[sname]
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ",")
														 << "\t" << currentExpectedSeqsFrac[missing]
														 << "\t" << expectedToMajorClass[missing];

//				haplotypesClassified
//								<< "\t" << "NA"
//								<< "\t" << "NA";
				if (nullptr != metaGroupData) {
					for (const auto &meta: metalevels) {
						haplotypesClassified << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
					}
				}
				if (!newColumnToks.empty()) {
					for (const auto &ele: newColumnEleToks) {
						haplotypesClassified << "\t" << ele;
					}
				}
				haplotypesClassified << std::endl;
			}
		}
	}

	OutputStream qmdOut(njh::files::make_path(setUp.pars_.directoryName_, "process.qmd"));
	qmdOut << plotSingleTargetResultsQmd << std::endl;
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, false);

	return 0;
}


}  // namespace njhseq

