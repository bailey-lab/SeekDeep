/*
 * SeekDeepUtilsRunner_benchmarkControlMixtures.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/ControlBenchmarking.h"
#include "SeekDeep/objects/ControlBenchmarking/ExpectedControls.hpp"
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/pmo.h>

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
  geom_text(aes(x = label, y = strain, label = paste0(relative_abundance *100, "%"), color = ifelse(relative_abundance *100 > 50, "black", "white")), angle = "-90", fontsize = 12) +
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

	std::set<std::string> selectTargets;
	std::set<std::string> excludeTargets;


	comparison allowableError;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Benchmark the results of many targets, as a results table or a PMO, against known control mixtures: for each control sample and target, classify its haplotypes as expected or false "
			"and report how many expected haplotypes were recovered and how far off their frequencies are, writes classifiedHaplotypes.tsv, performancePerTarget.tsv and the false haplotype tables";
	bfs::path pmoFnp;
	setUp.setOption(resultsFnp, "--resultsFnp",
									"results tab delimited file, each row is a haplotype, should have at least 5 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), 5)target name column (--targetNameColName), and a 6th col with the hap sequence (--popHapSeqColName) unless read in from --popSeqsDirFnp, this or --pmoFnp is required");
	setUp.setOption(pmoFnp, "--pmoFnp",
									"results as a PMO (portable microhaplotype object) json file, can be gzipped, samples are the library sample names, haplotypes are the representative microhaplotypes (named by microhaplotype_name if present) and within sample freq is the haplotype's umis (or reads, see --useReads) over the sample's total for that target, use instead of --resultsFnp");
	bool useReads = false;
	setUp.setOption(useReads, "--useReads",
									"with --pmoFnp, use the detected microhaplotypes' reads rather than their umis for the counts and within sample freqs, umis are used by default when the PMO has them");
	setUp.setOption(selectTargets, "--selectTargets", "only analyze these select targets", false);
	setUp.setOption(excludeTargets, "--excludeTargets", "exclude these Targets", false);

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
									"Expected Seqs directory of fasta files, in this directory should be a fasta file with the name of each target, needed unless using --expectedPmoFnp");
	bfs::path strainAllelesFnp;
	setUp.setOption(strainAllelesFnp, "--strainAllelesFnp",
									"Optional, for strains with more than one copy or allele at a target (e.g. a copy number variant), tab delimited with columns target, strain, allele and optionally copies (copies of that allele per genome, default 1), the alleles are named in the --expectedSeqsDirFnp fasta, each allele's share of the template is its strain's relative abundance times its copies, strains not listed have one copy named after the strain");
	uint32_t trimFront = 0;
	uint32_t trimBack = 0;
	setUp.setOption(trimFront, "--trimFront", "trim this many bases from the front of every expected seq from --expectedSeqsDirFnp, e.g. to match a pipeline that trims the edges of the insert");
	setUp.setOption(trimBack, "--trimBack", "trim this many bases from the back of every expected seq from --expectedSeqsDirFnp, e.g. to match a pipeline that trims the edges of the insert");
	bfs::path expectedPmoFnp;
	setUp.setOption(expectedPmoFnp, "--expectedPmoFnp",
									"What each control sample is expected to have as a PMO (can be gzipped), use instead of --sampleToMixture, --mixtureSetUp and --expectedSeqsDirFnp, library samples are the control samples, their specimen names the mixture names and their detected microhaplotypes the expected alleles with reads as pseudo counts for the within sample frequencies, make one with mixturesToExpectedPmo");
	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");
	setUp.setOption(popSeqsRegexPatRemoval, "--popSeqsRegexPatRemoval",
									"optional regex pattern to process the input pop sequences from --popSeqsDirFnp to make it match up with the input results folder");
	if (resultsFnp.empty() == pmoFnp.empty()) {
		setUp.failed_ = true;
		setUp.addWarning("need exactly one of --resultsFnp or --pmoFnp");
	}
	if (!pmoFnp.empty() && !popSeqsDirFnp.empty()) {
		setUp.failed_ = true;
		setUp.addWarning("--popSeqsDirFnp can't be used with --pmoFnp, the PMO already has the haplotype sequences");
	}
	if (useReads && pmoFnp.empty()) {
		setUp.failed_ = true;
		setUp.addWarning("--useReads only applies to --pmoFnp, a results table's counts come from --withinSampleReadCntColName");
	}

	setUp.setOption(conBenchPars.samplesToMixFnp_, "--sampleToMixture",
									"Sample To Mixture, 2 columns 1)sample, 2)MixName, needed unless using --expectedPmoFnp");
	setUp.setOption(conBenchPars.mixSetUpFnp_, "--mixtureSetUp",
									"Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance, needed unless using --expectedPmoFnp");
	if (expectedPmoFnp.empty()) {
		if (expectedSeqsDirFnp.empty() || conBenchPars.samplesToMixFnp_.empty() || conBenchPars.mixSetUpFnp_.empty()) {
			setUp.failed_ = true;
			setUp.addWarning("need --expectedPmoFnp or all of --sampleToMixture, --mixtureSetUp and --expectedSeqsDirFnp");
		}
	} else if (!expectedSeqsDirFnp.empty() || !conBenchPars.samplesToMixFnp_.empty() ||
						 !conBenchPars.mixSetUpFnp_.empty() || !strainAllelesFnp.empty() || trimFront > 0 || trimBack > 0) {
		setUp.failed_ = true;
		setUp.addWarning("--expectedPmoFnp can't be used with --sampleToMixture, --mixtureSetUp, --expectedSeqsDirFnp, --strainAllelesFnp, --trimFront or --trimBack, trim when making the expected PMO with mixturesToExpectedPmo");
	}
	setUp.setOption(skipMissingSamples, "--skipMissingSamples", "Skip control samples that have no results for a target instead of stopping, they're listed in missingSamples.txt");
	setUp.setOption(fillInMissingSamples, "--fillInMissingSamples", "Add control samples that have no results for a target to the output with 0 recovered haplotypes, needs --skipMissingSamples");

	setUp.setOption(columnName, "--newColumnName",
									"Name of a new Column to add to table, can be comma delimited to add multiple");
	setUp.setOption(elementStr, "--newColumnElement",
									"What to Add to the Table Under new Column, can be comma delimited when adding multiple new columns");

	setUp.processComparison(allowableError);
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName("benchmarkMultiTarAmpControlMixtures_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	const bfs::path inputFnp = pmoFnp.empty() ? resultsFnp : pmoFnp;
	if (expectedPmoFnp.empty()) {
		njh::files::checkExistenceThrow({inputFnp, expectedSeqsDirFnp,
																		 conBenchPars.samplesToMixFnp_, conBenchPars.mixSetUpFnp_},
																		__PRETTY_FUNCTION__);
		if (!strainAllelesFnp.empty()) {
			njh::files::checkExistenceThrow(std::vector<bfs::path>{strainAllelesFnp}, __PRETTY_FUNCTION__);
		}
	} else {
		njh::files::checkExistenceThrow(std::vector<bfs::path>{inputFnp, expectedPmoFnp}, __PRETTY_FUNCTION__);
	}


	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double>>> readCountsPerHapPerSample;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> cNameToPopUID;
	std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> hPopUIDPopSamps;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hPopUID_to_hConsensus;

	std::set<std::string> targetNamesSet;

	uint64_t maxLen = 0;
	//population seqs;
	std::unordered_map<std::string, std::vector<seqInfo>> popSeqs;
	std::unordered_map<std::string, std::unordered_set<std::string>> allSamplesInOutput;
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> allResultSeqs;

	auto passTargetFilter = [&selectTargets, &excludeTargets](const std::string& target) {
		return (selectTargets.empty() && excludeTargets.empty()) ||
		       (!selectTargets.empty() && njh::in(target, selectTargets)) ||
		       (!excludeTargets.empty() && njh::notIn(target, excludeTargets));
	};

	//add one haplotype call for a sample, calls with the same sequence in the same sample are combined
	auto addHap = [&](const std::string & target, const std::string & sample, const std::string & h_popUID,
										const std::string & hapSeq, double readCnt, double frac) {
		targetNamesSet.emplace(target);
		auto hapName = njh::pasteAsStr(sample, "__", h_popUID);
		seqInfo clus(hapName, hapSeq);
		clus.cnt_ = readCnt;
		clus.frac_ = frac;
		bool add = true;
		for (auto &seq: allResultSeqs[target][sample]) {
			if (seq.seq_ == clus.seq_) {
				add = false;
				seq.cnt_ += readCnt;
				seq.frac_ += frac;
				readCountsPerHapPerSample[target][sample][seq.name_] += readCnt;
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
	};

	if (!pmoFnp.empty()) {
		auto pmoObj = pmo::readPmo(pmoFnp);
		//a target can be in representative_microhaplotypes more than once, mhaps_target_id is the position in that list and mhap_id the position within that entry
		std::vector<std::string> targetNameByMhapsId;
		for (const auto &representatives: pmoObj.representative_microhaplotypes_.targets_) {
			const auto &target = pmoObj.target_info_.at(representatives.target_id_).target_name_;
			targetNameByMhapsId.emplace_back(target);
			if (!passTargetFilter(target)) {
				continue;
			}
			for (const auto mhapPos: iter::range(representatives.microhaplotypes_.size())) {
				const auto &mhap = representatives.microhaplotypes_[mhapPos];
				auto name = mhap.microhaplotype_name_.has_value()
										? mhap.microhaplotype_name_.value()
										: njh::pasteAsStr(target, ".", targetNameByMhapsId.size() - 1, ".", mhapPos);
				hPopUID_to_hConsensus[target][name] = mhap.seq_;
			}
		}
		auto mhapName = [&](uint32_t mhapsTargetId, uint32_t mhapId) {
			const auto &mhap = pmoObj.representative_microhaplotypes_.targets_.at(mhapsTargetId).microhaplotypes_.at(mhapId);
			return mhap.microhaplotype_name_.has_value()
						 ? mhap.microhaplotype_name_.value()
						 : njh::pasteAsStr(targetNameByMhapsId.at(mhapsTargetId), ".", mhapsTargetId, ".", mhapId);
		};
		//counts are umis if any detected microhaplotype has them (then all of them need to) unless --useReads, otherwise reads
		auto hasUmis = [](const pmo::MicrohaplotypeForTarget & mhap) {
			return mhap.umis_.has_value() && mhap.umis_.value() != std::numeric_limits<uint32_t>::max();
		};
		bool useUmis = false;
		if (!useReads) {
			uint64_t withUmis = 0;
			uint64_t withoutUmis = 0;
			for (const auto &run: pmoObj.detected_microhaplotypes_) {
				for (const auto &sampleResults: run.library_samples_) {
					for (const auto &targetResults: sampleResults.target_results_) {
						if (!passTargetFilter(targetNameByMhapsId.at(targetResults.mhaps_target_id_))) {
							continue;
						}
						for (const auto &mhap: targetResults.mhaps_) {
							if (hasUmis(mhap)) {
								++withUmis;
							} else {
								++withoutUmis;
							}
						}
					}
				}
			}
			if (withUmis > 0 && withoutUmis > 0) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << withUmis << " detected microhaplotypes in " << pmoFnp << " have umis but "
					 << withoutUmis << " don't (missing or NA), can't mix umis and reads, use --useReads to use reads for all of them" << "\n";
				throw std::runtime_error{ss.str()};
			}
			useUmis = withUmis > 0;
		}
		auto mhapCount = [&useUmis](const pmo::MicrohaplotypeForTarget & mhap) -> double {
			return useUmis ? mhap.umis_.value() : mhap.reads_;
		};
		std::cout << "using " << (useUmis ? "umis" : "reads") << " from " << pmoFnp << " for haplotype counts and within sample freqs" << std::endl;
		//only the first set of results for a sample and target is used, later ones are counted and reported as a warning
		std::map<std::string, std::set<std::string>> targetsDoneBySample;
		std::map<std::string, uint32_t> repeatedResults;
		for (const auto &run: pmoObj.detected_microhaplotypes_) {
			for (const auto &sampleResults: run.library_samples_) {
				const auto &sample = pmoObj.library_sample_info_.at(sampleResults.library_sample_id_).library_sample_name_;
				for (const auto &targetResults: sampleResults.target_results_) {
					const auto &target = targetNameByMhapsId.at(targetResults.mhaps_target_id_);
					if (!passTargetFilter(target)) {
						continue;
					}
					if (!targetsDoneBySample[sample].emplace(target).second) {
						++repeatedResults[sample];
						continue;
					}
					double total = 0;
					for (const auto &mhap: targetResults.mhaps_) {
						total += mhapCount(mhap);
					}
					for (const auto &mhap: targetResults.mhaps_) {
						auto name = mhapName(targetResults.mhaps_target_id_, mhap.mhap_id_);
						addHap(target, sample, name, njh::mapAt(njh::mapAt(hPopUID_to_hConsensus, target), name),
									 mhapCount(mhap), total > 0 ? mhapCount(mhap) / total : 0);
					}
				}
			}
		}
		for (const auto &repeated: repeatedResults) {
			std::cerr << __PRETTY_FUNCTION__ << ", warning, " << repeated.first << " has more than one set of results for "
								<< repeated.second << " targets in " << pmoFnp << ", only using the first of each" << std::endl;
		}
	} else {
		//population seqs from --popSeqsDirFnp are read the first time their target comes up in the table
		std::set<std::string> popSeqsLoaded;
		auto loadPopSeqs = [&](const std::string & target) {
			if (!popSeqsLoaded.emplace(target).second) {
				return;
			}
			auto popSeqsFnp = njh::files::make_path(popSeqsDirFnp, target + ".fasta");
			if (!exists(popSeqsFnp)) {
				popSeqsFnp = njh::files::make_path(popSeqsDirFnp, target + ".fasta.gz");
			}
			if (!exists(popSeqsFnp)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "error, could not find population seqs file for target: " << target
					 << ", should be " << njh::files::make_path(popSeqsDirFnp, target + ".fasta") << " or "
					 << njh::files::make_path(popSeqsDirFnp, target + ".fasta.gz") << "\n";
				throw std::runtime_error{ss.str()};
			}
			SeqInput reader{SeqIOOptions(popSeqsFnp, SeqIOOptions::getInFormatFromFnp(popSeqsFnp))};
			reader.openIn();
			seqInfo seq;
			while (reader.readNextRead(seq)) {
				readVec::getMaxLength(seq, maxLen);
				seq.name_ = std::regex_replace(seq.name_, std::regex{popSeqsRegexPatRemoval}, "");
				hPopUID_to_hConsensus[target][seq.name_] = seq.seq_;
			}
		};
		VecStr requiredColumns{sampleColName, withinSampleReadCntColName,
													 withinSampleFreqColName, withinSampleReadCntColName,
													 popHapIdColName,
													 targetNameColName};
		if (!exists(popSeqsDirFnp)) {
			requiredColumns.emplace_back(popHapSeqColName);
		}
		TableReader sampInfoReader(TableIOOpts::genTabFileIn(resultsFnp, true));
		sampInfoReader.header_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);
		VecStr row;
		while (sampInfoReader.getNextRow(row)) {
			auto target = row[sampInfoReader.header_.getColPos(targetNameColName)];
			//check if only including select targets or excluding targets
			if (!passTargetFilter(target)) {
				continue;
			}
			auto sample = row[sampInfoReader.header_.getColPos(sampleColName)];
			auto readCnt = njh::StrToNumConverter::stoToNum<double>(
							row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);
			auto frac = njh::StrToNumConverter::stoToNum<double>(
							row[sampInfoReader.header_.getColPos(withinSampleFreqColName)]);
			auto h_popUID = row[sampInfoReader.header_.getColPos(popHapIdColName)];
			std::string hapSeq;
			if (popSeqsDirFnp.empty()) {
				hapSeq = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
				hPopUID_to_hConsensus[target][h_popUID] = hapSeq;
			} else {
				loadPopSeqs(target);
				const auto &targetSeqs = hPopUID_to_hConsensus[target];
				auto search = targetSeqs.find(h_popUID);
				if (targetSeqs.end() == search) {
					VecStr foundNames;
					for (const auto &popSeq: targetSeqs) {
						foundNames.emplace_back(popSeq.first);
					}
					njh::naturalSortNameSet(foundNames);
					if (foundNames.size() > 10) {
						foundNames.resize(10);
						foundNames.emplace_back("...");
					}
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "couldn't find " << h_popUID << " for target " << target
						 << " in " << popSeqsDirFnp << " after removing --popSeqsRegexPatRemoval " << popSeqsRegexPatRemoval
						 << " from the names, found: " << njh::conToStr(foundNames, ", ") << "\n";
					throw std::runtime_error{ss.str()};
				}
				hapSeq = search->second;
			}
			addHap(target, sample, h_popUID, hapSeq, readCnt, frac);
		}
	}
	for (const auto &tarPopHaps: hPopUID_to_hConsensus) {
		for (const auto &popHap: tarPopHaps.second) {
			popSeqs[tarPopHaps.first].emplace_back(popHap.first, popHap.second);
		}
	}
	if(targetNamesSet.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "targetNamesSet is empty, if excluding or selecting only specific targets there may now be no targets" << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto targetNames = std::vector<std::string>(targetNamesSet.begin(), targetNamesSet.end());
	njh::naturalSortNameSet(targetNames);

	//what each control sample is expected to have
	ExpectedControls expected;
	if (!expectedPmoFnp.empty()) {
		expected = ExpectedControls::fromPmo(pmo::readPmo(expectedPmoFnp));
		VecStr expectedTargets;
		for (const auto &target: expected.targetNames()) {
			if (passTargetFilter(target)) {
				expectedTargets.emplace_back(target);
			} else {
				expected.targets_.erase(target);
			}
		}
		VecStr noExpectations;
		for (const auto &target: targetNames) {
			if (!njh::in(target, expectedTargets)) {
				noExpectations.emplace_back(target);
			}
		}
		if (!noExpectations.empty()) {
			std::cerr << __PRETTY_FUNCTION__ << ", warning, " << "no expected results in " << expectedPmoFnp
								<< " for the following targets so they won't be benchmarked: " << njh::conToStr(noExpectations, ",")
								<< std::endl;
		}
		targetNames = expectedTargets;
	} else {
		ExpectedControls::FromMixturesPars mixPars;
		mixPars.samplesToMixFnp_ = conBenchPars.samplesToMixFnp_;
		mixPars.mixSetUpFnp_ = conBenchPars.mixSetUpFnp_;
		mixPars.expectedSeqsDirFnp_ = expectedSeqsDirFnp;
		mixPars.strainAllelesFnp_ = strainAllelesFnp;
		mixPars.targetNameColName_ = targetNameColName;
		mixPars.trimFront_ = trimFront;
		mixPars.trimBack_ = trimBack;
		//if the mixture tables are split by target, use their targets
		auto tableTargets = ExpectedControls::mixtureTableTargets(mixPars);
		if (!tableTargets.empty()) {
			targetNames.clear();
			for (const auto &target: tableTargets) {
				if (passTargetFilter(target)) {
					targetNames.emplace_back(target);
				}
			}
		}
		expected = ExpectedControls::fromMixtures(mixPars, targetNames);
	}
	if (targetNames.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no targets to benchmark, if excluding or selecting only specific targets there may now be no targets" << "\n";
		throw std::runtime_error{ss.str()};
	}
	maxLen = std::max<uint64_t>(maxLen, expected.maxSeqLen());

	//check for samples in analysis
	std::map<std::string, std::set<std::string>> missingSamples;
	for (const auto &target: targetNames) {
		for (const auto &sname: expected.samples(target)) {
			if (!njh::in(sname, allSamplesInOutput[target])) {
				missingSamples[target].emplace(sname);
			}
		}
	}

	if (!skipMissingSamples && !missingSamples.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
			 << "missing the following samples in " << inputFnp << " for targets" << "\n";
		for (const auto &missing: missingSamples) {
			ss << "tar: " << missing.first << ", samples: " << njh::conToStrEndSpecial(missing.second, ", ", " and ") << "\n";
		}
		throw std::runtime_error{ss.str()};
	} else if (skipMissingSamples && !missingSamples.empty()) {
		OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples.txt"));
		OutputStream outMissing(outOptsMissing);
		outMissing << "target\tsample" << std::endl;
		for (const auto &target: targetNames) {
			if (!njh::in(target, missingSamples)) {
				continue;
			}
			for (const auto &samp: missingSamples.at(target)) {
				outMissing << target << "\t" << samp << std::endl;
			}
		}
	}

	if (expectedPmoFnp.empty()) {
		bfs::copy(njh::files::normalize(conBenchPars.samplesToMixFnp_), njh::files::make_path(setUp.pars_.directoryName_, "samplesToMix.tsv"));
	} else {
		OutputStream samplesToMixOut(njh::files::make_path(setUp.pars_.directoryName_, "samplesToMix.tsv"));
		expected.genSamplesToMixTab().outPutContents(samplesToMixOut, "\t");
	}
	{
		OutputStream mixOut(njh::files::make_path(setUp.pars_.directoryName_, "mixSetUps.tsv"));
		expected.genMixSetUpsTab(targetNameColName).outPutContents(mixOut, "\t");
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

	//the meta data and new columns that end every row
	auto addExtraColumns = [&](std::ostream &out, const std::string &sname) {
		if (nullptr != metaGroupData) {
			for (const auto &meta: metalevels) {
				out << "\t" << metaGroupData->groupData_[meta]->getGroupForSample(sname);
			}
		}
		for (const auto &ele: newColumnEleToks) {
			out << "\t" << ele;
		}
		out << std::endl;
	};

	//per sample expected haplotypes, the expected alleles with the same sequence combined
	struct CurrentExpected {
		std::vector<seqInfo> seqs_;
		std::map<std::string, double> fracs_;
		std::map<std::string, std::string> majorClass_;
		std::map<std::string, std::string> nameKey_;

		[[nodiscard]] double frac(const std::string &name) const {
			auto search = fracs_.find(name);
			return fracs_.end() == search ? 0 : search->second;
		}

		[[nodiscard]] std::string majorOrMinor(const std::string &name) const {
			auto search = majorClass_.find(name);
			return majorClass_.end() == search ? std::string("") : search->second;
		}
	};
	auto getCurrentExpected = [&expected](const std::string &target, const std::string &sname) {
		CurrentExpected ret;
		double maxExpFrac = 0;
		for (const auto &hap: expected.targets_.at(target).at(sname).collapsed()) {
			ret.seqs_.emplace_back(hap.name_, hap.seq_);
			ret.fracs_[hap.name_] = hap.frac_;
			ret.nameKey_[hap.name_] = hap.name_;
			maxExpFrac = std::max(maxExpFrac, hap.frac_);
		}
		for (const auto &expFrac: ret.fracs_) {
			ret.majorClass_[expFrac.first] = expFrac.second == maxExpFrac ? "major" : "minor";
		}
		return ret;
	};

	aligner alignerObj(maxLen, gapScoringParameters(5, 1, 0, 0, 0, 0));
	alignerObj.weighHomopolymers_ = setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_;
	alignerObj.countEndGaps_ = setUp.pars_.colOpts_.alignOpts_.countEndGaps_;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, false);
	for (const auto &target: targetNames) {
		if (setUp.pars_.verbose_) {
			std::cout << "target: " << target << std::endl;
		}

		for (const auto &sname: expected.samples(target)) {
			//missing samples are only here with --skipMissingSamples, --fillInMissingSamples adds them below
			if (njh::in(sname, missingSamples[target])) {
				continue;
			}
			const auto &mix = expected.targets_.at(target).at(sname).mix_;
			const std::vector<seqInfo> &resultSeqs = allResultSeqs[target][sname];

			std::map<std::string, uint32_t> resSeqToPos;
			double maxResFrac = 0;
			for (const auto pos: iter::range(resultSeqs.size())) {
				resSeqToPos[resultSeqs[pos].name_] = pos;
				if (resultSeqs[pos].frac_ > maxResFrac) {
					maxResFrac = resultSeqs[pos].frac_;
				}
			}
			std::map<std::string, std::string> resultsToMajorClass;
			for (const auto &res: resultSeqs) {
				resultsToMajorClass[res.name_] = res.frac_ == maxResFrac ? "major" : "minor";
			}

			auto current = getCurrentExpected(target, sname);
			auto res = njhseq::ControlBencher::benchmark(resultSeqs, current.seqs_, current.fracs_, current.nameKey_,
																									 alignerObj, allowableError);
			double total = 0;

			//haplotype classification
			for (const auto &seq: resultSeqs) {
				total += readCountsPerHapPerSample[target][sname][seq.name_];
				const auto &matchingExpected = res.resSeqToExpSeq_[seq.name_];
				haplotypesClassified << target
														 << "\t" << sname
														 << "\t" << mix
														 << "\t" << seq.name_
														 << "\t" << cNameToPopUID[target][seq.name_]
														 << "\t" << hPopUIDPopSamps[target][cNameToPopUID[target][seq.name_]].size()
														 << "\t" << readCountsPerHapPerSample[target][sname][seq.name_]
														 << "\t" << seq.frac_
														 << "\t" << (matchingExpected.empty() ? "FALSE" : "TRUE")
														 << "\t" << matchingExpected
														 << "\t" << current.frac(matchingExpected)
														 << "\t" << current.majorOrMinor(matchingExpected);
				addExtraColumns(haplotypesClassified, sname);
			}
			for (const auto &missing: res.missingExpecteds_) {
				haplotypesClassified << target
														 << "\t" << sname
														 << "\t" << mix
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << "NA"
														 << "\t" << missing
														 << "\t" << current.frac(missing)
														 << "\t" << current.majorOrMinor(missing);
				addExtraColumns(haplotypesClassified, sname);
			}
			//performance
			performanceOut << target
										 << "\t" << sname
										 << "\t" << mix
										 << "\t" << total
										 << "\t" << res.recoveredHaps_
										 << "\t" << res.falseHaps_
										 << "\t" << res.totalHaps()
										 << "\t" << res.expectedHapCnt_
										 << "\t" << res.hapRecoveryRate()
										 << "\t" << res.falseHapRate()
										 << "\t" << njh::conToStr(res.missingExpecteds_, ";")
										 << "\t" << res.RMSE();
			addExtraColumns(performanceOut, sname);

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
																				 << "\t" << mix
																				 << "\t" << falseHap.first
																				 << "\t" << readCountsPerHapPerSample[target][sname][falseHap.first]
																				 << "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
																				 << "\t" << exp.first
																				 << "\t" << current.frac(exp.first)
																				 << "\t" << current.majorOrMinor(exp.first)
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
					addExtraColumns(falseHaplotypesToExpClassified, sname);
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
																									<< "\t" << mix
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
					addExtraColumns(falseHaplotypesToOtherResultsClassified, sname);
				}
			}
		}

		if (fillInMissingSamples) {
			for (const auto &sname: missingSamples[target]) {
				const auto &mix = expected.targets_.at(target).at(sname).mix_;
				auto current = getCurrentExpected(target, sname);
				auto expectedNames = readVec::getNames(current.seqs_);
				//performance
				performanceOut << target
											 << "\t" << sname
											 << "\t" << mix
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << current.seqs_.size()
											 << "\t" << 0
											 << "\t" << 0
											 << "\t" << njh::conToStr(expectedNames, ";")
											 << "\t" << 0;
				addExtraColumns(performanceOut, sname);
				for (const auto &missing: expectedNames) {
					haplotypesClassified << target
															 << "\t" << sname
															 << "\t" << mix
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << "NA"
															 << "\t" << missing
															 << "\t" << current.frac(missing)
															 << "\t" << current.majorOrMinor(missing);
					addExtraColumns(haplotypesClassified, sname);
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
  geom_text(aes(x = label, y = strain, label = paste0(relative_abundance *100, "%"), color = ifelse(relative_abundance *100 > 50, "black", "white")), angle = "-90", fontsize = 12) +
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
	setUp.description_ = "Benchmark one target's results table (one row per sample haplotype, e.g. processClusters selectedClustersInfo) against known control mixtures, "
			"the single target version of benchmarkMultiTarAmpControlMixtures, writes the same classified haplotype and performance tables";
	setUp.setOption(resultsFnp, "--resultsFnp",
									"results tab delimited file, each row is a haplotype, should have at least 4 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), and a 5th col with the hap sequence (--popHapSeqColName) unless read in from --popSeqsFnp",
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

