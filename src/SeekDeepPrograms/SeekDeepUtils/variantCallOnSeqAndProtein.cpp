//
// Created by Nicholas Hathaway on 12/16/23.
//
#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {

int SeekDeepUtilsRunner::variantCallOnSeqAndProtein(
				const njh::progutils::CmdArgs &inputCommands) {
	bfs::path resultsFnp;
	std::string sampleColName = "s_Sample";
	std::string withinSampleReadCntColName = "c_ReadCnt";
	std::string popHapIdColName = "h_popUID";
	std::string popHapSeqColName = "h_Consensus";

	std::string targetNameColName = "p_name";

	bfs::path popSeqsDirFnp = "";

	bfs::path metaFnp = "";
	VCFOutput::comnbineVCFsPars combiningVcfPars;
	// bool doNotRescueVariantCallsAcrossTargets = false;
	std::string popSeqsRegexPatRemoval = R"(_([tf])?\d+(\.\d+)?$)";
	uint32_t numThreads = 1;

	CollapseAndCallVariantsPars collapseVarCallPars;
	std::set<std::string> selectTargets;
	std::set<std::string> selectSamples;
	std::set<std::string> excludeTargets;
	std::set<std::string> excludeSamples;

	bfs::path bedLocs;
	bool genomicLocsChangePeriodToDash = false;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(collapseVarCallPars.variantCallerRunPars.ploidy, "--ploidy", "Ploidy to force for the sample for the vcf files");
	combiningVcfPars.ploidy = collapseVarCallPars.variantCallerRunPars.ploidy;
	setUp.setOption(combiningVcfPars.doNotRescueVariantCallsAcrossTargets, "--doNotRescueVariantCallsAcrossTargets", "do Not Rescue Variant Calls Across Targets");
	setUp.setOption(combiningVcfPars.combinedOverlappingCallsAcrossTargets, "--combineOverlappingCallsAcrossTargets", "Rather than taking the best variant call for overlapping targets, sum them instead");

	setUp.setOption(excludeTargets, "--excludeTargets", "exclude these select targets from analyze");
	setUp.setOption(excludeSamples, "--excludeSamples", "exclude these select samples from analyze");
	setUp.setOption(selectTargets, "--selectTargets", "Only analyze these select targets");
	setUp.setOption(selectSamples, "--selectSamples", "Only analyze these select samples");
	setUp.setOption(bedLocs, "--genomicLocations", "a bed file with specific genomic locations to align to, location name needs to match target name");
	setUp.setOption(genomicLocsChangePeriodToDash, "--genomicLocationsChangePeriodToDash", "when supplying a location name, change periods to dashes in the name");

	setUp.setOption(collapseVarCallPars.exportLabIsolateSeqs, "--exportLabIsolateSeqs", "Export Lab Isolate Seqs, will export seqs with meta of site==LabIsolate");
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.occurrenceCutOff, "--variantOccurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.lowVariantCutOff, "--variantFrequencyCutOff", "Low Variant Cut Off, don't report variants below this frequency");
	collapseVarCallPars.variantCallerRunPars.totalReadDepthCutOff = 10;
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.totalReadDepthCutOff, "--variantTotalReadDepthCutOff", "Low Variant Total Read Depth Cut Off, don't report variants that have a summed total read depth less than this across all samples");
	collapseVarCallPars.calcPopMeasuresPars.lowVarFreq = collapseVarCallPars.variantCallerRunPars.lowVariantCutOff;
	collapseVarCallPars.transPars.setOptions(setUp, true);
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(collapseVarCallPars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	collapseVarCallPars.calcPopMeasuresPars.diagAlnPairwiseComps = !collapseVarCallPars.noDiagAlnPairwiseComps;
	//setOption(collapseVarCallPars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.numThreads, "--mappingNumThreads", "Number of threads to use for the alignment portion");
	setUp.setOption(collapseVarCallPars.metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "meta Fields To Calc Pop Diffs");
	setUp.setOption(collapseVarCallPars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

	setUp.setOption(collapseVarCallPars.variantCallerRunPars.complexVarPars.withinDist, "--complexVariantWithinDist", "The distance within to link complex variants");
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.complexVarPars.fractionOfCoveredSamples, "--complexVariantFractionOfCoveredSamples", "complex Variant Fraction Of Covered Samples");


	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");


	setUp.setOption(resultsFnp, "--resultsFnp",
									"results tab delimited file, each row is a haplotype, should have at least 5 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), 5)target name column (--targetNameColName), optionally 4th col with hap sequence (--popHapSeqColName) or read in from --popSeqsDirFnp",
									true);

	setUp.setOption(sampleColName, "--sampleColName", "sample Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleReadCntColName, "--withinSampleReadCntColName", "within Sample Read Cnt Col Column Name", false, "Results Column Names");
	setUp.setOption(popHapIdColName, "--popHapIdColName", "popHapIdColName", false, "Results Column Names");
	setUp.setOption(popHapSeqColName, "--popHapSeqColName",
									"population Haplotype Sequence Column Name, the seq to call variants on", false, "Results Column Names");
	setUp.setOption(targetNameColName, "--targetNameColName",
									"target Name Column Name, the column name in the table which indicates the different targets", false, "Results Column Names");

	setUp.setOption(popSeqsDirFnp, "--popSeqsDirFnp",
									"Population Sequences, in this directory should be a fasta file with the name of each target");

	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");
	setUp.setOption(popSeqsRegexPatRemoval, "--popSeqsRegexPatRemoval",
									"optional regex pattern to process the input pop sequences from --popSeqsFnp to make it match up with the input results folder");

	setUp.processDirectoryOutputName("variantCalling_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow({resultsFnp},
																	__PRETTY_FUNCTION__);

	Json::Value runLog;
	std::mutex runLogMut;
	njh::stopWatch fullWatch;
	fullWatch.setLapName("initial set up");

	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hPopUID_to_hConsensus;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hConsensus_to_hPopUID;

	std::set<std::string> targetNamesSet;
	std::set<std::string> inputSampleNamesSet;
	VecStr requiredColumns{sampleColName, withinSampleReadCntColName,
													withinSampleReadCntColName,
												 popHapIdColName,
												 targetNameColName};
	if (!exists(popSeqsDirFnp)) {
		requiredColumns.emplace_back(popHapSeqColName);
	}
	uint64_t maxLen = 0;
	//population seqs;
	auto sampInfoFnp = resultsFnp;

	//key1 == target, key2 == sample
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> allResultSeqs;
	//key1 == target, key2 == seq, value = popUID
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> seqToPopHapUID;

	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
	{
		TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));
		sampInfoReader.header_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);
		VecStr row;
		while (sampInfoReader.getNextRow(row)) {
			auto target = row[sampInfoReader.header_.getColPos(targetNameColName)];
			//filter to just the select targets if filtering for that
			if((!selectTargets.empty() && njh::notIn(target, selectTargets)) ||
				 (!excludeTargets.empty() && njh::in(target, excludeTargets))) {
				continue;
			}

			const auto& sample = row[sampInfoReader.header_.getColPos(sampleColName)];
			//filter to just the select samples if filtering for that
			if((!selectSamples.empty() && njh::notIn(sample, selectSamples)) ||
				 (!excludeSamples.empty() && njh::in(sample, excludeSamples))) {
				continue;
			}
			targetNamesSet.emplace(target);
		}
	}
	if (exists(popSeqsDirFnp)) {
		seqInfo seq;
		for (const auto &target: targetNamesSet) {
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
			if (!bfs::exists(popSeqsFnp)) {
				popSeqsFnp = njh::files::make_path(popSeqsDirFnp, target + ".fasta.gz");
			}
			SeqInput reader{SeqIOOptions(popSeqsFnp, SeqIOOptions::getInFormatFromFnp(popSeqsFnp))};
			reader.openIn();
			while (reader.readNextRead(seq)) {
				readVec::getMaxLength(seq, maxLen);
				seq.name_ = std::regex_replace(seq.name_, std::regex{popSeqsRegexPatRemoval}, "");
				hPopUID_to_hConsensus[target][seq.name_] = seq.seq_;
				hConsensus_to_hPopUID[target][seq.seq_] = seq.name_;
			}
		}
	}

	{
		std::unordered_map<std::string, std::unordered_set<std::string>> allSamplesInOutput;
		std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> hPopUIDPopSamps;
		std::unordered_map<std::string, std::unordered_map<std::string, std::string>> cNameToPopUID;
		std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double>>> readCountsPerHapPerSample;
		TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));
		sampInfoReader.header_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);
		VecStr row;
		while (sampInfoReader.getNextRow(row)) {
			auto target = row[sampInfoReader.header_.getColPos(targetNameColName)];
			//filter to just the select targets if filtering for that
			if ((!selectTargets.empty() && njh::notIn(target, selectTargets)) ||
			    (!excludeTargets.empty() && njh::in(target, excludeTargets))) {
				continue;
			}

			auto sample = row[sampInfoReader.header_.getColPos(sampleColName)];
			//filter to just the select samples if filtering for that
			if ((!selectSamples.empty() && njh::notIn(sample, selectSamples)) ||
			    (!excludeSamples.empty() && njh::in(sample, excludeSamples))) {
				continue;
			}
			targetNamesSet.emplace(target);
			if(!isDoubleStr(row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)])) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " for row: "<< "\n";
				ss << njh::conToStr(row, "\t") << "\n";
				ss << "Within Sample Read Count column, " << withinSampleReadCntColName << ":" << row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)] << ", does not look like a number" << "\n";
				if(allWhiteSpaceStr(row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)])) {
					ss << "value can't be empty, has to be a number" << "\n";
				}
				throw std::runtime_error{ss.str()};
			}

			auto readCnt = njh::StrToNumConverter::stoToNum<double>(
							row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);
			const auto& h_popUID = row[sampInfoReader.header_.getColPos(popHapIdColName)];
			auto hapName = njh::pasteAsStr(sample, "__", h_popUID);
			std::string hapSeq;
			if (popSeqsDirFnp.empty()) {
				hPopUID_to_hConsensus[target][h_popUID] = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
				hapSeq = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
				hConsensus_to_hPopUID[target][hapSeq] = h_popUID;
			} else {
				hapSeq = njh::mapAt(njh::mapAt(hPopUID_to_hConsensus, target), h_popUID);
			}
			seqInfo clus(hapName, hapSeq);
			clus.cnt_ = readCnt;
			bool add = true;
			for (auto &seq: allResultSeqs[target][sample]) {
				if (seq.seq_ == clus.seq_) {
					add = false;
					seq.cnt_ += readCnt;
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
			inputSampleNamesSet.emplace(sample);
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	if (!exists(popSeqsDirFnp)) {
		std::unordered_map<std::string, std::vector<seqInfo>> popSeqs;
		for (const auto &tarPopHaps: hPopUID_to_hConsensus) {
			for (const auto &popHap: tarPopHaps.second) {
				popSeqs[tarPopHaps.first].emplace_back(popHap.first, popHap.second);
			}
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	//write out seqs
	for(const auto & tar : allResultSeqs) {
		auto tarDir = njh::files::make_path(setUp.pars_.directoryName_, tar.first);
		njh::files::makeDir(njh::files::MkdirPar{tarDir});
		auto outSeqOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, tar.first, "inputSeqs.fasta.gz"));
		SeqOutput writer(outSeqOpts);
		writer.openOut();
		for(const auto & sampSeqs : tar.second) {
			for( auto  seq : sampSeqs.second) {
				MetaDataInName meta;
				if(nullptr != metaGroupData) {
					meta = metaGroupData->getMetaForSample(sampSeqs.first);
				}
				meta.addMeta("sample", sampSeqs.first, true);
				meta.addMeta("readCount", std::max(1.0, std::round(seq.cnt_)), true);
				meta.addMeta("originalIdentifier", hConsensus_to_hPopUID[tar.first][seq.seq_]);
				meta.resetMetaInName(seq.name_);
				writer.write(seq);
			}
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> genomicLocs;
	if(!bedLocs.empty()) {
		auto locs = getBeds(bedLocs);
		for(auto & loc : locs) {
			if(genomicLocsChangePeriodToDash) {
				loc->name_ = njh::replaceString(loc->name_, ".", "-");
			}
			if(njh::in(loc->name_, genomicLocs)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have location: " << loc->name_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			genomicLocs[loc->name_] = loc;
		}
		VecStr missingTargets;
		for(const auto & tar : targetNamesSet) {
			if(njh::notIn(tar, genomicLocs)) {
				missingTargets.emplace_back(tar);
			}
		}
		if(!missingTargets.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "missing the following targets from " << bedLocs << "\n";
			ss << njh::conToStr(missingTargets, "\n") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	runLog["initial set up time"] = fullWatch.totalTime();
	fullWatch.startNewLap("run variant calling on each target");
	auto & runLogTargetTimes = runLog["targets"];
	njh::concurrent::LockableQueue<std::string> targetNamesQueue(targetNamesSet);
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	// sleep(100000);

	{
		std::function<void()> callVariants = [&collapseVarCallPars,&targetNamesQueue, &setUp,
					&genomicLocs, &runLogMut, &runLogTargetTimes]() {
			std::string target;
			while(targetNamesQueue.getVal(target)) {
				njh::stopWatch watch;
				Json::Value currentLog;
				currentLog["target"] = target;
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				auto inputSeqsOpts = SeqIOOptions::genFastaInGz(njh::files::make_path(setUp.pars_.directoryName_, target, "inputSeqs.fasta.gz"));
				auto inputSeqs = SeqInput::getSeqVec<seqInfo>(inputSeqsOpts);
				const auto varCallDirPath = njh::files::make_path(setUp.pars_.directoryName_, target,  "variantCalling");
				auto collapseVarCallParsForTar = collapseVarCallPars;
				collapseVarCallParsForTar.identifier = target;
				collapseVarCallParsForTar.outputDirectory = varCallDirPath;
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				if(njh::in(target, genomicLocs)) {
					collapseVarCallParsForTar.refSeqRegion = GenomicRegion(*njh::mapAt(genomicLocs, target));
				}
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				collapseVarCallParsForTar.calcPopMeasuresPars.seqCountCutOffPloidyCalc_ = 2000;
				collapseAndCallVariants(collapseVarCallParsForTar, inputSeqs);
				currentLog["totalTime"] = watch.totalTime();
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				{
					std::lock_guard<std::mutex> lock(runLogMut);
					runLogTargetTimes.append(currentLog);
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(callVariants, numThreads);
	}


	auto reportsDir = njh::files::make_path(setUp.pars_.directoryName_, "reports");
	njh::files::makeDir(njh::files::MkdirPar{reportsDir});

	//vcfs
	auto reportsVcfsDir = njh::files::make_path(reportsDir, "vcfs");
	njh::files::makeDir(njh::files::MkdirPar{reportsVcfsDir});
	//info
	auto reportsInfoDir = njh::files::make_path(reportsDir, "info");
	njh::files::makeDir(njh::files::MkdirPar{reportsInfoDir});
	//summary
	auto reportsSummaryDir = njh::files::make_path(reportsDir, "summary");
	njh::files::makeDir(njh::files::MkdirPar{reportsSummaryDir});
	//popGenetics
	auto reportsPopGeneticsDir = njh::files::make_path(reportsDir, "popGenetics");
	njh::files::makeDir(njh::files::MkdirPar{reportsPopGeneticsDir});

	auto targetNamesVec = std::vector<std::string>(targetNamesSet.begin(), targetNamesSet.end());
	njh::naturalSortNameSet(targetNamesVec);
	fullWatch.startNewLap("get out the actual samples output were");
	//get out the actual samples output were
	std::set<std::string> sampleNamesSet;
	for (const auto& tar: targetNamesVec) {
		auto outputMetaFnp = njh::files::make_path(setUp.pars_.directoryName_, "/", tar, "/variantCalling/uniqueSeqs_meta.tab.txt.gz");
		TableReader tabReader(TableIOOpts::genTabFileIn(outputMetaFnp));
		VecStr row;
		while(tabReader.getNextRow(row)) {
			sampleNamesSet.emplace(row[tabReader.header_.getColPos("sample")]);
		}
	}

	auto collectAndWriteMasterFile = [&targetNamesVec,&setUp](const std::string & pathUnderTarDirectory, const bfs::path & outputFnp) {
		std::vector<bfs::path> divMeasuresFnps;
		for (const auto& tar: targetNamesVec) {
			auto divFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/", pathUnderTarDirectory);
			if (bfs::exists(divFnp) && 0 != njh::files::bfs::file_size(divFnp)) {
				divMeasuresFnps.emplace_back(divFnp);
			}
		}
		if (!divMeasuresFnps.empty()) {
			njh::files::bfs::path firstFileFnp = divMeasuresFnps.front();
			TableReader firstTable(TableIOOpts(InOptions(firstFileFnp), "\t", true));
			OutputStream out(outputFnp);
			out << njh::conToStr(firstTable.header_.columnNames_, "\t") << '\n'; {
				VecStr firstTableRow;
				while (firstTable.getNextRow(firstTableRow)) {
					out << njh::conToStr(firstTableRow, "\t") << '\n';
				}
			}
			for (const auto& file: divMeasuresFnps) {
				if (file != firstFileFnp) {
					TableReader currentTable(TableIOOpts(InOptions(file), "\t", true));
					VecStr currentRow;
					if (!std::equal(firstTable.header_.columnNames_.begin(), firstTable.header_.columnNames_.end(),
													currentTable.header_.columnNames_.begin(), currentTable.header_.columnNames_.end())) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "header for " << file << " doesn't match other columns" << "\n";
						ss << "expected header: " << njh::conToStr(firstTable.header_.columnNames_) << "\n";
						ss << "found    header: " << njh::conToStr(currentTable.header_.columnNames_) << '\n';
						throw std::runtime_error{ss.str()};
													}
					while (currentTable.getNextRow(currentRow)) {
						out << njh::conToStr(currentRow, "\t") << '\n';
					}
				}
			}
		}
	};

	//concatenate diversity calls
	{
		fullWatch.startNewLap("gather sequence diversity");
		//sequence diversity
		collectAndWriteMasterFile(njh::pasteAsStr("/variantCalling/divMeasures.tab.txt"),
			njh::files::make_path(reportsPopGeneticsDir, "allDiversityMeasures.tsv.gz"));
	}
	{
		fullWatch.startNewLap("gather translated diversity");
		//translated diversity
		collectAndWriteMasterFile(njh::pasteAsStr("/variantCalling/variantCalls/translatedDivMeasures.tab.txt"),
	njh::files::make_path(reportsPopGeneticsDir, "allTranslatedDivMeasures.tsv.gz"));
	}

	if(!collapseVarCallPars.metaFieldsToCalcPopDiffs.empty()){
		fullWatch.startNewLap("gather sub meta fields diversity info");
		for(const auto & metaField : collapseVarCallPars.metaFieldsToCalcPopDiffs) {
			{
				collectAndWriteMasterFile(
					njh::pasteAsStr("/variantCalling/perMetaFields/", metaField, "_divMeasures.tab.txt.gz"),
					njh::files::make_path(reportsPopGeneticsDir, njh::pasteAsStr("all_", metaField, "_divMeasures.tsv.gz" )));

				collectAndWriteMasterFile(
					njh::pasteAsStr("/variantCalling/perMetaFields/", metaField, "_diffMeasures.tab.txt.gz"),
					njh::files::make_path(reportsPopGeneticsDir, njh::pasteAsStr("all_", metaField, "_diffMeasures.tsv.gz")));

				collectAndWriteMasterFile(
					njh::pasteAsStr("/variantCalling/perMetaFields/", metaField, "_pairwiseDiffMeasures.tab.txt.gz"),
					njh::files::make_path(reportsPopGeneticsDir, njh::pasteAsStr("all_", metaField, "_pairwiseDiffMeasures.tsv.gz")));
			}
		}
	}

	if(collapseVarCallPars.exportLabIsolateSeqs) {
		//refSeqs.fasta.gz
		auto reportsRefSeqsDir = njh::files::make_path(reportsDir, "refSeqs");
		njh::files::makeDir(njh::files::MkdirPar{reportsRefSeqsDir});
		for (const auto& tar: targetNamesVec) {
			auto refSeqsFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/", "variantCalling/refSeqs.fasta.gz");
			if (bfs::exists(refSeqsFnp) && 0 != njh::files::bfs::file_size(refSeqsFnp)) {
				auto outputFnp = njh::files::make_path(reportsRefSeqsDir, tar + ".fasta.gz");
				bfs::copy_file(refSeqsFnp, outputFnp);
			}
		}
	}

	//all samples and target coverage counts
	{
		fullWatch.startNewLap("all samples and target coverage counts");
		OutputStream readCountsOut(njh::files::make_path(reportsSummaryDir, "readCountsPerSamplePerTarget.tsv.gz"));

		OutputStream hapCountsOut(njh::files::make_path(reportsSummaryDir, "hapCountsPerSamplePerTarget.tsv.gz"));
		hapCountsOut << "target\t" << njh::conToStr(sampleNamesSet, "\t") << std::endl;
		readCountsOut << "target\t" << njh::conToStr(sampleNamesSet, "\t") << std::endl;

		for(const auto & tar : targetNamesVec) {
			hapCountsOut << tar;
			readCountsOut << tar;
			for(const auto & samp : sampleNamesSet) {
				hapCountsOut << "\t" << allResultSeqs[tar][samp].size();
				uint32_t currentReadCount = 0;
				for(const auto & r : allResultSeqs[tar][samp]) {
					currentReadCount += static_cast<uint32_t>(std::round(r.cnt_));
				}
				readCountsOut << "\t" << currentReadCount;
			}
			hapCountsOut << std::endl;
			readCountsOut << std::endl;
		}
	}

	//copy in meta
	if(!metaFnp.empty() && exists(metaFnp) && metaGroupData) {
		metaGroupData->writeOutMetaFile(njh::files::make_path(reportsInfoDir, "meta.tsv"),	inputSampleNamesSet);
	}
	TranslatorByAlignment::GetGenomicLocationsForAminoAcidPositionsRet locs;
	if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
		fullWatch.startNewLap("add known mutations to bed files");
		//add known to bed files
		TranslatorByAlignment::GetGenomicLocationsForAminoAcidPositionsPars parsForBedFileGen;
		parsForBedFileGen.gffFnp = collapseVarCallPars.transPars.gffFnp_;
		parsForBedFileGen.outOpts = OutOptions(njh::files::make_path(reportsInfoDir, "genomicLocsForKnownAAChanges.bed"));
		auto gprefix = bfs::path(collapseVarCallPars.transPars.lzPars_.genomeFnp).replace_extension("");
		auto twoBitFnp = gprefix.string() + ".2bit";
		parsForBedFileGen.twoBitFnp = twoBitFnp;
		parsForBedFileGen.proteinMutantTypingFnp = collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_;

		locs = TranslatorByAlignment::getGenomicLocationsForAminoAcidPositions(parsForBedFileGen);

		OutputStream transcriptOut(njh::files::make_path(reportsInfoDir, "transcriptLocsForKnownAAChanges.bed"));
		for(const auto & t : locs.transcriptLocs) {
			transcriptOut << t.toDelimStrWithExtra() << std::endl;
		}
	}
	fullWatch.startNewLap("add what genes are intersecting with the");
	//add what genes are intersecting with the
	{
		table overlappingGeneInfo(VecStr{"target", "geneInfo"});
		for(const auto & target : targetNamesVec) {
			auto geneBedFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/geneInfos"),false,
				std::vector<std::regex>{std::regex{".*.bed"}},
				std::vector<std::regex>{std::regex{".*_exonIntronPositions.bed"}});
			for(const auto  & f : geneBedFiles) {
				auto beds = getBeds(f.first);
				std::string geneInfo;
				for(const auto & b : beds) {
					if(!b->extraFields_.empty()) {
						geneInfo = b->extraFields_[0];
					}
				}
				overlappingGeneInfo.addRow(target,geneInfo);
			}
		}
		table::splitColWithMetaPars splitPars;
		splitPars.column_ = "geneInfo";
		splitPars.removeEmptyColumn_ = true;
		auto splitTab = table::splitColWithMeta(overlappingGeneInfo, splitPars);
		OutputStream geneInfoTabout(njh::files::make_path(reportsInfoDir, "targetsIntersectingWithGenesInfo.tsv"));
		splitTab.outPutContents(geneInfoTabout, "\t");
	}
	if(!genomicLocs.empty() && !locs.genomicLocs.empty()) {
		OutputStream intersectionWithKnownLocsOut(njh::files::make_path(reportsInfoDir, "targetsIntersectingWithLocsForKnownAAChanges.bed"));
		std::vector<std::shared_ptr<Bed6RecordCore>> allLocs;
		allLocs.reserve(genomicLocs.size());
		for (const auto& g: genomicLocs) {
			if(njh::in(g.first, targetNamesSet)) {
				allLocs.emplace_back(g.second);
			}
		}
		BedUtility::coordSort(allLocs);
		for (const auto& loc: allLocs) {
			VecStr intersected;
			for (const auto& knownLoc: locs.genomicLocs) {
				if (knownLoc.overlaps(*loc, 1)) {
					intersected.emplace_back(knownLoc.name_);
				}
			}
			if (!intersected.empty()) {
				intersectionWithKnownLocsOut << loc->toDelimStrWithExtra() << "\t" << njh::conToStr(intersected, ",") << std::endl;
			}
		}
	}
	fullWatch.startNewLap("gather unmapped read counts and gather translation filter counts");
	//gather unmapped read counts and gather translation filter counts
	std::unordered_map<std::string, uint32_t> unmmappedHapCounts;
	std::unordered_map<std::string, uint32_t> untranslatableHapCounts;
	for(const auto & target : targetNamesVec) {
		{
			auto unmappableFnp = njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/seqsUnableToBeMapped.txt");
			auto allLines = njh::files::getAllLines(unmappableFnp);
			uint32_t count = 0;
			for(const auto & l : allLines) {
				if(!l.empty()) {
					++count;
				}
			}
			unmmappedHapCounts[target] = count;
		}
		{
			auto untranslatable = njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/seqsTranslationFiltered.txt");
			auto allLines = njh::files::getAllLines(untranslatable);
			uint32_t count = 0;
			for(const auto & l : allLines) {
				if(!l.empty()) {
					++count;
				}
			}
			untranslatableHapCounts[target] = count;
		}
	}

	{
		OutputStream unmmappedHapCountsOut(njh::files::make_path(reportsSummaryDir, "unmappedUntranslatableHapCounts.tsv.gz"));
		unmmappedHapCountsOut << "target\tunmappedHaps\tmappedButUntranslatable" << std::endl;
		for (const auto& target: targetNamesVec) {
			unmmappedHapCountsOut << target
					<< "\t" << unmmappedHapCounts[target]
					<< "\t" << untranslatableHapCounts[target]
					<< std::endl;
		}
	}

	fullWatch.startNewLap("gather vcfs");
	//gather vcfs
	std::vector<bfs::path> proteinVcfs;
	std::vector<bfs::path> genomicVcfs;
	std::vector<bfs::path> complexGenomicVcfs;
	for (const auto& target: targetNamesVec) {
		auto proteinVcfFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"),false,
				std::vector<std::regex>{std::regex{".*-protein.vcf.gz"}});
		auto genomicVcfFiles = njh::files::listAllFiles(
			njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"), false,
			std::vector<std::regex>{std::regex{".*-genomic.vcf.gz"}},
			std::vector<std::regex>{
				std::regex{".*-complex-genomic.vcf.gz"}
			});
		auto complexGenomicVcfFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"),false,
		                                                       std::vector<std::regex>{
			                                                       std::regex{".*-complex-genomic.vcf.gz"}
		                                                       });
		for (const auto& pvcf: proteinVcfFiles) {
			proteinVcfs.emplace_back(pvcf.first);
		}
		for (const auto& gvcf: genomicVcfFiles) {
			genomicVcfs.emplace_back(gvcf.first);
		}
		for (const auto& gvcf: complexGenomicVcfFiles) {
			complexGenomicVcfs.emplace_back(gvcf.first);
		}
	}

	//combine vcfs with handling of overlapping variant calls
	//processing protein vcfs;
	std::vector<GenomicRegion> knownAAVariantRegions;
	knownAAVariantRegions.reserve(locs.transcriptLocs.size());
	for (const auto& b: locs.transcriptLocs) {
		knownAAVariantRegions.emplace_back(b);
	}
	if(!proteinVcfs.empty()){
		fullWatch.startNewLap("combine protein vcfs");
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		auto firstPVcf = VCFOutput::comnbineVCFs(proteinVcfs, sampleNamesSet, combiningVcfPars);
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		{
			OutputStream pvcf(njh::files::make_path(reportsVcfsDir, "allProteinVariantCalls.vcf.gz"));
			firstPVcf.writeOutFixedAndSampleMeta(pvcf);
		}
		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream pvcfOutFile(njh::files::make_path(reportsVcfsDir, "knownAAChangesProteinVariantCalls.vcf.gz"));

			firstPVcf.writeOutFixedAndSampleMeta(pvcfOutFile, knownAAVariantRegions);
		}
		{
			// aminoAcidInfo::infos::allInfo
			OutputStream aminoAcidChangesTable(njh::files::make_path(reportsSummaryDir, "AAChangesInfo.tsv.gz"));
			aminoAcidChangesTable << "Gene_ID"
					<< "\t" << "Gene_Transcript_ID"
					<< "\t" << "Gene"
					<< "\t" << "Mutation_Name"
					<< "\t" << "ExonicFunc"
					<< "\t" << "AA_Change"
					<< "\t" << "Targeted"
			    << "\t" << "CoveredBy"
					<< "\t" << "sample";
					aminoAcidChangesTable << "\t" << "AA_Position";
					aminoAcidChangesTable << "\t" << "reference_AA_cnt"
			    << "\t" << "reference_AA_freq"
					<< "\t" << "alternate_AA_cnt"
					<< "\t" << "alternate_AA_freq"
					<< "\t" << "coverage_AA_cnt"
					<< "\t" << "AlleleCount"
					<< "\t" << "AlleleFrequency"
					<< "\t" << "SampleCount"
					<< "\t" << "SamplePrevalence";
			aminoAcidChangesTable << std::endl;
			for(const auto & rec : firstPVcf.records_) {
				std::string targeted = "No";

				auto genomicRegion = rec.genRegion();
				for(const auto & knownRegion : knownAAVariantRegions) {
					if(genomicRegion.overlaps(knownRegion)) {
						targeted = "Yes";
						break;
					}
				}
				for(const auto & sample : rec.sampleFormatInfos_) {

					auto TYPE = tokenizeString(rec.info_.getMeta("TYPE"), ",");

					auto DP = sample.second.getMeta("DP");
					auto sample_ADs = tokenizeString(sample.second.getMeta("AD"), ",");
					auto sample_AFs = tokenizeString(sample.second.getMeta("AF"), ",");
					auto ACs = tokenizeString(rec.info_.getMeta("AC"), ",");
					auto AFs = tokenizeString(rec.info_.getMeta("AF"), ",");
					auto SCs = tokenizeString(rec.info_.getMeta("SC"), ",");
					auto PREVs = tokenizeString(rec.info_.getMeta("PREV"), ",");
					for(const auto & altEnum : iter::enumerate(rec.alts_)) {
						auto ref = rec.ref_;
						std::string refTriCodeName;
						for(const auto c : ref) {
							auto currentTriCode = aminoAcidInfo::infos::allInfo.at(c).triCode_;
							currentTriCode[0] = static_cast<char>(toupper(currentTriCode[0]));
							refTriCodeName+= currentTriCode;
						}
						auto alt = altEnum.element;
						std::string altTriCodeName;
						for(const auto c : alt) {
							if(c != 'X' && c != 'x') {
								auto currentTriCode = aminoAcidInfo::infos::allInfo.at(c).triCode_;
								currentTriCode[0] = static_cast<char>(toupper(currentTriCode[0]));
								altTriCodeName+= currentTriCode;
							} else {
								std::string currentTriCode = "XXX";
								altTriCodeName+= currentTriCode;
							}
						}

						std::string ExonicFunc;
						if (TYPE[altEnum.index] == "snp") {
							ExonicFunc = "missense_variant";
						} else if (TYPE[altEnum.index] == "del") {
							ExonicFunc = "conservative_inframe_deletion";
						} else if (TYPE[altEnum.index] == "ins") {
							ExonicFunc = "conservative_inframe_insertion";
						}


						std::string geneName = rec.info_.getMeta("GeneName");

						aminoAcidChangesTable << rec.info_.getMeta("GeneID")
							<< "\t" << rec.chrom_
							<< "\t" << geneName
							<< "\t" << njh::pasteAsStr(geneName, "-", refTriCodeName, rec.pos_, altTriCodeName)
							<< "\t" << ExonicFunc
							<< "\t" << njh::pasteAsStr(refTriCodeName, rec.pos_, altTriCodeName)
							<< "\t" << targeted
						  << "\t" << rec.info_.getMeta("TARGET")
							<< "\t" << sample.first;

						aminoAcidChangesTable << "\t" << rec.pos_;

						if("." == DP) {
							//no coverage
							aminoAcidChangesTable << "\t" << "0"
									<< "\t" << "0"
									<< "\t" << "0"
									<< "\t" << "0"
									<< "\t" << "0";
						} else {
							aminoAcidChangesTable << "\t" << sample_ADs[0]
							<< "\t" << sample_AFs[0]
							<< "\t" << sample_ADs[altEnum.index + 1]
							<< "\t" << sample_AFs[altEnum.index + 1]
							<< "\t" << DP;
						}
						aminoAcidChangesTable << "\t" << ACs[altEnum.index]
							<< "\t" << AFs[altEnum.index]
							<< "\t" << SCs[altEnum.index]
							<< "\t" << PREVs[altEnum.index];
						aminoAcidChangesTable << std::endl;
					}
				}
			}
		}
	}

	//process genomic
	if(!genomicVcfs.empty()) {
		fullWatch.startNewLap("combine genomic vcfs");
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		auto firstGVcf = VCFOutput::comnbineVCFs(genomicVcfs, sampleNamesSet, combiningVcfPars);
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		{
			OutputStream gvcfOutFile(njh::files::make_path(reportsVcfsDir, "allGenomicVariantCalls.vcf.gz"));
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile);
		}
		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream gvcfOutFile(njh::files::make_path(reportsVcfsDir, "knownAAChangesGenomicVariantCalls.vcf.gz"));
			std::vector<GenomicRegion> knownSnpVariantRegions;
			knownSnpVariantRegions.reserve(locs.genomicLocs.size());
			for (const auto& b: locs.genomicLocs) {
				knownSnpVariantRegions.emplace_back(b);
			}
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile, knownSnpVariantRegions);
		}
	}

	//process complex genomic
	if(!complexGenomicVcfs.empty()) {
		fullWatch.startNewLap("combine genomic vcfs");
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		{
			auto firstGVcf = VCFOutput::comnbineVCFs(complexGenomicVcfs, sampleNamesSet, combiningVcfPars);
			OutputStream gvcfOutFile(njh::files::make_path(reportsVcfsDir, "allComplexGenomicVariantCalls.vcf.gz"));
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile);
		}
	}

	//combine summary tables
	{
		fullWatch.startNewLap("gather summary tables");
		//translated diversity
		std::vector<bfs::path> summaryFnps;
		for (const auto& tar: targetNamesVec) {
			auto summaryFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/variantCalling/summaryTable.tab.txt.gz");
			if (bfs::exists(summaryFnp) && 0 != njh::files::bfs::file_size(summaryFnp)) {
				summaryFnps.emplace_back(summaryFnp);
			}
		}
		if (!summaryFnps.empty()) {
			njh::files::bfs::path firstFileFnp = summaryFnps.front();
			TableReader firstTable(TableIOOpts(InOptions(firstFileFnp), "\t", true));
			OutputStream out(njh::files::make_path(reportsSummaryDir, "allSummaryTables.tab.txt.gz"));
			out << njh::conToStr(firstTable.header_.columnNames_, "\t") << '\n'; {
				VecStr firstTableRow;
				while (firstTable.getNextRow(firstTableRow)) {
					out << njh::conToStr(firstTableRow, "\t") << '\n';
				}
			}
			for (const auto& file: summaryFnps) {
				if (file != firstFileFnp) {
					TableReader currentTable(TableIOOpts(InOptions(file), "\t", true));
					VecStr currentRow;
					if (!std::equal(firstTable.header_.columnNames_.begin(), firstTable.header_.columnNames_.end(),
													currentTable.header_.columnNames_.begin(), currentTable.header_.columnNames_.end())) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "header for " << file << " doesn't match other columns" << "\n";
						ss << "expected header: " << njh::conToStr(firstTable.header_.columnNames_) << "\n";
						ss << "found    header: " << njh::conToStr(currentTable.header_.columnNames_) << '\n';
						throw std::runtime_error{ss.str()};
													}
					while (currentTable.getNextRow(currentRow)) {
						out << njh::conToStr(currentRow, "\t") << '\n';
					}
				}
			}
		}
	}
	//create counts of summary table
	{
		table allSummaryTable(TableIOOpts::genTabFileIn(njh::files::make_path(reportsSummaryDir, "allSummaryTables.tab.txt.gz")));
		{
			auto genomicLocCounts = allSummaryTable.countGroupColumns(VecStr{ "chrom", "0based_start", "0based_end","target", "length", "strand"});
			VecStr isMax;
			std::unordered_map<std::string, uint32_t> maxCounts;
			for(const auto & row : genomicLocCounts) {
				auto n = njh::StrToNumConverter::stoToNum<uint32_t>(row[genomicLocCounts.getColPos("n")]);
				if(n > maxCounts[row[genomicLocCounts.getColPos("target")]]) {
					maxCounts[row[genomicLocCounts.getColPos("target")]] = n;
				}
			}
			for(const auto & row : genomicLocCounts) {
				auto n = njh::StrToNumConverter::stoToNum<uint32_t>(row[genomicLocCounts.getColPos("n")]);
				if(n == maxCounts[row[genomicLocCounts.getColPos("target")]]) {
					isMax.emplace_back("true");
				} else {
					isMax.emplace_back("false");
				}
			}
			genomicLocCounts.addColumn(isMax, "isMaxCount");
			genomicLocCounts.naturlSortTable("target", false);
			genomicLocCounts.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(reportsInfoDir, "genomicLocPerTargetsCounts.tab.txt.gz")));

		}
		{
			auto proteinLocCounts = allSummaryTable.countGroupColumns(VecStr{"transcript", "transcript_1based_start", "transcript_1based_end", "target",  "transcript_length"});
			proteinLocCounts.addColumn(VecStr{"+"}, "transcript_strand");
			proteinLocCounts = proteinLocCounts.getColumns(VecStr{"transcript", "transcript_1based_start", "transcript_1based_end", "target",  "transcript_length", "transcript_strand", "n"});
			VecStr isMax;
			std::unordered_map<std::string, uint32_t> maxCounts;
			for(const auto & row : proteinLocCounts) {
				auto n = njh::StrToNumConverter::stoToNum<uint32_t>(row[proteinLocCounts.getColPos("n")]);
				if(n > maxCounts[row[proteinLocCounts.getColPos("target")]]) {
					maxCounts[row[proteinLocCounts.getColPos("target")]] = n;
				}
			}
			for(const auto & row : proteinLocCounts) {
				auto n = njh::StrToNumConverter::stoToNum<uint32_t>(row[proteinLocCounts.getColPos("n")]);
				if(n == maxCounts[row[proteinLocCounts.getColPos("target")]]) {
					isMax.emplace_back("true");
				} else {
					isMax.emplace_back("false");
				}
			}
			proteinLocCounts.addColumn(isMax, "isMaxCount");
			proteinLocCounts.naturlSortTable("target", false);
			proteinLocCounts.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(reportsInfoDir, "proteinLocPerTargetsCounts.tab.txt.gz")));
		}
	}

	//run log
	{
		fullWatch.startNewLap("end");
		runLog["run_times"] = fullWatch.toJson();
		OutputStream logOut(njh::files::make_path(reportsDir, "log.txt"));
		logOut << runLog << std::endl;
	}

	return 0;

}

}  // namespace njhseq