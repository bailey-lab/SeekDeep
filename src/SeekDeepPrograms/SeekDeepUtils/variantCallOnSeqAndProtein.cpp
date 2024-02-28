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
	bool doNotRescueVariantCallsAccrossTargets = false;
	std::string popSeqsRegexPatRemoval = R"(_([tf])?\d+(\.\d+)?$)";
	uint32_t numThreads = 1;

	CollapseAndCallVariantsPars collapseVarCallPars;
	std::set<std::string> selectTargets;
	std::set<std::string> selectSamples;

	bfs::path bedLocs;
	bool genomicLocsChangePeriodToDash = false;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(doNotRescueVariantCallsAccrossTargets, "--doNotRescueVariantCallsAccrossTargets", "do Not Rescue Variant Calls Accross Targets");

	setUp.setOption(selectTargets, "--selectTargets", "Only analzye these select targets");
	setUp.setOption(selectSamples, "--selectSamples", "Only analzye these select samples");
	setUp.setOption(bedLocs, "--genomicLocations", "a bed file with specific genomic locations to align to, location name needs to match target name");
	setUp.setOption(genomicLocsChangePeriodToDash, "--genomicLocationsChangePeriodToDash", "when supplying a location name, change periods to dashes in the name");


	setUp.setOption(collapseVarCallPars.variantCallerRunPars.occurrenceCutOff, "--variantOccurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.lowVariantCutOff, "--variantFrequencyCutOff", "Low Variant Cut Off, don't report variants below this frequency");
	collapseVarCallPars.calcPopMeasuresPars.lowVarFreq = collapseVarCallPars.variantCallerRunPars.lowVariantCutOff;
	collapseVarCallPars.transPars.setOptions(setUp, true);
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(collapseVarCallPars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	collapseVarCallPars.calcPopMeasuresPars.diagAlnPairwiseComps = !collapseVarCallPars.noDiagAlnPairwiseComps;
	//setOption(collapseVarCallPars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.numThreads, "--mappingNumThreads", "Number of threads to use for the alignment portion");
	setUp.setOption(collapseVarCallPars.metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "meta Fields To Calc Pop Diffs");


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

	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double>>> readCountsPerHapPerSample;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> cNameToPopUID;
	std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> hPopUIDPopSamps;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hPopUID_to_hConsensus;

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
	std::unordered_map<std::string, std::vector<seqInfo>> popSeqs;
	auto sampInfoFnp = resultsFnp;

	std::unordered_map<std::string, std::unordered_set<std::string>> allSamplesInOutput;
	//key1 == target, key2 == sample
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> allResultSeqs;
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
			if(!selectTargets.empty() && !njh::in(target, selectTargets)) {
				continue;
			}

			auto sample = row[sampInfoReader.header_.getColPos(sampleColName)];
			//filter to just the select samples if filtering for that
			if(!selectSamples.empty() && njh::notIn(sample, selectSamples)) {
				continue;;
			}
			targetNamesSet.emplace(target);

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

	if (!exists(popSeqsDirFnp)) {
		for (const auto &tarPopHaps: hPopUID_to_hConsensus) {
			for (const auto &popHap: tarPopHaps.second) {
				popSeqs[tarPopHaps.first].emplace_back(popHap.first, popHap.second);
			}
		}
	}

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
				meta.resetMetaInName(seq.name_);
				writer.write(seq);
			}
		}
	}

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
	runLog["initial set up time"] = fullWatch.totalTime();
	fullWatch.startNewLap("run variant calling on each target");
	auto & runLogTargetTimes = runLog["targets"];
	njh::concurrent::LockableQueue<std::string> targetNamesQueue(targetNamesSet);

	std::function<void()> callVariants = [&collapseVarCallPars,&targetNamesQueue, &setUp,
		&genomicLocs, &runLogMut, &runLogTargetTimes]() {
		std::string target;
		while(targetNamesQueue.getVal(target)) {
			njh::stopWatch watch;
			Json::Value currentLog;
			currentLog["target"] = target;
			auto inputSeqsOpts = SeqIOOptions::genFastaInGz(njh::files::make_path(setUp.pars_.directoryName_, target, "inputSeqs.fasta.gz"));
			auto inputSeqs = SeqInput::getSeqVec<seqInfo>(inputSeqsOpts);
			const auto varCallDirPath = njh::files::make_path(setUp.pars_.directoryName_, target,  "variantCalling");
			auto collapseVarCallParsForTar = collapseVarCallPars;
			collapseVarCallParsForTar.identifier = target;
			collapseVarCallParsForTar.outputDirectory = varCallDirPath;
			if(njh::in(target, genomicLocs)) {
				collapseVarCallParsForTar.refSeqRegion = GenomicRegion(*njh::mapAt(genomicLocs, target));
			}
			collapseVarCallParsForTar.calcPopMeasuresPars.seqCountCutOffPloidyCalc_ = 2000;
			collapseAndCallVariants(collapseVarCallParsForTar, inputSeqs);
			currentLog["totalTime"] = watch.totalTime();
			{
				std::lock_guard<std::mutex> lock(runLogMut);
				runLogTargetTimes.append(currentLog);
			}
		}
	};

	njh::concurrent::runVoidFunctionThreaded(callVariants, numThreads);

	auto reportsDir = njh::files::make_path(setUp.pars_.directoryName_, "reports");
	njh::files::makeDir(njh::files::MkdirPar{reportsDir});




	//concatenate diversity calls

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

	{
		fullWatch.startNewLap("gather sequence diversity");
		//sequence diversity
		std::vector<bfs::path> seqDivMeasuresFnps;
		for (const auto& tar: targetNamesVec) {
			auto divFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/variantCalling/divMeasures.tab.txt");
			if (bfs::exists(divFnp) && 0 != njh::files::bfs::file_size(divFnp)) {
				seqDivMeasuresFnps.emplace_back(divFnp);
			}
		}
		if (!seqDivMeasuresFnps.empty()) {
			njh::files::bfs::path firstFileFnp = seqDivMeasuresFnps.front();
			TableReader firstTable(TableIOOpts(InOptions(firstFileFnp), "\t", true));
			OutputStream out(njh::files::make_path(reportsDir, "allDiversityMeasures.tsv.gz"));
			out << njh::conToStr(firstTable.header_.columnNames_, "\t") << '\n'; {
				VecStr firstTableRow;
				while (firstTable.getNextRow(firstTableRow)) {
					out << njh::conToStr(firstTableRow, "\t") << '\n';
				}
			}
			for (const auto& file: seqDivMeasuresFnps) {
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
	{
		fullWatch.startNewLap("gather translated diversity");
		//translated diversity
		std::vector<bfs::path> transDivMeasuresFnps;
		for (const auto& tar: targetNamesVec) {
			auto divFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/variantCalling/variantCalls/translatedDivMeasures.tab.txt");
			if (bfs::exists(divFnp) && 0 != njh::files::bfs::file_size(divFnp)) {
				transDivMeasuresFnps.emplace_back(divFnp);
			}
		}
		if (!transDivMeasuresFnps.empty()) {
			njh::files::bfs::path firstFileFnp = transDivMeasuresFnps.front();
			TableReader firstTable(TableIOOpts(InOptions(firstFileFnp), "\t", true));
			OutputStream out(njh::files::make_path(reportsDir, "allTranslatedDivMeasures.tsv.gz"));
			out << njh::conToStr(firstTable.header_.columnNames_, "\t") << '\n'; {
				VecStr firstTableRow;
				while (firstTable.getNextRow(firstTableRow)) {
					out << njh::conToStr(firstTableRow, "\t") << '\n';
				}
			}
			for (const auto& file: transDivMeasuresFnps) {
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

	//all samples and target coverage counts
	{
		fullWatch.startNewLap("all samples and target coverage counts");
		OutputStream readCountsOut(njh::files::make_path(reportsDir, "readCountsPerSamplePerTarget.tsv.gz"));

		OutputStream hapCountsOut(njh::files::make_path(reportsDir, "hapCountsPerSamplePerTarget.tsv.gz"));
		hapCountsOut << "target\t" << njh::conToStr(sampleNamesSet, "\t") << std::endl;
		readCountsOut << "target\t" << njh::conToStr(sampleNamesSet, "\t") << std::endl;

		for(const auto & tar : targetNamesVec) {
			hapCountsOut << tar;
			readCountsOut << tar;
			for(const auto & samp : sampleNamesSet) {
				hapCountsOut << "\t" << allResultSeqs[tar][samp].size();
				uint32_t currentReadCount = 0;
				for(const auto & r : allResultSeqs[tar][samp]) {
					currentReadCount += r.cnt_;
				}
				readCountsOut << "\t" << currentReadCount;
			}
			hapCountsOut << std::endl;
			readCountsOut << std::endl;
		}
	}

	//copy in meta
	if(!metaFnp.empty() && exists(metaFnp) && metaGroupData) {

		// bfs::copy_file(metaFnp, njh::files::make_path(reportsDir, "meta.tsv"));
		metaGroupData->writeOutMetaFile(njh::files::make_path(reportsDir, "meta.tsv"),	inputSampleNamesSet);
	}
	TranslatorByAlignment::GetGenomicLocationsForAminoAcidPositionsRet locs;
	if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
		fullWatch.startNewLap("add known mutations to bed files");
		//add known to bed files
		TranslatorByAlignment::GetGenomicLocationsForAminoAcidPositionsPars parsForBedFileGen;
		parsForBedFileGen.gffFnp = collapseVarCallPars.transPars.gffFnp_;
		parsForBedFileGen.outOpts = OutOptions(njh::files::make_path(reportsDir, "genomicLocsForKnownAAChanges.bed"));
		auto gprefix = bfs::path(collapseVarCallPars.transPars.lzPars_.genomeFnp).replace_extension("");
		auto twoBitFnp = gprefix.string() + ".2bit";
		parsForBedFileGen.twoBitFnp = twoBitFnp;
		parsForBedFileGen.proteinMutantTypingFnp = collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_;

		 locs = TranslatorByAlignment::getGenomicLocationsForAminoAcidPositions(parsForBedFileGen);

		OutputStream transcriptOut(njh::files::make_path(reportsDir, "transcriptLocsForKnownAAChanges.bed"));
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
		OutputStream geneInfoTabout(njh::files::make_path(reportsDir, "targetsIntersectingWithGenesInfo.tsv"));
		splitTab.outPutContents(geneInfoTabout, "\t");
	}
	if(!genomicLocs.empty() && !locs.genomicLocs.empty()) {
		OutputStream intersectionWithKnownLocsOut(njh::files::make_path(reportsDir, "targetsIntersectingWithLocsForKnownAAChanges.bed"));
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
			intersectionWithKnownLocsOut << loc->toDelimStrWithExtra() << "\t" << njh::conToStr(intersected) << std::endl;
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
			auto untranslateable = njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/seqsTranslationFiltered.txt");
			auto allLines = njh::files::getAllLines(untranslateable);
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
		OutputStream unmmappedHapCountsOut(njh::files::make_path(reportsDir, "unmmappedUntranslatableHapCounts.tsv.gz"));
		unmmappedHapCountsOut << "target\tunmmapedHaps\tmappedButUntranslatable" << std::endl;
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
	for (const auto& target: targetNamesVec) {
		auto proteinVcfFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"),false,
				std::vector<std::regex>{std::regex{".*-protein.vcf.gz"}});
		auto genomicVcfFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"),false,
			std::vector<std::regex>{std::regex{".*-genomic.vcf.gz"}});
		for(const auto & pvcf : proteinVcfFiles) {
			proteinVcfs.emplace_back(pvcf.first);
		}
		for(const auto & gvcf : genomicVcfFiles) {
			genomicVcfs.emplace_back(gvcf.first);
		}
	}

	//combine vcfs with handling of overlapping variant calls
	//processing protein vcfs;

	if(!proteinVcfs.empty()){
		fullWatch.startNewLap("combine protein vcfs");
		auto firstPVcf = VCFOutput::comnbineVCFs(proteinVcfs, sampleNamesSet, doNotRescueVariantCallsAccrossTargets);
		{
			OutputStream pvcf(njh::files::make_path(reportsDir, "allProteinVariantCalls.vcf.gz"));
			firstPVcf.writeOutFixedAndSampleMeta(pvcf);
		}

		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream pvcfOutFile(njh::files::make_path(reportsDir, "knownAAChangesProteinVariantCalls.vcf.gz"));
			std::vector<GenomicRegion> knownAAVariantRegions;
			knownAAVariantRegions.reserve(locs.transcriptLocs.size());
			for (const auto& b: locs.transcriptLocs) {
				knownAAVariantRegions.emplace_back(b);
			}
			firstPVcf.writeOutFixedAndSampleMeta(pvcfOutFile, knownAAVariantRegions);
		}
	}

	//process genomic
	if(!genomicVcfs.empty()) {
		fullWatch.startNewLap("combine genomic vcfs");
		auto firstGVcf = VCFOutput::comnbineVCFs(genomicVcfs, sampleNamesSet, doNotRescueVariantCallsAccrossTargets);
		{
			OutputStream gvcfOutFile(njh::files::make_path(reportsDir, "allGenomicVariantCalls.vcf.gz"));
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile);
		}
		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream gvcfOutFile(njh::files::make_path(reportsDir, "knownAAChangesGenomicVariantCalls.vcf.gz"));
			std::vector<GenomicRegion> knownSnpVariantRegions;
			knownSnpVariantRegions.reserve(locs.genomicLocs.size());
			for (const auto& b: locs.genomicLocs) {
				knownSnpVariantRegions.emplace_back(b);
			}
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile, knownSnpVariantRegions);
		}
	}

	//run log
	{
		fullWatch.startNewLap("end");
		runLog["run_times"] = fullWatch.toJson();
		OutputStream logOut(njh::files::make_path(reportsDir, "log.txt"));
		logOut << runLog << std::endl;
	}

	//add flag to count singler changes per certain meta and the haplotype known aa typed
	return 0;

}

}  // namespace njhseq