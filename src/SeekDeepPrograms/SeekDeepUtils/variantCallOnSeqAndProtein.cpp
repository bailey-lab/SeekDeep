//
// Created by Nicholas Hathaway on 12/16/23.
//
#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {


int SeekDeepUtilsRunner::variantCallOnSeqAndProteinPost(
				const njh::progutils::CmdArgs &inputCommands) {


	bfs::path resultDirs;

	bfs::path metaFnp = "";
	bool doNotRescueVariantCallsAccrossTargets = false;
	uint32_t numThreads = 1;

	CollapseAndCallVariantsPars collapseVarCallPars;
	std::set<std::string> targetNamesSet;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(doNotRescueVariantCallsAccrossTargets, "--doNotRescueVariantCallsAccrossTargets", "do Not Rescue Variant Calls Accross Targets");

	setUp.setOption(targetNamesSet, "--selectTargets", "Only analzye these select targets", true);

	setUp.setOption(collapseVarCallPars.variantCallerRunPars.occurrenceCutOff, "--variantOccurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.lowVariantCutOff, "--variantFrequencyCutOff", "Low Variant Cut Off, don't report variants below this frequency");
	collapseVarCallPars.calcPopMeasuresPars.lowVarFreq = collapseVarCallPars.variantCallerRunPars.lowVariantCutOff;
	collapseVarCallPars.transPars.setOptions(setUp, true);
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(collapseVarCallPars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	collapseVarCallPars.calcPopMeasuresPars.diagAlnPairwiseComps = !collapseVarCallPars.noDiagAlnPairwiseComps;
	//setOption(collapseVarCallPars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");
	setUp.setOption(collapseVarCallPars.metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "meta Fields To Calc Pop Diffs");


	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");



	setUp.setOption(resultDirs, "--resultDirs", "resultDirs", true);

	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");

	setUp.processDirectoryOutputName("variantCalling_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow({resultDirs},
																	__PRETTY_FUNCTION__);



	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}


	//concatenate diversity calls
	auto reportsDir = njh::files::make_path(setUp.pars_.directoryName_, "reports");
	njh::files::makeDir(njh::files::MkdirPar{reportsDir});
	auto targetNamesVec = std::vector<std::string>(targetNamesSet.begin(), targetNamesSet.end());
	njh::naturalSortNameSet(targetNamesVec);

	//get out the actual samples output were
	std::set<std::string> sampleNamesSet;
	for (const auto& tar: targetNamesVec) {
		auto outputMetaFnp = njh::files::make_path(resultDirs.string(), "/", tar, "/variantCalling/uniqueSeqs_meta.tab.txt.gz");
		TableReader tabReader(TableIOOpts::genTabFileIn(outputMetaFnp));
		VecStr row;
		while(tabReader.getNextRow(row)) {
			sampleNamesSet.emplace(row[tabReader.header_.getColPos("sample")]);
		}
	}

	{
		//sequence diversity
		std::vector<bfs::path> seqDivMeasuresFnps;
		for (const auto& tar: targetNamesVec) {
			auto divFnp = njh::pasteAsStr(resultDirs.string(), "/", tar, "/variantCalling/divMeasures.tab.txt");
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
		//translated diversity
		std::vector<bfs::path> transDivMeasuresFnps;
		for (const auto& tar: targetNamesVec) {
			auto divFnp = njh::pasteAsStr(resultDirs.string(), "/", tar, "/variantCalling/variantCalls/translatedDivMeasures.tab.txt");
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



	//copy in meta
	if(!metaFnp.empty() && exists(metaFnp)) {
		bfs::copy_file(metaFnp, njh::files::make_path(reportsDir, "meta.tsv"));
	}
	TranslatorByAlignment::GetGenomicLocationsForAminoAcidPositionsRet locs;
	if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
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

	//add what genes are intersecting with the
	{
		table overlappingGeneInfo(VecStr{"target", "geneInfo"});
		for(const auto & target : targetNamesVec) {
			auto geneBedFiles = njh::files::listAllFiles(njh::files::make_path(resultDirs.string(), "/", target, "/variantCalling/variantCalls/geneInfos"),false,
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
	//gather unmapped read counts and gather translation filter counts
	std::unordered_map<std::string, uint32_t> unmmappedHapCounts;
	std::unordered_map<std::string, uint32_t> untranslatableHapCounts;
	for(const auto & target : targetNamesVec) {
		{
			auto unmappableFnp = njh::files::make_path(resultDirs.string(), "/", target, "/variantCalling/variantCalls/seqsUnableToBeMapped.txt");
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
			auto untranslateable = njh::files::make_path(resultDirs.string(), "/", target, "/variantCalling/variantCalls/seqsTranslationFiltered.txt");
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


	//gather vcfs
	std::unordered_map<std::string, std::vector<bfs::path>> proteinVcfs;
	std::unordered_map<std::string, std::vector<bfs::path>> genomicVcfs;

	for (const auto& target: targetNamesVec) {
		auto proteinVcfFiles = njh::files::listAllFiles(njh::files::make_path(resultDirs.string(), "/", target, "/variantCalling/variantCalls/"),false,
				std::vector<std::regex>{std::regex{".*-protein.vcf"}});
		auto genomicVcfFiles = njh::files::listAllFiles(njh::files::make_path(resultDirs.string(), "/", target, "/variantCalling/variantCalls/"),false,
			std::vector<std::regex>{std::regex{".*-genomic.vcf"}});
		for(const auto & pvcf : proteinVcfFiles) {
			proteinVcfs[target].emplace_back(pvcf.first);
		}
		for(const auto & gvcf : genomicVcfFiles) {
			genomicVcfs[target].emplace_back(gvcf.first);
		}
	}

	//processing protein vcfs;
	{
		auto firstPVcfFnp = proteinVcfs.begin()->second.front();
		auto firstPVcf = VCFOutput::readInHeader(firstPVcfFnp);
		{
			InputStream in(firstPVcfFnp);
			firstPVcf.addInRecordsFromFile(in);
		}
		//add in blanks for any missing samples for any of the records
		firstPVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		if(firstPVcf.records_.empty()) {
			firstPVcf.samples_ = VecStr{sampleNamesSet.begin(), sampleNamesSet.end()};
		}
		for(const auto & target : proteinVcfs) {
			for(const auto & pvcfFnp : target.second) {
				if(firstPVcfFnp == pvcfFnp) {
					continue;
				}
				auto currentPvcf = VCFOutput::readInHeader(pvcfFnp);

				{
					InputStream in(pvcfFnp);
					currentPvcf.addInRecordsFromFile(in);
					// firstPVcf.addInRecordsFromFile(in);
				}
				//add in header
				//given we know what generated these vcf files there's not need to check the FORMAT, FILTER, or INFO tags

				//add in new contigs if any
				for(const auto & contig : currentPvcf.contigEntries_) {
					if(!njh::in(contig.first, firstPVcf.contigEntries_)) {
						firstPVcf.contigEntries_.emplace(contig);
					} else if(!(firstPVcf.contigEntries_.at(contig.first) == contig.second)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding contig: " << contig.first << " but already have an entry for this but it doesn't match" << "\n";
						ss << "contig in master header: " << njh::json::writeAsOneLine(njh::json::toJson(firstPVcf.contigEntries_.at(contig.first))) << "\n";
						ss << "adding contig: " << njh::json::writeAsOneLine(njh::json::toJson(contig.second)) << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
				currentPvcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
				firstPVcf.records_.reserve(firstPVcf.records_.size() + currentPvcf.records_.size());
				for(auto & record : currentPvcf.records_) {
					firstPVcf.records_.emplace_back(std::move(record));
				}
				//std::move(currentPvcf.records_.begin(), currentPvcf.records_.back(). std::back_inserter(firstPVcf.records_));
			}
		}
		firstPVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		firstPVcf.sortRecords();
		//eliminate duplicate variant calls and optional merge variants calls across amplicons that overlap that may have been missing in one target for a sample but present others

		{

			if(firstPVcf.records_.size() > 1) {
				std::set<uint32_t> positionsToErase;
				std::set<uint32_t> currentPositionsToComp;
				auto compSamePositions = [&currentPositionsToComp,&firstPVcf,&positionsToErase,
				&doNotRescueVariantCallsAccrossTargets]() {
					if(currentPositionsToComp.size() > 1) {
						uint32_t bestPos = std::numeric_limits<uint32_t>::max();
						uint32_t bestNS = 0;
						for(const auto & checkPos : currentPositionsToComp) {
							auto currentNS = firstPVcf.records_[checkPos].info_.getMeta<uint32_t>("NS");
							if(currentNS > bestNS) {
								bestNS = currentNS;
								bestPos = checkPos;
							}
						}
						for(const auto & checkPos : currentPositionsToComp) {
							if(bestPos != checkPos) {
								positionsToErase.emplace(checkPos);
							}
						}
						if(!doNotRescueVariantCallsAccrossTargets) {

							std::set<std::string> blankSamples;
							std::regex blankDataPattern("\\.(,\\.)*");
							for(const auto & sampInfo : firstPVcf.records_[bestPos].sampleFormatInfos_) {
								//check to see if sample has all blank meta fields (e.g. . .,. .,.,. etc)
								if(std::all_of(sampInfo.second.meta_.begin(), sampInfo.second.meta_.end(),
									[&blankDataPattern](const auto & meta) {
										return std::regex_match(meta.second, blankDataPattern);
									})) {
									blankSamples.emplace(sampInfo.first);
								}
							}

							for(const auto & samp : blankSamples) {
								std::vector<uint32_t> rowPositionsWithNonBlankSamples;
								for(const auto & checkPos : currentPositionsToComp) {
									if (bestPos != checkPos) {
										//check to see if this blank sample has data for the overlapping variant call
										if (!std::all_of(firstPVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.begin(),
										                firstPVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.end(),
										                [&blankDataPattern](const auto& meta) {
											                return std::regex_match(meta.second, blankDataPattern);
										                }) &&
										                "." != firstPVcf.records_[checkPos].sampleFormatInfos_[samp].getMeta("DP")) {
											rowPositionsWithNonBlankSamples.emplace_back(checkPos);
										}
									}
								}

								if(!rowPositionsWithNonBlankSamples.empty()) {
									//replace with the row with the best read depth
									uint32_t bestReadDepth = 0;
									uint32_t bestRowPositionWithData = std::numeric_limits<uint32_t>::max();
									for(const auto currentPos : rowPositionsWithNonBlankSamples) {

										auto readDepth = firstPVcf.records_[currentPos].sampleFormatInfos_[samp].getMeta<uint32_t>("DP");
										if(readDepth > bestReadDepth) {
											bestReadDepth = readDepth;
											bestRowPositionWithData = currentPos;
										}
									}
									firstPVcf.records_[bestPos].sampleFormatInfos_[samp] = firstPVcf.records_[bestRowPositionWithData].sampleFormatInfos_[samp];
									//add to bestRefPos the NS, AN, AC
									//// NS will increase by 1 (though have to investigate that it's possible a sample is blank from the first initial call becauase it might already be counted)
									firstPVcf.records_[bestPos].info_.addMeta("NS",firstPVcf.records_[bestPos].info_.getMeta<uint32_t>("NS") + 1,true);
									//// AN will increase by the number non-zero allele calls for the sample
									auto ADs = tokenizeString(firstPVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD"), ",");
									auto ANcount = std::count_if(ADs.begin(), ADs.end(), [](const std::string & str){ return "0" != str;});
									firstPVcf.records_[bestPos].info_.addMeta("AN",firstPVcf.records_[bestPos].info_.getMeta<uint32_t>("AN") + ANcount,true);
									//// AC will increase the counts of allele that arent' 0 (and since we are adding samples with data there shouldn't be .)
									//// subsequently AF will also have to re-calculated
									////// will have to reconstruct the AC/AF for the best ref pos

									////// AC
									auto AC = firstPVcf.records_[bestPos].info_.getMeta("AC");
									auto AC_toks = tokenizeString(AC, ",");
									//the ADs go ref, variant1, variant2 etc
									if(AC_toks.size() + 1 != ADs.size()) {
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error " << "AC_toks.size() + 1 should equal ADs.size()" << "\n";
										ss << "AC_toks.size() + 1: " << AC_toks.size() + 1 << ", " << "ADs.size(): " << ADs.size() << "\n";
										throw std::runtime_error{ss.str()};
									}
									for (const auto& ADs_tok: iter::enumerate(ADs)) {
										//skip reference
										if (ADs_tok.index != 0) {
											//if ADs_tok.element does not equal 0 then add 1
											if (ADs_tok.element != "0") {
												AC_toks[ADs_tok.index - 1] = estd::to_string(
													njh::StrToNumConverter::stoToNum<uint32_t>(AC_toks[ADs_tok.index - 1]) + 1);
											}
										}
									}
									//replace the updated AC
									firstPVcf.records_[bestPos].info_.addMeta("AC", njh::conToStr(AC_toks, ","), true);
									//////// AF
									// re-calculate the AFs based on the new AC and AN
									double AN = firstPVcf.records_[bestPos].info_.getMeta<uint32_t>("AN");
									VecStr AFs;
									for(const auto & AC_tok : AC_toks) {
										AFs.emplace_back(estd::to_string(njh::StrToNumConverter::stoToNum<uint32_t>(AC_tok)/AN));
									}
									//replace the updated AF
									firstPVcf.records_[bestPos].info_.addMeta("AF", njh::conToStr(AFs, ","), true);
								}
							}
						}
					}

					currentPositionsToComp.clear();
				};
				for(uint32_t pos = 1; pos < firstPVcf.records_.size(); ++pos) {
					//if same exact position and refrence info than this needs to be compared
					if(firstPVcf.records_[pos].pos_ == firstPVcf.records_[pos-1].pos_ &&
						firstPVcf.records_[pos].ref_ == firstPVcf.records_[pos-1].ref_) {
						currentPositionsToComp.emplace(pos);
						currentPositionsToComp.emplace(pos-1);
					} else {
						compSamePositions();
					}
				}
				compSamePositions();

				for(const auto posToErase : iter::reversed(positionsToErase)) {
					firstPVcf.records_.erase(firstPVcf.records_.begin() + posToErase);
				}
			}
		}
		{
			OutputStream pvcf(njh::files::make_path(reportsDir, "allProteinVariantCalls.vcf"));
			firstPVcf.writeOutFixedAndSampleMeta(pvcf);
		}

		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream pvcfOutFile(njh::files::make_path(reportsDir, "knownAAChangesProteinVariantCalls.vcf"));
			std::vector<GenomicRegion> knownAAVariantRegions;
			knownAAVariantRegions.reserve(locs.transcriptLocs.size());
			for (const auto& b: locs.transcriptLocs) {
				knownAAVariantRegions.emplace_back(b);
			}
			firstPVcf.writeOutFixedAndSampleMeta(pvcfOutFile, knownAAVariantRegions);
		}
	}
	//process genomic
	{
		auto firstGVcfFnp = genomicVcfs.begin()->second.front();
		auto firstGVcf = VCFOutput::readInHeader(firstGVcfFnp);
		{
			InputStream in(firstGVcfFnp);
			firstGVcf.addInRecordsFromFile(in);
		}
		//add in blanks for any missing samples for any of the records
		firstGVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		if(firstGVcf.records_.empty()) {
			firstGVcf.samples_ = VecStr{sampleNamesSet.begin(), sampleNamesSet.end()};
		}
		for(const auto & target : genomicVcfs) {
			for(const auto & gvcfFnp : target.second) {
				if(firstGVcfFnp == gvcfFnp) {
					continue;
				}
				auto currentGvcf = VCFOutput::readInHeader(gvcfFnp);

				{
					InputStream in(gvcfFnp);
					currentGvcf.addInRecordsFromFile(in);
					// firstGVcf.addInRecordsFromFile(in);
				}
				//add in header
				//given we know what generated these vcf files there's not need to check the FORMAT, FILTER, or INFO tags

				//add in new contigs if any
				for(const auto & contig : currentGvcf.contigEntries_) {
					if(!njh::in(contig.first, firstGVcf.contigEntries_)) {
						firstGVcf.contigEntries_.emplace(contig);
					} else if(!(firstGVcf.contigEntries_.at(contig.first) == contig.second)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding contig: " << contig.first << " but already have an entry for this but it doesn't match" << "\n";
						ss << "contig in master header: " << njh::json::writeAsOneLine(njh::json::toJson(firstGVcf.contigEntries_.at(contig.first))) << "\n";
						ss << "adding contig: " << njh::json::writeAsOneLine(njh::json::toJson(contig.second)) << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
				currentGvcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
				firstGVcf.records_.reserve(firstGVcf.records_.size() + currentGvcf.records_.size());
				for(auto & record : currentGvcf.records_) {
					firstGVcf.records_.emplace_back(std::move(record));
				}
				//std::move(currentPvcf.records_.begin(), currentPvcf.records_.back(). std::back_inserter(firstPVcf.records_));
			}
		}
		firstGVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		firstGVcf.sortRecords();
		//eliminate duplicate variant calls
		{

			if(firstGVcf.records_.size() > 1) {
				std::set<uint32_t> positionsToErase;
				std::set<uint32_t> currentPositionsToComp;
				auto compSamePositions = [&currentPositionsToComp,&firstGVcf,&positionsToErase,
					&doNotRescueVariantCallsAccrossTargets]() {
					if(currentPositionsToComp.size() > 1) {
						uint32_t bestPos = std::numeric_limits<uint32_t>::max();
						uint32_t bestNS = 0;
						for(const auto & checkPos : currentPositionsToComp) {
							auto currentNS = firstGVcf.records_[checkPos].info_.getMeta<uint32_t>("NS");
							if(currentNS > bestNS) {
								bestNS = currentNS;
								bestPos = checkPos;
							}
						}
						for(const auto & checkPos : currentPositionsToComp) {
							if(bestPos != checkPos) {
								positionsToErase.emplace(checkPos);
							}
						}
						if(!doNotRescueVariantCallsAccrossTargets) {

							std::set<std::string> blankSamples;
							std::regex blankDataPattern("\\.(,\\.)*");
							for(const auto & sampInfo : firstGVcf.records_[bestPos].sampleFormatInfos_) {
								//check to see if sample has all blank meta fields (e.g. . .,. .,.,. etc)
								if(std::all_of(sampInfo.second.meta_.begin(), sampInfo.second.meta_.end(),
									[&blankDataPattern](const auto & meta) {
										return std::regex_match(meta.second, blankDataPattern);
									})) {
									blankSamples.emplace(sampInfo.first);
								}
							}

							for(const auto & samp : blankSamples) {
								std::vector<uint32_t> rowPositionsWithNonBlankSamples;
								for(const auto & checkPos : currentPositionsToComp) {
									if (bestPos != checkPos) {
										//check to see if this blank sample has data for the overlapping variant call
										if (!std::all_of(firstGVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.begin(),
										                firstGVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.end(),
										                [&blankDataPattern](const auto& meta) {
											                return std::regex_match(meta.second, blankDataPattern);
										                }) &&
										                "." != firstGVcf.records_[checkPos].sampleFormatInfos_[samp].getMeta("DP")) {
											rowPositionsWithNonBlankSamples.emplace_back(checkPos);
										}
									}
								}

								if(!rowPositionsWithNonBlankSamples.empty()) {
									//replace with the row with the best read depth
									uint32_t bestReadDepth = 0;
									uint32_t bestRowPositionWithData = std::numeric_limits<uint32_t>::max();
									for(const auto currentPos : rowPositionsWithNonBlankSamples) {

										auto readDepth = firstGVcf.records_[currentPos].sampleFormatInfos_[samp].getMeta<uint32_t>("DP");
										if(readDepth > bestReadDepth) {
											bestReadDepth = readDepth;
											bestRowPositionWithData = currentPos;
										}
									}
									firstGVcf.records_[bestPos].sampleFormatInfos_[samp] = firstGVcf.records_[bestRowPositionWithData].sampleFormatInfos_[samp];
									//add to bestRefPos the NS, AN, AC
									//// NS will increase by 1 (though have to investigate that it's possible a sample is blank from the first initial call becauase it might already be counted)
									firstGVcf.records_[bestPos].info_.addMeta("NS",firstGVcf.records_[bestPos].info_.getMeta<uint32_t>("NS") + 1,true);
									//// AN will increase by the number non-zero allele calls for the sample
									auto ADs = tokenizeString(firstGVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD"), ",");
									auto ANcount = std::count_if(ADs.begin(), ADs.end(), [](const std::string & str){ return "0" != str;});
									firstGVcf.records_[bestPos].info_.addMeta("AN",firstGVcf.records_[bestPos].info_.getMeta<uint32_t>("AN") + ANcount,true);
									//// AC will increase the counts of allele that arent' 0 (and since we are adding samples with data there shouldn't be .)
									//// subsequently AF will also have to re-calculated
									////// will have to reconstruct the AC/AF for the best ref pos

									////// AC
									auto AC = firstGVcf.records_[bestPos].info_.getMeta("AC");
									auto AC_toks = tokenizeString(AC, ",");
									//the ADs go ref, variant1, variant2 etc
									if(AC_toks.size() + 1 != ADs.size()) {
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error " << "AC_toks.size() + 1 should equal ADs.size()" << "\n";
										ss << "AC_toks.size() + 1: " << AC_toks.size() + 1 << ", " << "ADs.size(): " << ADs.size() << "\n";
										throw std::runtime_error{ss.str()};
									}
									for (const auto& ADs_tok: iter::enumerate(ADs)) {
										//skip reference
										if (ADs_tok.index != 0) {
											//if ADs_tok.element does not equal 0 then add 1
											if (ADs_tok.element != "0") {
												AC_toks[ADs_tok.index - 1] = estd::to_string(
													njh::StrToNumConverter::stoToNum<uint32_t>(AC_toks[ADs_tok.index - 1]) + 1);
											}
										}
									}
									//replace the updated AC
									firstGVcf.records_[bestPos].info_.addMeta("AC", njh::conToStr(AC_toks, ","), true);
									//////// AF
									// re-calculate the AFs based on the new AC and AN
									double AN = firstGVcf.records_[bestPos].info_.getMeta<uint32_t>("AN");
									VecStr AFs;
									for(const auto & AC_tok : AC_toks) {
										AFs.emplace_back(estd::to_string(njh::StrToNumConverter::stoToNum<uint32_t>(AC_tok)/AN));
									}
									//replace the updated AF
									firstGVcf.records_[bestPos].info_.addMeta("AF", njh::conToStr(AFs, ","), true);
								}
							}
						}
					}

					currentPositionsToComp.clear();
				};
				for(uint32_t pos = 1; pos < firstGVcf.records_.size(); ++pos) {
					if(firstGVcf.records_[pos].pos_ == firstGVcf.records_[pos-1].pos_ &&
						firstGVcf.records_[pos].ref_ == firstGVcf.records_[pos-1].ref_) {
						currentPositionsToComp.emplace(pos);
						currentPositionsToComp.emplace(pos-1);
						} else {
							compSamePositions();
						}
				}
				compSamePositions();
				// std::cout << njh::conToStr(positionsToErase, ",") << std::endl;
				for(const auto posToErase : iter::reversed(positionsToErase)) {
					firstGVcf.records_.erase(firstGVcf.records_.begin() + posToErase);
				}
			}
		}
		{
			OutputStream gvcfOutFile(njh::files::make_path(reportsDir, "allGenomicVariantCalls.vcf"));
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile);
		}
		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream gvcfOutFile(njh::files::make_path(reportsDir, "knownAAChangesGenomicVariantCalls.vcf"));
			std::vector<GenomicRegion> knownSnpVariantRegions;
			knownSnpVariantRegions.reserve(locs.genomicLocs.size());
			for (const auto& b: locs.genomicLocs) {
				knownSnpVariantRegions.emplace_back(b);
			}
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile, knownSnpVariantRegions);
		}
	}

	//filter vcfs to just known amino acid changes

	//add flag to count singler changes per certain meta and the haplotype known aa typed


	return 0;

}

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

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(doNotRescueVariantCallsAccrossTargets, "--doNotRescueVariantCallsAccrossTargets", "do Not Rescue Variant Calls Accross Targets");

	setUp.setOption(selectTargets, "--selectTargets", "Only analzye these select targets");
	setUp.setOption(selectSamples, "--selectSamples", "Only analzye these select samples");

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



	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
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
			auto c_ReadCnt = njh::StrToNumConverter::stoToNum<double>(
							row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);
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
			bool add = true;
			for (auto &seq: allResultSeqs[target][sample]) {
				if (seq.seq_ == clus.seq_) {
					add = false;
					seq.cnt_ += c_ReadCnt;
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


	njh::concurrent::LockableQueue<std::string> targetNamesQueue(targetNamesSet);

	std::function<void()> callVariants = [&collapseVarCallPars,&targetNamesQueue, &setUp]() {
		std::string target;
		while(targetNamesQueue.getVal(target)) {
			auto inputSeqsOpts = SeqIOOptions::genFastaInGz(njh::files::make_path(setUp.pars_.directoryName_, target, "inputSeqs.fasta.gz"));
			auto inputSeqs = SeqInput::getSeqVec<seqInfo>(inputSeqsOpts);
			const auto varCallDirPath = njh::files::make_path(setUp.pars_.directoryName_, target,  "variantCalling");
			auto collapseVarCallParsForTar = collapseVarCallPars;
			collapseVarCallParsForTar.identifier = target;
			collapseVarCallParsForTar.outputDirectory = varCallDirPath;
			collapseAndCallVariants(collapseVarCallParsForTar, inputSeqs);
		}
	};

	njh::concurrent::runVoidFunctionThreaded(callVariants, numThreads	);

	//concatenate diversity calls
	auto reportsDir = njh::files::make_path(setUp.pars_.directoryName_, "reports");
	njh::files::makeDir(njh::files::MkdirPar{reportsDir});
	auto targetNamesVec = std::vector<std::string>(targetNamesSet.begin(), targetNamesSet.end());
	njh::naturalSortNameSet(targetNamesVec);

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
	if(!metaFnp.empty() && exists(metaFnp)) {
		bfs::copy_file(metaFnp, njh::files::make_path(reportsDir, "meta.tsv"));
	}
	TranslatorByAlignment::GetGenomicLocationsForAminoAcidPositionsRet locs;
	if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
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


	//gather vcfs
	std::unordered_map<std::string, std::vector<bfs::path>> proteinVcfs;
	std::unordered_map<std::string, std::vector<bfs::path>> genomicVcfs;

	for (const auto& target: targetNamesVec) {
		auto proteinVcfFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"),false,
				std::vector<std::regex>{std::regex{".*-protein.vcf"}});
		auto genomicVcfFiles = njh::files::listAllFiles(njh::files::make_path(setUp.pars_.directoryName_, "/", target, "/variantCalling/variantCalls/"),false,
			std::vector<std::regex>{std::regex{".*-genomic.vcf"}});
		for(const auto & pvcf : proteinVcfFiles) {
			proteinVcfs[target].emplace_back(pvcf.first);
		}
		for(const auto & gvcf : genomicVcfFiles) {
			genomicVcfs[target].emplace_back(gvcf.first);
		}
	}

	//processing protein vcfs;
	{
		auto firstPVcfFnp = proteinVcfs.begin()->second.front();
		auto firstPVcf = VCFOutput::readInHeader(firstPVcfFnp);
		{
			InputStream in(firstPVcfFnp);
			firstPVcf.addInRecordsFromFile(in);
		}
		//add in blanks for any missing samples for any of the records
		firstPVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		if(firstPVcf.records_.empty()) {
			firstPVcf.samples_ = VecStr{sampleNamesSet.begin(), sampleNamesSet.end()};
		}
		for(const auto & target : proteinVcfs) {
			for(const auto & pvcfFnp : target.second) {
				if(firstPVcfFnp == pvcfFnp) {
					continue;
				}
				auto currentPvcf = VCFOutput::readInHeader(pvcfFnp);

				{
					InputStream in(pvcfFnp);
					currentPvcf.addInRecordsFromFile(in);
					// firstPVcf.addInRecordsFromFile(in);
				}
				//add in header
				//given we know what generated these vcf files there's not need to check the FORMAT, FILTER, or INFO tags

				//add in new contigs if any
				for(const auto & contig : currentPvcf.contigEntries_) {
					if(!njh::in(contig.first, firstPVcf.contigEntries_)) {
						firstPVcf.contigEntries_.emplace(contig);
					} else if(!(firstPVcf.contigEntries_.at(contig.first) == contig.second)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding contig: " << contig.first << " but already have an entry for this but it doesn't match" << "\n";
						ss << "contig in master header: " << njh::json::writeAsOneLine(njh::json::toJson(firstPVcf.contigEntries_.at(contig.first))) << "\n";
						ss << "adding contig: " << njh::json::writeAsOneLine(njh::json::toJson(contig.second)) << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
				currentPvcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
				firstPVcf.records_.reserve(firstPVcf.records_.size() + currentPvcf.records_.size());
				for(auto & record : currentPvcf.records_) {
					firstPVcf.records_.emplace_back(std::move(record));
				}
				//std::move(currentPvcf.records_.begin(), currentPvcf.records_.back(). std::back_inserter(firstPVcf.records_));
			}
		}
		firstPVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		firstPVcf.sortRecords();
		//eliminate duplicate variant calls and optional merge variants calls across amplicons that overlap that may have been missing in one target for a sample but present others

		{

			if(firstPVcf.records_.size() > 1) {
				std::set<uint32_t> positionsToErase;
				std::set<uint32_t> currentPositionsToComp;
				auto compSamePositions = [&currentPositionsToComp,&firstPVcf,&positionsToErase,
				&doNotRescueVariantCallsAccrossTargets]() {
					if(currentPositionsToComp.size() > 1) {
						uint32_t bestPos = std::numeric_limits<uint32_t>::max();
						uint32_t bestNS = 0;
						for(const auto & checkPos : currentPositionsToComp) {
							auto currentNS = firstPVcf.records_[checkPos].info_.getMeta<uint32_t>("NS");
							if(currentNS > bestNS) {
								bestNS = currentNS;
								bestPos = checkPos;
							}
						}
						for(const auto & checkPos : currentPositionsToComp) {
							if(bestPos != checkPos) {
								positionsToErase.emplace(checkPos);
							}
						}
						if(!doNotRescueVariantCallsAccrossTargets) {

							std::set<std::string> blankSamples;
							std::regex blankDataPattern("\\.(,\\.)*");
							for(const auto & sampInfo : firstPVcf.records_[bestPos].sampleFormatInfos_) {
								//check to see if sample has all blank meta fields (e.g. . .,. .,.,. etc)
								if(std::all_of(sampInfo.second.meta_.begin(), sampInfo.second.meta_.end(),
									[&blankDataPattern](const auto & meta) {
										return std::regex_match(meta.second, blankDataPattern);
									})) {
									blankSamples.emplace(sampInfo.first);
								}
							}

							for(const auto & samp : blankSamples) {
								std::vector<uint32_t> rowPositionsWithNonBlankSamples;
								for(const auto & checkPos : currentPositionsToComp) {
									if (bestPos != checkPos) {
										//check to see if this blank sample has data for the overlapping variant call
										if (!std::all_of(firstPVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.begin(),
										                firstPVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.end(),
										                [&blankDataPattern](const auto& meta) {
											                return std::regex_match(meta.second, blankDataPattern);
										                }) &&
										                "." != firstPVcf.records_[checkPos].sampleFormatInfos_[samp].getMeta("DP")) {
											rowPositionsWithNonBlankSamples.emplace_back(checkPos);
										}
									}
								}

								if(!rowPositionsWithNonBlankSamples.empty()) {
									//replace with the row with the best read depth
									uint32_t bestReadDepth = 0;
									uint32_t bestRowPositionWithData = std::numeric_limits<uint32_t>::max();
									for(const auto currentPos : rowPositionsWithNonBlankSamples) {

										auto readDepth = firstPVcf.records_[currentPos].sampleFormatInfos_[samp].getMeta<uint32_t>("DP");
										if(readDepth > bestReadDepth) {
											bestReadDepth = readDepth;
											bestRowPositionWithData = currentPos;
										}
									}
									firstPVcf.records_[bestPos].sampleFormatInfos_[samp] = firstPVcf.records_[bestRowPositionWithData].sampleFormatInfos_[samp];
									//add to bestRefPos the NS, AN, AC
									//// NS will increase by 1 (though have to investigate that it's possible a sample is blank from the first initial call becauase it might already be counted)
									firstPVcf.records_[bestPos].info_.addMeta("NS",firstPVcf.records_[bestPos].info_.getMeta<uint32_t>("NS") + 1,true);
									//// AN will increase by the number non-zero allele calls for the sample
									auto ADs = tokenizeString(firstPVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD"), ",");
									auto ANcount = std::count_if(ADs.begin(), ADs.end(), [](const std::string & str){ return "0" != str;});
									firstPVcf.records_[bestPos].info_.addMeta("AN",firstPVcf.records_[bestPos].info_.getMeta<uint32_t>("AN") + ANcount,true);
									//// AC will increase the counts of allele that arent' 0 (and since we are adding samples with data there shouldn't be .)
									//// subsequently AF will also have to re-calculated
									////// will have to reconstruct the AC/AF for the best ref pos

									////// AC
									auto AC = firstPVcf.records_[bestPos].info_.getMeta("AC");
									auto AC_toks = tokenizeString(AC, ",");
									//the ADs go ref, variant1, variant2 etc
									if(AC_toks.size() + 1 != ADs.size()) {
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error " << "AC_toks.size() + 1 should equal ADs.size()" << "\n";
										ss << "AC_toks.size() + 1: " << AC_toks.size() + 1 << ", " << "ADs.size(): " << ADs.size() << "\n";
										throw std::runtime_error{ss.str()};
									}
									for (const auto& ADs_tok: iter::enumerate(ADs)) {
										//skip reference
										if (ADs_tok.index != 0) {
											//if ADs_tok.element does not equal 0 then add 1
											if (ADs_tok.element != "0") {
												AC_toks[ADs_tok.index - 1] = estd::to_string(
													njh::StrToNumConverter::stoToNum<uint32_t>(AC_toks[ADs_tok.index - 1]) + 1);
											}
										}
									}
									//replace the updated AC
									firstPVcf.records_[bestPos].info_.addMeta("AC", njh::conToStr(AC_toks, ","), true);
									//////// AF
									// re-calculate the AFs based on the new AC and AN
									double AN = firstPVcf.records_[bestPos].info_.getMeta<uint32_t>("AN");
									VecStr AFs;
									for(const auto & AC_tok : AC_toks) {
										AFs.emplace_back(estd::to_string(njh::StrToNumConverter::stoToNum<uint32_t>(AC_tok)/AN));
									}
									//replace the updated AF
									firstPVcf.records_[bestPos].info_.addMeta("AF", njh::conToStr(AFs, ","), true);
								}
							}
						}
					}

					currentPositionsToComp.clear();
				};
				for(uint32_t pos = 1; pos < firstPVcf.records_.size(); ++pos) {
					//if same exact position and refrence info than this needs to be compared
					if(firstPVcf.records_[pos].pos_ == firstPVcf.records_[pos-1].pos_ &&
						firstPVcf.records_[pos].ref_ == firstPVcf.records_[pos-1].ref_) {
						currentPositionsToComp.emplace(pos);
						currentPositionsToComp.emplace(pos-1);
					} else {
						compSamePositions();
					}
				}
				compSamePositions();

				for(const auto posToErase : iter::reversed(positionsToErase)) {
					firstPVcf.records_.erase(firstPVcf.records_.begin() + posToErase);
				}
			}
		}
		{
			OutputStream pvcf(njh::files::make_path(reportsDir, "allProteinVariantCalls.vcf"));
			firstPVcf.writeOutFixedAndSampleMeta(pvcf);
		}

		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream pvcfOutFile(njh::files::make_path(reportsDir, "knownAAChangesProteinVariantCalls.vcf"));
			std::vector<GenomicRegion> knownAAVariantRegions;
			knownAAVariantRegions.reserve(locs.transcriptLocs.size());
			for (const auto& b: locs.transcriptLocs) {
				knownAAVariantRegions.emplace_back(b);
			}
			firstPVcf.writeOutFixedAndSampleMeta(pvcfOutFile, knownAAVariantRegions);
		}
	}
	//process genomic
	{
		auto firstGVcfFnp = genomicVcfs.begin()->second.front();
		auto firstGVcf = VCFOutput::readInHeader(firstGVcfFnp);
		{
			InputStream in(firstGVcfFnp);
			firstGVcf.addInRecordsFromFile(in);
		}
		//add in blanks for any missing samples for any of the records
		firstGVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		if(firstGVcf.records_.empty()) {
			firstGVcf.samples_ = VecStr{sampleNamesSet.begin(), sampleNamesSet.end()};
		}
		for(const auto & target : genomicVcfs) {
			for(const auto & gvcfFnp : target.second) {
				if(firstGVcfFnp == gvcfFnp) {
					continue;
				}
				auto currentGvcf = VCFOutput::readInHeader(gvcfFnp);

				{
					InputStream in(gvcfFnp);
					currentGvcf.addInRecordsFromFile(in);
					// firstGVcf.addInRecordsFromFile(in);
				}
				//add in header
				//given we know what generated these vcf files there's not need to check the FORMAT, FILTER, or INFO tags

				//add in new contigs if any
				for(const auto & contig : currentGvcf.contigEntries_) {
					if(!njh::in(contig.first, firstGVcf.contigEntries_)) {
						firstGVcf.contigEntries_.emplace(contig);
					} else if(!(firstGVcf.contigEntries_.at(contig.first) == contig.second)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding contig: " << contig.first << " but already have an entry for this but it doesn't match" << "\n";
						ss << "contig in master header: " << njh::json::writeAsOneLine(njh::json::toJson(firstGVcf.contigEntries_.at(contig.first))) << "\n";
						ss << "adding contig: " << njh::json::writeAsOneLine(njh::json::toJson(contig.second)) << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
				currentGvcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
				firstGVcf.records_.reserve(firstGVcf.records_.size() + currentGvcf.records_.size());
				for(auto & record : currentGvcf.records_) {
					firstGVcf.records_.emplace_back(std::move(record));
				}
				//std::move(currentPvcf.records_.begin(), currentPvcf.records_.back(). std::back_inserter(firstPVcf.records_));
			}
		}
		firstGVcf.addInBlnaksForAnyMissingSamples(sampleNamesSet);
		firstGVcf.sortRecords();
		//eliminate duplicate variant calls
		{

			if(firstGVcf.records_.size() > 1) {
				std::set<uint32_t> positionsToErase;
				std::set<uint32_t> currentPositionsToComp;
				auto compSamePositions = [&currentPositionsToComp,&firstGVcf,&positionsToErase,
					&doNotRescueVariantCallsAccrossTargets]() {
					if(currentPositionsToComp.size() > 1) {
						uint32_t bestPos = std::numeric_limits<uint32_t>::max();
						uint32_t bestNS = 0;
						for(const auto & checkPos : currentPositionsToComp) {
							auto currentNS = firstGVcf.records_[checkPos].info_.getMeta<uint32_t>("NS");
							if(currentNS > bestNS) {
								bestNS = currentNS;
								bestPos = checkPos;
							}
						}
						for(const auto & checkPos : currentPositionsToComp) {
							if(bestPos != checkPos) {
								positionsToErase.emplace(checkPos);
							}
						}
						if(!doNotRescueVariantCallsAccrossTargets) {

							std::set<std::string> blankSamples;
							std::regex blankDataPattern("\\.(,\\.)*");
							for(const auto & sampInfo : firstGVcf.records_[bestPos].sampleFormatInfos_) {
								//check to see if sample has all blank meta fields (e.g. . .,. .,.,. etc)
								if(std::all_of(sampInfo.second.meta_.begin(), sampInfo.second.meta_.end(),
									[&blankDataPattern](const auto & meta) {
										return std::regex_match(meta.second, blankDataPattern);
									})) {
									blankSamples.emplace(sampInfo.first);
								}
							}

							for(const auto & samp : blankSamples) {
								std::vector<uint32_t> rowPositionsWithNonBlankSamples;
								for(const auto & checkPos : currentPositionsToComp) {
									if (bestPos != checkPos) {
										//check to see if this blank sample has data for the overlapping variant call
										if (!std::all_of(firstGVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.begin(),
										                firstGVcf.records_[checkPos].sampleFormatInfos_[samp].meta_.end(),
										                [&blankDataPattern](const auto& meta) {
											                return std::regex_match(meta.second, blankDataPattern);
										                }) &&
										                "." != firstGVcf.records_[checkPos].sampleFormatInfos_[samp].getMeta("DP")) {
											rowPositionsWithNonBlankSamples.emplace_back(checkPos);
										}
									}
								}

								if(!rowPositionsWithNonBlankSamples.empty()) {
									//replace with the row with the best read depth
									uint32_t bestReadDepth = 0;
									uint32_t bestRowPositionWithData = std::numeric_limits<uint32_t>::max();
									for(const auto currentPos : rowPositionsWithNonBlankSamples) {

										auto readDepth = firstGVcf.records_[currentPos].sampleFormatInfos_[samp].getMeta<uint32_t>("DP");
										if(readDepth > bestReadDepth) {
											bestReadDepth = readDepth;
											bestRowPositionWithData = currentPos;
										}
									}
									firstGVcf.records_[bestPos].sampleFormatInfos_[samp] = firstGVcf.records_[bestRowPositionWithData].sampleFormatInfos_[samp];
									//add to bestRefPos the NS, AN, AC
									//// NS will increase by 1 (though have to investigate that it's possible a sample is blank from the first initial call becauase it might already be counted)
									firstGVcf.records_[bestPos].info_.addMeta("NS",firstGVcf.records_[bestPos].info_.getMeta<uint32_t>("NS") + 1,true);
									//// AN will increase by the number non-zero allele calls for the sample
									auto ADs = tokenizeString(firstGVcf.records_[bestPos].sampleFormatInfos_[samp].getMeta("AD"), ",");
									auto ANcount = std::count_if(ADs.begin(), ADs.end(), [](const std::string & str){ return "0" != str;});
									firstGVcf.records_[bestPos].info_.addMeta("AN",firstGVcf.records_[bestPos].info_.getMeta<uint32_t>("AN") + ANcount,true);
									//// AC will increase the counts of allele that arent' 0 (and since we are adding samples with data there shouldn't be .)
									//// subsequently AF will also have to re-calculated
									////// will have to reconstruct the AC/AF for the best ref pos

									////// AC
									auto AC = firstGVcf.records_[bestPos].info_.getMeta("AC");
									auto AC_toks = tokenizeString(AC, ",");
									//the ADs go ref, variant1, variant2 etc
									if(AC_toks.size() + 1 != ADs.size()) {
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error " << "AC_toks.size() + 1 should equal ADs.size()" << "\n";
										ss << "AC_toks.size() + 1: " << AC_toks.size() + 1 << ", " << "ADs.size(): " << ADs.size() << "\n";
										throw std::runtime_error{ss.str()};
									}
									for (const auto& ADs_tok: iter::enumerate(ADs)) {
										//skip reference
										if (ADs_tok.index != 0) {
											//if ADs_tok.element does not equal 0 then add 1
											if (ADs_tok.element != "0") {
												AC_toks[ADs_tok.index - 1] = estd::to_string(
													njh::StrToNumConverter::stoToNum<uint32_t>(AC_toks[ADs_tok.index - 1]) + 1);
											}
										}
									}
									//replace the updated AC
									firstGVcf.records_[bestPos].info_.addMeta("AC", njh::conToStr(AC_toks, ","), true);
									//////// AF
									// re-calculate the AFs based on the new AC and AN
									double AN = firstGVcf.records_[bestPos].info_.getMeta<uint32_t>("AN");
									VecStr AFs;
									for(const auto & AC_tok : AC_toks) {
										AFs.emplace_back(estd::to_string(njh::StrToNumConverter::stoToNum<uint32_t>(AC_tok)/AN));
									}
									//replace the updated AF
									firstGVcf.records_[bestPos].info_.addMeta("AF", njh::conToStr(AFs, ","), true);
								}
							}
						}
					}

					currentPositionsToComp.clear();
				};
				for(uint32_t pos = 1; pos < firstGVcf.records_.size(); ++pos) {
					if(firstGVcf.records_[pos].pos_ == firstGVcf.records_[pos-1].pos_ &&
						firstGVcf.records_[pos].ref_ == firstGVcf.records_[pos-1].ref_) {
						currentPositionsToComp.emplace(pos);
						currentPositionsToComp.emplace(pos-1);
						} else {
							compSamePositions();
						}
				}
				compSamePositions();
				// std::cout << njh::conToStr(positionsToErase, ",") << std::endl;
				for(const auto posToErase : iter::reversed(positionsToErase)) {
					firstGVcf.records_.erase(firstGVcf.records_.begin() + posToErase);
				}
			}
		}
		{
			OutputStream gvcfOutFile(njh::files::make_path(reportsDir, "allGenomicVariantCalls.vcf"));
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile);
		}
		if(!collapseVarCallPars.transPars.knownAminoAcidMutationsFnp_.empty()) {
			OutputStream gvcfOutFile(njh::files::make_path(reportsDir, "knownAAChangesGenomicVariantCalls.vcf"));
			std::vector<GenomicRegion> knownSnpVariantRegions;
			knownSnpVariantRegions.reserve(locs.genomicLocs.size());
			for (const auto& b: locs.genomicLocs) {
				knownSnpVariantRegions.emplace_back(b);
			}
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile, knownSnpVariantRegions);
		}
	}

	//filter vcfs to just known amino acid changes

	//add flag to count singler changes per certain meta and the haplotype known aa typed


	return 0;

}

}  // namespace njhseq