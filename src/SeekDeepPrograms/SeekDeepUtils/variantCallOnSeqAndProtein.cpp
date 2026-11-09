//
// Created by Nicholas Hathaway on 12/16/23.
//
#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/pmo/pmoUtils.hpp>
#include "SeekDeep/utils/VariantCallInput.hpp"

namespace njhseq {

namespace {

/**
 * \brief the amino acid three letter codes for a sequence with the first letter capitalized, e.g. "Asn" for N
 * \param aas the amino acids
 * \param allowX whether X/x is allowed and written as XXX
 * \return the concatenated codes
 */
std::string aaTriCodeName(const std::string & aas, const bool allowX) {
	std::string ret;
	for (const auto c : aas) {
		if (allowX && (c == 'X' || c == 'x')) {
			ret += "XXX";
		} else {
			auto currentTriCode = aminoAcidInfo::infos::allInfo.at(c).triCode_;
			currentTriCode[0] = static_cast<char>(toupper(currentTriCode[0]));
			ret += currentTriCode;
		}
	}
	return ret;
}

/**
 * \brief write the amino acid changes table, a row for every sample for every allele of every protein variant record
 * \param pvcf the combined protein vcf
 * \param knownAAVariantRegions the known amino acid change locations, used for the Targeted column
 * \param refInfoNextToAlt place the reference info next to each alternate's info rather than in its own row
 * \param out the stream to write to
 */
void writeAAChangesInfoTable(const VCFOutput & pvcf,
                             const std::vector<GenomicRegion> & knownAAVariantRegions,
                             const bool refInfoNextToAlt,
                             std::ostream & out) {
	out << "Gene_ID"
			<< "\t" << "Gene_Transcript_ID"
			<< "\t" << "Gene"
			<< "\t" << "Mutation_Name"
			<< "\t" << "ExonicFunc"
			<< "\t" << "AA_Change"
			<< "\t" << "Targeted"
			<< "\t" << "CoveredBy"
			<< "\t" << "sample";
	out << "\t" << "reference_AA_pos"
	    << "\t" << "reference_AA";
	if (refInfoNextToAlt) {
		out << "\t" << "alternate_AA";
		out << "\t" << "reference_AA_cnt"
		    << "\t" << "reference_AA_freq"
		    << "\t" << "alternate_AA_cnt"
		    << "\t" << "alternate_AA_freq"
		    << "\t" << "AA_withinSampleCoverage"
		    << "\t" << "AlleleCount"
		    << "\t" << "AlleleFrequency"
		    << "\t" << "SampleCount"
		    << "\t" << "SamplePrevalence";
	} else {
		out << "\t" << "AA"
		    << "\t" << "is_reference"
		    << "\t" << "AA_cnt"
		    << "\t" << "AA_freq"
		    << "\t" << "AA_withinSampleCoverage"
		    << "\t" << "AlleleCount"
		    << "\t" << "AlleleFrequency"
		    << "\t" << "AlleleTotal"
		    << "\t" << "SampleCount"
		    << "\t" << "SamplePrevalence"
		    << "\t" << "SampleTotal";
	}
	out << '\n';

	//the position of a FORMAT field in a record's sample values, will throw if the record doesn't have it
	auto getFieldIndex = [](const VCFOutput::VCFRecord & rec, const std::string & field) {
		const auto index = rec.formatKeyIndex(field);
		if (std::string::npos == index) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "no FORMAT field " << field << " for " << rec.chrom_ << ":" << rec.pos_ << ", current fields are: " << njh::conToStr(rec.formatKeys_, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
		return index;
	};
	//the rows for a record are built here and written once per sample
	std::string rows;
	for (const auto & rec : pvcf.records_) {
		std::string targeted = "No";
		const auto genomicRegion = rec.genRegion();
		for (const auto & knownRegion : knownAAVariantRegions) {
			if (genomicRegion.overlaps(knownRegion)) {
				targeted = "Yes";
				break;
			}
		}
		const auto TYPE = tokenizeString(rec.info_.getMeta("TYPE"), ",");
		const auto ACs = tokenizeString(rec.info_.getMeta("AC_REAL"), ",");
		const auto AFs = tokenizeString(rec.info_.getMeta("AF_REAL"), ",");
		const auto SCs = tokenizeString(rec.info_.getMeta("SC"), ",");
		const auto PREVs = tokenizeString(rec.info_.getMeta("PREV"), ",");
		const auto geneID = rec.info_.getMeta("GeneID");
		const auto geneName = rec.info_.getMeta("GeneName");
		const auto coveredBy = rec.info_.getMeta("TARGET");
		const auto refTriCodeName = aaTriCodeName(rec.ref_, false);
		const auto pos = estd::to_string(rec.pos_);
		//the columns before the sample name are the same for every sample
		const std::string recordStart = geneID + "\t" + rec.chrom_ + "\t" + geneName + "\t";
		const std::string recordEnd = "\t" + targeted + "\t" + coveredBy + "\t";

		//per alt row parts
		VecStr altRowStarts;
		VecStr altRowMids;
		VecStr altRowEnds;
		for (const auto & altEnum : iter::enumerate(rec.alts_)) {
			const auto & alt = altEnum.element;
			const auto altTriCodeName = aaTriCodeName(alt, true);
			std::string ExonicFunc;
			if (TYPE[altEnum.index] == "snp") {
				ExonicFunc = "missense_variant";
			} else if (TYPE[altEnum.index] == "del") {
				ExonicFunc = "conservative_inframe_deletion";
			} else if (TYPE[altEnum.index] == "ins") {
				ExonicFunc = "conservative_inframe_insertion";
			}
			altRowStarts.emplace_back(recordStart + njh::pasteAsStr(geneName, "-", refTriCodeName, rec.pos_, altTriCodeName)
			                          + "\t" + ExonicFunc
			                          + "\t" + njh::pasteAsStr(refTriCodeName, rec.pos_, altTriCodeName)
			                          + recordEnd);
			if (refInfoNextToAlt) {
				altRowMids.emplace_back("\t" + pos + "\t" + rec.ref_ + "\t" + alt);
				altRowEnds.emplace_back("\t" + ACs[altEnum.index] + "\t" + AFs[altEnum.index] + "\t" + SCs[altEnum.index] + "\t" + PREVs[altEnum.index] + '\n');
			} else {
				altRowMids.emplace_back("\t" + pos + "\t" + rec.ref_ + "\t" + alt + "\t" + "FALSE");
			}
		}

		const bool hasSamples = !rec.sampleValues_.empty();
		const auto DPIndex = hasSamples ? getFieldIndex(rec, "DP") : 0;
		const auto ADIndex = hasSamples ? getFieldIndex(rec, "AD") : 0;
		const auto AFIndex = hasSamples ? getFieldIndex(rec, "AF") : 0;
		if (refInfoNextToAlt) {
			for (const auto & sampleEnum : iter::enumerate(rec.sampleValues_)) {
				rows.clear();
				const auto & sampleName = pvcf.samples_.at(sampleEnum.index);
				const auto vals = VCFOutput::VCFRecord::splitSampleValues(sampleEnum.element);
				const std::string DP(vals.at(DPIndex));
				const auto sample_ADs = tokenizeString(std::string(vals.at(ADIndex)), ",");
				const auto sample_AFs = tokenizeString(std::string(vals.at(AFIndex)), ",");
				for (const auto altPos : iter::range(rec.alts_.size())) {
					rows += altRowStarts[altPos];
					rows += sampleName;
					rows += altRowMids[altPos];
					if ("." == DP) {
						//no coverage
						rows += "\t0\t0\t0\t0\t0";
					} else {
						rows += "\t" + sample_ADs[0] + "\t" + sample_AFs[0] + "\t" + sample_ADs[altPos + 1] + "\t" + sample_AFs[altPos + 1] + "\t" + DP;
					}
					rows += altRowEnds[altPos];
				}
				out << rows;
			}
			continue;
		}

		//reference counts
		const auto AN = njh::StrToNumConverter::stoToNum<uint32_t>(rec.info_.getMeta("AN_REAL"));
		const auto ref_AC = AN - vectorSum(vecStrToVecNum<uint32_t>(ACs));
		const auto ref_AF = 1 - vectorSum(vecStrToVecNum<double>(AFs));
		const auto NS = njh::StrToNumConverter::stoToNum<uint32_t>(rec.info_.getMeta("NS"));
		uint32_t ref_SC = 0;
		for (const auto & sampleValues : rec.sampleValues_) {
			const auto vals = VCFOutput::VCFRecord::splitSampleValues(sampleValues);
			if ("." != vals.at(DPIndex)) {
				const auto sample_ADs = tokenizeString(std::string(vals.at(ADIndex)), ",");
				if (njh::StrToNumConverter::stoToNum<uint32_t>(sample_ADs[0]) > 0) {
					++ref_SC;
				}
			}
		}
		const double ref_PREV = ref_SC / static_cast<double>(NS);
		const auto ANStr = estd::to_string(AN);
		const auto NSStr = estd::to_string(NS);
		for (const auto altPos : iter::range(rec.alts_.size())) {
			altRowEnds.emplace_back("\t" + ACs[altPos] + "\t" + AFs[altPos] + "\t" + ANStr + "\t" + SCs[altPos] + "\t" + PREVs[altPos] + "\t" + NSStr + '\n');
		}
		const std::string refRowStart = recordStart + njh::pasteAsStr(geneName, "-", refTriCodeName, rec.pos_)
		                                + "\t" + "NA"
		                                + "\t" + njh::pasteAsStr(refTriCodeName, rec.pos_)
		                                + recordEnd;
		const std::string refRowMid = "\t" + pos + "\t" + rec.ref_ + "\t" + rec.ref_ + "\t" + "TRUE";
		const std::string refRowEnd = "\t" + estd::to_string(ref_AC)
		                              + "\t" + estd::to_string(ref_AF)
		                              + "\t" + ANStr
		                              + "\t" + estd::to_string(ref_SC)
		                              + "\t" + estd::to_string(ref_PREV)
		                              + "\t" + NSStr + '\n';
		for (const auto & sampleEnum : iter::enumerate(rec.sampleValues_)) {
			rows.clear();
			const auto & sampleName = pvcf.samples_.at(sampleEnum.index);
			const auto vals = VCFOutput::VCFRecord::splitSampleValues(sampleEnum.element);
			const std::string DP(vals.at(DPIndex));
			const auto sample_ADs = tokenizeString(std::string(vals.at(ADIndex)), ",");
			const auto sample_AFs = tokenizeString(std::string(vals.at(AFIndex)), ",");
			for (const auto altPos : iter::range(rec.alts_.size())) {
				rows += altRowStarts[altPos];
				rows += sampleName;
				rows += altRowMids[altPos];
				if ("." == DP) {
					//no coverage
					rows += "\t0\t0\t0";
				} else {
					rows += "\t" + sample_ADs[altPos + 1] + "\t" + sample_AFs[altPos + 1] + "\t" + DP;
				}
				rows += altRowEnds[altPos];
			}
			//adding reference
			rows += refRowStart;
			rows += sampleName;
			rows += refRowMid;
			if ("." == DP) {
				//no coverage
				rows += "\t0\t0\t0";
			} else {
				rows += "\t" + sample_ADs[0] + "\t" + sample_AFs[0] + "\t" + DP;
			}
			rows += refRowEnd;
			out << rows;
		}
	}
}

} // namespace

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
	VCFOutput::combineVCFsPars combiningVcfPars;
	// bool doNotRescueVariantCallsAcrossTargets = false;
	std::string popSeqsRegexPatRemoval = R"(_([tf])?\d+(\.\d+)?$)";
	uint32_t numThreads = 1;
	bool aaSummaryRefInfoNextToAlt = false;
	bool skipMissingInputBedRegions = false;
	CollapseAndCallVariantsPars collapseVarCallPars;
	std::set<std::string> selectTargets;
	std::set<std::string> selectSamples;
	std::set<std::string> excludeTargets;
	std::set<std::string> excludeSamples;

	bfs::path bedLocs;
	bool genomicLocsChangePeriodToDash = false;
	bool requireFullKnownCoverage = true;
	bool doNotRequireFullKnownCoverage = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(doNotRequireFullKnownCoverage, "--doNotRequireFullKnownCoverage", "do Not Require Full Known Coverage e.g. will report covered by target if there is partial coverage");
	requireFullKnownCoverage = !doNotRequireFullKnownCoverage;
	setUp.setOption(aaSummaryRefInfoNextToAlt, "--aaSummaryRefInfoNextToAlt", "In AA Summary table place the Reference Info Next To Alternate info, rather than stacked");

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
	setUp.setOption(skipMissingInputBedRegions, "--skipMissingInputBedRegions", "when supplying a location skip warning on missing Input Bed Regions");


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
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.complexVarPars.fractionOfSamplesCoveringSpan, "--complexVariantFractionOfSamplesCoveringSpan", "Fraction of samples that must cover a complex variant's full span for it to be reported");


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
	bfs::path pmoFnp;
	PmoGenomePars pmoGenomePars;
	bfs::path perTargetVariantCallsDir;
	std::string perTargetVariantCallsSubDir;
	setUp.setOption(perTargetVariantCallsDir, "--perTargetVariantCallsDir", "Make the reports from variant calls already made per target instead of calling them, a target's calls (what collapseAndCallVariants writes) are in <this>/<target>/<--perTargetVariantCallsSubDir>, e.g. the variantCalling directories of processClusters and PathWeaver runProcessClustersOnRecon");
	setUp.setOption(perTargetVariantCallsSubDir, "--perTargetVariantCallsSubDir", "The directory within each target's directory of --perTargetVariantCallsDir with its variant calls, e.g. variantCalling for PathWeaver runProcessClustersOnRecon output");
	setUp.setOption(pmoFnp, "--pmo", "A PMO with the haplotypes in --resultsFnp, the variant calls are added to its haplotypes by sequence and it's written to reports/allPMO.json.gz", false, "PMO");
	pmoGenomePars.setOptions(setUp);
	setUp.setOption(popSeqsRegexPatRemoval, "--popSeqsRegexPatRemoval",
									"optional regex pattern to process the input pop sequences from --popSeqsFnp to make it match up with the input results folder");

	setUp.processDirectoryOutputName("variantCalling_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow({resultsFnp},
																	__PRETTY_FUNCTION__);

	Json::Value runLog;
	std::mutex runLogMut;

	if (!perTargetVariantCallsDir.empty() && !pmoFnp.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "--pmo can't be used with --perTargetVariantCallsDir, the calls aren't made here\n";
		throw std::runtime_error{ss.str()};
	}
	std::unique_ptr<pmo::PortableMicrohaplotypeObject> pmoObj;
	std::set<std::string> pmoTargets;
	std::mutex pmoMut;
	if (!pmoFnp.empty()) {
		pmoObj = std::make_unique<pmo::PortableMicrohaplotypeObject>(pmo::readPmo(pmoFnp));
		for (const auto & target : pmoObj->target_info_) {
			pmoTargets.emplace(target.target_name_);
		}
	}
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

	//the final haplotypes per target and sample
	VariantCallInput input;

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
		const std::regex popSeqsRegexPatRemovalPat{popSeqsRegexPatRemoval};
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
				seq.name_ = std::regex_replace(seq.name_, popSeqsRegexPatRemovalPat, "");
				hPopUID_to_hConsensus[target][seq.name_] = seq.seq_;
				hConsensus_to_hPopUID[target][seq.seq_] = seq.name_;
			}
		}
	}

	{
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
			const auto& h_popUID = row[sampInfoReader.header_.getColPos(popHapIdColName)];
			const std::string hapSeq = popSeqsDirFnp.empty()
					? row[sampInfoReader.header_.getColPos(popHapSeqColName)]
					: njh::mapAt(njh::mapAt(hPopUID_to_hConsensus, target), h_popUID);
			input.addHap(target, sample, row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)], h_popUID, hapSeq, row);
			inputSampleNamesSet.emplace(sample);
		}
	}
	//population haplotype ids come from the population seqs file names when given
	if (!popSeqsDirFnp.empty()) {
		for (const auto & target : hConsensus_to_hPopUID) {
			input.popHapIdBySeq_[target.first] = target.second;
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	//write out seqs
	for(const auto & tar : input.seqs_) {
		auto tarDir = njh::files::make_path(setUp.pars_.directoryName_, tar.first);
		njh::files::makeDir(njh::files::MkdirPar{tarDir});
		input.writeInputSeqs(tar.first, njh::files::make_path(setUp.pars_.directoryName_, tar.first, "inputSeqs.fasta.gz"), metaGroupData.get());
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> genomicLocs;
	if(!bedLocs.empty()) {
		genomicLocs = readVariantCallLocations(bedLocs, genomicLocsChangePeriodToDash);
		VecStr missingTargets;
		for(const auto & tar : targetNamesSet) {
			if(njh::notIn(tar, genomicLocs)) {
				missingTargets.emplace_back(tar);
			}
		}
		if(!missingTargets.empty() && !skipMissingInputBedRegions) {
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
					&genomicLocs, &runLogMut, &runLogTargetTimes, &pmoObj, &pmoTargets, &pmoMut, &pmoGenomePars,
					&perTargetVariantCallsDir, &perTargetVariantCallsSubDir]() {
			std::string target;
			while(targetNamesQueue.getVal(target)) {
				njh::stopWatch watch;
				Json::Value currentLog;
				currentLog["target"] = target;
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				const auto varCallDirPath = njh::files::make_path(setUp.pars_.directoryName_, target,  "variantCalling");
				if (!perTargetVariantCallsDir.empty()) {
					//link in the calls already made
					auto existingCalls = bfs::absolute(njh::files::make_path(perTargetVariantCallsDir, target, perTargetVariantCallsSubDir));
					if (!bfs::exists(njh::files::make_path(existingCalls, "uniqueSeqs.fasta.gz"))) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "no variant calls for target " << target << " in " << existingCalls << "\n";
						throw std::runtime_error{ss.str()};
					}
					bfs::create_directory_symlink(existingCalls, varCallDirPath);
					currentLog["totalTime"] = watch.totalTime();
					std::lock_guard<std::mutex> lock(runLogMut);
					runLogTargetTimes.append(currentLog);
					continue;
				}
				auto inputSeqsOpts = SeqIOOptions::genFastaInGz(njh::files::make_path(setUp.pars_.directoryName_, target, "inputSeqs.fasta.gz"));
				auto inputSeqs = SeqInput::getSeqVec<seqInfo>(inputSeqsOpts);
				auto collapseVarCallParsForTar = collapseVarCallPars;
				collapseVarCallParsForTar.identifier = target;
				collapseVarCallParsForTar.outputDirectory = varCallDirPath;
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				if(njh::in(target, genomicLocs)) {
					collapseVarCallParsForTar.refSeqRegion = GenomicRegion(*njh::mapAt(genomicLocs, target));
				}
				//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

				auto calls = collapseAndCallVariants(collapseVarCallParsForTar, inputSeqs);
				if (pmoObj && njh::in(target, pmoTargets)) {
					std::lock_guard<std::mutex> lock(pmoMut);
					auto added = addVariantCallsToPmo(*pmoObj, target, collapseVarCallParsForTar, calls, pmoGenomePars);
					currentLog["pmoHapsWithCalls"] = added.mhapsWithCalls_;
					currentLog["pmoHapsWithoutCalls"] = added.mhapsWithoutCalls_;
				}
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
	if (pmoObj) {
		pmoObj->validate();
		OutOptions pmoOutOpts(njh::files::make_path(reportsDir, "allPMO.json.gz"), ".gz");
		pmoOutOpts.overWriteFile_ = true;
		pmo::writePmo(*pmoObj, pmoOutOpts);
	}

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

	{//seqs_prev_freq.tsv.gz
		fullWatch.startNewLap("gather sequence freqs and prevs");
		//sequence diversity
		collectAndWriteMasterFile(njh::pasteAsStr("/variantCalling/seqs_prev_freq.tsv.gz"),
			njh::files::make_path(reportsPopGeneticsDir, "all_seqs_prev_freq.tsv.gz"));
	}
	{//translated_prev_freq.tsv.gz
		fullWatch.startNewLap("gather translated sequence freqs and prevs");
		//translated diversity
		collectAndWriteMasterFile(njh::pasteAsStr("/variantCalling/variantCalls/translated_prev_freq.tsv.gz"),
	njh::files::make_path(reportsPopGeneticsDir, "all_translated_prev_freq.tsv.gz"));
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
				hapCountsOut << "\t" << input.seqs_[tar][samp].size();
				uint32_t currentReadCount = 0;
				for(const auto & r : input.seqs_[tar][samp]) {
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
				std::vector<std::regex>{std::regex{".*_exonIntronPositions.bed"}, std::regex{".*_withUTR.bed"}});
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
		auto splitTab = table::splitColWithMeta(overlappingGeneInfo.getUniqueRows(), splitPars);
		OutputStream geneInfoTabout(njh::files::make_path(reportsInfoDir, "targetsIntersectingWithGenesInfo.tsv"));
		splitTab.outPutContents(geneInfoTabout, "\t");
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
		auto firstPVcf = VCFOutput::combineVCFs(proteinVcfs, sampleNamesSet, combiningVcfPars);
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
			OutputStream aminoAcidChangesTable(njh::files::make_path(reportsSummaryDir, "AAChangesInfo.tsv.gz"));
			writeAAChangesInfoTable(firstPVcf, knownAAVariantRegions, aaSummaryRefInfoNextToAlt, aminoAcidChangesTable);
		}
	}
	//process genomic
	if(!genomicVcfs.empty()) {
		fullWatch.startNewLap("combine genomic vcfs");
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		auto firstGVcf = VCFOutput::combineVCFs(genomicVcfs, sampleNamesSet, combiningVcfPars);
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
		fullWatch.startNewLap("combine complex genomic vcfs");
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		// std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		{
			auto firstGVcf = VCFOutput::combineVCFs(complexGenomicVcfs, sampleNamesSet, combiningVcfPars);
			OutputStream gvcfOutFile(njh::files::make_path(reportsVcfsDir, "allComplexGenomicVariantCalls.vcf.gz"));
			firstGVcf.writeOutFixedAndSampleMeta(gvcfOutFile);
		}
	}
	//combine summary tables, while going through the rows also gather the unique haplotypes and count the locations per target
	const VecStr perHaplotypeColumns{"target","seq","transcript","translatedSeq","transcript_fullAATyped","transcript_knownAATyped"};
	std::unique_ptr<table> summaryPerHaplotype;
	std::unique_ptr<table::GroupColumnsCounter> genomicLocCounter;
	std::unique_ptr<table::GroupColumnsCounter> proteinLocCounter;
	{
		fullWatch.startNewLap("gather summary tables");
		//summary tables
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
			firstTable.header_.checkForColumnsThrow(perHaplotypeColumns, __PRETTY_FUNCTION__);
			const auto perHaplotypePositions = getPositionsMultipleTargets(firstTable.header_.columnNames_, perHaplotypeColumns);
			summaryPerHaplotype = std::make_unique<table>(getTargetsAtPositions(firstTable.header_.columnNames_, perHaplotypePositions));
			//unique rows in the order they are first seen
			std::unordered_set<std::string> perHaplotypeSeen;
			genomicLocCounter = std::make_unique<table::GroupColumnsCounter>(firstTable.header_, VecStr{ "chrom", "0based_start", "0based_end","target", "length", "strand"});
			proteinLocCounter = std::make_unique<table::GroupColumnsCounter>(firstTable.header_, VecStr{"transcript", "transcript_1based_start", "transcript_1based_end", "target",  "transcript_length"});
			auto addSummaryRow = [&](const VecStr & row) {
				auto perHaplotypeRow = getTargetsAtPositions(row, perHaplotypePositions);
				if (perHaplotypeSeen.emplace(njh::conToStr(perHaplotypeRow, "\t")).second) {
					summaryPerHaplotype->addRow(perHaplotypeRow);
				}
				genomicLocCounter->addRow(row);
				proteinLocCounter->addRow(row);
			};
			OutputStream out(njh::files::make_path(reportsSummaryDir, "allSummaryTables.tab.txt.gz"));
			out << njh::conToStr(firstTable.header_.columnNames_, "\t") << '\n'; {
				VecStr firstTableRow;
				while (firstTable.getNextRow(firstTableRow)) {
					out << njh::conToStr(firstTableRow, "\t") << '\n';
					addSummaryRow(firstTableRow);
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
						addSummaryRow(currentRow);
					}
				}
			}
		}
	}
	//combine unique seqs summaries
	{
		fullWatch.startNewLap("gather unique seqs");
		//unique seqs fnps
		std::vector<bfs::path> uniqueSeqsFnps;
		std::vector<bfs::path> uniqueSeqs_labIsolateNamesFnps;
		// std::vector<bfs::path> uniqueSeqs_metaFnps; don't need to add the meta data cause this will be within the allSummary table
		for (const auto& tar: targetNamesVec) {
			auto uniqueSeqsFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/variantCalling/uniqueSeqs.fasta.gz");
			auto uniqueSeqs_labIsolateNamesFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/variantCalling/uniqueSeqs_labIsolateNames.tab.txt.gz");
			// auto uniqueSeqs_metaFnp = njh::pasteAsStr(setUp.pars_.directoryName_, "/", tar, "/variantCalling/uniqueSeqs_meta.tab.txt.gz");
			if (bfs::exists(uniqueSeqsFnp) && 0 != njh::files::bfs::file_size(uniqueSeqsFnp)) {
				uniqueSeqsFnps.emplace_back(uniqueSeqsFnp);
			}
			if (bfs::exists(uniqueSeqs_labIsolateNamesFnp) && 0 != njh::files::bfs::file_size(uniqueSeqs_labIsolateNamesFnp)) {
				uniqueSeqs_labIsolateNamesFnps.emplace_back(uniqueSeqs_labIsolateNamesFnp);
			}
			// if (bfs::exists(uniqueSeqs_metaFnp) && 0 != njh::files::bfs::file_size(uniqueSeqs_metaFnp)) {
			// 	uniqueSeqs_metaFnps.emplace_back(uniqueSeqs_metaFnp);
			// }
		}
		if (!uniqueSeqsFnps.empty()) {
			concatenateFiles(uniqueSeqsFnps, OutOptions(njh::files::make_path(reportsSummaryDir, "all_uniqueSeqs.fasta.gz")));
		}
		if (!uniqueSeqs_labIsolateNamesFnps.empty()) {
			njh::files::bfs::path firstFileFnp = uniqueSeqs_labIsolateNamesFnps.front();
			TableReader firstTable(TableIOOpts(InOptions(firstFileFnp), "\t", true));
			OutputStream out(njh::files::make_path(reportsSummaryDir, "all_uniqueSeqs_labIsolateNames.tab.txt.gz"));
			out << njh::conToStr(firstTable.header_.columnNames_, "\t") << '\n'; {
				VecStr firstTableRow;
				while (firstTable.getNextRow(firstTableRow)) {
					out << njh::conToStr(firstTableRow, "\t") << '\n';
				}
			}
			for (const auto& file: uniqueSeqs_labIsolateNamesFnps) {
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
		if(nullptr == summaryPerHaplotype) {
			//no summary tables were combined, reading the combined table throws the missing file error
			table allSummaryTable(TableIOOpts::genTabFileIn(njh::files::make_path(reportsSummaryDir, "allSummaryTables.tab.txt.gz")));
		}
		{
			//summary of per haplotype
			summaryPerHaplotype->outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(reportsSummaryDir, "allSummaryPerMicrohaplotype.tab.txt.gz")));
		}
		{
			auto genomicLocCounts = genomicLocCounter->genTable();
			VecStr isMax;
			std::unordered_map<std::string, uint32_t> maxCounts;
			for(const auto & row : genomicLocCounts) {
				auto n = njh::StrToNumConverter::stoToNum<uint32_t>(row[genomicLocCounts.getColPos("n")]);
				if(n > maxCounts[row[genomicLocCounts.getColPos("target")]] && "NA" != row[genomicLocCounts.getColPos("chrom")]) {
					maxCounts[row[genomicLocCounts.getColPos("target")]] = n;
				}
			}
			for(const auto & row : genomicLocCounts) {
				auto n = njh::StrToNumConverter::stoToNum<uint32_t>(row[genomicLocCounts.getColPos("n")]);
				if(n == maxCounts[row[genomicLocCounts.getColPos("target")]] && "NA" != row[genomicLocCounts.getColPos("chrom")]) {
					isMax.emplace_back("true");
					//if more than two maxes, then the first max will be added but not the others
					if (njh::notIn(row[genomicLocCounts.getColPos("target")], genomicLocs)) {
						genomicLocs[row[genomicLocCounts.getColPos("target")]] = std::make_shared<Bed6RecordCore>(njh::conToStr(row, "\t"));
					}
				} else {
					isMax.emplace_back("false");
				}
			}
			genomicLocCounts.addColumn(isMax, "isMaxCount");
			genomicLocCounts.naturlSortTable("target", false);
			genomicLocCounts.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(reportsInfoDir, "genomicLocPerTargetsCounts.tab.txt.gz")));
		}
		{
			auto proteinLocCounts = proteinLocCounter->genTable();
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
	if(!genomicLocs.empty() && !locs.genomicLocs.empty()) {
		fullWatch.startNewLap("add what known locs are intersecting with the targets");
		std::unordered_map<std::string, VecStr> coveredBy;
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
				if ((requireFullKnownCoverage && knownLoc.overlaps(*loc, knownLoc.length())) || (!requireFullKnownCoverage && knownLoc.overlaps(*loc, 1))) {
					intersected.emplace_back(knownLoc.name_);
					coveredBy[knownLoc.name_].emplace_back(loc->name_);
				}
			}
			if (!intersected.empty()) {
				intersectionWithKnownLocsOut << loc->toDelimStrWithExtra() << "\t" << njh::conToStr(intersected, ",") << std::endl;
			}
		}
		OutputStream knownLocCoverageOut(njh::files::make_path(reportsInfoDir, "knownLocCoverage.tsv"));
		knownLocCoverageOut << "knownLoc\tcovered_by" << std::endl;
		for (const auto& knownLoc: locs.genomicLocs) {
			if (njh::in(knownLoc.name_, coveredBy)) {
				knownLocCoverageOut << knownLoc.name_ << "\t" << njh::conToStr(coveredBy[knownLoc.name_], ",") << std::endl;
			} else {
				knownLocCoverageOut << knownLoc.name_ << "\t" << "not_covered" << std::endl;
			}
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