//
// The haplotypes variant calling is run on per target, built from final results the same way by variantCallOnSeqAndProtein and the population clustering programs
//

#include "VariantCallInput.hpp"
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>

namespace njhseq {

void VariantCallInput::addHap(const std::string & target, const std::string & sample, const std::string & readCnt,
		const std::string & popHapId, const std::string & seq, const VecStr & row) {
	if (!isDoubleStr(readCnt)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " for row: " << "\n";
		ss << njh::conToStr(row, "\t") << "\n";
		ss << "Within Sample Read Count column: " << readCnt << ", does not look like a number" << "\n";
		if (allWhiteSpaceStr(readCnt)) {
			ss << "value can't be empty, has to be a number" << "\n";
		}
		throw std::runtime_error{ss.str()};
	}
	targets_.emplace(target);
	samples_.emplace(sample);
	const auto cnt = njh::StrToNumConverter::stoToNum<double>(readCnt);
	popHapIdBySeq_[target][seq] = popHapId;
	seqInfo clus(njh::pasteAsStr(sample, "__", popHapId), seq);
	clus.cnt_ = cnt;
	bool add = true;
	for (auto & current : seqs_[target][sample]) {
		if (current.seq_ == clus.seq_) {
			add = false;
			current.cnt_ += cnt;
			break;
		}
	}
	if (add) {
		seqs_[target][sample].emplace_back(clus);
	}
	readVec::getMaxLength(clus, maxLen_);
}

void VariantCallInput::addTargetRows(const table & tab, const std::string & target, const VariantCallInputColumns & cols) {
	tab.checkForColumnsThrow(VecStr{cols.sample_, cols.readCnt_, cols.popHapId_, cols.popHapSeq_, cols.target_}, __PRETTY_FUNCTION__);
	const auto targetCol = tab.getColPos(cols.target_);
	const auto sampleCol = tab.getColPos(cols.sample_);
	const auto readCntCol = tab.getColPos(cols.readCnt_);
	const auto popHapIdCol = tab.getColPos(cols.popHapId_);
	const auto popHapSeqCol = tab.getColPos(cols.popHapSeq_);
	for (const auto & row : tab.content_) {
		if (row[targetCol] == target) {
			addHap(target, row[sampleCol], row[readCntCol], row[popHapIdCol], row[popHapSeqCol], row);
		}
	}
}

void VariantCallInput::writeInputSeqs(const std::string & target, const bfs::path & fnp, const MultipleGroupMetaData * meta) const {
	SeqOutput writer(SeqIOOptions::genFastaOutGz(fnp));
	writer.openOut();
	auto targetSeqs = seqs_.find(target);
	if (seqs_.end() == targetSeqs) {
		return;
	}
	const auto & popHapIds = popHapIdBySeq_.at(target);
	for (const auto & sampSeqs : targetSeqs->second) {
		for (auto seq : sampSeqs.second) {
			MetaDataInName seqMeta;
			if (nullptr != meta) {
				seqMeta = meta->getMetaForSample(sampSeqs.first);
			}
			seqMeta.addMeta("sample", sampSeqs.first, true);
			seqMeta.addMeta("readCount", std::max(1.0, std::round(seq.cnt_)), true);
			seqMeta.addMeta("originalIdentifier", popHapIds.at(seq.seq_));
			seqMeta.resetMetaInName(seq.name_);
			writer.write(seq);
		}
	}
}

std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> readVariantCallLocations(const bfs::path & bedFnp, bool changePeriodToDash) {
	std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> ret;
	for (auto & loc : getBeds(bedFnp)) {
		if (changePeriodToDash) {
			loc->name_ = njh::replaceString(loc->name_, ".", "-");
		}
		if (njh::in(loc->name_, ret)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have location: " << loc->name_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		ret[loc->name_] = loc;
	}
	return ret;
}

void FinalHapsVariantCallPars::setOptions(seqSetUp & setUp, CollapseAndCallVariantsPars & varCallPars) {
	setUp.setOption(genomicLocationsFnp, "--genomicLocations", "A bed file with the genomic location of targets to align to when calling variants, the name is the target, targets not in it aren't variant called", false, "Variant Calling");
	setUp.setOption(genomicLocationsChangePeriodToDash, "--genomicLocationsChangePeriodToDash", "When supplying --genomicLocations, change periods to dashes in the names", false, "Variant Calling");
	varCallPars.variantCallerRunPars.totalReadDepthCutOff = 10;
	setUp.setOption(varCallPars.variantCallerRunPars.totalReadDepthCutOff, "--variantTotalReadDepthCutOff", "Don't report variants that have a summed total read depth less than this across all samples", false, "Variant Calling");
	setUp.setOption(varCallPars.exportLabIsolateSeqs, "--exportLabIsolateSeqs", "Export the haplotypes of samples with the meta site==LabIsolate when calling variants", false, "Variant Calling");
	setUp.setOption(varCallPars.metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta fields to calculate population differences on when calling variants", false, "Variant Calling");
}

TranslatorByAlignment::TranslatorByAlignmentResult callVariantsOnFinalHaps(const table & finalResults,
		const std::string & target,
		CollapseAndCallVariantsPars varCallPars,
		const std::shared_ptr<Bed6RecordCore> & location,
		const bfs::path & metaFnp,
		const bfs::path & inputSeqsFnp) {
	std::unique_ptr<MultipleGroupMetaData> meta;
	if (!metaFnp.empty() && bfs::exists(metaFnp)) {
		meta = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
	VariantCallInput input;
	input.addTargetRows(finalResults, target);
	input.writeInputSeqs(target, inputSeqsFnp, meta.get());
	auto inputSeqs = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaInGz(inputSeqsFnp));
	varCallPars.identifier = target;
	if (nullptr != location) {
		varCallPars.refSeqRegion = GenomicRegion(*location);
	}
	return collapseAndCallVariants(varCallPars, inputSeqs);
}

}  // namespace njhseq
