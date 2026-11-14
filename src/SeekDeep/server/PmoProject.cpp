//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2026 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "PmoProject.hpp"
#include <njhseq/objects/BioDataObject/pmo.h>

namespace njhseq {

namespace {
//the width needed so names made from positions sort in order
std::string paddedPos(size_t pos, size_t total) {
	auto ret = estd::to_string(pos);
	auto width = estd::to_string(total > 0 ? total - 1 : 0).size();
	return std::string(width > ret.size() ? width - ret.size() : 0, '0') + ret;
}
}  // namespace

PmoProject::PmoProject(const bfs::path & pmoFnp) {
	auto pmoObj = pmo::readPmo(pmoFnp);
	header_ = njh::json::parse(pmoObj.pmo_header_.to_json().dump());

	for (const auto & library : pmoObj.library_sample_info_) {
		librarySamples_.emplace_back(library.library_sample_name_);
		const auto & specimen = pmoObj.specimen_info_.at(library.specimen_id_);
		specimenNames_[library.library_sample_name_] = specimen.specimen_name_;
		auto specimenJson = specimen.to_json();
		for (const auto & field : specimenJson.items()) {
			if ("specimen_name" != field.key() && field.value().is_primitive() && !field.value().is_null()) {
				auto value = field.value().is_string() ? field.value().get<std::string>() : field.value().dump();
				//values become part of the sub group urls
				for (auto & letter : value) {
					if ('/' == letter || '\\' == letter || '#' == letter || '?' == letter || '%' == letter) {
						letter = '-';
					}
				}
				meta_[library.library_sample_name_][field.key()] = value;
			}
		}
	}
	njh::sort(librarySamples_);

	//meta fields that split samples up, leaving out fields with one value or a different value for every sample
	std::map<std::string, std::set<std::string>> levels;
	for (const auto & sampleMeta : meta_) {
		for (const auto & field : sampleMeta.second) {
			levels[field.first].emplace(field.second);
		}
	}
	for (const auto & level : levels) {
		std::set<std::string> values = level.second;
		for (const auto & sample : librarySamples_) {
			values.emplace(metaValue(sample, level.first));
		}
		if (values.size() > 1 && values.size() < librarySamples_.size()) {
			metaGroups_.emplace_back(level.first);
		}
	}

	//a target can be in representative_microhaplotypes several times and mhap_id is the position within that entry
	std::vector<std::string> targetNameByMhapsId;
	std::vector<std::vector<uint32_t>> hapIdsByMhapsId;
	std::map<std::string, std::unordered_map<std::string, uint32_t>> hapIdBySeq;
	for (const auto & representatives : pmoObj.representative_microhaplotypes_.targets_) {
		const auto & info = pmoObj.target_info_.at(representatives.target_id_);
		targetNameByMhapsId.emplace_back(info.target_name_);
		hapIdsByMhapsId.emplace_back();
		auto & target = targets_[info.target_name_];
		target.name_ = info.target_name_;
		target.forwardPrimer_ = info.forward_primer_.seq_;
		target.reversePrimer_ = info.reverse_primer_.seq_;
		target.fakePrimers_ = info.target_attributes_.has_value() && njh::in(std::string(pmo::PMO_FAKE_PRIMERS_ATTRIBUTE), info.target_attributes_.value());
		auto & targetHapIds = hapIdBySeq[info.target_name_];
		for (const auto & mhap : representatives.microhaplotypes_) {
			auto search = targetHapIds.find(mhap.seq_);
			if (targetHapIds.end() == search) {
				search = targetHapIds.emplace(mhap.seq_, target.hapSeqs_.size()).first;
				target.hapSeqs_.emplace_back(mhap.seq_);
				target.hapNames_.emplace_back(mhap.microhaplotype_name_.has_value() ? mhap.microhaplotype_name_.value() : "");
				auto & aminoAcids = target.hapAminoAcids_.emplace_back();
				if (mhap.associated_protein_variants_.has_value()) {
					for (const auto & variant : mhap.associated_protein_variants_.value()) {
						const auto & loc = variant.protein_location_;
						aminoAcids.emplace_back(AminoAcidCall{loc.chrom_, variant.gene_name_.value_or(""), variant.alternative_gene_name_.value_or(""), loc.start_,
								loc.ref_seq_.value_or(""), loc.alt_seq_.value_or("")});
						hasGeneNames_ = hasGeneNames_ || variant.alternative_gene_name_.has_value();
					}
				}
				hasAminoAcidCalls_ = hasAminoAcidCalls_ || !aminoAcids.empty();
			}
			hapIdsByMhapsId.back().emplace_back(search->second);
		}
	}
	for (auto & target : targets_) {
		for (const auto hapPos : iter::range(target.second.hapNames_.size())) {
			if (target.second.hapNames_[hapPos].empty()) {
				target.second.hapNames_[hapPos] = target.first + "." + paddedPos(hapPos, target.second.hapNames_.size());
			}
		}
	}

	std::map<std::string, uint32_t> repeatedResults;
	for (const auto & run : pmoObj.detected_microhaplotypes_) {
		for (const auto & sampleResults : run.library_samples_) {
			const auto & sample = pmoObj.library_sample_info_.at(sampleResults.library_sample_id_).library_sample_name_;
			for (const auto & targetResults : sampleResults.target_results_) {
				auto & target = targets_.at(targetNameByMhapsId.at(targetResults.mhaps_target_id_));
				if (njh::in(sample, target.detected_)) {
					++repeatedResults[sample];
					continue;
				}
				std::map<uint32_t, uint32_t> readsByHapId;
				for (const auto & mhap : targetResults.mhaps_) {
					readsByHapId[hapIdsByMhapsId.at(targetResults.mhaps_target_id_).at(mhap.mhap_id_)] += mhap.reads_;
				}
				auto & detected = target.detected_[sample];
				for (const auto & hapReads : readsByHapId) {
					detected.emplace_back(DetectedHap{hapReads.first, hapReads.second});
				}
				std::stable_sort(detected.begin(), detected.end(), [](const DetectedHap & first, const DetectedHap & second) {
					return first.reads_ > second.reads_;
				});
			}
		}
	}
	for (const auto & repeated : repeatedResults) {
		std::cerr << __PRETTY_FUNCTION__ << ", warning, " << repeated.first << " has more than one set of results for " << repeated.second << " targets, only using the first of each" << std::endl;
	}
}

VecStr PmoProject::targetNames() const {
	return njh::getVecOfMapKeys(targets_);
}

bool PmoProject::hasTarget(const std::string & target) const {
	return njh::in(target, targets_);
}

VecStr PmoProject::sampleNames(const std::string & target) const {
	return njh::getVecOfMapKeys(targets_.at(target).detected_);
}

bool PmoProject::hasSample(const std::string & target, const std::string & sample) const {
	return hasTarget(target) && njh::in(sample, targets_.at(target).detected_);
}

std::string PmoProject::metaValue(const std::string & sample, const std::string & group) const {
	auto sampleMeta = meta_.find(sample);
	if (meta_.end() != sampleMeta) {
		auto value = sampleMeta->second.find(group);
		if (sampleMeta->second.end() != value) {
			return value->second;
		}
	}
	return "NA";
}

VecStr PmoProject::groupNames(const std::string & target) const {
	VecStr ret;
	for (const auto & group : metaGroups_) {
		if (subGroupNames(target, group).size() > 1) {
			ret.emplace_back(group);
		}
	}
	return ret;
}

bool PmoProject::hasGroup(const std::string & target, const std::string & group) const {
	return hasTarget(target) && njh::in(group, metaGroups_);
}

VecStr PmoProject::subGroupNames(const std::string & target, const std::string & group) const {
	std::set<std::string> ret;
	for (const auto & sample : targets_.at(target).detected_) {
		ret.emplace(metaValue(sample.first, group));
	}
	return VecStr(ret.begin(), ret.end());
}

bool PmoProject::hasSubGroup(const std::string & target, const std::string & group, const std::string & subGroup) const {
	return hasGroup(target, group) && njh::in(subGroup, subGroupNames(target, group));
}

VecStr PmoProject::subGroupSampleNames(const std::string & target, const std::string & group, const std::string & subGroup) const {
	VecStr ret;
	for (const auto & sample : targets_.at(target).detected_) {
		if (subGroup == metaValue(sample.first, group)) {
			ret.emplace_back(sample.first);
		}
	}
	return ret;
}

PmoProject::TargetTables PmoProject::genTables(const Target & target, const VecStr & samples) const {
	struct HapStats {
		double sumOfFracs_{0};
		uint64_t reads_{0};
		VecStr samples_;
		VecStr clusterNames_;
	};
	//the fraction of each haplotype in each sample
	std::map<uint32_t, HapStats> hapStats;
	std::map<std::string, uint64_t> sampleReads;
	std::vector<uint32_t> cois;
	uint64_t totalReads = 0;
	uint64_t totalClusters = 0;
	for (const auto & sample : samples) {
		const auto & detected = target.detected_.at(sample);
		uint64_t reads = 0;
		for (const auto & hap : detected) {
			reads += hap.reads_;
		}
		sampleReads[sample] = reads;
		totalReads += reads;
		totalClusters += detected.size();
		cois.emplace_back(detected.size());
		for (const auto clusPos : iter::range(detected.size())) {
			auto & stats = hapStats[detected[clusPos].hapId_];
			stats.sumOfFracs_ += 0 == reads ? 0 : detected[clusPos].reads_ / static_cast<double>(reads);
			stats.reads_ += detected[clusPos].reads_;
			stats.samples_.emplace_back(sample);
			stats.clusterNames_.emplace_back(sample + "." + estd::to_string(clusPos));
		}
	}
	const double sampleCnt = samples.size();
	auto popInfoCols = [&target, &hapStats, &samples, &sampleCnt, &totalReads](uint32_t hapId) {
		const auto & stats = hapStats.at(hapId);
		const double hapSampleCnt = stats.samples_.size();
		return toVecStr(target.hapNames_.at(hapId), samples.size(), "", stats.sumOfFracs_ / sampleCnt, stats.sumOfFracs_,
				stats.sumOfFracs_ / hapSampleCnt, 0 == totalReads ? 0.0 : stats.reads_ / static_cast<double>(totalReads),
				stats.samples_.size(), hapSampleCnt / sampleCnt, stats.reads_, stats.samples_.size(), target.hapSeqs_.at(hapId));
	};

	TargetTables ret;

	//sample table
	VecStr sampCols{"s_Sample", "p_name", "h_popUID", "p_TotalPopulationSampCnt", "h_AATyped", "h_PopFrac", "h_SumOfAllFracs",
		"h_AvgFracFoundAt", "h_ReadFrac", "h_SampCnt", "h_SampFrac", "h_ReadCnt", "h_ClusterCnt", "h_Consensus",
		"s_Name", "s_ReadCntTotUsed", "s_FinalClusterCnt", "s_COI"};
	for (const auto & group : metaGroups_) {
		sampCols.emplace_back("s_" + group);
	}
	addOtherVec(sampCols, VecStr{"c_clusterID", "c_name", "c_AveragedFrac", "c_ReadFrac", "c_ReadCnt", "c_RepCnt"});
	ret.sampInfo_ = table(sampCols);
	for (const auto & sample : samples) {
		const auto & detected = target.detected_.at(sample);
		const auto reads = sampleReads.at(sample);
		for (const auto clusPos : iter::range(detected.size())) {
			const double frac = 0 == reads ? 0 : detected[clusPos].reads_ / static_cast<double>(reads);
			auto row = toVecStr(sample, target.name_);
			addOtherVec(row, popInfoCols(detected[clusPos].hapId_));
			addOtherVec(row, toVecStr(sample, reads, detected.size(), detected.size()));
			for (const auto & group : metaGroups_) {
				row.emplace_back(metaValue(sample, group));
			}
			addOtherVec(row, toVecStr(clusPos, njh::pasteAsStr(sample, ".", clusPos, "_f", frac), frac, frac, detected[clusPos].reads_, 1));
			ret.sampInfo_.content_.emplace_back(row);
		}
	}

	//population table
	double coiMean = 0;
	double coiMedian = 0;
	uint32_t coiMin = 0;
	uint32_t coiMax = 0;
	if (!cois.empty()) {
		std::sort(cois.begin(), cois.end());
		coiMean = std::accumulate(cois.begin(), cois.end(), 0.0) / cois.size();
		coiMedian = 0 == cois.size() % 2 ? (cois[cois.size() / 2 - 1] + cois[cois.size() / 2]) / 2.0 : cois[cois.size() / 2];
		coiMin = cois.front();
		coiMax = cois.back();
	}
	//he = [n/(n - 1)][1 - sum(pi^2)]
	double sumOfSquares = 0;
	for (const auto & stats : hapStats) {
		sumOfSquares += std::pow(stats.second.samples_.size() / static_cast<double>(totalClusters), 2);
	}
	const double he = totalClusters > 1 ? (totalClusters / (totalClusters - 1.0)) * (1 - sumOfSquares) : 0;

	VecStr popCols{"h_popUID", "p_Name", "p_TotalInputReadCnt", "p_TotalInputClusterCnt", "p_TotalPopulationSampCnt", "p_TotalHaplotypes",
		"p_ExpectedHeterozygosity", "p_meanCoi", "p_medianCoi", "p_minCoi", "p_maxCoi"};
	for (const auto & group : metaGroups_) {
		popCols.emplace_back("g_" + group + "_subGroupCounts");
		popCols.emplace_back("g_" + group + "_subGroupFracs");
	}
	addOtherVec(popCols, VecStr{"h_PopFrac", "h_AATyped", "h_SumOfAllFracs", "h_AvgFracFoundAt", "h_ReadFrac", "h_SampCnt", "h_SampFrac", "h_ReadCnt", "h_ClusterCnt", "h_clusterNames", "h_Consensus"});
	ret.popInfo_ = table(popCols);
	VecStr hapNames;
	for (const auto & stats : hapStats) {
		const auto hapId = stats.first;
		const auto & hap = stats.second;
		const double hapSampleCnt = hap.samples_.size();
		hapNames.emplace_back(target.hapNames_.at(hapId));
		auto row = toVecStr(target.hapNames_.at(hapId), target.name_, totalReads, totalClusters, samples.size(), hapStats.size(), he, coiMean, coiMedian, coiMin, coiMax);
		for (const auto & group : metaGroups_) {
			std::map<std::string, uint32_t> counts;
			for (const auto & sample : hap.samples_) {
				++counts[metaValue(sample, group)];
			}
			std::string countsStr;
			std::string fracsStr;
			for (const auto & count : counts) {
				countsStr += njh::pasteAsStr(count.first, ":", count.second, ";");
				fracsStr += njh::pasteAsStr(count.first, ":", count.second / hapSampleCnt, ";");
			}
			addOtherVec(row, VecStr{countsStr, fracsStr});
		}
		addOtherVec(row, toVecStr(hap.sumOfFracs_ / sampleCnt, "", hap.sumOfFracs_, hap.sumOfFracs_ / hapSampleCnt,
				0 == totalReads ? 0.0 : hap.reads_ / static_cast<double>(totalReads), hap.samples_.size(), hapSampleCnt / sampleCnt,
				hap.reads_, hap.samples_.size(), njh::conToStr(hap.clusterNames_, ","), target.hapSeqs_.at(hapId)));
		ret.popInfo_.content_.emplace_back(row);
	}

	//haplotype by sample table
	njh::sort(hapNames);
	VecStr hapIdCols{"#PopUID"};
	addOtherVec(hapIdCols, samples);
	ret.hapIdTab_ = table(hapIdCols);
	std::unordered_map<std::string, uint32_t> hapRowPos;
	for (const auto & hapName : hapNames) {
		hapRowPos[hapName] = ret.hapIdTab_.content_.size();
		VecStr row(hapIdCols.size(), "0");
		row[0] = hapName;
		ret.hapIdTab_.content_.emplace_back(row);
	}
	for (const auto samplePos : iter::range(samples.size())) {
		const auto reads = sampleReads.at(samples[samplePos]);
		for (const auto & hap : target.detected_.at(samples[samplePos])) {
			ret.hapIdTab_.content_[hapRowPos.at(target.hapNames_.at(hap.hapId_))][samplePos + 1] =
					estd::to_string(0 == reads ? 0.0 : hap.reads_ / static_cast<double>(reads));
		}
	}
	return ret;
}

const PmoProject::TargetTables & PmoProject::targetTables(const std::string & target) {
	std::lock_guard<std::mutex> lock(tablesMut_);
	auto search = tablesByTarget_.find(target);
	if (tablesByTarget_.end() != search) {
		return search->second;
	}
	return tablesByTarget_.emplace(target, genTables(targets_.at(target), sampleNames(target))).first->second;
}

void PmoProject::computeGroup(const std::string & target, const std::string & group) {
	auto & subGroupTabs = tablesBySubGroup_[target][group];
	std::map<std::string, std::set<std::string>> popUids;
	for (const auto & subGroup : subGroupNames(target, group)) {
		auto & tabs = subGroupTabs[subGroup];
		tabs = genTables(targets_.at(target), subGroupSampleNames(target, group, subGroup));
		auto uids = tabs.popInfo_.getColumn("h_popUID");
		popUids[subGroup] = std::set<std::string>(uids.begin(), uids.end());
	}
	VecStr groupInfoCols{"g_GroupName", "p_TotalInputReadCnt", "p_TotalInputClusterCnt", "p_TotalPopulationSampCnt",
		"p_TotalHaplotypes", "p_TotalUniqueHaplotypes", "p_ExpectedHeterozygosity", "p_meanCoi", "p_medianCoi", "p_minCoi", "p_maxCoi", "g_hapsFoundOnlyInThisGroup"};
	table groupTab(groupInfoCols);
	for (auto & subGroup : subGroupTabs) {
		VecStr unique;
		for (const auto & uid : popUids[subGroup.first]) {
			bool found = false;
			for (const auto & other : popUids) {
				if (other.first != subGroup.first && njh::in(uid, other.second)) {
					found = true;
					break;
				}
			}
			if (!found) {
				unique.emplace_back(uid);
			}
		}
		const auto groupName = group + ":" + subGroup.first;
		auto & popInfo = subGroup.second.popInfo_;
		popInfo.addColumn({groupName}, "g_GroupName");
		popInfo.addColumn({njh::conToStr(unique)}, "g_hapsFoundOnlyInThisGroup");
		popInfo.addColumn({estd::to_string(unique.size())}, "p_TotalUniqueHaplotypes");
		groupTab.rbind(popInfo.getColumns(groupInfoCols), false);

		auto & sampInfo = subGroup.second.sampInfo_;
		auto sampCols = sampInfo.columnNames_;
		removeElement(sampCols, std::string("s_Sample"));
		VecStr outCols{"s_Sample", "g_GroupName"};
		addOtherVec(outCols, sampCols);
		sampInfo.addColumn({groupName}, "g_GroupName");
		sampInfo = sampInfo.getColumns(outCols);
	}
	groupTab = groupTab.getUniqueRows();
	groupTab.sortTable("g_GroupName", false);
	groupInfoByTarget_[target][group] = groupTab;
}

const PmoProject::TargetTables & PmoProject::subGroupTables(const std::string & target, const std::string & group, const std::string & subGroup) {
	std::lock_guard<std::mutex> lock(tablesMut_);
	if (!njh::in(group, groupInfoByTarget_[target])) {
		computeGroup(target, group);
	}
	return tablesBySubGroup_.at(target).at(group).at(subGroup);
}

const table & PmoProject::groupInfo(const std::string & target, const std::string & group) {
	std::lock_guard<std::mutex> lock(tablesMut_);
	if (!njh::in(group, groupInfoByTarget_[target])) {
		computeGroup(target, group);
	}
	return groupInfoByTarget_.at(target).at(group);
}

table PmoProject::targetsInfo() {
	table ret(VecStr{"target", "forward_primer", "reverse_primer", "fake_primers", "samples", "haplotypes", "reads", "expected_heterozygosity", "mean_coi", "max_coi"});
	for (const auto & target : targets_) {
		const auto & popInfo = targetTables(target.first).popInfo_;
		const bool hasData = !popInfo.content_.empty();
		ret.content_.emplace_back(toVecStr(target.first, target.second.forwardPrimer_, target.second.reversePrimer_, target.second.fakePrimers_ ? "yes" : "no", target.second.detected_.size(), popInfo.content_.size(),
				hasData ? popInfo.getColumn("p_TotalInputReadCnt").front() : "0",
				hasData ? popInfo.getColumn("p_ExpectedHeterozygosity").front() : "0",
				hasData ? popInfo.getColumn("p_meanCoi").front() : "0",
				hasData ? popInfo.getColumn("p_maxCoi").front() : "0"));
	}
	return ret;
}

table PmoProject::librarySamplesInfo() const {
	VecStr cols{"library_sample", "specimen", "targets_detected", "total_reads"};
	for (const auto & group : metaGroups_) {
		cols.emplace_back(group);
	}
	table ret(cols);
	for (const auto & sample : librarySamples_) {
		uint32_t targetsDetected = 0;
		uint64_t reads = 0;
		for (const auto & target : targets_) {
			auto detected = target.second.detected_.find(sample);
			if (target.second.detected_.end() != detected) {
				++targetsDetected;
				for (const auto & hap : detected->second) {
					reads += hap.reads_;
				}
			}
		}
		auto row = toVecStr(sample, specimenNames_.at(sample), targetsDetected, reads);
		for (const auto & group : metaGroups_) {
			row.emplace_back(metaValue(sample, group));
		}
		ret.content_.emplace_back(row);
	}
	return ret;
}

std::string PmoProject::popSeqsFasta(const std::string & target) {
	const auto & popInfo = targetTables(target).popInfo_;
	const auto uidCol = popInfo.getColPos("h_popUID");
	const auto fracCol = popInfo.getColPos("h_AvgFracFoundAt");
	const auto seqCol = popInfo.getColPos("h_Consensus");
	std::stringstream ret;
	for (const auto & row : popInfo.content_) {
		ret << ">" << row[uidCol] << "_f" << row[fracCol] << "\n" << row[seqCol] << "\n";
	}
	return ret.str();
}

std::string PmoProject::sampleSeqsFasta(const std::string & target, const std::string & sample) const {
	const auto & tar = targets_.at(target);
	const auto & detected = tar.detected_.at(sample);
	uint64_t reads = 0;
	for (const auto & hap : detected) {
		reads += hap.reads_;
	}
	std::stringstream ret;
	for (const auto clusPos : iter::range(detected.size())) {
		ret << ">" << sample << "." << clusPos << "_f" << (0 == reads ? 0.0 : detected[clusPos].reads_ / static_cast<double>(reads)) << "\n"
				<< tar.hapSeqs_.at(detected[clusPos].hapId_) << "\n";
	}
	return ret.str();
}

VecStr PmoProject::librarySampleNames() const {
	return librarySamples_;
}

bool PmoProject::hasLibrarySample(const std::string & sample) const {
	return njh::in(sample, specimenNames_);
}

table PmoProject::sampleAcrossTargets(const std::string & sample) {
	const VecStr cols{"s_Sample", "p_name", "h_popUID", "c_clusterID", "c_AveragedFrac", "c_ReadCnt", "s_ReadCntTotUsed", "s_COI",
		"h_SampCnt", "h_SampFrac", "h_PopFrac", "p_TotalPopulationSampCnt", "h_Consensus"};
	table ret(cols);
	for (const auto & target : targets_) {
		if (!njh::in(sample, target.second.detected_)) {
			continue;
		}
		const auto & sampInfo = targetTables(target.first).sampInfo_;
		const auto sampleCol = sampInfo.getColPos("s_Sample");
		std::vector<uint32_t> colPositions;
		for (const auto & col : cols) {
			colPositions.emplace_back(sampInfo.getColPos(col));
		}
		for (const auto & row : sampInfo.content_) {
			if (row[sampleCol] == sample) {
				VecStr outRow;
				for (const auto colPos : colPositions) {
					outRow.emplace_back(row[colPos]);
				}
				ret.content_.emplace_back(outRow);
			}
		}
	}
	return ret;
}

VecStr PmoProject::targetsNotDetected(const std::string & sample) const {
	VecStr ret;
	for (const auto & target : targets_) {
		if (!njh::in(sample, target.second.detected_)) {
			ret.emplace_back(target.first);
		}
	}
	return ret;
}

std::map<std::string, std::string> PmoProject::sampleMeta(const std::string & sample) const {
	std::map<std::string, std::string> ret;
	auto search = meta_.find(sample);
	if (meta_.end() != search) {
		ret = search->second;
	}
	ret["specimen_name"] = specimenNames_.at(sample);
	return ret;
}

std::string PmoProject::sampleAcrossTargetsFasta(const std::string & sample) const {
	std::stringstream ret;
	for (const auto & target : targets_) {
		auto detected = target.second.detected_.find(sample);
		if (target.second.detected_.end() == detected) {
			continue;
		}
		uint64_t reads = 0;
		for (const auto & hap : detected->second) {
			reads += hap.reads_;
		}
		for (const auto clusPos : iter::range(detected->second.size())) {
			ret << ">" << target.first << "." << clusPos << "_f" << (0 == reads ? 0.0 : detected->second[clusPos].reads_ / static_cast<double>(reads)) << "\n"
					<< target.second.hapSeqs_.at(detected->second[clusPos].hapId_) << "\n";
		}
	}
	return ret.str();
}

VecStr PmoProject::metaGroupNames() const {
	return metaGroups_;
}

VecStr PmoProject::metaSubGroupNames(const std::string & group) const {
	std::set<std::string> ret;
	for (const auto & sample : librarySamples_) {
		ret.emplace(metaValue(sample, group));
	}
	return VecStr(ret.begin(), ret.end());
}

bool PmoProject::hasMetaSubGroup(const std::string & group, const std::string & subGroup) const {
	return njh::in(group, metaGroups_) && njh::in(subGroup, metaSubGroupNames(group));
}

VecStr PmoProject::metaSubGroupSampleNames(const std::string & group, const std::string & subGroup) const {
	VecStr ret;
	for (const auto & sample : librarySamples_) {
		if (subGroup == metaValue(sample, group)) {
			ret.emplace_back(sample);
		}
	}
	return ret;
}

Json::Value PmoProject::genHaplotypeRainbow(const VecStr & samples, const std::set<std::string> & targets) const {
	const auto skipTarget = [&targets](const std::string & target) {
		return !targets.empty() && !njh::in(target, targets);
	};
	VecStr withData;
	VecStr withoutData;
	for (const auto & sample : samples) {
		bool detected = false;
		for (const auto & target : targets_) {
			if (!skipTarget(target.first) && njh::in(sample, target.second.detected_)) {
				detected = true;
				break;
			}
		}
		(detected ? withData : withoutData).emplace_back(sample);
	}
	std::unordered_map<std::string, uint32_t> samplePositions;
	for (const auto pos : iter::range(withData.size())) {
		samplePositions[withData[pos]] = pos;
	}

	Json::Value ret;
	ret["samples"] = njh::json::toJson(withData);
	ret["samplesWithoutData"] = njh::json::toJson(withoutData);
	auto & targetsOut = ret["targets"];
	targetsOut = Json::arrayValue;
	auto & calls = ret["calls"];
	calls = Json::arrayValue;
	uint32_t targetPos = 0;
	for (const auto & target : targets_) {
		if (skipTarget(target.first)) {
			continue;
		}
		std::map<uint32_t, uint32_t> sampleCounts;
		for (const auto & detected : target.second.detected_) {
			if (njh::in(detected.first, samplePositions)) {
				for (const auto & hap : detected.second) {
					++sampleCounts[hap.hapId_];
				}
			}
		}
		auto hapIds = njh::getVecOfMapKeys(sampleCounts);
		std::sort(hapIds.begin(), hapIds.end(), [&sampleCounts, &target](uint32_t first, uint32_t second) {
			if (sampleCounts.at(first) != sampleCounts.at(second)) {
				return sampleCounts.at(first) > sampleCounts.at(second);
			}
			return target.second.hapNames_.at(first) < target.second.hapNames_.at(second);
		});
		std::unordered_map<uint32_t, uint32_t> rankByHapId;
		Json::Value targetOut;
		targetOut["name"] = target.first;
		targetOut["haps"] = Json::arrayValue;
		targetOut["hapSampleCounts"] = Json::arrayValue;
		for (const auto rank : iter::range(hapIds.size())) {
			rankByHapId[hapIds[rank]] = rank;
			targetOut["haps"].append(target.second.hapNames_.at(hapIds[rank]));
			targetOut["hapSampleCounts"].append(sampleCounts.at(hapIds[rank]));
		}
		uint32_t typedSamples = 0;
		for (const auto & detected : target.second.detected_) {
			auto samplePos = samplePositions.find(detected.first);
			if (samplePositions.end() == samplePos) {
				continue;
			}
			++typedSamples;
			for (const auto & hap : detected.second) {
				calls.append(samplePos->second);
				calls.append(targetPos);
				calls.append(rankByHapId.at(hap.hapId_));
				calls.append(hap.reads_);
			}
		}
		targetOut["samples"] = typedSamples;
		targetsOut.append(targetOut);
		++targetPos;
	}

	ret["hasAminoAcidCalls"] = hasAminoAcidCalls_;
	ret["metaFields"] = njh::json::toJson(metaGroups_);
	auto & sampleMetaOut = ret["sampleMeta"];
	sampleMetaOut = Json::arrayValue;
	for (const auto & sample : withData) {
		Json::Value values = Json::arrayValue;
		for (const auto & group : metaGroups_) {
			values.append(metaValue(sample, group));
		}
		sampleMetaOut.append(values);
	}
	return ret;
}

const Json::Value & PmoProject::haplotypeRainbow(const std::string & group, const std::string & subGroup) {
	std::lock_guard<std::mutex> lock(tablesMut_);
	const auto key = std::make_pair(group, subGroup);
	auto search = rainbows_.find(key);
	if (rainbows_.end() != search) {
		return search->second;
	}
	return rainbows_.emplace(key, genHaplotypeRainbow(group.empty() ? librarySamples_ : metaSubGroupSampleNames(group, subGroup))).first->second;
}

VecStr PmoProject::knownLibrarySamples(const VecStr & samples) const {
	std::unordered_set<std::string> wanted(samples.begin(), samples.end());
	VecStr ret;
	for (const auto & sample : librarySamples_) {
		if (njh::in(sample, wanted)) {
			ret.emplace_back(sample);
		}
	}
	return ret;
}

std::set<std::string> PmoProject::knownTargets(const VecStr & targets) const {
	std::set<std::string> ret;
	for (const auto & target : targets) {
		if (njh::in(target, targets_)) {
			ret.emplace(target);
		}
	}
	return ret;
}

Json::Value PmoProject::haplotypeRainbowFor(const VecStr & samples, const VecStr & targets) const {
	return genHaplotypeRainbow(knownLibrarySamples(samples), knownTargets(targets));
}

bool PmoProject::hasAminoAcidCalls() const {
	return hasAminoAcidCalls_;
}

table PmoProject::genAminoAcidChanges(const VecStr & samples, const std::set<std::string> & targets) const {
	//a position is a transcript and 1-based amino acid position
	using Position = std::pair<std::string, uint32_t>;
	struct AlleleCounts {
		uint32_t haps_{0};
		uint32_t samples_{0};
		double fracSum_{0};
	};
	struct PositionCounts {
		std::string gene_;
		std::string geneName_;
		std::string refAA_;
		std::set<std::string> targets_;
		uint32_t samples_{0};
		uint32_t haps_{0};
		std::map<std::string, AlleleCounts> alleles_;
	};
	struct SampleAllele {
		uint64_t reads_{0};
		uint32_t haps_{0};
	};
	std::map<Position, PositionCounts> positions;
	//key1 = position, key2 = target, key3 = amino acid
	using SampleCalls = std::map<Position, std::map<std::string, std::map<std::string, SampleAllele>>>;
	auto genSampleCalls = [this, &positions, &targets](const std::string & sample) {
		SampleCalls ret;
		for (const auto & target : targets_) {
			if (!targets.empty() && !njh::in(target.first, targets)) {
				continue;
			}
			auto detected = target.second.detected_.find(sample);
			if (target.second.detected_.end() == detected) {
				continue;
			}
			for (const auto & hap : detected->second) {
				for (const auto & call : target.second.hapAminoAcids_.at(hap.hapId_)) {
					auto & allele = ret[Position{call.transcript_, call.pos_}][target.first][call.aa_];
					allele.reads_ += hap.reads_;
					++allele.haps_;
					auto & posCounts = positions[Position{call.transcript_, call.pos_}];
					posCounts.gene_ = call.gene_;
					posCounts.geneName_ = call.geneName_;
					posCounts.refAA_ = call.refAA_;
				}
			}
		}
		return ret;
	};
	//like combining the per target vcfs, a position uses the target typed in the most samples (first by name on ties)
	//and samples without that target use the target with the most reads for them
	std::map<Position, std::map<std::string, uint32_t>> samplesPerTarget;
	for (const auto & sample : samples) {
		for (const auto & position : genSampleCalls(sample)) {
			for (const auto & target : position.second) {
				++samplesPerTarget[position.first][target.first];
			}
		}
	}
	std::map<Position, std::string> mainTarget;
	for (const auto & position : samplesPerTarget) {
		uint32_t most = 0;
		for (const auto & target : position.second) {
			if (target.second > most) {
				most = target.second;
				mainTarget[position.first] = target.first;
			}
		}
	}
	for (const auto & sample : samples) {
		for (const auto & position : genSampleCalls(sample)) {
			auto bestSearch = position.second.find(mainTarget.at(position.first));
			if (position.second.end() == bestSearch) {
				uint64_t bestReads = 0;
				for (auto targetSearch = position.second.begin(); targetSearch != position.second.end(); ++targetSearch) {
					uint64_t reads = 0;
					for (const auto & allele : targetSearch->second) {
						reads += allele.second.reads_;
					}
					if (position.second.end() == bestSearch || reads > bestReads) {
						bestSearch = targetSearch;
						bestReads = reads;
					}
				}
			}
			const auto & best = bestSearch->second;
			uint64_t bestReads = 0;
			uint32_t bestHaps = 0;
			for (const auto & allele : best) {
				bestReads += allele.second.reads_;
				bestHaps += allele.second.haps_;
			}
			auto & posCounts = positions[position.first];
			++posCounts.samples_;
			posCounts.haps_ += bestHaps;
			posCounts.targets_.emplace(bestSearch->first);
			for (const auto & allele : best) {
				auto & counts = posCounts.alleles_[allele.first];
				counts.haps_ += allele.second.haps_;
				++counts.samples_;
				//haplotypes without reads count equally
				counts.fracSum_ += bestReads > 0 ? allele.second.reads_ / static_cast<double>(bestReads) : allele.second.haps_ / static_cast<double>(bestHaps);
			}
		}
	}

	table ret(VecStr{"Gene_ID", "Gene_Transcript_ID", "AA_Change", "reference_AA_pos", "reference_AA", "AA", "is_reference",
			"AlleleCount", "AlleleFrequency", "AlleleTotal", "SampleCount", "SamplePrevalence", "SampleTotal", "CoveredBy"});
	if (hasGeneNames_) {
		ret.columnNames_.insert(ret.columnNames_.begin() + 1, "Gene");
	}
	for (auto & position : positions) {
		auto & posCounts = position.second;
		if (0 == posCounts.samples_) {
			continue;
		}
		posCounts.alleles_[posCounts.refAA_];
		//the reference first then the most frequent
		auto aas = njh::getVecOfMapKeys(posCounts.alleles_);
		std::stable_sort(aas.begin(), aas.end(), [&posCounts](const std::string & first, const std::string & second) {
			if ((first == posCounts.refAA_) != (second == posCounts.refAA_)) {
				return first == posCounts.refAA_;
			}
			return posCounts.alleles_.at(first).fracSum_ > posCounts.alleles_.at(second).fracSum_;
		});
		for (const auto & aa : aas) {
			const auto & counts = posCounts.alleles_.at(aa);
			const bool isRef = aa == posCounts.refAA_;
			ret.content_.emplace_back(toVecStr(posCounts.gene_, position.first.first,
					njh::pasteAsStr(posCounts.refAA_, position.first.second, isRef ? "" : aa),
					position.first.second, posCounts.refAA_, aa, isRef ? "TRUE" : "FALSE",
					counts.haps_, counts.fracSum_ / posCounts.samples_, posCounts.haps_,
					counts.samples_, counts.samples_ / static_cast<double>(posCounts.samples_), posCounts.samples_,
					njh::conToStr(posCounts.targets_, "::")));
			if (hasGeneNames_) {
				ret.content_.back().insert(ret.content_.back().begin() + 1, posCounts.geneName_.empty() ? "NA" : posCounts.geneName_);
			}
		}
	}
	return ret;
}

const table & PmoProject::aminoAcidChanges(const std::string & group, const std::string & subGroup) {
	std::lock_guard<std::mutex> lock(tablesMut_);
	const auto key = std::make_pair(group, subGroup);
	auto search = aminoAcidChanges_.find(key);
	if (aminoAcidChanges_.end() != search) {
		return search->second;
	}
	return aminoAcidChanges_.emplace(key, genAminoAcidChanges(group.empty() ? librarySamples_ : metaSubGroupSampleNames(group, subGroup))).first->second;
}

table PmoProject::aminoAcidChangesFor(const VecStr & samples, const VecStr & targets) const {
	return genAminoAcidChanges(knownLibrarySamples(samples), knownTargets(targets));
}

}  // namespace njhseq
