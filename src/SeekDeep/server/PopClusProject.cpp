/*
 * PopClusProject.cpp
 *
 *  Created on: Sep 19, 2016
 *      Author: nick
 */

//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2019 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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


#include "PopClusProject.hpp"
#include <unordered_map>
#include <unistd.h>

namespace njhseq {

PopClusProject::PopClusProject(const Json::Value & configJson) :
		config_(configJson) {
	njh::json::MemberChecker checker(configJson);
	checker.failMemberCheckThrow( { "shortName", "projectName", "mainDir" },
			__PRETTY_FUNCTION__);
	auto coreJsonFnp = njh::files::make_path(configJson["mainDir"],
			"coreInfo.json");
	if (!bfs::exists(coreJsonFnp)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, project "
				<< configJson["projectName"] << "configuration "
				<< " main directory doesn't contain coreInfo.json file" << "\n";
		ss << coreJsonFnp << " doesn't exist" << "\n";
		throw std::runtime_error { ss.str() };
	}
	collection_ = std::make_unique<collapse::SampleCollapseCollection>(
			njh::json::parseFile(coreJsonFnp.string()));
	shortName_ = config_["shortName"].asString();
	projectName_ = config_["projectName"].asString();

	tabs_.popInfo_ = std::make_unique<TableCache>(TableIOOpts(InOptions(collection_->getPopInfoPath()), "\t", true));
	tabs_.sampInfo_ = std::make_unique<TableCache>(TableIOOpts(InOptions(collection_->getSampInfoPath()), "\t", true));
	tabs_.hapIdTab_ = std::make_unique<TableCache>(TableIOOpts(InOptions(collection_->getHapIdTabPath()), "\t", true));

	//set up group meta data
	if(nullptr != collection_->groupDataPaths_){
		for(const auto & group : collection_->groupDataPaths_->allGroupPaths_){
			//groups without files are computed when requested
			if(!bfs::exists(group.second.groupInfoFnp_)){
				continue;
			}
			topGroupTabs_[group.first] = std::make_unique<TableCache>(TableIOOpts(InOptions(group.second.groupInfoFnp_), "\t", true));
			for(const auto & subGroup : group.second.groupPaths_){
				subGroupTabs_[group.first][subGroup.first].popInfo_ = std::make_unique<TableCache>(TableIOOpts(InOptions(subGroup.second.popFileFnp_), "\t", true));
				subGroupTabs_[group.first][subGroup.first].sampInfo_ = std::make_unique<TableCache>(TableIOOpts(InOptions(subGroup.second.sampFileFnp_), "\t", true));
				subGroupTabs_[group.first][subGroup.first].hapIdTab_ = std::make_unique<TableCache>(TableIOOpts(InOptions(subGroup.second.hapIdTabFnp_), "\t", true));
			}
		}
	}

	sampleClustersInfoFnp_ = njh::files::make_path(collection_->masterOutputDir_, "sampleClustersInfo.tab.txt.gz");
	if(bfs::exists(sampleClustersInfoFnp_) && bfs::exists(njh::BlockGzipOutputStream::indexPath(sampleClustersInfoFnp_))){
		sampleClustersIndex_ = std::make_unique<njh::BlockGzipIndex>(sampleClustersInfoFnp_);
	}

	auto extractionProfileFnp = njh::files::make_path(configJson["mainDir"].asString(), "extractionInfo", "extractionProfile.tab.txt");
	auto extractionStatsFnp = njh::files::make_path(configJson["mainDir"].asString(), "extractionInfo", "extractionStats.tab.txt");

	if(bfs::exists(extractionProfileFnp)){
		extractionProfileTab_ = std::make_unique<TableCache>(TableIOOpts(InOptions(extractionProfileFnp), "\t", true));
	}

	if(bfs::exists(extractionStatsFnp)){
		extractionStatsTab_ = std::make_unique<TableCache>(TableIOOpts(InOptions(extractionStatsFnp), "\t", true));
	}

}


void PopClusProject::registerSeqFiles(SeqCache & cache) {
	auto getSeqFormat = [](const bfs::path & seqFnp){
		auto ret = njh::files::getExtension(seqFnp);
		if("gz" == ret){
			if (njh::endsWith(seqFnp.string(), ".fasta.gz")) {
				ret = ".fasta.gz";
			} else if (njh::endsWith(seqFnp.string(), ".fastq.gz")) {
				ret = ".fastq.gz";
			}
		}
		return ret;
	};
	auto popHapFile = collection_->getPopFinalHapsPath().string();
	cache.updateAddCache(shortName_,
			SeqIOOptions(popHapFile,
					SeqIOOptions::getInFormat(getSeqFormat(popHapFile)),
					true));
	//without sample directories samples are added by registerSampleSeqs when requested
	auto samplesDir = njh::files::make_path(collection_->masterOutputDir_, "samplesOutput");
	if(!bfs::exists(samplesDir) || bfs::is_empty(samplesDir)){
		return;
	}
	for(const auto & samp : collection_->popNames_.samples_){
		auto sampHapFile = collection_->getSampleFinalHapsPath(samp).string();
		if(bfs::exists(sampHapFile)){
			cache.updateAddCache(shortName_ + "_" + samp,
					SeqIOOptions(sampHapFile,
							SeqIOOptions::getInFormat(getSeqFormat(sampHapFile)),
							true));
		}
	}
}

void PopClusProject::registerSampleSeqs(const std::string & sample, SeqCache & cache) {
	auto uid = shortName_ + "_" + sample;
	if(cache.containsRecord(uid) || nullptr == sampleClustersIndex_ || !sampleClustersIndex_->has(sample)){
		return;
	}
	std::string line;
	InputStream headerIn{InOptions(sampleClustersInfoFnp_)};
	njh::files::crossPlatGetline(headerIn, line);
	auto header = tokenizeString(line, "\t", true);
	auto colPos = [this, &header](const std::string & colName) {
		auto search = std::find(header.begin(), header.end(), colName);
		if(header.end() == search){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, " << sampleClustersInfoFnp_ << " doesn't have column " << colName << "\n";
			throw std::runtime_error{ss.str()};
		}
		return static_cast<size_t>(std::distance(header.begin(), search));
	};
	const auto statusCol = colPos("c_status");
	const auto nameCol = colPos("c_name");
	const auto seqCol = colPos("c_seq");
	const auto qualCol = colPos("c_qual");

	std::stringstream finalSeqs;
	std::stringstream block(sampleClustersIndex_->readBlock(sample));
	while(njh::files::crossPlatGetline(block, line)){
		auto toks = tokenizeString(line, "\t", true);
		if("excluded" != toks[statusCol]){
			finalSeqs << "@" << toks[nameCol] << "\n" << toks[seqCol] << "\n+\n" << toks[qualCol] << "\n";
		}
	}
	//samples with only excluded clusters have no final haplotypes to serve
	if(finalSeqs.str().empty()){
		return;
	}
	//per process since the working directory defaults to the shared /tmp/
	auto seqDir = njh::files::make_path(cache.workingDir_, njh::pasteAsStr("sampleSeqs_", getpid()), shortName_);
	njh::files::makeDirP(njh::files::MkdirPar(seqDir.string()));
	auto seqFnp = njh::files::make_path(seqDir, sample + ".fastq");
	OutOptions seqOutOpts(seqFnp);
	seqOutOpts.overWriteFile_ = true;
	{
		OutputStream seqOut(seqOutOpts);
		seqOut << finalSeqs.str();
	}
	cache.updateAddCache(uid, SeqIOOptions::genFastqIn(seqFnp, true));
}

const table & PopClusProject::topGroupInfo(const std::string & group) {
	if(njh::in(group, topGroupTabs_)){
		return topGroupTabs_.at(group)->get();
	}
	return computedGroup(group).groupInfo_;
}

const table & PopClusProject::subGroupPopInfo(const std::string & group, const std::string & subGroup) {
	if(njh::in(group, topGroupTabs_)){
		return subGroupTabs_.at(group).at(subGroup).popInfo_->get();
	}
	return computedGroup(group).subGroups_.at(subGroup).popInfo_;
}

const table & PopClusProject::subGroupSampInfo(const std::string & group, const std::string & subGroup) {
	if(njh::in(group, topGroupTabs_)){
		return subGroupTabs_.at(group).at(subGroup).sampInfo_->get();
	}
	return computedGroup(group).subGroups_.at(subGroup).sampInfo_;
}

const table & PopClusProject::subGroupHapIdTab(const std::string & group, const std::string & subGroup) {
	if(njh::in(group, topGroupTabs_)){
		return subGroupTabs_.at(group).at(subGroup).hapIdTab_->get();
	}
	return computedGroup(group).subGroups_.at(subGroup).hapIdTab_;
}

VecStr PopClusProject::subGroupSampleNames(const std::string & group, const std::string & subGroup) {
	if(njh::in(group, topGroupTabs_)){
		return collection_->groupDataPaths_->allGroupPaths_.at(group).groupPaths_.at(subGroup).readInSampNames();
	}
	return computedGroup(group).subGroups_.at(subGroup).sampNames_;
}

PopClusProject::ComputedSubGroup PopClusProject::computeSubGroup(const std::set<std::string> & samples) {
	ComputedSubGroup ret;
	collection_->loadInPreviousPop(samples);
	ret.popInfo_ = collection_->genPopulationCollapseInfo();

	//population columns recalculated for just these samples
	const auto & popCollapsed = collection_->popCollapse_->collapsed_;
	std::unordered_map<std::string, VecStr> popInfoByUid;
	std::set<std::string> popUids;
	for(const auto & clus : popCollapsed.clusters_){
		popInfoByUid[clus.getStubName(false)] = clus.getPopInfoVec(
				popCollapsed.info_.totalReadCount_,
				popCollapsed.info_.numberOfClusters_,
				popCollapsed.info_.infos_.size());
		popUids.emplace(clus.getStubName(true));
	}
	collection_->popCollapse_ = nullptr;

	const auto & allSampInfo = tabs_.sampInfo_->get();
	const auto sampleCol = allSampInfo.getColPos("s_Sample");
	const auto popUidCol = allSampInfo.getColPos("h_popUID");
	std::vector<uint32_t> popInfoCols;
	for(const auto & colName : sampleCluster::getPopInfoHeaderVec()){
		popInfoCols.emplace_back(allSampInfo.getColPos(colName));
	}
	table sampInfo(allSampInfo.columnNames_);
	for(const auto & row : allSampInfo.content_){
		if(njh::in(row[sampleCol], samples)){
			auto subGroupRow = row;
			const auto & popInfo = popInfoByUid.at(row[popUidCol]);
			for(size_t pos = 0; pos < popInfoCols.size(); ++pos){
				subGroupRow[popInfoCols[pos]] = popInfo[pos];
			}
			sampInfo.content_.emplace_back(subGroupRow);
		}
	}

	//only as many replicate columns as these samples have
	uint32_t maxRunCount = 1;
	for(uint32_t run = 2; njh::in(njh::pasteAsStr("R", run, "_name"), allSampInfo.columnNames_); ++run){
		const auto runNameCol = allSampInfo.getColPos(njh::pasteAsStr("R", run, "_name"));
		if(std::any_of(sampInfo.content_.begin(), sampInfo.content_.end(), [&runNameCol](const VecStr & row){ return !row[runNameCol].empty(); })){
			maxRunCount = run;
		}
	}
	VecStr sampInfoCols(allSampInfo.columnNames_.begin(), allSampInfo.columnNames_.begin() + allSampInfo.getColPos("R1_name"));
	addOtherVec(sampInfoCols, tokenizeString(sampleCluster::getRepsInfoHeader(maxRunCount, true, "\t"), "\t"));
	ret.sampInfo_ = sampInfo.getColumns(sampInfoCols);
	ret.sampNames_ = ret.sampInfo_.getColumnLevels("s_Name");

	const auto & allHapIds = tabs_.hapIdTab_->get();
	VecStr hapIdCols{"#PopUID"};
	for(const auto & samp : samples){
		if(njh::in(samp, allHapIds.columnNames_)){
			hapIdCols.emplace_back(samp);
		}
	}
	ret.hapIdTab_ = allHapIds.getColumns(hapIdCols).extractByComp("#PopUID",
			[&popUids](const std::string & str) {return njh::in(str, popUids);});
	return ret;
}

const PopClusProject::ComputedGroup & PopClusProject::computedGroup(const std::string & group) {
	std::lock_guard<std::mutex> lock(computedGroupsMut_);
	auto search = computedGroups_.find(group);
	if(computedGroups_.end() != search){
		return search->second;
	}
	ComputedGroup ret;
	for(const auto & subGroup : collection_->groupMetaData_->groupData_.at(group)->subGroupToSamples_){
		ret.subGroups_[subGroup.first] = computeSubGroup(subGroup.second);
	}

	std::unordered_map<std::string, VecStr> popUids;
	for(const auto & subGroup : ret.subGroups_){
		popUids[subGroup.first] = subGroup.second.popInfo_.getColumnLevels("h_popUID");
	}
	std::unordered_map<std::string, VecStr> uniquePopUids;
	for(const auto & subGroupUids : popUids){
		for(const auto & uid : subGroupUids.second){
			bool unique = true;
			for(const auto & otherUids : popUids){
				if(subGroupUids.first != otherUids.first && njh::in(uid, otherUids.second)){
					unique = false;
					break;
				}
			}
			if(unique){
				uniquePopUids[subGroupUids.first].emplace_back(uid);
			}
		}
	}

	VecStr groupInfoColNames { "g_GroupName", "p_TotalInputReadCnt",
			"p_TotalInputClusterCnt", "p_TotalPopulationSampCnt",
			"p_TotalHaplotypes", "p_TotalUniqueHaplotypes", "p_ExpectedHeterozygosity",
			"p_meanCoi", "p_medianCoi", "p_minCoi",
			"p_maxCoi", "g_hapsFoundOnlyInThisGroup"};
	table groupInfo(groupInfoColNames);
	for(auto & subGroup : ret.subGroups_){
		const auto groupName = group + ":" + subGroup.first;
		auto & popInfo = subGroup.second.popInfo_;
		popInfo.addColumn({groupName}, "g_GroupName");
		popInfo.addColumn({njh::conToStr(uniquePopUids[subGroup.first])}, "g_hapsFoundOnlyInThisGroup");
		popInfo.addColumn({estd::to_string(uniquePopUids[subGroup.first].size())}, "p_TotalUniqueHaplotypes");
		groupInfo.rbind(popInfo.getColumns(groupInfoColNames), false);

		auto & sampInfo = subGroup.second.sampInfo_;
		auto sampInfoCols = sampInfo.columnNames_;
		removeElement(sampInfoCols, std::string("s_Sample"));
		VecStr outCols{"s_Sample", "g_GroupName"};
		addOtherVec(outCols, sampInfoCols);
		sampInfo.addColumn({groupName}, "g_GroupName");
		sampInfo = sampInfo.getColumns(outCols);
	}
	ret.groupInfo_ = groupInfo.getUniqueRows();
	ret.groupInfo_.sortTable("g_GroupName", false);
	return computedGroups_.emplace(group, std::move(ret)).first->second;
}

}  // namespace njhseq
