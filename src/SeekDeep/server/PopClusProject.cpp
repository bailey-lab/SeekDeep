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
			topGroupTabs_[group.first] = std::make_unique<TableCache>(TableIOOpts(InOptions(group.second.groupInfoFnp_), "\t", true));
			for(const auto & subGroup : group.second.groupPaths_){
				subGroupTabs_[group.first][subGroup.first].popInfo_ = std::make_unique<TableCache>(TableIOOpts(InOptions(subGroup.second.popFileFnp_), "\t", true));
				subGroupTabs_[group.first][subGroup.first].sampInfo_ = std::make_unique<TableCache>(TableIOOpts(InOptions(subGroup.second.sampFileFnp_), "\t", true));
				subGroupTabs_[group.first][subGroup.first].hapIdTab_ = std::make_unique<TableCache>(TableIOOpts(InOptions(subGroup.second.hapIdTabFnp_), "\t", true));
			}
		}
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

}  // namespace njhseq
