//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
/*
 * popClusterViewer.cpp
 *
 *  Created on: Jan 13, 2015
 *      Author: nickhathaway
 */

#include "popClusterViewer.hpp"

namespace bibseq {



pcv::pcv(cppcms::service& srv, std::map<std::string, std::string> config) :
		bibseq::seqApp(srv, config), config_(config){
	bool pass = configTest(config, requiredOptions(), "pcv");
	if(pass){
		std::cout << "Passed config test" << std::endl;
	}else{
		std::cout << "Didn't pass config test " << std::endl;
	}

	rootName_ =  config["name"];
	configDir_ = config["configDir"];
	debug_ = config["debug"] == "true";
	//add html
	pages_.emplace("mainPageHtml",
			make_path(config["resources"] + "pcv/mainPage.html"));
	pages_.emplace("mainProjectPageHtml",
			make_path(config["resources"] + "pcv/mainProjectPage.html"));
	pages_.emplace("redirectPage",
			make_path(config["resources"] + "pcv/redirectPage.html"));
	pages_.emplace("individualSample",
			make_path(config["resources"] + "pcv/individualSample.html"));
	pages_.emplace("extractionStats",
			make_path(config["resources"] + "pcv/extractionStats.html"));
	pages_.emplace("subGroupsPage",
			make_path(config["resources"] + "pcv/subGroupsPage.html"));
	pages_.emplace("groupMainPage",
			make_path(config["resources"] + "pcv/groupMainPage.html"));
	pages_.emplace("minTree",
			make_path(config["resources"] + "pcv/minTree.html"));
	for(auto & fCache : pages_){
		fCache.second.replaceStr("/pcv", rootName_);
	}

	//add pcv javascript and css
	jsAndCss_.emplace("jsPcv", getVectorOfMapKeys(bib::files::listAllFiles(config["resources"] + "pcv/js",false,VecStr{})));
	jsAndCss_.emplace("cssPcv", getVectorOfMapKeys(bib::files::listAllFiles(config["resources"] + "pcv/css",false,VecStr{})));


	//main page
	dispMapRoot(&pcv::mainPage, this);
	dispMap(&pcv::projectsNames, this, "projectsNames");
	dispMap_1word(&pcv::mainProjectPage, this, "mainProjectPage");
	dispMap_1word(&pcv::projectName, this, "projectName");

	dispMap_1word(&pcv::getPopInfo, this, "popInfo");
	dispMap_1word(&pcv::getPosSeqData, this, "popSeqData");
	dispMap_1word(&pcv::getPopInfoForSamps, this, "popInfoForSamps");
	dispMap_1word(&pcv::getPosSeqDataForSamps, this, "popSeqDataForSamps");

	dispMap_1word(&pcv::getSampInfo, this, "sampInfo");
	dispMap_1word(&pcv::getSampleNames, this, "sampleNames");

	//group info
	dispMap_1word(&pcv::getGroupNames, this, "getGroupNames");
	dispMap_3word(&pcv::getGroupPopInfo, this, "groupPopInfo");
	dispMap_3word(&pcv::getGroupPopSeqData, this, "groupPopSeqData");
	dispMap_3word(&pcv::getGroupPopInfoForSamps, this, "groupPopInfoForSamps");
	dispMap_3word(&pcv::getGroupPopSeqDataForSamps, this, "groupPopSeqDataForSamps");

	dispMap_3word(&pcv::getGroupSampInfo, this, "groupSampInfo");
	dispMap_3word(&pcv::getGroupSampleNames, this, "groupSampleNames");

	dispMap_3word(&pcv::showGroupMainPage, this, "showGroupMainPage");
	dispMap_2word(&pcv::showSubGroupsPage, this, "showSubGroupsPage");
	dispMap_2word(&pcv::getSubGroupsForGroup, this, "getSubGroupsForGroup");
	dispMap_2word(&pcv::getGroupPopInfos, this, "getGroupPopInfos");


	dispMap_1word(&pcv::getFracCutOff, this, "getFracCutOff");
	dispMap_1word(&pcv::setFracCutOff, this, "setFracCutOff");

	dispMap_1word(&pcv::getAllSampleNames, this, "allSampleNames");

	dispMap_1word(&pcv::getSampleNamesEncoding, this, "sampleNamesEncoding");
	dispMap_1word(&pcv::getEncodingForSampleNames, this, "encodingForSampleNames");


	dispMap_1word(&pcv::showMinTree, this, "showMinTree");
	dispMap_1word(&pcv::getMinTreeData, this, "minTreeData");


	dispMap_1word(&pcv::showExtractionInfo, this, "showExtractionInfo");
	dispMap_1word(&pcv::getIndexExtractionInfo, this, "getIndexExtractionInfo");
	dispMap_1word(&pcv::getSampleExtractionInfo, this, "getSampleExtractionInfo");


	dispMap_2word(&pcv::individualSamplePage, this, "individualSamplePage");
	dispMap_2word(&pcv::getSeqData, this, "seqData");
	dispMap_2word(&pcv::showMinTreeForSample, this, "showMinTreeForSample");
	dispMap_2word(&pcv::getMinTreeDataForSample, this, "minTreeDataForSample");

	dispMap(&pcv::jsPcv, this, "jsPcv");
	dispMap(&pcv::cssPcv, this, "cssPcv");

	mapper().root(rootName_);
	std::cout << "Finished set up" << std::endl;
}

void pcv::main(std::string url) {
	if (!dispatcher().dispatch(url)) {
		// Set the 404 status.
		response().status(cppcms::http::response::not_found);
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::jsPcv(){
	ret_js();
	auto search = jsAndCss_.find("jsPcv");
	response().out() << search->second.get("/pcv", rootName_);
}

void pcv::cssPcv(){
	ret_css();
	auto search = jsAndCss_.find("cssPcv");
	response().out() << search->second.get("/pcv", rootName_);
}

bool pcv::hasProject(const std::string shortProjName){
	return projectModels_.find(shortProjName) != projectModels_.end();
}

void pcv::projectsNames(){
	loadInProjects();
	ret_json();
	response().out() << bib::json::toJson(getVectorOfMapKeys(projectModels_));
}

void pcv::loadInProjects(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	auto files = bib::files::listAllFiles(configDir_, false, VecStr{});
	for(const auto & f : files){
		if(!f.second){
			Json::Reader jReader;
			Json::Value configJson;
			jReader.parse(bib::files::get_file_contents(f.first, debug_), configJson);
			if(!bib::in(configJson["shortName"].asString(), projectModels_)){
				projectModels_.emplace(configJson["shortName"].asString(), pcm(configJson, seqs_));
			}
		}
	}
}

void pcv::mainPage(){
	auto search = pages_.find("mainPageHtml");
	response().out() << search->second.get("/pcv", rootName_);
}



void pcv::mainProjectPage(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto search = pages_.find("mainProjectPageHtml");
		response().out() << search->second.get("/pcv", rootName_);
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::showMinTree(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto treeHtml = pages_.find("minTree")->second.get("/pcv", rootName_);
		treeHtml = replaceString(treeHtml, "MIN_TREE_JSON_LINK", rootName_  + "/minTreeData" + "/" + shortProjName);
		response().out() << treeHtml;
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::individualSamplePage(std::string shortProjName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"sampName",sampName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto proj = projectModels_.find(shortProjName);
		sampName = proj->second.decodeSampEncoding(sampName);
		if(bib::in(sampName, proj->second.clusteredSampleNames_)){
			auto search = pages_.find("individualSample");
			response().out() << search->second.get("/pcv", rootName_);
		}else{
			auto search = pages_.find("redirectPage");
			std::cout << "SampName: " << sampName << " not found, redirecting" << std::endl;
			response().out() << search->second.get("/pcv", rootName_);
		}
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::showMinTreeForSample(std::string shortProjName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"sampName", sampName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto proj = projectModels_.find(shortProjName);
		std::string encodedName = sampName;
		sampName = proj->second.decodeSampEncoding(sampName);
		if(bib::in(sampName, proj->second.clusteredSampleNames_)){
			auto treeHtml = pages_.find("minTree")->second.get("/pcv", rootName_);
			treeHtml = replaceString(treeHtml, "MIN_TREE_JSON_LINK", rootName_ + "/minTreeDataForSample/" + shortProjName + "/" + encodedName);
			response().out() << treeHtml;
		}else{
			auto search = pages_.find("redirectPage");
			std::cout << "SampName: " << sampName << " not found, redirecting" << std::endl;
			response().out() << search->second.get("/pcv", rootName_);
		}
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}


void pcv::showExtractionInfo(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto search = pages_.find("extractionStats");
		response().out() << search->second.get("/pcv", rootName_);
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::showGroupMainPage(std::string shortProjName, std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}, {"subGroup", subGroup}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto proj = projectModels_.find(shortProjName);
		auto search = proj->second.groupInfosDirNames_.find(group);
		if(search == proj->second.groupInfosDirNames_.end()){
			std::cout << "No group: " << group << "\n";
			std::cout << "options: " << vectorToString(getVectorOfMapKeys(proj->second.groupInfosDirNames_), ", ") << "\n";
			auto search = pages_.find("redirectPage");
			std::cout << "group: " << group << " not found, redirecting" << std::endl;
			response().out() << search->second.get("/pcv", rootName_);
		}else{
			auto subSearch = search->second.find(subGroup);
			if(subSearch == search->second.end()){
				std::cout << "No subgroup: " << subGroup << " in group: " << group << "\n";
				std::cout << "options: " << vectorToString(getVectorOfMapKeys(search->second), ", ") << "\n";
				auto search = pages_.find("redirectPage");
				std::cout << "group: " << group << " not found, redirecting" << std::endl;
				response().out() << search->second.get("/pcv", rootName_);
			}else{
				auto search = pages_.find("groupMainPage");
				response().out() << search->second.get("/pcv", rootName_);
			}
		}
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::showSubGroupsPage(std::string shortProjName, std::string group){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		auto proj = projectModels_.find(shortProjName);
		auto search = proj->second.groupInfosDirNames_.find(group);
		if(search == proj->second.groupInfosDirNames_.end()){
			std::cout << "No group: " << group << "\n";
			std::cout << "options: " << vectorToString(getVectorOfMapKeys(proj->second.groupInfosDirNames_), ", ") << "\n";
			auto search = pages_.find("redirectPage");
			std::cout << "group: " << group << " not found, redirecting" << std::endl;
			response().out() << search->second.get("/ssv", rootName_);
		}else{
			auto search = pages_.find("subGroupsPage");
			response().out() << search->second.get("/ssv", rootName_);
		}
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get("/pcv", rootName_);
	}
}

void pcv::setFracCutOff(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
	  auto postData = request().post();
	  auto postJson = bib::json::toJson(postData);
	  std::cout << postJson << std::endl;
	  projectModels_.find(shortProjName)->second.setFracCutOff(bib::lexical_cast<double>(postJson["input"].asString()));
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getFracCutOff(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		ret_json();
		response().out() <<projectModels_.find(shortProjName)->second.getFracCutOff();
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}


void pcv::getGroupNames(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if(hasProject(shortProjName)){
		ret_json();
		auto ret = projectModels_.find(shortProjName)->second.getGroupNames();
		response().out() << ret ;
	}else{
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::projectName(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getProjectName();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getSampleNames(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getSampleNames();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getAllSampleNames(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getAllSampleNames();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getSampleNamesEncoding(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getSampleNamesEncoding();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getEncodingForSampleNames(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getEncodingForSampleNames();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getSampInfo(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec{};
		for(const auto & kv : postData){
			if(kv.first == "sampNames[]"){
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out() <<  projectModels_.find(shortProjName)->second.getSampInfo(sampNamesVec);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getPosSeqData(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getPosSeqData();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getPopInfo(std::string shortProjName) {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getPopInfo();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getPosSeqDataForSamps(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec{};
		for(const auto & kv : postData){
			if(kv.first == "sampNames[]"){
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out() <<  projectModels_.find(shortProjName)->second.getPopSeqDataForSamps(sampNamesVec);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getPopInfoForSamps(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec{};
		for(const auto & kv : postData){
			if(kv.first == "sampNames[]"){
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out() <<  projectModels_.find(shortProjName)->second.getPopInfoForSamps(sampNamesVec);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}






void pcv::getMinTreeData(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getMinTreeData();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getSeqData(std::string shortProjName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"sampName", sampName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getSeqData(sampName);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}



void pcv::getMinTreeDataForSample(std::string shortProjName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"sampName", sampName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getMinTreeDataForSample(sampName);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getIndexExtractionInfo(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getIndexExtractionInfo();
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getSampleExtractionInfo(std::string shortProjName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec{};
		for(const auto & kv : postData){
			if(kv.first == "sampNames[]"){
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out() <<  projectModels_.find(shortProjName)->second.getSampleExtractionInfo(sampNamesVec);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getGroupPopInfo(std::string shortProjName, std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}, {"subGroup", subGroup}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getGroupPopInfo(group,subGroup );
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getGroupPopSeqData(std::string shortProjName, std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}, {"subGroup", subGroup}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getGroupPopSeqData(group,subGroup );
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getGroupPopInfoForSamps(std::string shortProjName, std::string group,
		std::string subGroup) {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, { {
			"shortProjName", shortProjName }, { "group", group }, { "subGroup",
			subGroup } }), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec { };
		for (const auto & kv : postData) {
			if (kv.first == "sampNames[]") {
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out()
				<< projectModels_.find(shortProjName)->second.getGroupPopInfoForSamps(
						group, subGroup, sampNamesVec);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getGroupPopSeqDataForSamps(std::string shortProjName, std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}, {"subGroup", subGroup}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec { };
		for (const auto & kv : postData) {
			if (kv.first == "sampNames[]") {
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out() <<  projectModels_.find(shortProjName)->second.getGroupPopSeqDataForSamps(group,subGroup, sampNamesVec);
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}


void pcv::getGroupSampInfo(std::string shortProjName, std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}, {"subGroup", subGroup}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		auto postData = request().post();
		VecStr sampNamesVec{};
		for(const auto & kv : postData){
			if(kv.first == "sampNames[]"){
				sampNamesVec.emplace_back(kv.second);
			}
		}
		response().out() <<  projectModels_.find(shortProjName)->second.getGroupSampInfo(group,subGroup,sampNamesVec );
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}

void pcv::getGroupSampleNames(std::string shortProjName, std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}, {"subGroup", subGroup}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getGroupSampleNames(group,subGroup );
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}
void pcv::getSubGroupsForGroup(std::string shortProjName, std::string group){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getSubGroupsForGroup(group );
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}


void pcv::getGroupPopInfos(std::string shortProjName, std::string group){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"shortProjName", shortProjName}, {"group", group}}), std::cout, debug_);
	if (hasProject(shortProjName)) {
		ret_json();
		response().out() <<  projectModels_.find(shortProjName)->second.getGroupPopInfos(group );
	} else {
		std::cout << "Doesn't contain project " << shortProjName << std::endl;
	}
}


int genProjectConfig(std::map<std::string, std::string> inputCommands) {
	bibseq::seqSetUp setUp(inputCommands);
	std::string mainDir = "";
	std::string extractionDir = "";
	std::string indexToDir = "";
	std::string sampNames = "";
	std::string projectName = "";
	std::string shortName = "";
	bool debug = false;
	setUp.setOption(debug, "-debug", "Run In Debug Mode");
	setUp.setOption(shortName, "-shortName", "Short Name for url for project",
			true);
	setUp.setOption(mainDir, "-mainDir", "Name of the Master Result Directory",
			true);
	if (mainDir.back() != '/') {
		mainDir.append("/");
	}
	setUp.setOption(extractionDir, "-extractionDir",
			"Name of the directory where extraction was done");
	setUp.setOption(indexToDir, "-indexToDir",
			"File, first column is index name, second is the name of the file extraction was done on",
			extractionDir != "");
	setUp.setOption(sampNames, "-sampNames",
			"A file, first column is index name, second is sample name, the rest of the columns are MID names",
			extractionDir != "");
	setUp.setOption(projectName, "-projectName", "Name of the Project", true);
	setUp.finishSetUp(std::cout);
	//
	Json::Value config;
	config["mainDir"] = bib::files::appendAsNeededRet(bib::files::normalize(mainDir).string(), "/");
	config["shortName"] = shortName;
	config["projectName"] = projectName;
	config["debug"] = debug;
	if(extractionDir != ""){
		config["extractionDir"] = bib::files::appendAsNeededRet(bib::files::normalize(extractionDir).string(), "/");
		config["indexToDir"] = bib::files::normalize(indexToDir).string();
		config["sampNames"] = bib::files::normalize(sampNames).string();
	}else{
		config["extractionDir"] = "";
		config["indexToDir"] = "";
		config["sampNames"] = "";
	}

	std::cout << config << std::endl;
	return 0;
}

int popClusteringViewer(std::map<std::string, std::string> inputCommands){
	bibseq::seqSetUp setUp(inputCommands);
	std::string configDir = "";

	uint32_t port = 9881;
	std::string name = "pcv";
	std::string resourceDirName = "";

	setUp.setOption(resourceDirName, "-resourceDirName", "Name of the resource Directory where the js and hmtl is located", true);
	if(resourceDirName.back() != '/'){
		resourceDirName.append("/");
	}
	setUp.setOption(configDir, "-configDir", "Name of the Master Result Directory", true);
	if(configDir.back() != '/'){
		configDir.append("/");
	}
	setUp.setOption(port, "-port", "Port Number to Serve On");
	setUp.setOption(name, "-name", "Name of root of the server");
	setUp.processDebug();
	setUp.finishSetUp(std::cout);
	name = "/" + name;
  auto config = server_config(name, port);
  //
  std::map<std::string, std::string> appConfig;
  appConfig["name"] = name;
  appConfig["configDir"] = configDir;
  appConfig["resources"] = resourceDirName;
  appConfig["js"] = resourceDirName + "js/";
  appConfig["css"] = resourceDirName + "css/";
  appConfig["debug"] = convertBoolToString(setUp.debug_);
  std::cout << "localhost:"  << port << name << std::endl;
	try {
		cppcms::service app(config);
		app.applications_pool().mount(
				cppcms::applications_factory<pcv>(appConfig));
		app.run();
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
	return 0;
}
} /* namespace bibseq */
