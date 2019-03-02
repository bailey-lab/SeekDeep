/*
 * pcv.cpp
 *
 *  Created on: Nov 25, 2016
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

#include "pcv.hpp"

namespace njhseq {



pcv::pcv(const Json::Value & config) :
		njhseq::SeqApp(config) {
	configDir_ = config["configDir"].asString();
	resourceDir_ = config["resources"].asString();

	jsFiles_->addFiles(
			njh::files::gatherFiles(njh::files::make_path(resourceDir_, "pcv/js"),
					".js"));
	cssFiles_->addFiles(
			njh::files::gatherFiles(njh::files::make_path(resourceDir_, "pcv/css"),
					".css"));
	loadInCollections();

	addScripts(njh::files::make_path(resourceDir_, "pcv"));
}


void pcv::loadInCollections(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);

	auto files = njh::files::listAllFiles(configDir_.string(), false, {std::regex{".*.config$"}});
	for(const auto & f : files){
		if(njh::beginsWith(f.first.filename().string(), ".")){
			continue;
		}
		if(!f.second){
			Json::Value configJson = njh::json::parseFile(f.first.string());
			njh::json::MemberChecker checker(configJson);
			checker.failMemberCheckThrow( { "shortName", "projectName", "mainDir" },
					__PRETTY_FUNCTION__);
			if (!njh::in(configJson["shortName"].asString(), collections_)) {
				auto coreJsonFnp = njh::files::make_path(configJson["mainDir"],
						"coreInfo.json");
				if (!bfs::exists(coreJsonFnp)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ": Error, project " << njh::bashCT::bold
							<< configJson["projectName"].asString() << njh::bashCT::reset
							<< " configuration from " << f.first
							<< " main directory doesn't contain coreInfo.json file" << "\n";
					ss << coreJsonFnp << " doesn't exist" << "\n";
					throw std::runtime_error { ss.str() };
				}
				Json::Value coreJson = njh::json::parseFile(coreJsonFnp.string());
				if (0 == coreJson["popNames_"]["samples_"].size()) {
					std::cerr << __PRETTY_FUNCTION__ << " folder "
							<< njh::bashCT::boldRed(coreJson["masterOutputDir_"].asString())
							<< " contains no data, not adding" << std::endl;
				} else {
					collections_.emplace(configJson["shortName"].asString(),
							std::make_unique<PopClusProject>(configJson));
					collections_[configJson["shortName"].asString()]->registerSeqFiles(
							*seqs_);
				}
			}
		}
	}
}



std::shared_ptr<restbed::Resource> pcv::extractionPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "extractionInfo" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						extractionPageHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getExtractionProfileData() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getExtractionProfileInfo" }, { "projectName",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getExtractionProfileDataHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getExtractionStatsData() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getExtractionStatsInfo" }, { "projectName",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getExtractionStatsDataHanlder(session);
					}));
	return resource;
}

void pcv::extractionPageHanlder(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	if (njh::in(projectName, collections_)) {
		auto body = genHtmlDoc(rootName_, pages_.at("extractionStats.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as " << projectName
				<< ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}

void pcv::getExtractionProfileDataHanlder(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	if (njh::in(projectName, collections_)) {
		Json::Value ret;
		if(nullptr != collections_.at(projectName)->extractionProfileTab_){
			auto tab = collections_.at(projectName)->extractionProfileTab_->get();
			tab.trimElementsAtFirstOccurenceOf("(");
			for(auto & row : tab.content_){
				row[tab.getColPos("name")] = njh::pasteAsStr(row[tab.getColPos("extractionDir")], "_", row[tab.getColPos("name")]);
			}
			ret = tableToJsonByRow(tab,"name");
		}
		auto body = njh::json::writeAsOneLine(ret);
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateAppJsonHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}

void pcv::getExtractionStatsDataHanlder(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	if (njh::in(projectName, collections_)) {
		Json::Value ret;
		if(nullptr != collections_.at(projectName)->extractionStatsTab_){
			auto tab = collections_.at(projectName)->extractionStatsTab_->get();
			tab.trimElementsAtFirstOccurenceOf("(");
			ret = tableToJsonByRow(tab,"extractionDir");
		}
		auto body = njh::json::writeAsOneLine(ret);
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateAppJsonHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}

void pcv::projectNamesHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	Json::Value ret;
	ret["projects"]= njh::json::toJson(getVectorOfMapKeys(collections_));
	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

void pcv::mainPageHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto body = genHtmlDoc(rootName_, pages_.at("mainPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}


void pcv::redirect(std::shared_ptr<restbed::Session> session,
		std::string errorMessage) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	std::cerr << errorMessage << std::endl;
	auto body = genHtmlDoc(rootName_, pages_.at("redirectPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}

void pcv::mainProjectPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	if (njh::in(projectName, collections_)) {
		auto body = genHtmlDoc(rootName_, pages_.at("mainProjectPage.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}

void pcv::samplePageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string sampleName = request->get_path_parameter("sampleName");

	if (njh::in(projectName, collections_)) {
		if(collections_[projectName]->collection_->hasSample(sampleName)){
			auto body = genHtmlDoc(rootName_, pages_.at("sampleMainPage.js"));
			const std::multimap<std::string, std::string> headers =
					HeaderFactory::initiateTxtHtmlHeader(body);
			session->close(restbed::OK, body, headers);
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, no such sample as "
					<< sampleName << " " << "in project " << projectName << ", options are "
					<< njh::conToStr(collections_[projectName]->collection_->passingSamples_, ", ") << "\n";
			ss << "Redirecting..." << "\n";
			redirect(session, ss.str());
		}
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}

void pcv::groupInfoPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");

	if (njh::in(projectName, collections_)) {
		if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
			if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
				auto body = genHtmlDoc(rootName_, pages_.at("groupInfoPage.js"));
				const std::multimap<std::string, std::string> headers =
						HeaderFactory::initiateTxtHtmlHeader(body);
				session->close(restbed::OK, body, headers);
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
						<< " " << "in project " << projectName << ", options are "
						<< njh::conToStr(
								getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
								", ") << "\n";
				ss << "Redirecting..." << "\n";
				redirect(session, ss.str());
			}
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
			ss << "Redirecting..." << "\n";
			redirect(session, ss.str());
		}
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}


void pcv::getGroupsPopInfosHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getGroupsPopInfosPostHandler(ses, body);
					}));
}

void pcv::getGroupsPopInfosPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	Json::Value ret;
	if (njh::in(projectName, collections_)) {
		if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
			if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
				ret["groupNames"] = njh::json::toJson(collections_[projectName]->collection_->groupMetaData_->groupData_.at(groupName)->subGroupsLevels_);
				ret["popInfo"] = tableToJsonByRow(collections_[projectName]->topGroupTabs_.at(groupName)->get(), "g_GroupName", VecStr{});
			} else {
				std::cerr << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
						<< " " << "in project " << projectName << ", options are "
						<< njh::conToStr(
								getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
								", ") << "\n";
			}
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
		}
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}


void pcv::getProjectNameHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto projectName = request->get_path_parameter("projectName");
	Json::Value ret;
	ret["projectName"] = "";
	if (njh::in(projectName, collections_)) {
		ret["projectName"] = collections_[projectName]->projectName_;
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
	}
	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

void pcv::getSampleNamesHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto projectName = request->get_path_parameter("projectName");
	Json::Value ret;
	ret["samples"] = "";
	if (njh::in(projectName, collections_)) {

		ret["samples"] = njh::json::toJson(collections_[projectName]->collection_->passingSamples_);
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
	}
	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

void pcv::getGroupNamesHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto projectName = request->get_path_parameter("projectName");
	Json::Value ret;
	ret["groups"] = "";
	if (njh::in(projectName, collections_)) {
		if (nullptr != collections_[projectName]->collection_->groupMetaData_) {
			ret["groups"] =
					njh::json::toJson(
							getVectorOfMapKeys(
									collections_[projectName]->collection_->groupMetaData_->groupData_));
		}
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
	}
	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}


VecStr genColorsStrsForNamesMultiples(VecStr popNames){
	removeDuplicates(popNames);
	auto popColors = getColorsForNames(popNames);
	VecStr popColorsStrs;
	auto keys = getVectorOfMapKeys(popColors);
	njh::sort(keys);
	for(const auto & pColorKey : keys){
		popColorsStrs.emplace_back(popColors[pColorKey].getHexStr());
	}
	return popColorsStrs;
}

void pcv::getSampleInfoTabPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "sampNames" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto request = session->get_request();
		std::string projectName = request->get_path_parameter("projectName");
		if (njh::in(projectName, collections_)) {
			auto sampNames = njh::json::jsonArrayToVec<std::string>(postData["sampNames"], [](const Json::Value & val){ return val.asString();});
			auto sampTable_ = collections_[projectName]->tabs_.sampInfo_->get();
			auto containsSampName = [&sampNames](const std::string & str) {
				return njh::in(str, sampNames);
			};
			auto trimedTab = sampTable_.extractByComp("s_Name", containsSampName);
			std::string coiColName = "s_FinalClusterCnt";
			if(njh::in(std::string("s_COI"), sampTable_.columnNames_)){
				coiColName = "s_COI";
			}
			VecStr visibleColumns = VecStr { "s_Sample",
				"h_popUID", "h_SampCnt", "h_SampFrac", "s_ReadCntTotUsed",
				coiColName, "c_clusterID", "c_AveragedFrac", "c_ReadCnt",
				"c_RepCnt" };
			if(njh::in(std::string("bestExpected"), sampTable_.columnNames_)){
				visibleColumns.emplace_back("bestExpected");
			}
			ret = tableToJsonByRow(trimedTab, "s_Name", visibleColumns);

			auto popNames = trimedTab.getColumn("h_popUID");
			removeDuplicates(popNames);
			njh::sort(popNames);
			auto popColorsStrs = genColorsStrsForNamesMultiples(popNames);
			ret["popColors"] = njh::json::toJson(popColorsStrs);
			ret["popUIDs"] = njh::json::toJson(popNames);
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
		}
	}



	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::getSampleInfoTabHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getSampleInfoTabPostHandler(ses, body);
					}));
}


void pcv::getPopSeqsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value seqData;
	if (checker.failMemberCheck( { "popUIDs" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto request = session->get_request();
		std::string projectName = request->get_path_parameter("projectName");
		if (njh::in(projectName, collections_)) {
			uint32_t sesUid = std::numeric_limits<uint32_t>::max();
			//check to see if there is a session already started associated with this seq
			if(!postData.isMember("sessionUID")){
				sesUid = startSeqCacheSession();
			}else{
				sesUid = postData["sessionUID"].asUInt();
			}
			auto popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
					[](const Json::Value & val) {return val.asString();});
			seqsBySession_[sesUid]->cache_.at(projectName).reload();
			seqsBySession_[sesUid]->cache_.at(projectName).toggleSeqs(
					[&popUIDs](const readObject & seq) {
						return njh::in(seq.seqBase_.getStubName(false), popUIDs);
			});
			//make sure seqs aren't empty, viewer doesn't know how to handle that, if it isn't make sure to remove the placeholder seq if it is there
			seqsBySession_[sesUid]->cache_.at(projectName).ensureNonEmptyReads();
			seqData = seqsBySession_[sesUid]->getJson(projectName);
			seqData["sessionUID"] = njh::json::toJson(sesUid);

		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
		}
	}
	/**@todo check other headers for connection close
	 *
	 */
	auto retBody = njh::json::writeAsOneLine(seqData);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::getPopSeqsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getPopSeqsPostHandler(ses, body);
					})); //s_FinalClusterCnt
}

void pcv::getSampSeqsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value seqData;
	if (checker.failMemberCheck( { "sampleName" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto request = session->get_request();
		std::string projectName = request->get_path_parameter("projectName");
		if (njh::in(projectName, collections_)) {
			auto sampleName = postData["sampleName"].asString();
			if (collections_[projectName]->collection_->hasSample(sampleName)) {
				uint32_t sesUid = startSeqCacheSession();
				seqData = seqsBySession_[sesUid]->getJson(
						projectName + "_" + sampleName);
				seqData["sessionUID"] = njh::json::toJson(sesUid);
			} else {
				std::cerr << __PRETTY_FUNCTION__ << ": error, no such sample as "
						<< sampleName << " " << "in project " << projectName << ", options are "
						<< njh::conToStr(collections_[projectName]->collection_->passingSamples_, ", ") << std::endl;
			}
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << std::endl;
		}
	}
	/**@todo check other headers for connection close
	 *
	 */
	auto retBody = njh::json::writeAsOneLine(seqData);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}



void pcv::getSampSeqsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getSampSeqsPostHandler(ses, body);
					}));
}

void pcv::getPopInfoPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value popInfo;
	if (checker.failMemberCheck( { "popUIDs" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto request = session->get_request();
		std::string projectName = request->get_path_parameter("projectName");
		if (njh::in(projectName, collections_)) {

			auto popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
					[](const Json::Value & val) {return val.asString();});
			auto trimedPopTab =
					collections_[projectName]->tabs_.popInfo_->get().extractByComp(
							"h_popUID",
							[&popUIDs](const std::string & str) {return njh::in(str, popUIDs);});
			popInfo = tableToJsonByRow(trimedPopTab, "h_popUID", VecStr { }, VecStr {
					"p_TotalInputReadCnt", "p_TotalInputClusterCnt",
					"p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanCoi",
					"p_medianCoi", "p_minCoi", "p_maxCoi" });
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
		}
	}
	/**@todo check other headers for connection close
	 *
	 */
	auto retBody = njh::json::writeAsOneLine(popInfo);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::getHapIdTablePostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "popUIDs", "samples" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto request = session->get_request();
		std::string projectName = request->get_path_parameter("projectName");
		if (njh::in(projectName, collections_)) {
			auto popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
					[](const Json::Value & val) {return val.asString();});
			auto samples = njh::json::jsonArrayToVec<std::string>(postData["samples"],
					[](const Json::Value & val) {return val.asString();});
			auto hapIdTab = collections_[projectName]->tabs_.hapIdTab_->get().getColumns(concatVecs(VecStr{"#PopUID"}, samples)).extractByComp(
					"#PopUID",
					[&popUIDs](const std::string & str) {return njh::in(str, popUIDs);});
			auto trimedHapIdTab =
					hapIdTab.extractByComp(
							"#PopUID",
							[&popUIDs](const std::string & str) {return njh::in(str, popUIDs);});
			ret = tableToJsonByRow(trimedHapIdTab, "#PopUID");

		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_)) << "\n";
		}
	}
	/**@todo check other headers for connection close
	 *
	 */
	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}



void pcv::getPopInfoHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
		getPopInfoPostHandler(ses, body);
					}));
}


//


std::shared_ptr<restbed::Resource> pcv::projectNames(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {"projectsNames"} }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						projectNamesHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::mainPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						mainPageHandler(session);
					}));
	return resource;
}


std::shared_ptr<restbed::Resource> pcv::groupInfoPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupInfoPage" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupInfoPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getGroupsPopInfos() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupsPopInfos" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getGroupsPopInfosHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::samplePage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "sampleMainPage" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"sampleName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						samplePageHandler(session);
					}));
	return resource;
}


std::shared_ptr<restbed::Resource> pcv::mainProjectPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "mainProjectPage" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						mainProjectPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getProjectName() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "projectName" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getProjectNameHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getSampleNames() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "sampleNames" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getSampleNamesHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getGroupNames() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "getGroupNames" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getGroupNamesHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getSampleInfoTab() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, { "sampInfo" },
			{ "projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getSampleInfoTabHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getPopSeqs() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, { "popSeqData" },
			{ "projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getPopSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getSampSeqs() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, { "sampSeqData" },
			{ "projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getSampSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getPopInfo() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "popInfo" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getPopInfoHandler(session);
					}));
	return resource;
}



void pcv::groupMainProjectPageHanlder(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	std::string subGroupName = request->get_path_parameter("subGroupName");

	if (njh::in(projectName, collections_)) {
		if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
			if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
				if(njh::in(subGroupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_)){
					auto body = genHtmlDoc(rootName_, pages_.at("groupMainPage.js"));
					const std::multimap<std::string, std::string> headers =
							HeaderFactory::initiateTxtHtmlHeader(body);
					session->close(restbed::OK, body, headers);
				}else{
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ": error, no such sub group as " << subGroupName
							<< " in group " << groupName << " in project " << projectName << ", options are "
							<< njh::conToStr(
									getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_),
									", ") << "\n";
					ss << "Redirecting..." << "\n";
					redirect(session, ss.str());
				}
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
						<< " " << "in project " << projectName << ", options are "
						<< njh::conToStr(
								getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
								", ") << "\n";
				ss << "Redirecting..." << "\n";
				redirect(session, ss.str());
			}
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
			ss << "Redirecting..." << "\n";
			redirect(session, ss.str());
		}
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		ss << "Redirecting..." << "\n";
		redirect(session, ss.str());
	}
}

void pcv::groupGetSampleNamesHanlder(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	std::string subGroupName = request->get_path_parameter("subGroupName");

	Json::Value ret;
	if (njh::in(projectName, collections_)) {
		if (nullptr != collections_[projectName]->collection_->groupDataPaths_) {
			if (njh::in(groupName,
					collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
				if (njh::in(subGroupName,
						collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(
								groupName).groupPaths_)) {
					ret["groupSamples"] =
							njh::json::toJson(
									collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_.at(subGroupName).readInSampNames());
				} else {
					std::cerr << __PRETTY_FUNCTION__ << ": error, no such sub group as "
							<< subGroupName << " in group " << groupName << " in project "
							<< projectName << ", options are "
							<< njh::conToStr(
									getVectorOfMapKeys(
											collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(
													groupName).groupPaths_), ", ") << "\n";
				}
			} else {
				std::cerr << __PRETTY_FUNCTION__ << ": error, no such group as "
						<< groupName << " " << "in project " << projectName
						<< ", options are "
						<< njh::conToStr(
								getVectorOfMapKeys(
										collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
								", ") << "\n";
			}
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no group data for "
					<< projectName << "\n";
		}
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
				<< projectName << ", options are "
				<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::groupGetSampleInfoTabPostHanlder(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	std::string subGroupName = request->get_path_parameter("subGroupName");

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "sampNames" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		if (njh::in(projectName, collections_)) {
			if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
				if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
					if(njh::in(subGroupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_)){
						auto sampNames = njh::json::jsonArrayToVec<std::string>(postData["sampNames"], [](const Json::Value & val){ return val.asString();});
						auto sampTable_ = collections_[projectName]->subGroupTabs_.at(groupName).at(subGroupName).sampInfo_->get();
						auto containsSampName = [&sampNames](const std::string & str) {
							return njh::in(str, sampNames);
						};
						auto trimedTab = sampTable_.extractByComp("s_Name", containsSampName);
						std::string coiColName = "s_FinalClusterCnt";
						if(njh::in(std::string("s_COI"), sampTable_.columnNames_)){
							coiColName = "s_COI";
						}
						VecStr visibleColumns = VecStr { "s_Sample", "g_GroupName",
							"h_popUID", "h_SampCnt", "h_SampFrac", "s_ReadCntTotUsed",
							coiColName, "c_clusterID", "c_AveragedFrac", "c_ReadCnt",
							"c_RepCnt" };
						if(njh::in(std::string("bestExpected"), sampTable_.columnNames_)){
							visibleColumns.emplace_back("bestExpected");
						}
						ret = tableToJsonByRow(trimedTab, "s_Name", visibleColumns);

						auto popNames = trimedTab.getColumn("h_popUID");
						removeDuplicates(popNames);
						njh::sort(popNames);
						auto popColorsStrs = genColorsStrsForNamesMultiples(popNames);
						ret["popColors"] = njh::json::toJson(popColorsStrs);
						ret["popUIDs"] = njh::json::toJson(popNames);
					}else{
						std::cerr << __PRETTY_FUNCTION__ << ": error, no such sub group as " << subGroupName
								<< " in group " << groupName << " in project " << projectName << ", options are "
								<< njh::conToStr(
										getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_),
										", ") << "\n";
					}
				} else {
					std::cerr << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
							<< " " << "in project " << projectName << ", options are "
							<< njh::conToStr(
									getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
									", ") << "\n";
				}
			}else{
				std::cerr << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
			}
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		}
	}


	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::groupGetSampleInfoTabHanlder(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						groupGetSampleInfoTabPostHanlder(ses, body);
					}));
}

void pcv::groupGetPopSeqsPostHanlder(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	std::string subGroupName = request->get_path_parameter("subGroupName");

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "popUIDs" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		if (njh::in(projectName, collections_)) {
			if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
				if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
					if(njh::in(subGroupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_)){
						uint32_t sesUid;
						//check to see if there is a session already started associated with this seq
						if(!postData.isMember("sessionUID")){
							sesUid = startSeqCacheSession();
						}else{
							sesUid = postData["sessionUID"].asUInt();
						}
						auto popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
								[](const Json::Value & val) {return val.asString();});
						seqsBySession_[sesUid]->cache_.at(projectName).reload();
						if (seqsBySession_[sesUid]->cache_.at(projectName).reads_->size()
								!= popUIDs.size()) {
							for (auto & seq : *seqsBySession_[sesUid]->cache_.at(projectName).reads_) {
								if (njh::in(seq.seqBase_.getStubName(false), popUIDs)) {
									seq.seqBase_.on_ = true;
								} else {
									seq.seqBase_.on_ = false;
								}
							}
						}
						//make sure seqs aren't empty, viewer doesn't know how to handle that, if it isn't make sure to remove the placeholder seq if it is there
						seqsBySession_[sesUid]->cache_.at(projectName).ensureNonEmptyReads();
						ret = seqsBySession_[sesUid]->getJson(projectName);
						ret["sessionUID"] = njh::json::toJson(sesUid);
					}else{
						std::cerr << __PRETTY_FUNCTION__ << ": error, no such sub group as " << subGroupName
								<< " in group " << groupName << " in project " << projectName << ", options are "
								<< njh::conToStr(
										getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_),
										", ") << "\n";
					}
				} else {
					std::cerr << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
							<< " " << "in project " << projectName << ", options are "
							<< njh::conToStr(
									getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
									", ") << "\n";
				}
			}else{
				std::cerr << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
			}
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		}
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::groupGetPopSeqsHanlder(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						groupGetPopSeqsPostHanlder(ses, body);
					}));
}

void pcv::groupGetPopInfoPostHanlder(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	std::string subGroupName = request->get_path_parameter("subGroupName");

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "popUIDs" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		if (njh::in(projectName, collections_)) {
			if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
				if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
					if(njh::in(subGroupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_)){
						auto popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
								[](const Json::Value & val) {return val.asString();});
						auto trimedPopTab =
								collections_[projectName]->subGroupTabs_.at(groupName).at(subGroupName).popInfo_->get().extractByComp(
										"h_popUID",
										[&popUIDs](const std::string & str) {return njh::in(str, popUIDs);});
						ret = tableToJsonByRow(trimedPopTab, "h_popUID", VecStr { }, VecStr {
								"p_TotalInputReadCnt", "g_GroupName","g_hapsFoundOnlyInThisGroup",
								"p_TotalUniqueHaplotypes", "p_TotalInputClusterCnt",
								"p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanCoi",
								"p_medianCoi", "p_minCoi", "p_maxCoi"});
					}else{
						std::cerr << __PRETTY_FUNCTION__ << ": error, no such sub group as " << subGroupName
								<< " in group " << groupName << " in project " << projectName << ", options are "
								<< njh::conToStr(
										getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_),
										", ") << "\n";
					}
				} else {
					std::cerr << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
							<< " " << "in project " << projectName << ", options are "
							<< njh::conToStr(
									getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
									", ") << "\n";
				}
			}else{
				std::cerr << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
			}
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		}
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::groupGetHapIdTablePostHanlder(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	std::string projectName = request->get_path_parameter("projectName");
	std::string groupName = request->get_path_parameter("groupName");
	std::string subGroupName = request->get_path_parameter("subGroupName");

	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "popUIDs", "samples" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		if (njh::in(projectName, collections_)) {
			if(nullptr != collections_[projectName]->collection_->groupDataPaths_){
				if (njh::in(groupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_)) {
					if(njh::in(subGroupName, collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_)){
						auto popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
								[](const Json::Value & val) {return val.asString();});
						auto samples = njh::json::jsonArrayToVec<std::string>(postData["samples"],
								[](const Json::Value & val) {return val.asString();});
						auto hapIdTab = collections_[projectName]->subGroupTabs_.at(groupName).at(subGroupName).hapIdTab_->get().getColumns(concatVecs(VecStr{"#PopUID"}, samples)).extractByComp(
								"#PopUID",
								[&popUIDs](const std::string & str) {return njh::in(str, popUIDs);});

						auto trimedHapIdTab = hapIdTab.extractByComp(
										"#PopUID",
										[&popUIDs](const std::string & str) {return njh::in(str, popUIDs);});
						ret = tableToJsonByRow(trimedHapIdTab, "#PopUID");
					}else{
						std::cerr << __PRETTY_FUNCTION__ << ": error, no such sub group as " << subGroupName
								<< " in group " << groupName << " in project " << projectName << ", options are "
								<< njh::conToStr(
										getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_.at(groupName).groupPaths_),
										", ") << "\n";
					}
				} else {
					std::cerr << __PRETTY_FUNCTION__ << ": error, no such group as " << groupName
							<< " " << "in project " << projectName << ", options are "
							<< njh::conToStr(
									getVectorOfMapKeys(collections_[projectName]->collection_->groupDataPaths_->allGroupPaths_),
									", ") << "\n";
				}
			}else{
				std::cerr << __PRETTY_FUNCTION__ << ": error, no group data for "  << projectName << "\n";
			}
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": error, no such project as "
					<< projectName << ", options are "
					<< njh::conToStr(njh::getVecOfMapKeys(collections_), ", ") << "\n";
		}
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void pcv::groupGetPopInfoHanlder(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						groupGetPopInfoPostHanlder(ses, body);
					}));
}




void pcv::groupGetHapIdTableHanlder(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
		groupGetHapIdTablePostHanlder(ses, body);
					}));
}

void pcv::getHapIdTableHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
		getHapIdTablePostHandler(ses, body);
					}));
}


std::shared_ptr<restbed::Resource> pcv::groupMainProjectPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupMainPage" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ },
		{"subGroupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupMainProjectPageHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::groupGetSampleNames(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupSampleNames" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ },
		{"subGroupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupGetSampleNamesHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::groupGetSampleInfoTab(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupSampInfo" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ },
		{"subGroupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupGetSampleInfoTabHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::groupGetPopSeqs(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupPopSeqData" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ },
		{"subGroupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupGetPopSeqsHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::groupGetPopInfo(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupPopInfo" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ },
		{"subGroupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupGetPopInfoHanlder(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::getHapIdTable(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "hapIdTable" }, {
					"projectName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getHapIdTableHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> pcv::groupGetHapIdTable(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "groupHapIdTable" },
		{"projectName", UrlPathFactory::pat_wordNumsDash_ },
		{"groupName", UrlPathFactory::pat_wordNumsDash_ },
		{"subGroupName", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						groupGetHapIdTableHanlder(session);
					}));
	return resource;
}





std::vector<std::shared_ptr<restbed::Resource>> pcv::getAllResources() {
	auto ret = super::getAllResources();
	ret.emplace_back(mainPage());
	ret.emplace_back(projectNames());

	//project
	ret.emplace_back(mainProjectPage());
	ret.emplace_back(getProjectName());
	ret.emplace_back(getSampleNames());
	ret.emplace_back(getGroupNames());
	ret.emplace_back(getSampleInfoTab());
	ret.emplace_back(getPopSeqs());
	ret.emplace_back(getPopInfo());
	ret.emplace_back(getHapIdTable());

	//sample
	ret.emplace_back(samplePage());
	ret.emplace_back(getSampSeqs());

	//top group
	ret.emplace_back(groupInfoPage());
	ret.emplace_back(getGroupsPopInfos());

	//sub group
	ret.emplace_back(groupMainProjectPage());
	ret.emplace_back(groupGetSampleNames());
	ret.emplace_back(groupGetSampleInfoTab());
	ret.emplace_back(groupGetPopSeqs());
	ret.emplace_back(groupGetPopInfo());
	ret.emplace_back(groupGetHapIdTable());

	//extraction
	ret.emplace_back(extractionPage());
	ret.emplace_back(getExtractionProfileData());
	ret.emplace_back(getExtractionStatsData());
	return ret;
}



VecStr pcv::requiredOptions() const {
	return concatVecs(super::requiredOptions(), VecStr { "configDir",
			"resources" });
}

}  // namespace njhseq

