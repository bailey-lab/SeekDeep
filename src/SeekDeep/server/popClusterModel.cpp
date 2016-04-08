/*
 * popClusterModel.cpp
 *
 *  Created on: Aug 26, 2015
 *      Author: Nick Hathaway
 */
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
#include "popClusterModel.hpp"

namespace bibseq {
VecStr getColorsStrsForNamesMultiples(VecStr popNames){
	removeDuplicates(popNames);
	auto popColors = getColorsForNames(popNames);
	VecStr popColorsStrs;
	auto keys = getVectorOfMapKeys(popColors);
	bib::sort(keys);
	for(const auto & pColorKey : keys){
		popColorsStrs.emplace_back("#" + popColors[pColorKey].hexStr_);
	}
	return popColorsStrs;
}

pcm::pcm(Json::Value config, std::shared_ptr<seqCache> seqs) :
		config_(config), seqs_(seqs){
	bool pass = true;
	if(!config_.isMember("shortName") ||
			!config_.isMember("mainDir") ||
			!config_.isMember("projectName") ){
		pass = false;
	}
	if(config_.isMember("debug")){
		debug_ = config_["debug"].asBool();
	}
	mainDir_ = config_["mainDir"].asString();
	projectName_ = config_["projectName"].asString();

	//printVector( config_.getMemberNames());
	if(!pass){
		Json::FastWriter jwriter;
		throw std::runtime_error{"Error in config of pcm\n" + jwriter.write(config) };
	}

	if (config_.isMember("extractionDir") && config_["extractionDir"].asString() != "") {
		extractionDir_ = config_["extractionDir"].asString();
		extractInfo_ = collectExtractionInfoDirectName(
				config_["extractionDir"].asString(), config_["indexToDir"].asString(),
				config_["sampNames"].asString());
	}

	sampTable_ = table(mainDir_ + "selectedClustersInfo.tab.txt", "\t", true);

	auto fracs = vecStrToVecNum<double>(sampTable_.getColumn("c_AveragedFrac"));
	fracCutOff_ = vectorMinimum(fracs);
	//get samp names
	auto sampCounts = countVec(sampTable_.getColumn("s_Name"));
	clusteredSampleNames_ = getVectorOfMapKeys(sampCounts);
	//read pop table
	popTable_ = table(mainDir_ + "population/populationCluster.tab.txt", "\t", true);
	popTable_.sortTable("h_popUID", false);

	//encode sample names
	std::set<std::string> allUniSampNames(clusteredSampleNames_.begin(), clusteredSampleNames_.end());
	if(!extractInfo_.profileBySampTab_.content_.empty()){
		for(const auto & samp : extractInfo_.profileBySampTab_.getColumn("Sample")){
			allUniSampNames.emplace(samp);
		}
	}
	allSampleNames_ = std::vector<std::string>(allUniSampNames.begin(), allUniSampNames.end());
	for(const auto & pos : iter::range(allSampleNames_.size())){
		codedNumToSampName_[pos] = allSampleNames_[pos];
		sampNameToCodedNum_[allSampleNames_[pos]] = pos;
	}

	if(bfs::exists(mainDir_ + "/groups")){
		auto sampFiles = bib::files::listAllFiles(mainDir_ + "/groups", true, VecStr{"sampFile.tab.txt"});
		for(const auto & sampF : sampFiles){
			auto toks = tokenizeString(sampF.first.string(), "/");
			groupInfosDirNames_[toks[toks.size() - 3]][toks[toks.size() - 2]] = sampF.first.parent_path().string();
		}
	}
	if(debug_){
		std::cout << "Finished set up" << std::endl;
	}
}

std::string pcm::messStrFactory(const std::string & funcName)const{
	return bib::err::F() << "[" << getCurrentDate() << "] " << funcName;
}

std::string pcm::messStrFactory(const std::string & funcName, const MapStrStr & args){
	VecStr argsVec;
	for(const auto & kv : args){
		argsVec.emplace_back(kv.first + " = " + kv.second);
	}
	std::string argStrs = messStrFactory(funcName) + " [" + vectorToString(argsVec, ", ") + "]";
	return argStrs;
}

void pcm::setFracCutOff(double fracCutOff){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
  fracCutOff_ = fracCutOff;
  fracCutOff_ /=100;
}

double pcm::getFracCutOff(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	return roundDecPlaces((fracCutOff_ * 100),2);
}


Json::Value pcm::getGroupNames(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	auto gNames = getVectorOfMapKeys(groupInfosDirNames_);
	Json::Value ret = bib::json::toJson(getVectorOfMapKeys(groupInfosDirNames_));
	return ret;
}

//json
Json::Value pcm::getProjectName(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	Json::Value ret = projectName_;
	return ret;
}

Json::Value pcm::getSampleNames(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	Json::Value ret = bib::json::toJson(clusteredSampleNames_);
	return ret;
}

Json::Value pcm::getAllSampleNames(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	Json::Value ret = bib::json::toJson(allSampleNames_);
	return ret;
}

Json::Value pcm::getSampleNamesEncoding(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	Json::Value ret;
	for(const auto & sampName : sampNameToCodedNum_){
		ret[sampName.first] = sampName.second;
	}
	return ret;
}

Json::Value pcm::getEncodingForSampleNames(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	Json::Value ret;
	for(const auto & sampName : sampNameToCodedNum_){
		ret[sampName.second] = sampName.first;
	}
	return ret;
}

Json::Value pcm::getSampInfo(const VecStr & sampNames){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"sampNames", vectorToString(sampNames, ",")}}), std::cout, debug_);
	Json::Value ret;

	auto sampNameToCodedNum = sampNameToCodedNum_;
	auto containsSampName = [&sampNames, &sampNameToCodedNum](const std::string & str) {
		return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
	};
	auto trimedTab = sampTable_.extractByComp("s_Name", containsSampName);
	ret = cppcmsJsonToJson(tableToJsonRowWise(trimedTab,"s_Name", VecStr{}));
	auto popNames = trimedTab.getColumn("h_popUID");
	removeDuplicates(popNames);
	auto popColorsStrs = getColorsStrsForNamesMultiples(popNames);
	ret["popColors"] = bib::json::toJson(popColorsStrs);
	return ret;
}

void pcm::loadInPopSeqs(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	std::string popSeqDataUid = projectName_ + "_seqs";
	if(!seqs_->containsRecord(popSeqDataUid) || !seqs_->recordValid(popSeqDataUid)){
		/**todo make safter for non fastq pop clustering */
		//read pop seqs
		auto files = bib::files::listAllFiles(mainDir_ + "population/", false, {std::regex{".*\\.fastq"}});
		SeqIOOptions popOpts;
		popOpts.firstName_ = files.begin()->first.string();
		popOpts.inFormat_ = SeqIOOptions::getInFormat(bib::files::getExtension(files.begin()->first.string()));
		popOpts.processed_ = true;
		SeqInput reader(popOpts);
		reader.openIn();
		auto reads = reader.readAllReads<readObject>();
		for(auto & read : reads){
			read.seqBase_.name_ = read.seqBase_.name_.substr(0, read.seqBase_.name_.rfind("_"));
		}
		seqs_->updateAddCache(popSeqDataUid, std::make_shared<std::vector<readObject>>(std::vector<readObject>(reads)));
	}
}
Json::Value pcm::getPosSeqData(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	loadInPopSeqs();
	std::string popSeqDataUid = projectName_ + "_seqs";
	return cppcmsJsonToJson(seqs_->getJson(popSeqDataUid));
}

Json::Value pcm::getPopInfo() {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);

	auto ret = cppcmsJsonToJson(tableToJsonRowWise(popTable_, "h_popUID", VecStr { "h_ReadCnt" },
			VecStr { "p_TotalInputReadCnt", "p_TotalInputClusterCnt",
					"p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanMoi",
					"p_medianMoi", "p_minMoi", "p_maxMoi" }));
	return ret;
}
Json::Value pcm::getPopSeqDataForSamps(const VecStr & sampNames){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"sampNames", vectorToString(sampNames, ",")}}), std::cout, debug_);
	loadInPopSeqs();
	std::string popSeqDataUid = projectName_ + "_seqs";
	auto sampNameToCodedNum = sampNameToCodedNum_;
	auto containsSampName = [&sampNames, &sampNameToCodedNum](const std::string & str) {
		return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
	};
	auto trimedTab = sampTable_.extractByComp("s_Name", containsSampName);
	Json::Value ret = cppcmsJsonToJson(tableToJsonRowWise(trimedTab,"s_Name", VecStr{}));
	auto popNames = trimedTab.getColumn("h_popUID");
	removeDuplicates(popNames);
	auto currentReads = readVec::getSeqsWithNames(*seqs_->getRecord(popSeqDataUid),popNames);
	std::string popSeqDataUidCurrent = projectName_ + "_seqs_selected";
	seqs_->updateAddCache(popSeqDataUidCurrent, std::make_shared<std::vector<readObject>>(std::vector<readObject>(currentReads)));
	return cppcmsJsonToJson(seqs_->getJson(popSeqDataUidCurrent));
}

Json::Value pcm::getPopInfoForSamps(const VecStr & sampNames){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"sampNames", vectorToString(sampNames, ",")}}), std::cout, debug_);
	auto sampNameToCodedNum = sampNameToCodedNum_;
	auto containsSampName = [&sampNames, &sampNameToCodedNum](const std::string & str) {
		return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
	};
	auto trimedTab = sampTable_.extractByComp("s_Name", containsSampName);
	Json::Value ret = cppcmsJsonToJson(tableToJsonRowWise(trimedTab,"s_Name", VecStr{}));
	auto popNames = trimedTab.getColumn("h_popUID");
	removeDuplicates(popNames);
	auto trimedPopTab = popTable_.extractByComp("h_popUID", [&popNames](const std::string & str){return bib::in(str, popNames);});
	return cppcmsJsonToJson(tableToJsonRowWise(trimedPopTab, "h_popUID", VecStr { "h_ReadCnt" },
			VecStr { "p_TotalInputReadCnt", "p_TotalInputClusterCnt",
					"p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanMoi",
					"p_medianMoi", "p_minMoi", "p_maxMoi" }));
}






std::string pcm::decodeSampEncoding(const std::string& sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"sampName", sampName}} ), std::cout, debug_);
	return codedNumToSampName_[bib::lexical_cast<uint32_t>(sampName)];
}



Json::Value pcm::getSeqData(std::string encodeSampName) {
	std::string deCodedSampName = decodeSampEncoding(encodeSampName);
	bib::scopedMessage mess(
			messStrFactory(__PRETTY_FUNCTION__, { { "encodeSampName",
					encodeSampName }, { "deCodedSampName", deCodedSampName } }),
			std::cout, debug_);
	if (bib::in(deCodedSampName, clusteredSampleNames_)) {
		std::string fileName = "";
		/*auto files = bib::files::listAllFiles(
				bib::appendAsNeededRet(mainDir_, "/") + "final/", false,
				{ std::regex(deCodedSampName + R"(\.fast.*)") });*/
		auto files = bib::files::listAllFiles(
				bib::appendAsNeededRet(mainDir_, "/") + "final/", false,
				VecStr { deCodedSampName + ".fast" });
		if (files.empty()) {
			std::cerr << "No files starting with "
					<< bib::appendAsNeededRet(mainDir_, "/") + "final/" + deCodedSampName
					<< " found in " << bib::appendAsNeededRet(mainDir_, "/") + "final/"
					<< std::endl;
		} else {
			fileName = files.begin()->first.string();
			std::string seqUid = projectName_ + "_" + deCodedSampName;
			if (!seqs_->containsRecord(seqUid) || !seqs_->recordValid(seqUid)) {
				SeqIOOptions sampOpts;
				sampOpts.inFormat_ = SeqIOOptions::getInFormat(bib::files::getExtension(fileName));
				sampOpts.firstName_ = fileName;
				sampOpts.processed_ = true;
				SeqInput reader(sampOpts);
				reader.openIn();
				seqs_->updateAddCache(seqUid,
											std::make_shared<std::vector<readObject>>(reader.readAllReads<readObject>()));
			}
			return cppcmsJsonToJson(seqs_->getJson(seqUid));
		}
		if (fexists(fileName)) {

		} else {
			std::cerr << "File " << fileName << " does not exist" << std::endl;
		}
	} else {
		std::cerr << "SampName: " << deCodedSampName << " not found" << std::endl;
	}
	return Json::Value();
}


Json::Value pcm::getIndexExtractionInfo(){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__), std::cout, debug_);
	if(extractInfo_.allStatsTab_.content_.empty()){
		if(debug_){
			std::cerr <<  __PRETTY_FUNCTION__ << ": Extraction Info not loaded" << std::endl;
		}
	}else{

		auto statsCopy = extractInfo_.allStatsTab_;
		statsCopy.trimElementsAtFirstOccurenceOf("(");
		auto ret = tableToJsonRowWise(statsCopy, "IndexName", VecStr{"SamllFragments(len<50)","failedForwardPrimer", "failedQualityFiltering", "contamination"});
		return cppcmsJsonToJson(ret);
	}
	return Json::Value();
}

Json::Value pcm::getSampleExtractionInfo(const VecStr & sampNames){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"sampNames", vectorToString(sampNames, ",")}}), std::cout, debug_);
	if(extractInfo_.allStatsTab_.content_.empty()){
		if(debug_){
			std::cerr <<  __PRETTY_FUNCTION__ << ": Extraction Info not loaded" << std::endl;
		}
	}else{
		auto sampTabCopy = extractInfo_.profileBySampTab_;
		sampTabCopy.trimElementsAtFirstOccurenceOf("(");
		auto sampNameToCodedNum = sampNameToCodedNum_;
		auto containsSampName = [&sampNames, &sampNameToCodedNum](const std::string & str) {
			return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
		};
		sampTabCopy = sampTabCopy.extractByComp("Sample", containsSampName);
		for(auto & row : sampTabCopy.content_){
			row[sampTabCopy.getColPos("Sample")] = row[sampTabCopy.getColPos("Sample")] + "_" +  row[sampTabCopy.getColPos("MidName")];
		}
		auto ret = tableToJsonRowWise(sampTabCopy, "Sample", VecStr{});
		return cppcmsJsonToJson(ret);
	}
	return Json::Value();
}


Json::Value pcm::getGroupPopInfoForSamps(std::string group,
		std::string subGroup, const VecStr & sampNames) {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {
			{ "group", group }, { "subGroup", subGroup }, { "sampNames",
					vectorToString(sampNames, ",") } }), std::cout, debug_);
	if (setUpGroup(group, subGroup)) {
		Json::Value ret;
		auto sampNameToCodedNum = sampNameToCodedNum_;
		auto containsSampName =
				[&sampNames, &sampNameToCodedNum](const std::string & str) {
					return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
				};
		auto trimedTab = groupInfos_[group][subGroup].sampTable_.extractByComp(
				"s_Name", containsSampName);
		trimedTab.sortTable("s_Sample", false);
		ret = cppcmsJsonToJson(
				tableToJsonRowWise(trimedTab, "s_Sample", VecStr { }));
		auto popNames = trimedTab.getColumn("h_popUID");
		removeDuplicates(popNames);
		auto trimedPopTab = groupInfos_[group][subGroup].popTable_.extractByComp(
				"h_popUID",
				[&popNames](const std::string & str) {return bib::in(str, popNames);});
		return cppcmsJsonToJson(
				tableToJsonRowWise(trimedPopTab, "h_popUID", VecStr { "h_ReadCnt" },
						VecStr { "p_TotalInputReadCnt", "p_TotalInputClusterCnt",
								"p_TotalPopulationSampCnt", "p_TotalHaplotypes",
								"p_TotalUniqueHaplotypes", "p_meanMoi", "p_medianMoi",
								"p_minMoi", "p_maxMoi" }));
	} else {
		std::cout << "group: " << group << ":" << subGroup
				<< " not found, redirecting" << std::endl;
	}
	return Json::Value();
}
Json::Value pcm::getGroupPopSeqDataForSamps(std::string group, std::string subGroup, const VecStr & sampNames){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}, {"subGroup", subGroup},{"sampNames", vectorToString(sampNames, ",")}	}), std::cout, debug_);
	if(setUpGroup(group, subGroup)){
		std::string seqUid = projectName_ + "_" + group + "_" + subGroup + "_seqs";
		if(!seqs_->containsRecord(seqUid) || ! seqs_->recordValid(seqUid)){
			seqs_->updateAddCache(seqUid, std::make_shared<std::vector<readObject>>(groupInfos_[group][subGroup].popReads_));
		}
		Json::Value ret;
		auto sampNameToCodedNum = sampNameToCodedNum_;
		auto containsSampName = [&sampNames, &sampNameToCodedNum](const std::string & str) {
			return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
		};
		auto trimedTab = groupInfos_[group][subGroup].sampTable_.extractByComp("s_Name", containsSampName);
		trimedTab.sortTable("s_Sample",false);
		ret = cppcmsJsonToJson(tableToJsonRowWise(trimedTab, "s_Sample", VecStr{}));
		auto popNames = trimedTab.getColumn("h_popUID");
		removeDuplicates(popNames);
		auto currentReads = readVec::getSeqsWithNames(*seqs_->getRecord(seqUid),popNames);
		std::string popSeqDataUidCurrent = seqUid + "_selected";
		seqs_->updateAddCache(popSeqDataUidCurrent, std::make_shared<std::vector<readObject>>(std::vector<readObject>(currentReads)));
		return cppcmsJsonToJson(seqs_->getJson(popSeqDataUidCurrent));
	}else{
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
	}
	return Json::Value();
}


Json::Value pcm::getGroupPopInfo(std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}, {"subGroup", subGroup}	}), std::cout, debug_);
	if(setUpGroup(group, subGroup)){
		auto ret = tableToJsonRowWise(groupInfos_[group][subGroup].popTable_, "h_popUID", VecStr{"h_ReadCnt"}, VecStr{"p_TotalInputReadCnt","p_TotalInputClusterCnt","p_TotalPopulationSampCnt","p_TotalHaplotypes","p_TotalUniqueHaplotypes","p_meanMoi","p_medianMoi","p_minMoi","p_maxMoi"});
		return cppcmsJsonToJson(ret);
	}else{
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
	}
	return Json::Value();
}
Json::Value pcm::getGroupPopSeqData(std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}, {"subGroup", subGroup}	}), std::cout, debug_);
	if(setUpGroup(group, subGroup)){
		std::string seqUid = projectName_ + "_" + group + "_" + subGroup;
		if(!seqs_->containsRecord(seqUid) || ! seqs_->recordValid(seqUid)){
			seqs_->updateAddCache(seqUid, std::make_shared<std::vector<readObject>>(groupInfos_[group][subGroup].popReads_));
		}
		return cppcmsJsonToJson(seqs_->getJson(seqUid));
	}else{
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
	}
	return Json::Value();
}

Json::Value pcm::getGroupSampInfo(std::string group, std::string subGroup, const VecStr & sampNames){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {
			{ "group", group }, { "subGroup", subGroup }, { "sampNames",
					vectorToString(sampNames) } }), std::cout, debug_);
	if (setUpGroup(group, subGroup)) {
		Json::Value ret;
		auto sampNameToCodedNum = sampNameToCodedNum_;
		auto containsSampName = [&sampNames, &sampNameToCodedNum](const std::string & str) {
			return bib::in(estd::to_string(sampNameToCodedNum[str]), sampNames);
		};
		auto trimedTab = groupInfos_[group][subGroup].sampTable_.extractByComp("s_Name", containsSampName);
		trimedTab.sortTable("s_Sample",false);
		ret = cppcmsJsonToJson(tableToJsonRowWise(trimedTab, "s_Sample", VecStr{}));
		auto popNames = trimedTab.getColumn("h_popUID");
		ret["popColors"] = bib::json::toJson(getColorsStrsForNamesMultiples(popNames));
		return ret;
	}else{
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
	}
	return Json::Value();
}

Json::Value pcm::getGroupSampleNames(std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}, {"subGroup", subGroup}	}), std::cout, debug_);
	if(setUpGroup(group, subGroup)){
		Json::Value ret = bib::json::toJson(groupInfos_[group][subGroup].clusteredSampleNames_);
		return ret;
	}else{
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
	}
	return Json::Value();
}
Json::Value pcm::getSubGroupsForGroup(std::string group){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}} ), std::cout, debug_);
	auto search = groupInfosDirNames_.find(group);
	if(search == groupInfosDirNames_.end()){
		std::cout << "No group: " << group << "\n";
		std::cout << "options: " << vectorToString(getVectorOfMapKeys(groupInfosDirNames_), ", ") << "\n";
		std::cout << "group: " << group << " not found" << std::endl;
	}else{
		auto keys = getVectorOfMapKeys(search->second);
		bib::sort(keys);
		Json::Value ret = bib::json::toJson(keys);
		return ret;
	}
	return Json::Value();
}



Json::Value pcm::getGroupPopInfos(std::string group){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}} ), std::cout, debug_);
	auto search = groupInfosDirNames_.find(group);
	Json::Value ret;

	if(search == groupInfosDirNames_.end()){
		std::cout << "No group: " << group << "\n";
		std::cout << "options: " << vectorToString(getVectorOfMapKeys(groupInfosDirNames_), ", ") << "\n";
		std::cout << "group: " << group << " not found, redirecting" << std::endl;
	}else{
		auto keys = getVectorOfMapKeys(search->second);
		bib::sort(keys);
		for(const auto & k : keys){
			setUpGroup(group, k);
		}
		VecStr groupInfoColNames { "g_GroupName", "p_TotalInputReadCnt",
				"p_TotalInputClusterCnt", "p_TotalPopulationSampCnt",
				"p_TotalHaplotypes", "p_TotalUniqueHaplotypes", "p_meanMoi", "p_medianMoi", "p_minMoi",
				"p_maxMoi", "g_hapsFoundOnlyInThisGroup"};
		table outTab(groupInfoColNames);
		for(const auto & k : keys){
			outTab.rbind(groupInfos_[group][k].popTable_.getColumns(groupInfoColNames), false);
		}
		outTab = outTab.getUniqueRows();
		outTab.sortTable("g_GroupName", false);
		ret = cppcmsJsonToJson(tableToJsonRowWise(outTab, "g_GroupName", VecStr{}));
	}
	return ret;
}

bool pcm::setUpGroup(std::string group, std::string subGroup){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__ ,{{"group", group}, {"subGroup", subGroup}	}), std::cout, debug_);
	auto search = groupInfosDirNames_.find(group);
	if(search == groupInfosDirNames_.end()){
		std::cout << "No group: " << group << "\n";
		std::cout << "options: " << vectorToString(getVectorOfMapKeys(groupInfosDirNames_), ", ") << "\n";
		return false;
	}else{
		auto subSearch = search->second.find(subGroup);
		if(subSearch == search->second.end()){
			std::cout << "No subgroup: " << subGroup << " in group: " << group << "\n";
			std::cout << "options: " << vectorToString(getVectorOfMapKeys(search->second), ", ") << "\n";
			return false;
		}else{
			if(groupInfos_[group].find(subGroup) == groupInfos_[group].end()){
				/**@todo should check to see if file exists*/
				table sampTab = table(subSearch->second + "/sampFile.tab.txt", "\t", true);
				table popTab = table(subSearch->second + "/popFile.tab.txt", "\t", true);
				auto popNames = popTab.getColumn("h_popUID");
				std::vector<readObject> currentPopReads;

				loadInPopSeqs();
				std::string popSeqDataUid = projectName_ + "_seqs";
				auto popReadsPtr = seqs_->getRecord(popSeqDataUid);
				std::vector<readObject> & popReads = *popReadsPtr;
				for(const auto & read : popReads){
					if(bib::in(read.seqBase_.name_, popNames)){
						currentPopReads.emplace_back(read);
					}
				}
				groupInfos_[group][subGroup] = popInfo(sampTab, popTab, currentPopReads);
			}
			return true;
		}
	}
}



} /* namespace bibseq */
