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
ExtractionInfo collectExtractionInfoDirectName(const std::string & dirName,
		const std::string & indexToDir, const std::string & sampNames) {
	std::string nameDelim = "_extractor_";
	VecStr dirs;
	table indexNamesTab(indexToDir, "\t", true);
	table mainTableExtractionProfile;
	table mainTableExtractionStats;
	std::unordered_map<std::string, std::string> nameToIndex;
	for (const auto & row : indexNamesTab.content_) {
		if (row.size() != 2) {
			std::cerr << "Error in parsing " << indexToDir
					<< ", should have two columns, not " << row.size() << std::endl;
			std::cerr << vectorToString(row, "\t") << std::endl;
			exit(1);
		}
		dirs.emplace_back(row[1]);
		nameToIndex[row[1]] = row[0];
	}
	uint32_t count = 0;
	VecStr oldProfileColNames;

	VecStr oldStatsColNames;
	table inTab(sampNames, "whitespace", false);
	//goes sample name -> pairs of midname - index name
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> sampleDirWithSubDirs;
	for(const auto & rowPos : iter::range(inTab.content_.size())){
		const auto & row = inTab.content_[rowPos];
		if(row.empty() || row[0].front() == '#'){
			continue;
		}
		if (row.size() < 3) {
			throw std::runtime_error { bib::err::F()
					<< "setUpSampleDirs: rows should have at least 3 columns, row: "
					<< rowPos << "has " << row.size() };
		}
		for(const auto & colPos : iter::range<uint32_t>(2,row.size())){
			if(row[colPos] == "" || allWhiteSpaceStr(row[colPos])){
				continue;
			}
			//Discrepancy  between people naming MID01 and MID1 so replacing MID0 with MID
			sampleDirWithSubDirs[row[1]].emplace_back(replaceString(row[colPos], "MID0", "MID"), row[0]);
		}
	}

	//extraction info by index by mid;
	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> extractionInfo;

	for (const auto &dir : dirs) {
		auto fullDirName = bib::files::appendAsNeededRet(dirName, "/") + dir;
		table inProfileTab(fullDirName + "/extractionProfile.tab.txt", "\t", true);
		//i have accidentally added one more tab than is needed to
		//extractionProfile which makes it so all rows have an empty
		//at the end that correspond to any column names
		for (auto & row : inProfileTab.content_) {
			if (row.back() == "") {
				row.erase(row.begin() + row.size() - 1);
			}
		}
		table inStatsTab(fullDirName + "/extractionStats.tab.txt", "\t", true);
		std::string indexName = nameToIndex[dir];
		if (count == 0) {
			oldStatsColNames = inStatsTab.columnNames_;
			oldProfileColNames = inProfileTab.columnNames_;
		}
		inStatsTab.addColumn(VecStr { indexName }, "IndexName");
		inProfileTab.addColumn(VecStr { indexName }, "IndexName");
		for(const auto & row : inProfileTab.content_){
			auto midName = row[inProfileTab.getColPos("name")];
			auto midPos = midName.rfind("MID");
			midName = replaceString(midName.substr(midPos), "MID0", "MID");
			extractionInfo[row[inProfileTab.getColPos("IndexName")]][midName] = row;
		}
		if (count == 0) {
			mainTableExtractionStats = inStatsTab;
			mainTableExtractionProfile = inProfileTab;
		} else {
			mainTableExtractionStats.rbind(inStatsTab);
			mainTableExtractionProfile.rbind(inProfileTab);
		}
		++count;
	}

	table outSampleInfo(catenateVectors(VecStr{"Sample"}, mainTableExtractionProfile.columnNames_));
	for(const auto & samp : sampleDirWithSubDirs){
		for(const auto & indexMidNames : samp.second){
			auto addingRow = extractionInfo[indexMidNames.second][indexMidNames.first];
			if(addingRow.empty()){
				addingRow = std::vector<std::string>(mainTableExtractionProfile.content_.front().size(), "0");
			}
			addingRow[mainTableExtractionProfile.getColPos("IndexName")] = indexMidNames.second;
			addingRow[mainTableExtractionProfile.getColPos("name")] = indexMidNames.first;
			outSampleInfo.content_.emplace_back(catenateVectors(VecStr{samp.first}, addingRow));
		}
	}

	auto outSampleColName = catenateVectors(VecStr{"Sample"}, catenateVectors(VecStr{"IndexName"}, oldProfileColNames));
	auto profileColNames = catenateVectors(VecStr{"IndexName"}, oldProfileColNames);
	auto statsColName = catenateVectors(VecStr{"IndexName"}, oldStatsColNames);


	outSampleInfo = outSampleInfo.getColumns(outSampleColName);
	outSampleInfo.sortTable("Sample", false);
	mainTableExtractionProfile = mainTableExtractionProfile.getColumns(profileColNames);
	mainTableExtractionStats = mainTableExtractionStats.getColumns(statsColName);
	outSampleInfo.columnNames_[outSampleInfo.getColPos("name")] = "MidName";
	outSampleInfo.setColNamePositions();
	return ExtractionInfo(mainTableExtractionStats, mainTableExtractionProfile, outSampleInfo);
}
ExtractionInfo collectExtractionInfo(const std::string & dirName, const std::string & indexToDir, const std::string & sampNames){
	std::string nameDelim = "_extractor_";
	auto dirs = getNewestDirs(dirName, nameDelim);
	table indexNamesTab(indexToDir,"\t", true);
  table mainTableExtractionProfile;
  table mainTableExtractionStats;
	std::unordered_map<std::string, std::string> nameToIndex;
	for (const auto & row : indexNamesTab.content_) {
		if (row.size() != 2) {
			std::cerr << "Error in parsing " << indexToDir
					<< ", should have two columns, not " << row.size() << std::endl;
			std::cerr << vectorToString(row, "\t") << std::endl;
			exit(1);
		}
		auto lastPeriod = row[1].rfind(".");
		nameToIndex[row[1].substr(0, lastPeriod)] = row[0];
	}
	uint32_t count = 0;
	VecStr oldProfileColNames;

	VecStr oldStatsColNames;
	table inTab(sampNames, "whitespace", false);
	//goes sample name -> pairs of midname - index name
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> sampleDirWithSubDirs;
	for(const auto & rowPos : iter::range(inTab.content_.size())){
		const auto & row = inTab.content_[rowPos];
		if(row.empty() || row[0].front() == '#'){
			continue;
		}
		if (row.size() < 3) {
			throw std::runtime_error { bib::err::F()
					<< "setUpSampleDirs: rows should have at least 3 columns, row: "
					<< rowPos << "has " << row.size() };
		}
		for(const auto & colPos : iter::range<uint32_t>(2,row.size())){
			if(row[colPos] == "" || allWhiteSpaceStr(row[colPos])){
				continue;
			}
			//Discrepancy  between people naming MID01 and MID1 so replacing MID0 with MID
			sampleDirWithSubDirs[row[1]].emplace_back(replaceString(row[colPos], "MID0", "MID"), row[0]);
		}
	}

	//extraction info by index by mid;
	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> extractionInfo;

	for (const auto &dir : dirs) {
		table inProfileTab(dir + "/extractionProfile.tab.txt", "\t", true);
		//i have accidentally added one more tab than is needed to
		//extractionProfile which makes it so all rows have an empty
		//at the end that correspond to any column names
		for (auto & row : inProfileTab.content_) {
			if (row.back() == "") {
				row.erase(row.begin() + row.size() - 1);
			}
		}
		table inStatsTab(dir + "/extractionStats.tab.txt", "\t", true);
		std::string dirName = bib::files::getFileName(dir);
		auto pos = dirName.find(nameDelim);
		std::string indexName = nameToIndex[dirName.substr(0, pos)];
		if (count == 0) {
			oldStatsColNames = inStatsTab.columnNames_;
			oldProfileColNames = inProfileTab.columnNames_;
		}
		inStatsTab.addColumn(VecStr { indexName }, "IndexName");
		inProfileTab.addColumn(VecStr { indexName }, "IndexName");
		for(const auto & row : inProfileTab.content_){
			auto midName = row[inProfileTab.getColPos("name")];
			auto midPos = midName.rfind("MID");
			midName = replaceString(midName.substr(midPos), "MID0", "MID");
			extractionInfo[row[inProfileTab.getColPos("IndexName")]][midName] = row;
		}
		if (count == 0) {
			mainTableExtractionStats = inStatsTab;
			mainTableExtractionProfile = inProfileTab;
		} else {
			mainTableExtractionStats.rbind(inStatsTab);
			mainTableExtractionProfile.rbind(inProfileTab);
		}
		++count;
	}

	table outSampleInfo(catenateVectors(VecStr{"Sample"}, mainTableExtractionProfile.columnNames_));
	for(const auto & samp : sampleDirWithSubDirs){
		for(const auto & indexMidNames : samp.second){
			auto addingRow = extractionInfo[indexMidNames.second][indexMidNames.first];
			if(addingRow.empty()){
				addingRow = std::vector<std::string>(mainTableExtractionProfile.content_.front().size(), "0");
			}
			addingRow[mainTableExtractionProfile.getColPos("IndexName")] = indexMidNames.second;
			addingRow[mainTableExtractionProfile.getColPos("name")] = indexMidNames.first;
			outSampleInfo.content_.emplace_back(catenateVectors(VecStr{samp.first}, addingRow));
		}
	}

	auto outSampleColName = catenateVectors(VecStr{"Sample"}, catenateVectors(VecStr{"IndexName"}, oldProfileColNames));
	auto profileColNames = catenateVectors(VecStr{"IndexName"}, oldProfileColNames);
	auto statsColName = catenateVectors(VecStr{"IndexName"}, oldStatsColNames);


	outSampleInfo = outSampleInfo.getColumns(outSampleColName);
	outSampleInfo.sortTable("Sample", false);
	mainTableExtractionProfile = mainTableExtractionProfile.getColumns(profileColNames);
	mainTableExtractionStats = mainTableExtractionStats.getColumns(statsColName);
	outSampleInfo.columnNames_[outSampleInfo.getColPos("name")] = "MidName";
	outSampleInfo.setColNamePositions();
	return ExtractionInfo(mainTableExtractionStats, mainTableExtractionProfile, outSampleInfo);
}


pcv::pcv(cppcms::service& srv, std::map<std::string, std::string> config) :
		bibseq::seqApp(srv, config), config_(config){
	bool pass = configTest(config, requiredOptions(), "pcv");
	if(pass){
		std::cout << "Passed config test" << std::endl;
	}else{
		std::cout << "Didn't pass config test " << std::endl;
	}

	rootName_ = config["name"];
	mainDir_ = config["mainDir"];
	projectName_ = config["projectName"];
	pages_.emplace("mainPageHtml",
			make_path(config["resources"] + "pcv/mainPage.html"));
	pages_.emplace("redirectPage",
			make_path(config["resources"] + "html/redirectPage.html"));
	pages_.emplace("individualSample",
			make_path(config["resources"] + "pcv/individualSample.html"));
	pages_.emplace("extractionStats",
			make_path(config["resources"] + "pcv/extractionStats.html"));
	pages_.emplace("subGroupsPage",
			make_path(config["resources"] + "pcv/subGroupsPage.html"));
	pages_.emplace("groupMainPage",
			make_path(config["resources"] + "pcv/groupMainPage.html"));
	for (auto & fCache : pages_) {
		fCache.second.replaceStr("/ssv", rootName_);
	}
	std::cout << "ExtractionDir: " << config["extractionDir"] << std::endl;
	if(config["extractionDir"] != ""){
		std::string nameDelim = "_extractor_";
		auto dirs = getNewestDirs(config["extractionDir"], nameDelim);
		for(const auto & dir : dirs){
			std::cout << "Dir: " << dir << std::endl;
		}
		extractInfo_ = collectExtractionInfoDirectName(config["extractionDir"],config["indexToDir"], config["sampNames"]);
	}



	//main page
	dispMapRoot(&pcv::mainPage, this);

	dispMap(&pcv::getProjectName, this, "projectName");

	dispMap(&pcv::getPopInfo, this, "popInfo");
	dispMap(&pcv::getPosSeqData, this, "popSeqData");
	dispMap(&pcv::getPopProtenData, this, "popProteinData");
	dispMap_1arg(&pcv::getSampInfo, this, "sampInfo", "(\\w+)");
	dispMap(&pcv::getSampleNames, this, "sampleNames");

	//group info
	dispMap(&pcv::getGroupNames, this, "getGroupNames");
	dispMap_2arg(&pcv::getGroupPopInfo, this, "groupPopInfo", "(\\w+)/(\\w+)");
	dispMap_2arg(&pcv::getGroupPopSeqData, this, "groupPopSeqData", "(\\w+)/(\\w+)");
	dispMap_2arg(&pcv::getGroupPopProtenData, this, "groupPopProteinData", "(\\w+)/(\\w+)");
	dispMap_3arg(&pcv::getGroupSampInfo, this, "groupSampInfo", "(\\w+)/(\\w+)/(\\w+)");
	dispMap_2arg(&pcv::getGroupSampleNames, this, "groupSampleNames", "(\\w+)/(\\w+)");

	dispMap_2arg(&pcv::showGroupMainPage, this, "showGroupMainPage", "(\\w+)/(\\w+)");
	dispMap_1arg(&pcv::showSubGroupsPage, this, "showSubGroupsPage", "(\\w+)");
	dispMap_1arg(&pcv::getSubGroupsForGroup, this, "getSubGroupsForGroup", "(\\w+)");
	dispMap_1arg(&pcv::getGroupPopInfos, this, "getGroupPopInfos", "(\\w+)");


	dispMap(&pcv::getFracCutOff, this, "getFracCutOff");
	dispMap(&pcv::setFracCutOff, this, "setFracCutOff");

	dispMap(&pcv::getAllSampleNames, this, "allSampleNames");

	dispMap(&pcv::getSampleNamesEncoding, this, "sampleNamesEncoding");
	dispMap(&pcv::getEncodingForSampleNames, this, "encodingForSampleNames");


	dispMap(&pcv::showMinTree, this, "showMinTree");
	dispMap(&pcv::getMinTreeData, this, "minTreeData");


	dispMap(&pcv::showExtractionInfo, this, "showExtractionInfo");
	dispMap(&pcv::getIndexExtractionInfo, this, "getIndexExtractionInfo");
	dispMap_1arg(&pcv::getSampleExtractionInfo, this, "getSampleExtractionInfo", "(\\w+)");


	dispMap_1arg(&pcv::individualSamplePage, this, "individualSamplePage", "(\\w+)");
	dispMap_1arg(&pcv::getSeqData, this, "seqData", "(\\w+)");
	dispMap_1arg(&pcv::getProteinData, this, "proteinData", "(\\w+)");
	dispMap_1arg(&pcv::showMinTreeForSample, this, "showMinTreeForSample", "(\\w+)");
	dispMap_1arg(&pcv::getMinTreeDataForSample, this, "minTreeDataForSample", "(\\w+)");


	mapper().root(rootName_);

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

	std::cout << "Finished set up" << std::endl;
}


void pcv::setFracCutOff(){
	bib::scopedMessage mess("setFracCutOff",std::cout, true);
  auto postData = request().post();
  auto postJson = bib::json::toJson(postData);
  std::cout << postJson << std::endl;
  fracCutOff_ = bib::lexical_cast<double>(postJson["input"].asString());
  fracCutOff_ /=100;
  std::cout  << fracCutOff_ << std::endl;
}

void pcv::getFracCutOff(){
	ret_json();
	response().out() << roundDecPlaces((fracCutOff_ * 100),2);
}


void pcv::getGroupNames(){
	ret_json();
	cppcms::json::value ret = getVectorOfMapKeys(groupInfosDirNames_);
	response().out() << ret;
}

void pcv::mainPage(){
	auto search = pages_.find("mainPageHtml");
	response().out() << search->second.get("/ssv", rootName_);
}

//json
void pcv::getProjectName(){
	ret_json();
	cppcms::json::value ret = projectName_;
	response().out() << ret;
}

void pcv::getSampleNames(){
	ret_json();
	cppcms::json::value ret = clusteredSampleNames_;
	response().out() << ret;
}

void pcv::getAllSampleNames(){
	ret_json();
	cppcms::json::value ret = allSampleNames_;
	response().out() << ret;
}

void pcv::getSampleNamesEncoding(){
	ret_json();
	cppcms::json::value ret;
	for(const auto & sampName : sampNameToCodedNum_){
		ret[sampName.first] = sampName.second;
	}
	response().out() << ret;
}

void pcv::getEncodingForSampleNames(){
	ret_json();
	cppcms::json::value ret;
	for(const auto & sampName : sampNameToCodedNum_){
		ret[sampName.second] = sampName.first;
	}
	response().out() << ret;
}

void pcv::getSampInfo(std::string sampNames){
	ret_json();
	cppcms::json::value ret;

	auto sampToks = bibseq::tokenizeString(sampNames, "DELIM");
	auto sampNameToCodedNum = sampNameToCodedNum_;
	auto containsSampName = [&sampToks, &sampNameToCodedNum](const std::string & str) {
		return bib::in(estd::to_string(sampNameToCodedNum[str]), sampToks);
	};
	auto trimedTab = sampTable_.extractByComp("s_Name", containsSampName);
	ret = tableToJsonRowWise(trimedTab,"s_Name", VecStr{});
	auto popCounts = bibseq::countVec(trimedTab.getColumn("h_popUID"));
	auto popColors = bib::njhColors(popCounts.size());
	bibseq::VecStr popColorsStrs(popColors.size());
	uint32_t count = 0;
	uint32_t halfCount = 0;
	for(const auto & cPos : iter::range(popColors.size())) {
		uint32_t pos = 0;
		if(cPos %2 == 0) {
			pos = popColors.size()/2 + halfCount;
			++halfCount;
		} else {
			pos = count;
			++count;
		}
		popColorsStrs[cPos] = "#" + popColors[pos].hexStr_;
	}
	ret["popColors"] = popColorsStrs;
	response().out() << ret;
}

void pcv::loadInPopSeqs(){
	bib::scopedMessage mess(bib::err::F() << "loadInPopSeqs", std::cout, true);
	std::string popSeqDataUid = projectName_ + "_seqs";
	if(!seqs_.containsRecord(popSeqDataUid) || !seqs_.recordValid(popSeqDataUid)){
		/**todo make safter for non fastq pop clustering */
		//read pop seqs
		auto files = bib::files::listAllFiles(mainDir_ + "population/", false, {std::regex{".*\\.fastq"}});
		readObjectIO reader;
		reader.read("fastq", files.begin()->first.string(), true);
		for(auto & read : reader.reads){
			read.seqBase_.name_ = read.seqBase_.name_.substr(0, read.seqBase_.name_.rfind("_"));
		}
		if(seqs_.containsRecord(popSeqDataUid)){
			seqs_.updateCache(popSeqDataUid, std::make_shared<std::vector<readObject>>(std::vector<readObject>(reader.reads)));
		}else{
			seqs_.addToCache(popSeqDataUid, std::make_shared<std::vector<readObject>>(std::vector<readObject>(reader.reads)));
		}
	}
}
void pcv::getPosSeqData(){
	ret_json();
	loadInPopSeqs();
	std::string popSeqDataUid = projectName_ + "_seqs";
	response().out() << seqs_.getJson(popSeqDataUid);;
}

void pcv::getPopInfo() {
	ret_json();
	auto ret = tableToJsonRowWise(popTable_, "h_popUID", VecStr { "h_ReadCnt" },
			VecStr { "p_TotalInputReadCnt", "p_TotalInputClusterCnt",
					"p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanMoi",
					"p_medianMoi", "p_minMoi", "p_maxMoi" });
	response().out() << ret;
}



void pcv::getPopProtenData(){
	ret_json();
	bib::scopedMessage mess(bib::err::F() << "getPopProtenData", std::cout, true);
	std::string proteinUid = projectName_ + "_protein";
	if(!seqs_.containsRecord(proteinUid) || !seqs_.recordValid(proteinUid)){
		auto proteins = popTable_.getColumn("h_Protein");
		auto names = popTable_.getColumn("h_popUID");
		std::vector<readObject> popReads;
		for(auto pos : iter::range(names.size())){
			popReads.emplace_back(seqInfo(names[pos], proteins[pos]));
		}
		std::unordered_map<std::string, std::vector<uint32_t>> proteinCheck;
		for(const auto & pEnum : iter::enumerate(popReads)){
			proteinCheck[pEnum.element.seqBase_.seq_].emplace_back(pEnum.index);
		}
		if(seqs_.containsRecord(proteinUid)){
			seqs_.updateCache(proteinUid, std::make_shared<std::vector<readObject>>(popReads));
		}else{
			seqs_.addToCache(proteinUid, std::make_shared<std::vector<readObject>>(popReads));
		}
		/*for(const auto & pc : proteinCheck){
			std::cout << pc.first << " : " << vectorToString(pc.second, ",") << std::endl;
		}*/
	}

	response().out() << seqs_.getJson(proteinUid);;
}

void pcv::getMinTreeData(){
	ret_json();
	if(!calculatedTreeData_){
		loadInPopSeqs();
		std::string popSeqDataUid = projectName_ + "_seqs";
		//get min tree data
		uint64_t maxLength = 0;
		auto popReadsPtr = seqs_.getRecord(popSeqDataUid);
		std::vector<readObject> & popReads = *popReadsPtr;
		readVec::getMaxLength(popReads, maxLength);
		aligner alignerObj(maxLength, gapScoringParameters(5,1), substituteMatrix::createDegenScoreMatrix(2, -2));
	  std::function<uint32_t(const readObject & ,
	  		const readObject &, aligner, bool)> misFun = getMismatches<readObject>;
		auto misDistances = getDistanceCopy(popReads, 2, misFun,
				alignerObj, true);
	  readDistGraph<uint32_t> graphMis(misDistances, popReads);
		std::vector<std::string> popNames;
	  for(const auto & n : graphMis.nodes_){
	  	popNames.emplace_back(n->name_);
	  }
		auto popColors = bib::njhColors(popReads.size());
		bibseq::VecStr popColorsStrs(popColors.size());
		uint32_t count = 0;
		uint32_t halfCount = 0;
		for(const auto & cPos : iter::range(popColors.size())) {
			uint32_t pos = 0;
			if(cPos %2 == 0) {
				pos = popColors.size()/2 + halfCount;
				++halfCount;
			} else {
				pos = count;
				++count;
			}
			popColorsStrs[cPos] = "#" + popColors[pos].hexStr_;
		}
		std::unordered_map<std::string,bib::color> nameColors;
		for(auto pos : iter::range(popNames.size())){
			nameColors[popNames[pos]] = popColorsStrs[pos];
		}
		minTreeData_ = graphMis.toJsonMismatchGraphAll(bib::color("#000000"), nameColors);
	}
	response().out() << minTreeData_;
}

void pcv::showMinTree(){
	response().out() << genHtmlStrForPsuedoMintree(rootName_ + "/minTreeData");
}

std::string pcv::decodeSampEncoding(const std::string& sampName){
	return codedNumToSampName_[bib::lexical_cast<uint32_t>(sampName)];
}

void pcv::individualSamplePage(std::string sampName){
	sampName = decodeSampEncoding(sampName);
	bib::scopedMessage mess(bib::err::F() << "individualSamplePage; " << "sampName: " << sampName, std::cout, true);
	if(bib::in(sampName, clusteredSampleNames_)){
		auto search = pages_.find("individualSample");
		response().out() << search->second.get("/ssv", rootName_);
	}else{
		auto search = pages_.find("redirectPage");
		std::cout << "SampName: " << sampName << " not found, redirecting" << std::endl;
		response().out() << search->second.get("/ssv", rootName_);
	}
}

void pcv::getSeqData(std::string sampName){
	sampName = decodeSampEncoding(sampName);
	bib::scopedMessage mess(bib::err::F()<< "getSeqData; " << "sampName: " << sampName, std::cout, true);
	if(bib::in(sampName, clusteredSampleNames_)){
		auto fileName = bib::files::appendAsNeededRet(mainDir_, "/") + "final/" + sampName + ".fastq";
		if(fexists(fileName)){
			ret_json();
			std::string seqUid = projectName_ + "_" + sampName ;
			if(!seqs_.containsRecord(seqUid) || ! seqs_.recordValid(seqUid)){
				readObjectIO reader;
				reader.read("fastq", fileName, true);
				if(seqs_.containsRecord(seqUid)){
					seqs_.updateCache(seqUid, std::make_shared<std::vector<readObject>>(reader.reads));
				}else{
					seqs_.addToCache(seqUid, std::make_shared<std::vector<readObject>>(reader.reads));
				}
			}
			response().out() << seqs_.getJson(seqUid);
		}else{
			std::cout << "File " << fileName << " does not exist" << std::endl;
		}
	}else{
		std::cout << "SampName: " << sampName << " not found" << std::endl;
	}
}

void pcv::getProteinData(std::string sampName){
	sampName = decodeSampEncoding(sampName);
  bool forceStartM = false;
  bool transcribeToRNAFirst = false;
  uint64_t start = 0;
	bib::scopedMessage mess(bib::err::F()<< "getProteinData; " << "sampName: " << sampName, std::cout, true);
	if(bib::in(sampName, clusteredSampleNames_)){
		auto fileName = bib::files::appendAsNeededRet(mainDir_, "/") + "final/" + sampName + ".fastq";
		if(fexists(fileName)){
			ret_json();
			std::string seqUid = projectName_ + "_" + sampName + "_protein";
			if(!seqs_.containsRecord(seqUid) || ! seqs_.recordValid(seqUid)){
				readObjectIO reader;
				reader.read("fastq", fileName, true);
				readVec::convertReadsToProteinFromcDNA(reader.reads, transcribeToRNAFirst,
				                                           start, forceStartM);
				if(seqs_.containsRecord(seqUid)){
					seqs_.updateCache(seqUid, std::make_shared<std::vector<readObject>>(reader.reads));
				}else{
					seqs_.addToCache(seqUid, std::make_shared<std::vector<readObject>>(reader.reads));
				}
			}
			response().out() << seqs_.getJson(seqUid);
		}else{
			std::cout << "File " << fileName << " does not exist" << std::endl;
		}
	}else{
		std::cout << "SampName: " << sampName << " not found" << std::endl;
	}
}

void pcv::getMinTreeDataForSample(std::string sampName){
	sampName = decodeSampEncoding(sampName);
	bib::scopedMessage mess(bib::err::F()<< "getMinTreeDataForSample; " << "sampName: " << sampName, std::cout, true);
	if(bib::in(sampName, clusteredSampleNames_)){
		auto fileName = bib::files::appendAsNeededRet(mainDir_, "/") + "final/" + sampName + ".fastq";
		if(fexists(fileName)){
			auto search = sampleMinTreeDataCache_.find(sampName);
			Json::Value graphJsonData;
			if(search == sampleMinTreeDataCache_.end()){
				ret_json();
				readObjectIO reader;
				reader.read("fastq", fileName, true);
				uint64_t maxLength = 0;
				readVec::getMaxLength(reader.reads, maxLength);
				aligner alignerObj(maxLength, gapScoringParameters(5,1), substituteMatrix::createDegenScoreMatrix(2, -2));
			  std::function<uint32_t(const readObject & ,
			  		const readObject &, aligner, bool)> misFun = getMismatches<readObject>;
				auto misDistances = getDistanceCopy(reader.reads, 2, misFun,
						alignerObj, true);
			  readDistGraph<uint32_t> graphMis(misDistances, reader.reads);
				std::vector<std::string> popNames;
			  for(const auto & n : graphMis.nodes_){
			  	popNames.emplace_back(n->name_);
			  }
				auto popColors = bib::njhColors(reader.reads.size());
				bibseq::VecStr popColorsStrs(popColors.size());
				uint32_t count = 0;
				uint32_t halfCount = 0;
				for(const auto & cPos : iter::range(popColors.size())) {
					uint32_t pos = 0;
					if(cPos %2 == 0) {
						pos = popColors.size()/2 + halfCount;
						++halfCount;
					} else {
						pos = count;
						++count;
					}
					popColorsStrs[cPos] = "#" + popColors[pos].hexStr_;
				}
				std::unordered_map<std::string,bib::color> nameColors;
				for(auto pos : iter::range(popNames.size())){
					nameColors[popNames[pos]] = popColorsStrs[pos];
				}
				graphJsonData = graphMis.toJsonMismatchGraphAll(bib::color("#000000"), nameColors);
				sampleMinTreeDataCache_[sampName] = graphJsonData;
			}else{
				graphJsonData = search->second;
			}
			response().out() << graphJsonData;
		}else{
			std::cout << "File " << fileName << " does not exist" << std::endl;
		}
	}else{
		std::cout << "SampName: " << sampName << " not found" << std::endl;
	}
}
//
void pcv::showMinTreeForSample(std::string sampName){
	std::string encodedName = sampName;
	sampName = decodeSampEncoding(sampName);
	std::cout << "showMinTreeForSample; " << "sampName: " << sampName << std::endl;
	if(bib::in(sampName, clusteredSampleNames_)){
		response().out() << genHtmlStrForPsuedoMintree(rootName_ + "/minTreeDataForSample/" + encodedName);
	}else{
		auto search = pages_.find("redirectPage");
		std::cout << "SampName: " << sampName << " not found, redirecting" << std::endl;
		response().out() << search->second.get("/ssv", rootName_);
	}
}


void pcv::showExtractionInfo(){
	bib::scopedMessage mess(bib::err::F()<< "showExtractionInfo", std::cout, true);
	auto search = pages_.find("extractionStats");
	response().out() << search->second.get("/ssv", rootName_);
}

void pcv::getIndexExtractionInfo(){
	bib::scopedMessage mess(bib::err::F()<< "getIndexExtractionInfo", std::cout, true);
	if(extractInfo_.allStatsTab_.content_.empty()){
		std::cerr << "getIndexExtractionInfo: Extraction Info not loaded" << std::endl;
	}else{
		ret_json();
		auto statsCopy = extractInfo_.allStatsTab_;
		statsCopy.trimElementsAtFirstOccurenceOf("(");
		auto ret = tableToJsonRowWise(statsCopy, "IndexName", VecStr{"SamllFragments(len<50)","failedForwardPrimer", "failedQualityFiltering", "contamination"});
		response().out() << ret;
	}
}

void pcv::getSampleExtractionInfo(std::string sampNames){
	bib::scopedMessage mess(bib::err::F()<< "getSampleExtractionInfo", std::cout, true);
	if(extractInfo_.allStatsTab_.content_.empty()){
		std::cerr << "getSampleExtractionInfo: Extraction Info not loaded" << std::endl;
	}else{
		ret_json();
		auto sampTabCopy = extractInfo_.profileBySampTab_;
		sampTabCopy.trimElementsAtFirstOccurenceOf("(");

		auto sampToks = bibseq::tokenizeString(sampNames, "DELIM");
		printVector(sampToks);
		auto sampNameToCodedNum = sampNameToCodedNum_;
		auto containsSampName = [&sampToks, &sampNameToCodedNum](const std::string & str) {
			return bib::in(estd::to_string(sampNameToCodedNum[str]), sampToks);
		};

		sampTabCopy = sampTabCopy.extractByComp("Sample", containsSampName);
		//sampTabCopy.outPutContentOrganized(std::cout);

		for(auto & row : sampTabCopy.content_){
			row[sampTabCopy.getColPos("Sample")] = row[sampTabCopy.getColPos("Sample")] + "_" +  row[sampTabCopy.getColPos("MidName")];
		}
		auto ret = tableToJsonRowWise(sampTabCopy, "Sample", VecStr{});
		response().out() << ret;
	}
}

void pcv::getGroupPopInfo(std::string group, std::string subGroup){
	bib::scopedMessage mess(bib::err::F()<< "getGroupPopInfo", std::cout, true);
	if(setUpGroup(group, subGroup)){
		ret_json();
		auto ret = tableToJsonRowWise(groupInfos_[group][subGroup].popTable_, "h_popUID", VecStr{"h_ReadCnt"}, VecStr{"p_TotalInputReadCnt","p_TotalInputClusterCnt","p_TotalPopulationSampCnt","p_TotalHaplotypes","p_TotalUniqueHaplotypes","p_meanMoi","p_medianMoi","p_minMoi","p_maxMoi"});
		response().out() << ret;
	}else{
		//auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
		//response().out() << search->second.get("/ssv", rootName_);
	}
}
void pcv::getGroupPopSeqData(std::string group, std::string subGroup){
	bib::scopedMessage mess(bib::err::F()<< "getGroupPopSeqData", std::cout, true);
	if(setUpGroup(group, subGroup)){
		ret_json();
		std::string seqUid = projectName_ + "_" + group + "_" + subGroup;
		if(!seqs_.containsRecord(seqUid) || ! seqs_.recordValid(seqUid)){
			if(seqs_.containsRecord(seqUid)){
				seqs_.updateCache(seqUid, std::make_shared<std::vector<readObject>>(groupInfos_[group][subGroup].popReads_));
			}else{
				seqs_.addToCache(seqUid, std::make_shared<std::vector<readObject>>(groupInfos_[group][subGroup].popReads_));
			}
		}
		response().out() << seqs_.getJson(seqUid);
	}else{
		//auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
		//response().out() << search->second.get("/ssv", rootName_);
	}
}
void pcv::getGroupPopProtenData(std::string group, std::string subGroup){
	bib::scopedMessage mess(bib::err::F()<< "getGroupPopProtenData", std::cout, true);
	if(setUpGroup(group, subGroup)){
		ret_json();
		std::string seqUid = projectName_ + "_" + group + "_" + subGroup + "_protein";
		if(!seqs_.containsRecord(seqUid) || ! seqs_.recordValid(seqUid)){
			if(seqs_.containsRecord(seqUid)){
				seqs_.updateCache(seqUid, std::make_shared<std::vector<readObject>>(groupInfos_[group][subGroup].popReadsTranslated_));
			}else{
				seqs_.addToCache(seqUid, std::make_shared<std::vector<readObject>>(groupInfos_[group][subGroup].popReadsTranslated_));
			}
		}
		response().out() << seqs_.getJson(seqUid);
	}else{
		//auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
		//response().out() << search->second.get("/ssv", rootName_);
	}

}
void pcv::getGroupSampInfo(std::string group, std::string subGroup, std::string sampNames){
	bib::scopedMessage mess(bib::err::F()<< "getGroupSampInfo", std::cout, true);
	if(setUpGroup(group, subGroup)){
		ret_json();
		cppcms::json::value ret;

		auto sampToks = bibseq::tokenizeString(sampNames, "DELIM");
		auto sampNameToCodedNum = sampNameToCodedNum_;
		auto containsSampName = [&sampToks, &sampNameToCodedNum](const std::string & str) {
			return bib::in(estd::to_string(sampNameToCodedNum[str]), sampToks);
		};
		auto trimedTab = groupInfos_[group][subGroup].sampTable_.extractByComp("s_Name", containsSampName);
		trimedTab.sortTable("s_Sample",false);
		ret = tableToJsonRowWise(trimedTab, "s_Sample", VecStr{});
		auto popCounts = bibseq::countVec(trimedTab.getColumn("h_popUID"));
		auto popColors = bib::njhColors(popCounts.size());
		bibseq::VecStr popColorsStrs(popColors.size());
		uint32_t count = 0;
		uint32_t halfCount = 0;
		for(const auto & cPos : iter::range(popColors.size())) {
			uint32_t pos = 0;
			if(cPos %2 == 0) {
				pos = popColors.size()/2 + halfCount;
				++halfCount;
			} else {
				pos = count;
				++count;
			}
			popColorsStrs[cPos] = "#" + popColors[pos].hexStr_;
		}
		ret["popColors"] = popColorsStrs;
		response().out() << ret;
	}else{
		//auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
		//response().out() << search->second.get("/ssv", rootName_);
	}

}

void pcv::getGroupSampleNames(std::string group, std::string subGroup){
	bib::scopedMessage mess(bib::err::F()<< "getGroupSampleNames", std::cout, true);
	if(setUpGroup(group, subGroup)){
		ret_json();
		cppcms::json::value ret = groupInfos_[group][subGroup].clusteredSampleNames_;

		response().out() << ret;
	}else{
		//auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << ":" << subGroup<< " not found, redirecting" << std::endl;
		//response().out() << search->second.get("/ssv", rootName_);
	}
}
void pcv::getSubGroupsForGroup(std::string group){
	bib::scopedMessage mess(bib::err::F()<< "getSubGroupsForGroup", std::cout, true);
	std::cout << "getSubGroupsForGroup: GroupName: " << group  << "\n";
	auto search = groupInfosDirNames_.find(group);
	if(search == groupInfosDirNames_.end()){
		std::cout << "No group: " << group << "\n";
		std::cout << "options: " << vectorToString(getVectorOfMapKeys(groupInfosDirNames_), ", ") << "\n";
		std::cout << "group: " << group << " not found" << std::endl;
	}else{
		auto keys = getVectorOfMapKeys(search->second);
		bib::sort(keys);
		ret_json();
		cppcms::json::value ret = keys;
		response().out() << ret;
	}
}
void pcv::showGroupMainPage(std::string group, std::string subGroup){
	bib::scopedMessage mess(bib::err::F()<< "showGroupMainPage", std::cout, true);
	std::cout << "showGroupMainPage: GroupName: " << group << " subGroup: " << subGroup << "\n";
	auto search = groupInfosDirNames_.find(group);
	if(search == groupInfosDirNames_.end()){
		std::cout << "No group: " << group << "\n";
		std::cout << "options: " << vectorToString(getVectorOfMapKeys(groupInfosDirNames_), ", ") << "\n";
		auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << " not found, redirecting" << std::endl;
		response().out() << search->second.get("/ssv", rootName_);
	}else{
		auto subSearch = search->second.find(subGroup);
		if(subSearch == search->second.end()){
			std::cout << "No subgroup: " << subGroup << " in group: " << group << "\n";
			std::cout << "options: " << vectorToString(getVectorOfMapKeys(search->second), ", ") << "\n";
			auto search = pages_.find("redirectPage");
			std::cout << "group: " << group << " not found, redirecting" << std::endl;
			response().out() << search->second.get("/ssv", rootName_);
		}else{
			auto search = pages_.find("groupMainPage");
			response().out() << search->second.get("/ssv", rootName_);
		}
	}
}
void pcv::showSubGroupsPage(std::string group){
	bib::scopedMessage mess(bib::err::F()<< "showSubGroupsPage", std::cout, true);
	std::cout << "showSubGroupsPage: GroupName: " << group << "\n";
	auto search = groupInfosDirNames_.find(group);
	if(search == groupInfosDirNames_.end()){
		std::cout << "No group: " << group << "\n";
		std::cout << "options: " << vectorToString(getVectorOfMapKeys(groupInfosDirNames_), ", ") << "\n";
		auto search = pages_.find("redirectPage");
		std::cout << "group: " << group << " not found, redirecting" << std::endl;
		response().out() << search->second.get("/ssv", rootName_);
	}else{
		auto search = pages_.find("subGroupsPage");
		response().out() << search->second.get("/ssv", rootName_);
	}
}

void pcv::getGroupPopInfos(std::string group){
	bib::scopedMessage mess(bib::err::F()<< "getGroupPopInfos", std::cout, true);
	std::cout << "getGroupPopInfos: GroupName: " << group << "\n";
	auto search = groupInfosDirNames_.find(group);
	cppcms::json::value ret;
	ret_json();
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
			outTab.rbind(groupInfos_[group][k].popTable_.getColumns(groupInfoColNames));
		}
		outTab = outTab.getUniqueRows();
		outTab.sortTable("g_GroupName", false);
		ret = tableToJsonRowWise(outTab, "g_GroupName", VecStr{});
	}
	response().out() << ret;
}

bool pcv::setUpGroup(std::string group, std::string subGroup){
	bib::scopedMessage mess(bib::err::F()<< "setUpGroup", std::cout, true);
	std::cout << "setUpGroup: GroupName: " << group << " subGroup: " << subGroup << "\n";
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
				//get min tree data
				//uint64_t maxLength = 0;
				std::cout << "recordValid: " << bib::colorBool(seqs_.recordValid(popSeqDataUid));
				auto popReadsPtr = seqs_.getRecord(popSeqDataUid);
				std::vector<readObject> & popReads = *popReadsPtr;
				for(const auto & read : popReads){
					if(bib::in(read.seqBase_.name_, popNames)){
						currentPopReads.emplace_back(read);
					}
				}
				std::vector<readObject> currentPopTranslatedReads = currentPopReads;
				readVec::convertReadsToProteinFromcDNA(currentPopTranslatedReads, false);
				groupInfos_[group][subGroup] = popInfo(sampTab, popTab, currentPopReads, currentPopTranslatedReads);
			}
			return true;
		}
	}
}


int popClusteringViewer(std::map<std::string, std::string> inputCommands){
	bibseq::seqSetUp setUp(inputCommands);
	std::string mainDir = "";
	std::string extractioinDir = "";
	uint32_t port = 9881;
	std::string name = "pcv";
	std::string resourceDirName = "";
	std::string projectName = "";
	setUp.setOption(resourceDirName, "-resourceDirName", "Name of the resource Directory where the js and hmtl is located", true);
	if(resourceDirName.back() != '/'){
		resourceDirName.append("/");
	}
	setUp.setOption(mainDir, "-mainDir", "Name of the Master Result Directory", true);
	if(mainDir.back() != '/'){
		mainDir.append("/");
	}
	setUp.setOption(port, "-port", "Port Number to Serve On");
	setUp.setOption(name, "-name", "Name of root of the server");
	std::string extractionDir = "";
	std::string indexToDir = "";
	std::string sampNames = "";

	setUp.setOption(indexToDir, "-indexToDir", "File, first column is index name, second is the name of the file extraction was done on");
	setUp.setOption(extractionDir, "-extractionDir", "Name of the directory where extraction was done");
	setUp.setOption(sampNames, "-sampNames", "A file, first column is index name, second is sample name, the rest of the columns are MID names");

	setUp.setOption(projectName, "-projectName", "Name of the Project", true);
	setUp.finishSetUp(std::cout);
	name = "/" + name;
  auto config = server_config(name, port);
  //
  std::map<std::string, std::string> appConfig;
  appConfig["name"] = name;
  appConfig["mainDir"] = mainDir;
  appConfig["resources"] = resourceDirName;
  appConfig["js"] = resourceDirName + "js/";
  appConfig["css"] = resourceDirName + "css/";
  appConfig["projectName"] = projectName;
  appConfig["extractionDir"] = extractionDir;
  appConfig["indexToDir"] = indexToDir;
  appConfig["sampNames"] = sampNames;
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
