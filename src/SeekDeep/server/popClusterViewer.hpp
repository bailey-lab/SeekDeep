#pragma once
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
 * popClusterViewer.hpp
 *
 *  Created on: Jan 13, 2015
 *      Author: nickhathaway
 */

#include <seqServer/apps/seqApp.hpp>
#include <seqServer/utils.h>
#include <bibcpp.h>



namespace bibseq {

template<typename T>
uint32_t getMismatches(const T & read1,
				const T & read2,
				aligner alignerObj, bool weightHomopolymers){
	alignerObj.alignVec(read1.seqBase_,read2.seqBase_, false);
	alignerObj.profilePrimerAlignment(read1.seqBase_, read2.seqBase_, weightHomopolymers);
	return alignerObj.comp_.hqMismatches_;
};

template<typename T>
Json::Value getMinMismatchTreeJson(const std::vector<T> & reads, aligner & alignerObj,
		uint32_t numThreads, bool weightHomopolymers,double hueStart, double hueStop,
    double lumStart, double lumStop,
    double satStart, double satStop, std::string backgroundColor){
	std::function<uint32_t(const T & ,
	  		const T &, aligner, bool)> misFun = getMismatches<T>;
		auto misDistances = getDistanceCopy(reads, numThreads, misFun,
				alignerObj, weightHomopolymers);
	  readDistGraph<uint32_t> graphMis(misDistances, reads);
		std::vector<std::string> names;
	  for(const auto & n : graphMis.nodes_){
	  	names.emplace_back(n->name_);
	  }
	  if(hueStop == 360 && hueStart == 0){
	  	hueStop = 360 - (360.0/names.size());
	  }
		auto nameColors = bib::getColorsForNames(names, hueStart, hueStop,
				lumStart, lumStop, satStop, satStart);
	  Json::Value graphJson = graphMis.toJsonMismatchGraphAll(bib::color(backgroundColor), nameColors);
	  return graphJson;
}

struct ExtractionInfo {
	ExtractionInfo(const table & allStatsTab, const table & allprofileTab,
			const table & profileBySampTab) :
			allStatsTab_(allStatsTab), allProfileTab_(allprofileTab), profileBySampTab_(
					profileBySampTab) {

	}
	ExtractionInfo(){}

	table allStatsTab_;
	table allProfileTab_;
	table profileBySampTab_;
};

ExtractionInfo collectExtractionInfo(const std::string & dirName, const std::string & indexToDir, const std::string & sampNames);


namespace bfs = boost::filesystem;

class pcv: public bibseq::seqApp {
private:


	struct popInfo {
		popInfo(){}
		popInfo(const table & sampTable, const table & popTable,
				const std::vector<readObject> & popReads,
				const std::vector<readObject> & popReadsTranslated) :
				sampTable_(sampTable), popTable_(popTable), popReads_(popReads), popReadsTranslated_(
						popReadsTranslated) {
			clusteredSampleNames_ = sampTable_.getColumnLevels("s_Name");
			bib::sort(clusteredSampleNames_);
		}
		table sampTable_;
		table popTable_;
		std::vector<readObject> popReads_;
		std::vector<readObject> popReadsTranslated_;
		VecStr clusteredSampleNames_;
	};

	table sampTable_;
	table popTable_;
	VecStr clusteredSampleNames_;
	double fracCutOff_ = 0;

	std::unordered_map<std::string, std::unordered_map<std::string, popInfo>> groupInfos_;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> groupInfosDirNames_;

	VecStr allSampleNames_;
	std::string projectName_;

	std::string rootName_;
	std::string mainDir_;


	ExtractionInfo extractInfo_;

	std::map<std::string, std::string> config_;


	std::unordered_map<std::string, Json::Value> sampleMinTreeDataCache_;
	Json::Value minTreeData_;
	bool calculatedTreeData_ = false;

	std::unordered_map<std::string, uint32_t> sampNameToCodedNum_;
	std::unordered_map<uint32_t, std::string> codedNumToSampName_;

	static bfs::path make_path(const bfs::path fn) {
		return fn;
	}

public:
	pcv(cppcms::service& srv, std::map<std::string, std::string> config);

	virtual VecStr requiredOptions() const{
		return VecStr{"mainDir", "resources", "projectName"};
	}

	void loadInPopSeqs();

	void setFracCutOff();
	void getFracCutOff();

	//html
	//main page
	void mainPage();
	void individualSamplePage(std::string sampName);

	//json
	void getProjectName();
	void getSampleNames();
	void getAllSampleNames();
	void getSampleNamesEncoding();

	void getEncodingForSampleNames();

	std::string decodeSampEncoding(const std::string& sampName);
	void getSampInfo(std::string sampNames);
	void getPosSeqData();
	void getPopInfo();

	void getPopProtenData();


	void showExtractionInfo();
	void getIndexExtractionInfo();
	void getSampleExtractionInfo(std::string sampNames);


	void showMinTree();
	void getMinTreeData();

	void getSeqData(std::string sampName);
	void getProteinData(std::string sampName);

	void showMinTreeForSample(std::string sampName);
	void getMinTreeDataForSample(std::string sampName);

	//group info
	void getGroupPopInfo(std::string group, std::string subGroup);
	void getGroupPopSeqData(std::string group, std::string subGroup);
	void getGroupPopProtenData(std::string group, std::string subGroup);
	void getGroupSampInfo(std::string group, std::string subGroup, std::string sampName);
	void getGroupSampleNames(std::string group, std::string subGroup);

	void showGroupMainPage(std::string group, std::string subGroup);
	void showSubGroupsPage(std::string group);
	void getSubGroupsForGroup(std::string group);

	void getGroupNames();
	bool setUpGroup(std::string group, std::string subGroup);


	void getGroupPopInfos(std::string group);

};

int popClusteringViewer(std::map<std::string, std::string> inputCommands);


} /* namespace bibseq */
