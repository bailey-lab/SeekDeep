#pragma once
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






namespace bfs = boost::filesystem;

class pcm {
private:


	struct popInfo {
		popInfo(){}
		popInfo(const table & sampTable, const table & popTable,
				const std::vector<readObject> & popReads) :
				sampTable_(sampTable), popTable_(popTable), popReads_(popReads){
			clusteredSampleNames_ = sampTable_.getColumnLevels("s_Name");
			bib::sort(clusteredSampleNames_);
		}
		table sampTable_;
		table popTable_;
		std::vector<readObject> popReads_;
		VecStr clusteredSampleNames_;
	};

	table sampTable_;
	table popTable_;
	VecStr clusteredSampleNames_;
	double fracCutOff_ = 0;

	std::unordered_map<std::string, std::unordered_map<std::string, popInfo>> groupInfos_;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> groupInfosDirNames_;

	VecStr allSampleNames_;
	std::unordered_map<std::string, uint32_t> sampNameToCodedNum_;
	std::unordered_map<uint32_t, std::string> codedNumToSampName_;

	std::string mainDir_;
	std::string projectName_;
	std::string shortName_;

	std::string extractionDir_;

	ExtractionInfo extractInfo_;
	Json::Value config_;
	std::shared_ptr<seqCache> seqs_;

	bool debug_ = false;



	friend class pcv;
public:
	pcm(Json::Value config, std::shared_ptr<seqCache> seqs);



	void loadInPopSeqs();

	void setFracCutOff(double fracCutOff);
	double getFracCutOff();



	//json
	Json::Value getProjectName();
	Json::Value getSampleNames();
	Json::Value getAllSampleNames();
	Json::Value getSampleNamesEncoding();
	Json::Value getEncodingForSampleNames();

	std::string decodeSampEncoding(const std::string& sampName);
	Json::Value getSampInfo(const VecStr & sampNames);
	Json::Value getPosSeqData();
	Json::Value getPopInfo();
	Json::Value getPopSeqDataForSamps(const VecStr & sampNames);
	Json::Value getPopInfoForSamps(const VecStr & sampNames);



	Json::Value getIndexExtractionInfo();
	Json::Value getSampleExtractionInfo(const VecStr & sampNames);

	Json::Value getSeqData(std::string sampName);


	//group info
	Json::Value getGroupPopInfo(std::string group, std::string subGroup);
	Json::Value getGroupPopSeqData(std::string group, std::string subGroup);
	Json::Value getGroupPopInfoForSamps(std::string group, std::string subGroup, const VecStr & sampNames);
	Json::Value getGroupPopSeqDataForSamps(std::string group, std::string subGroup, const VecStr & sampNames);

	Json::Value getGroupSampInfo(std::string group, std::string subGroup, const VecStr & sampNames);
	Json::Value getGroupSampleNames(std::string group, std::string subGroup);
	Json::Value getSubGroupsForGroup(std::string group);
	Json::Value getGroupNames();
	bool setUpGroup(std::string group, std::string subGroup);
	Json::Value getGroupPopInfos(std::string group);

	std::string messStrFactory(const std::string & funcName)const;
	std::string messStrFactory(const std::string & funcName, const MapStrStr & args);
};




} /* namespace bibseq */
