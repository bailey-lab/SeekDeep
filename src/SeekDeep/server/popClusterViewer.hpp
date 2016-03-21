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
#include "popClusterModel.hpp"


namespace bibseq {

namespace bfs = boost::filesystem;

class pcv: public bibseq::seqApp {
private:


	std::map<std::string, std::string> config_;
	std::string configDir_;
	bfs::path serverResourceDir_;
	std::string rootName_;
	std::map<std::string, pcm> projectModels_;



public:
	pcv(cppcms::service& srv, std::map<std::string, std::string> config);

	virtual VecStr requiredOptions() const{
		return VecStr{"configDir", "resources"};
	}

	void jsPcv();
	void cssPcv();

	void projectsNames();

	void loadInProjects();

	void setFracCutOff(std::string shortProjName);
	void getFracCutOff(std::string shortProjName);

	//html
	//main page
	void mainPage();
	void mainProjectPage(std::string shortProjName);
	void individualSamplePage(std::string shortProjName,std::string sampName);

	void projectName(std::string shortProjName);
	//json
	void getSampleNames(std::string shortProjName);
	void getAllSampleNames(std::string shortProjName);
	void getSampleNamesEncoding(std::string shortProjName);

	void getEncodingForSampleNames(std::string shortProjName);

	void getSampInfo  				(std::string shortProjName);
	void getPosSeqData				(std::string shortProjName);
	void getPopInfo   				(std::string shortProjName);
	void getPosSeqDataForSamps(std::string shortProjName);
	void getPopInfoForSamps   (std::string shortProjName);



	void showExtractionInfo     (std::string shortProjName);
	void getIndexExtractionInfo (std::string shortProjName);
	void getSampleExtractionInfo(std::string shortProjName);


	void getSeqData             (std::string shortProjName, std::string sampName);

	//group info
	void getGroupPopInfo     			 (std::string shortProjName, std::string group, std::string subGroup);
	void getGroupPopSeqData  			 (std::string shortProjName, std::string group, std::string subGroup);
	void getGroupPopInfoForSamps	 (std::string shortProjName, std::string group, std::string subGroup);
	void getGroupPopSeqDataForSamps(std::string shortProjName, std::string group, std::string subGroup);
	void getGroupSampInfo    			 (std::string shortProjName, std::string group, std::string subGroup);
	void getGroupSampleNames 			 (std::string shortProjName, std::string group, std::string subGroup);

	void showGroupMainPage   (std::string shortProjName, std::string group, std::string subGroup);
	void showSubGroupsPage   (std::string shortProjName, std::string group);
	void getSubGroupsForGroup(std::string shortProjName, std::string group);

	void getGroupNames(std::string shortProjName);


	void getGroupPopInfos(std::string shortProjName, std::string group);

	bool hasProject(const std::string shortProjName);

	virtual void main(std::string url);

};

int genProjectConfig(const bib::progutils::CmdArgs & inputCommands);
int popClusteringViewer(const bib::progutils::CmdArgs & inputCommands);


} /* namespace bibseq */
