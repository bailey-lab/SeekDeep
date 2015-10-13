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
#include "popClusterModel.hpp"


namespace bibseq {

namespace bfs = boost::filesystem;

class pcv: public bibseq::seqApp {
private:


	std::map<std::string, std::string> config_;
	std::string configDir_;
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


	void showMinTree            (std::string shortProjName);
	void getMinTreeData         (std::string shortProjName);

	void getSeqData             (std::string shortProjName, std::string sampName);

	void showMinTreeForSample   (std::string shortProjName, std::string sampName);
	void getMinTreeDataForSample(std::string shortProjName, std::string sampName);

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

int genProjectConfig(std::map<std::string, std::string> inputCommands);
int popClusteringViewer(std::map<std::string, std::string> inputCommands);


} /* namespace bibseq */
