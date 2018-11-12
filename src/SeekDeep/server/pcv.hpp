#pragma once

/*
 * pcv.hpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include <seqServer/apps/SeqApp.hpp>
#include <seqServer/utils.h>
#include <njhcpp.h>
#include "SeekDeep/server/PopClusProject.hpp"



namespace njhseq {


class pcv: public njhseq::SeqApp {
public:
	pcv(const Json::Value & config);

private:

	void projectNamesHandler(std::shared_ptr<restbed::Session> session);
	void mainPageHandler(std::shared_ptr<restbed::Session> session);


	////
	/// project
	////

	void mainProjectPageHandler(std::shared_ptr<restbed::Session> session);
	void getProjectNameHandler(std::shared_ptr<restbed::Session> session);
	void getSampleNamesHandler(std::shared_ptr<restbed::Session> session);
	void getGroupNamesHandler(std::shared_ptr<restbed::Session> session);

	void getSampleInfoTabPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getSampleInfoTabHandler(std::shared_ptr<restbed::Session> session);

	void getPopSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getPopSeqsHandler(std::shared_ptr<restbed::Session> session);

	void getPopInfoPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getPopInfoHandler(std::shared_ptr<restbed::Session> session);

	void getHapIdTablePostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getHapIdTableHandler(std::shared_ptr<restbed::Session> session);


	////
	/// sample
	////

	void samplePageHandler(std::shared_ptr<restbed::Session> session);
	void getSampSeqsHandler(std::shared_ptr<restbed::Session> session);
	void getSampSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);

	////
	/// top group
	////

	void groupInfoPageHandler(std::shared_ptr<restbed::Session> session);
	void getGroupsPopInfosHandler(std::shared_ptr<restbed::Session> session);
	void getGroupsPopInfosPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);


	////
	/// sub group
	////
	void groupMainProjectPageHanlder(std::shared_ptr<restbed::Session> session);
	void groupGetSampleNamesHanlder(std::shared_ptr<restbed::Session> session);

	void groupGetSampleInfoTabPostHanlder(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void groupGetSampleInfoTabHanlder(std::shared_ptr<restbed::Session> session);
	void groupGetPopSeqsPostHanlder(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void groupGetPopSeqsHanlder(std::shared_ptr<restbed::Session> session);
	void groupGetPopInfoPostHanlder(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void groupGetPopInfoHanlder(std::shared_ptr<restbed::Session> session);

	void groupGetHapIdTablePostHanlder(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void groupGetHapIdTableHanlder(std::shared_ptr<restbed::Session> session);




	///
	/// extraction
	///

	void extractionPageHanlder(std::shared_ptr<restbed::Session> session);
	void getExtractionProfileDataHanlder(std::shared_ptr<restbed::Session> session);
	void getExtractionStatsDataHanlder(std::shared_ptr<restbed::Session> session);

	typedef njhseq::SeqApp super;
	bfs::path configDir_;
	bfs::path resourceDir_;

	std::map<std::string, std::unique_ptr<PopClusProject>> collections_;


	void loadInCollections();

	void redirect(std::shared_ptr<restbed::Session> session, std::string errorMessage);


public:
	virtual std::vector<std::shared_ptr<restbed::Resource>> getAllResources();
	virtual VecStr requiredOptions() const;

	///
	/// main page
	///
	std::shared_ptr<restbed::Resource> projectNames();
	std::shared_ptr<restbed::Resource> mainPage();

	///
	/// project
	///
	std::shared_ptr<restbed::Resource> mainProjectPage();
	std::shared_ptr<restbed::Resource> getProjectName();
	std::shared_ptr<restbed::Resource> getSampleNames();
	std::shared_ptr<restbed::Resource> getGroupNames();
	std::shared_ptr<restbed::Resource> getSampleInfoTab();
	std::shared_ptr<restbed::Resource> getPopSeqs();
	std::shared_ptr<restbed::Resource> getPopInfo();

	std::shared_ptr<restbed::Resource> getHapIdTable();

	///
	/// sample
	///
	std::shared_ptr<restbed::Resource> samplePage();
	std::shared_ptr<restbed::Resource> getSampSeqs();

	///
	/// top group
	///

	std::shared_ptr<restbed::Resource> groupInfoPage();
	std::shared_ptr<restbed::Resource> getGroupsPopInfos();

	///
	/// sub group
	///
	std::shared_ptr<restbed::Resource> groupMainProjectPage();
	std::shared_ptr<restbed::Resource> groupGetSampleNames();
	std::shared_ptr<restbed::Resource> groupGetSampleInfoTab();
	std::shared_ptr<restbed::Resource> groupGetPopSeqs();
	std::shared_ptr<restbed::Resource> groupGetPopInfo();

	std::shared_ptr<restbed::Resource> groupGetHapIdTable();

	///
	/// extraction
	///
	std::shared_ptr<restbed::Resource> extractionPage();
	std::shared_ptr<restbed::Resource> getExtractionProfileData();
	std::shared_ptr<restbed::Resource> getExtractionStatsData();


};



}  // namespace njhseq



