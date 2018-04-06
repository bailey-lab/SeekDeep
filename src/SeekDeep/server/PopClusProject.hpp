#pragma once
/*
 * PopClusProject.hpp
 *
 *  Created on: Sep 19, 2016
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
#include <bibcpp.h>



namespace bibseq {

/**@brief struct to hold a pointer to table caches for sample and population information
 *
 */
struct ClusInfoTabs{
	std::unique_ptr<TableCache> sampInfo_;
	std::unique_ptr<TableCache> popInfo_;
	std::unique_ptr<TableCache> hapIdTab_;
};

/**@brief class to hold information on a clustering project to help in serving it's infomration
 *
 */
class PopClusProject {
public:

	/**@brief construct with json configuration file
	 *
	 * should have at least the following fields, "shortName", "projectName", "mainDir"
	 *
	 * @param configJson
	 */
	PopClusProject(const Json::Value & configJson);

	std::string shortName_;/**< short name for url viewer */
	std::string projectName_;/**< full length nmae of project*/

	Json::Value config_; /**< the configuration used to create project class*/
	std::unique_ptr<collapse::SampleCollapseCollection> collection_;/**< the collection class that holds info on the clustering */

	std::unique_ptr<TableCache> extractionProfileTab_;
	std::unique_ptr<TableCache> extractionStatsTab_;

	ClusInfoTabs tabs_; /**< holds pointers to caches for the pop and samp data*/

	std::unordered_map<std::string, std::unordered_map<std::string,ClusInfoTabs>> subGroupTabs_;/**< holds the sub groups population tables*/

	std::unordered_map<std::string, std::unique_ptr<TableCache>> topGroupTabs_;/**< holds the top group summarized info tabes*/

	/**@brief Register all the sequences paths to the SeqCache being used by the viewer for serving
	 *
	 * @param cache the cache to add seq files to
	 */
	void registerSeqFiles(SeqCache & cache);

};

}  // namespace bibseq



