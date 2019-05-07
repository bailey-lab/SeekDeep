
//  ServerRunner.cpp
//
//  Created by Nick Hathaway on 2015/06/24.
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
    
#include "SeekDeepServerRunner.hpp"

#include <seqServer/apps/pars.h>
    
namespace njhseq {

SeekDeepServerRunner::SeekDeepServerRunner()
    : njh::progutils::ProgramRunner(
    		{
				addFunc("genProjectConfig", genProjectConfig, false),
				addFunc("popClusteringViewer", popClusteringViewer, false)
    		},//
                    "SeekDeepServer") {}






void errorHandler(const int statusCode, const std::exception& exception,
		const std::shared_ptr<restbed::Session>& session) {
	std::cerr << "statusCode: " << statusCode << std::endl;
	std::cerr << exception.what() << std::endl;
	if(nullptr != session && session->is_open()){
		session->close(statusCode, exception.what(), { { "Server", "Restbed" } });
	}
}


int SeekDeepServerRunner::popClusteringViewer(const njh::progutils::CmdArgs & inputCommands){
	njhseq::seqSetUp setUp(inputCommands);
	SeqAppCorePars corePars;
	corePars.port_ = 9881;
	corePars.name_ = "pcv";
	bfs::path configDir = "";
	bfs::path resourceDirName = njh::files::make_path(SeekDeep_INSTALLDIR, "etc/serverResources").string();
	setUp.setOption(resourceDirName, "--resourceDirName",
			"Name of the resource Directory where the js and html is located",
			!bfs::exists(resourceDirName));

	setUp.setOption(configDir, "--configDir", "Name of the Master Result Directory", true);

	setUp.processDebug();
	setUp.processVerbose();
	corePars.setCoreOptions(setUp);
	setUp.finishSetUp(std::cout);
	resourceDirName = njh::appendAsNeededRet(resourceDirName.string(), "/");
	configDir = njh::appendAsNeededRet(configDir.string(), "/");

	VecStr warnings;
	if(!bfs::exists(resourceDirName)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, resource directory " << resourceDirName << " doesn't exist";
		warnings.emplace_back(ss.str());
	}
	if(!bfs::exists(configDir)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, configuration directory " << configDir << " doesn't exist";
		warnings.emplace_back(ss.str());
	}
	if(!warnings.empty()){
		throw std::runtime_error{njh::conToStr(warnings, "\n")};
	}
  //
  Json::Value appConfig;
  corePars.addCoreOpts(appConfig);
  appConfig["configDir"] = njh::json::toJson(configDir);
  appConfig["resources"] = njh::json::toJson(resourceDirName);
  if(setUp.pars_.verbose_){
  	std::cout << corePars.getAddress() << std::endl;
  }
	pcv pcViewer(appConfig);
	auto resources = pcViewer.getAllResources();
	auto settings = std::make_shared<restbed::Settings>();
	settings->set_port(corePars.port_);
	settings->set_default_header("Connection", "close");
	settings->set_bind_address(corePars.bindAddress_);
	restbed::Service service;
	service.set_error_handler(errorHandler);
	for(const auto & resource : resources){
		service.publish(resource);
	}
	try {
		service.start(settings);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
	}


	return 0;
}


int SeekDeepServerRunner::genProjectConfig(
		const njh::progutils::CmdArgs & inputCommands) {
	njhseq::seqSetUp setUp(inputCommands);
	std::string mainDir = "";
	std::string projectName = "";
	std::string shortName = "";
	bool debug = false;
	setUp.setOption(projectName, "--projectName", "Name of the Project", true);
	setUp.setOption(debug, "--debug", "Run In Debug Mode");
	setUp.setOption(shortName, "--shortName", "Short Name for url for project, so only characters allowed in urls should bed used");
	if ("" == shortName) {
		shortName = projectName;
	}
	setUp.setOption(mainDir, "--mainDir", "Name of the Master Result Directory created by SeekDeep processClusters",
			true);
	njh::appendAsNeeded(mainDir, "/");
	OutOptions outOpts(shortName, ".config");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	//
	Json::Value config;
	config["mainDir"] = njh::appendAsNeededRet(
			njh::files::normalize(mainDir).string(), "/");
	config["shortName"] = shortName;
	config["projectName"] = projectName;
	config["debug"] = debug;
	OutputStream out(outOpts);
	out << config << std::endl;
	return 0;
}


} // namespace njhseq
