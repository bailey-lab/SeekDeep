#pragma once
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2026 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "SeekDeep/server/PmoProject.hpp"

namespace njhseq {

/**@brief Viewer that serves directly from a PMO file, each target in the PMO gets its own pages
 *
 */
class pmv: public njhseq::SeqApp {
public:
	/**
	 * @param config needs the SeqApp options and "pmo" and "resources"
	 */
	explicit pmv(const Json::Value & config);

	std::vector<std::shared_ptr<restbed::Resource>> getAllResources() override;
	VecStr requiredOptions() const override;

private:
	typedef njhseq::SeqApp super;
	using Session = std::shared_ptr<restbed::Session>;
	using PathSegs = std::vector<UrlPathFactory::PathSeg>;

	std::unique_ptr<PmoProject> project_;
	bfs::path resourceDir_;
	bfs::path seqsDir_;

	static const std::string pat_name_;

	static std::string pathPar(const Session & session, const std::string & name);
	static void respondJson(const Session & session, const Json::Value & ret);
	static void respondNotFound(const Session & session, const std::string & message);
	void respondPage(const Session & session, const std::string & page);

	std::shared_ptr<restbed::Resource> getRoute(const PathSegs & segs, std::function<void(const Session &)> handler);
	std::shared_ptr<restbed::Resource> postRoute(const PathSegs & segs, std::function<void(const Session &, const Json::Value &)> handler);

	void registerSeqs(const std::string & uid, const std::string & fasta);

	Json::Value sampInfoJson(const table & sampTab, const VecStr & sampNames, bool hasGroupCol) const;
	static Json::Value popInfoJson(const table & popTab, const VecStr & popUIDs, const VecStr & excludeFromNum);
	static Json::Value hapIdTabJson(const table & hapIdTab, const VecStr & popUIDs, const VecStr & samples);
	Json::Value popSeqsJson(const std::string & target, const Json::Value & postData);
};

}  // namespace njhseq
