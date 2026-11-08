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
#include "pmv.hpp"
#include "pcv.hpp"
#include <njhseq/IO/OutputStream.hpp>
#include <njhseq/objects/collapseObjects/SampleCollapseCollection.hpp>
#include <unistd.h>

namespace njhseq {

//names in a PMO can have spaces and other characters in them, the requests have them url encoded
const std::string pmv::pat_name_{".+"};

pmv::pmv(const Json::Value & config) :
		njhseq::SeqApp(config) {
	resourceDir_ = config["resources"].asString();
	jsFiles_->addFiles(njh::files::gatherFiles(njh::files::make_path(resourceDir_, "pmv/js"), ".js"));
	cssFiles_->addFiles(njh::files::gatherFiles(njh::files::make_path(resourceDir_, "pmv/css"), ".css"));
	addScripts(njh::files::make_path(resourceDir_, "pmv"));
	//per process since the working directory defaults to the shared /tmp/
	seqsDir_ = njh::files::make_path(seqs_->workingDir_, njh::pasteAsStr("pmoSeqs_", getpid()));
	project_ = std::make_unique<PmoProject>(config["pmo"].asString());
}

VecStr pmv::requiredOptions() const {
	return concatVecs(super::requiredOptions(), VecStr{"pmo", "resources"});
}

std::string pmv::pathPar(const Session & session, const std::string & name) {
	return restbed::Uri::decode(session->get_request()->get_path_parameter(name));
}

void pmv::respondJson(const Session & session, const Json::Value & ret) {
	auto body = njh::json::writeAsOneLine(ret);
	auto headers = HeaderFactory::initiateAppJsonHeader(body);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, body, headers);
}

void pmv::respondNotFound(const Session & session, const std::string & message) {
	std::cerr << message << std::endl;
	session->close(restbed::NOT_FOUND, message, {{"Content-Length", estd::to_string(message.size())}, {"Connection", "close"}});
}

void pmv::respondPage(const Session & session, const std::string & page) {
	auto body = genHtmlDoc(rootName_, pages_.at(page));
	session->close(restbed::OK, body, HeaderFactory::initiateTxtHtmlHeader(body));
}

std::shared_ptr<restbed::Resource> pmv::getRoute(const PathSegs & segs, std::function<void(const Session &)> handler) {
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl(segs));
	resource->set_method_handler("GET", std::function<void(Session)>([this, handler](Session session) {
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		handler(session);
	}));
	return resource;
}

std::shared_ptr<restbed::Resource> pmv::postRoute(const PathSegs & segs, std::function<void(const Session &, const Json::Value &)> handler) {
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl(segs));
	resource->set_method_handler("POST", std::function<void(Session)>([this, handler](Session session) {
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		size_t contentLength = session->get_request()->get_header("Content-Length", 0);
		session->fetch(contentLength, std::function<void(Session, const restbed::Bytes &)>([handler](Session ses, const restbed::Bytes & body) {
			handler(ses, njh::json::parse(std::string(body.begin(), body.end())));
		}));
	}));
	return resource;
}

void pmv::registerSeqs(const std::string & uid, const std::string & fasta) {
	if (seqs_->containsRecord(uid)) {
		return;
	}
	njh::files::makeDirP(njh::files::MkdirPar(seqsDir_.string()));
	//targets and samples are numbered since the names might not be file safe
	auto seqFnp = njh::files::make_path(seqsDir_, njh::pasteAsStr("seqs_", seqs_->cache_.size(), ".fasta"));
	{
		OutOptions seqOutOpts(seqFnp);
		seqOutOpts.overWriteFile_ = true;
		OutputStream seqOut(seqOutOpts);
		seqOut << fasta;
	}
	seqs_->updateAddCache(uid, SeqIOOptions::genFastaIn(seqFnp, true));
}

Json::Value pmv::sampInfoJson(const table & sampTab, const VecStr & sampNames, bool hasGroupCol) const {
	const std::set<std::string> sampNamesSet(sampNames.begin(), sampNames.end());
	auto trimmedTab = sampTab.extractByComp("s_Name", [&sampNamesSet](const std::string & str) {
		return njh::in(str, sampNamesSet);
	});
	VecStr visibleColumns{"s_Sample"};
	if (hasGroupCol) {
		visibleColumns.emplace_back("g_GroupName");
	}
	addOtherVec(visibleColumns, VecStr{"h_popUID", "h_SampCnt", "h_SampFrac", "s_ReadCntTotUsed", "s_COI", "c_clusterID", "c_AveragedFrac", "c_ReadCnt"});
	auto ret = tableToJsonByRow(trimmedTab, "s_Name", visibleColumns);
	auto popNames = trimmedTab.getColumn("h_popUID");
	removeDuplicates(popNames);
	njh::sort(popNames);
	ret["popColors"] = njh::json::toJson(genColorsStrsForNamesMultiples(popNames));
	ret["popUIDs"] = njh::json::toJson(popNames);
	return ret;
}

Json::Value pmv::popInfoJson(const table & popTab, const VecStr & popUIDs, const VecStr & excludeFromNum) {
	const std::set<std::string> popUIDsSet(popUIDs.begin(), popUIDs.end());
	auto trimmedTab = popTab.extractByComp("h_popUID", [&popUIDsSet](const std::string & str) {
		return njh::in(str, popUIDsSet);
	});
	return tableToJsonByRow(trimmedTab, "h_popUID", VecStr{}, excludeFromNum);
}

Json::Value pmv::hapIdTabJson(const table & hapIdTab, const VecStr & popUIDs, const VecStr & samples) {
	const std::set<std::string> popUIDsSet(popUIDs.begin(), popUIDs.end());
	auto trimmedTab = hapIdTab.getColumns(concatVecs(VecStr{"#PopUID"}, samples)).extractByComp("#PopUID", [&popUIDsSet](const std::string & str) {
		return njh::in(str, popUIDsSet);
	});
	return tableToJsonByRow(trimmedTab, "#PopUID");
}

Json::Value pmv::popSeqsJson(const std::string & target, const Json::Value & postData) {
	registerSeqs(target, project_->popSeqsFasta(target));
	const uint32_t sesUid = postData.isMember("sessionUID") ? postData["sessionUID"].asUInt() : startSeqCacheSession();
	auto popUIDs = njh::json::jsonArrayToSet<std::string>(postData["popUIDs"], [](const Json::Value & val) {return val.asString();});
	auto & cache = seqsBySession_[sesUid]->cache_;
	//sessions started before this target was registered won't have it yet
	if (!njh::in(target, cache)) {
		seqsBySession_[sesUid]->updateAddCache(target, seqs_->cache_.at(target).ioOpts_.opts_);
	}
	auto & record = cache.at(target);
	record.reload();
	record.toggleSeqs([&popUIDs](const readObject & seq) {
		return njh::in(seq.seqBase_.getStubName(false), popUIDs);
	});
	//viewer doesn't know how to handle empty seqs
	record.ensureNonEmptyReads();
	auto ret = seqsBySession_[sesUid]->getJson(target);
	ret["sessionUID"] = njh::json::toJson(sesUid);
	return ret;
}

std::vector<std::shared_ptr<restbed::Resource>> pmv::getAllResources() {
	auto ret = super::getAllResources();
	const auto toStrVec = [](const Json::Value & vals) {
		return njh::json::jsonArrayToVec<std::string>(vals, [](const Json::Value & val) {return val.asString();});
	};
	auto seg = [this](const std::string & routeName, const VecStr & parNames) {
		PathSegs built{{rootName_}};
		if (!routeName.empty()) {
			built.emplace_back(UrlPathFactory::PathSeg{routeName});
		}
		for (const auto & par : parNames) {
			built.emplace_back(UrlPathFactory::PathSeg{par, pat_name_});
		}
		return built;
	};
	auto targetOr404 = [this](const Session & session, const std::function<void(const std::string &)> & func) {
		auto target = pathPar(session, "target");
		if (!project_->hasTarget(target)) {
			respondNotFound(session, "no such target as " + target);
			return;
		}
		func(target);
	};
	auto subGroupOr404 = [this](const Session & session, const std::function<void(const std::string &, const std::string &, const std::string &)> & func) {
		auto target = pathPar(session, "target");
		auto group = pathPar(session, "groupName");
		auto subGroup = pathPar(session, "subGroupName");
		if (!project_->hasSubGroup(target, group, subGroup)) {
			respondNotFound(session, "no such sub group as " + subGroup + " of " + group + " for " + target);
			return;
		}
		func(target, group, subGroup);
	};
	const VecStr groupPopInfoExclude{"p_TotalInputReadCnt", "g_GroupName", "g_hapsFoundOnlyInThisGroup", "p_TotalUniqueHaplotypes", "p_TotalInputClusterCnt",
		"p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanCoi", "p_medianCoi", "p_minCoi", "p_maxCoi"};
	const VecStr popInfoExclude{"p_TotalInputReadCnt", "p_TotalInputClusterCnt", "p_TotalPopulationSampCnt", "p_TotalHaplotypes", "p_meanCoi", "p_medianCoi", "p_minCoi", "p_maxCoi"};

	//top level
	ret.emplace_back(getRoute(seg("", {}), [this](const Session & session) {
		respondPage(session, "mainPage.js");
	}));
	ret.emplace_back(getRoute(seg("projectsNames", {}), [this](const Session & session) {
		Json::Value out;
		out["projects"] = njh::json::toJson(project_->targetNames());
		respondJson(session, out);
	}));
	ret.emplace_back(getRoute(seg("pmoHeader", {}), [this](const Session & session) {
		respondJson(session, project_->header_);
	}));
	ret.emplace_back(getRoute(seg("targetsInfo", {}), [this](const Session & session) {
		respondJson(session, tableToJsonByRow(project_->targetsInfo(), "target"));
	}));
	ret.emplace_back(getRoute(seg("librarySamplesInfo", {}), [this](const Session & session) {
		respondJson(session, tableToJsonByRow(project_->librarySamplesInfo(), "library_sample", VecStr{}, VecStr{"specimen", "library_sample"}));
	}));

	//target
	ret.emplace_back(getRoute(seg("mainProjectPage", {"target"}), [this, targetOr404](const Session & session) {
		targetOr404(session, [this, &session](const std::string &) {
			respondPage(session, "mainProjectPage.js");
		});
	}));
	ret.emplace_back(getRoute(seg("projectName", {"target"}), [targetOr404](const Session & session) {
		targetOr404(session, [&session](const std::string & target) {
			Json::Value out;
			out["projectName"] = target;
			respondJson(session, out);
		});
	}));
	ret.emplace_back(getRoute(seg("sampleNames", {"target"}), [this, targetOr404](const Session & session) {
		targetOr404(session, [this, &session](const std::string & target) {
			Json::Value out;
			out["samples"] = njh::json::toJson(project_->sampleNames(target));
			respondJson(session, out);
		});
	}));
	ret.emplace_back(getRoute(seg("getGroupNames", {"target"}), [this, targetOr404](const Session & session) {
		targetOr404(session, [this, &session](const std::string & target) {
			Json::Value out;
			out["groups"] = njh::json::toJson(project_->groupNames(target));
			respondJson(session, out);
		});
	}));
	ret.emplace_back(postRoute(seg("sampInfo", {"target"}), [this, targetOr404, toStrVec](const Session & session, const Json::Value & post) {
		targetOr404(session, [this, &session, &post, &toStrVec](const std::string & target) {
			respondJson(session, sampInfoJson(project_->targetTables(target).sampInfo_, toStrVec(post["sampNames"]), false));
		});
	}));
	ret.emplace_back(postRoute(seg("popInfo", {"target"}), [this, targetOr404, toStrVec, popInfoExclude](const Session & session, const Json::Value & post) {
		targetOr404(session, [this, &session, &post, &toStrVec, &popInfoExclude](const std::string & target) {
			respondJson(session, popInfoJson(project_->targetTables(target).popInfo_, toStrVec(post["popUIDs"]), popInfoExclude));
		});
	}));
	ret.emplace_back(postRoute(seg("hapIdTable", {"target"}), [this, targetOr404, toStrVec](const Session & session, const Json::Value & post) {
		targetOr404(session, [this, &session, &post, &toStrVec](const std::string & target) {
			respondJson(session, hapIdTabJson(project_->targetTables(target).hapIdTab_, toStrVec(post["popUIDs"]), toStrVec(post["samples"])));
		});
	}));
	ret.emplace_back(postRoute(seg("popSeqData", {"target"}), [this, targetOr404](const Session & session, const Json::Value & post) {
		targetOr404(session, [this, &session, &post](const std::string & target) {
			respondJson(session, popSeqsJson(target, post));
		});
	}));

	//one sample across all targets
	ret.emplace_back(getRoute(seg("librarySampleNames", {}), [this](const Session & session) {
		Json::Value out;
		out["samples"] = njh::json::toJson(project_->librarySampleNames());
		respondJson(session, out);
	}));
	ret.emplace_back(getRoute(seg("sampleAcrossTargetsPage", {"sampleName"}), [this](const Session & session) {
		auto sample = pathPar(session, "sampleName");
		if (!project_->hasLibrarySample(sample)) {
			respondNotFound(session, "no such library sample as " + sample);
			return;
		}
		respondPage(session, "sampleAcrossTargetsPage.js");
	}));
	ret.emplace_back(getRoute(seg("sampleAcrossTargets", {"sampleName"}), [this](const Session & session) {
		auto sample = pathPar(session, "sampleName");
		if (!project_->hasLibrarySample(sample)) {
			respondNotFound(session, "no such library sample as " + sample);
			return;
		}
		auto calls = project_->sampleAcrossTargets(sample);
		auto out = tableToJsonByRow(calls, "p_name", VecStr{"p_name", "h_popUID", "c_clusterID", "c_AveragedFrac", "c_ReadCnt", "s_COI", "h_SampCnt", "h_SampFrac"});
		//colored by rank within a target so the major haplotype is the same color in every target
		auto ranks = calls.getColumn("c_clusterID");
		removeDuplicates(ranks);
		njh::sort(ranks);
		out["popColors"] = njh::json::toJson(genColorsStrsForNamesMultiples(ranks));
		out["popUIDs"] = njh::json::toJson(ranks);
		out["targetsNotDetected"] = njh::json::toJson(project_->targetsNotDetected(sample));
		out["targetCount"] = njh::json::toJson(static_cast<uint32_t>(project_->targetNames().size()));
		for (const auto & field : project_->sampleMeta(sample)) {
			out["meta"][field.first] = field.second;
		}
		respondJson(session, out);
	}));
	ret.emplace_back(getRoute(seg("sampleAcrossTargetsSeqs", {"sampleName"}), [this](const Session & session) {
		auto sample = pathPar(session, "sampleName");
		if (!project_->hasLibrarySample(sample) || project_->targetsNotDetected(sample).size() == project_->targetNames().size()) {
			respondNotFound(session, "no haplotypes for " + sample);
			return;
		}
		auto uid = "pmoAllTargets_" + sample;
		registerSeqs(uid, project_->sampleAcrossTargetsFasta(sample));
		auto sesUid = startSeqCacheSession();
		auto out = seqsBySession_[sesUid]->getJson(uid);
		out["sessionUID"] = njh::json::toJson(sesUid);
		respondJson(session, out);
	}));

	//haplotype rainbow of all samples or of a sub group across all targets
	ret.emplace_back(getRoute(seg("metaGroups", {}), [this](const Session & session) {
		Json::Value out;
		out["groups"] = Json::objectValue;
		for (const auto & group : project_->metaGroupNames()) {
			out["groups"][group] = njh::json::toJson(project_->metaSubGroupNames(group));
		}
		respondJson(session, out);
	}));
	ret.emplace_back(getRoute(seg("haplotypeRainbowPage", {}), [this](const Session & session) {
		respondPage(session, "haplotypeRainbowPage.js");
	}));
	ret.emplace_back(getRoute(seg("haplotypeRainbow", {}), [this](const Session & session) {
		respondJson(session, project_->haplotypeRainbow());
	}));
	auto metaSubGroupOr404 = [this](const Session & session, const std::function<void(const std::string &, const std::string &)> & func) {
		auto group = pathPar(session, "groupName");
		auto subGroup = pathPar(session, "subGroupName");
		if (!project_->hasMetaSubGroup(group, subGroup)) {
			respondNotFound(session, "no such sub group as " + subGroup + " of " + group);
			return;
		}
		func(group, subGroup);
	};
	ret.emplace_back(getRoute(seg("haplotypeRainbowPage", {"groupName", "subGroupName"}), [this, metaSubGroupOr404](const Session & session) {
		metaSubGroupOr404(session, [this, &session](const std::string &, const std::string &) {
			respondPage(session, "haplotypeRainbowPage.js");
		});
	}));
	ret.emplace_back(getRoute(seg("haplotypeRainbow", {"groupName", "subGroupName"}), [this, metaSubGroupOr404](const Session & session) {
		metaSubGroupOr404(session, [this, &session](const std::string & group, const std::string & subGroup) {
			respondJson(session, project_->haplotypeRainbow(group, subGroup));
		});
	}));

	//sample
	ret.emplace_back(getRoute(seg("sampleMainPage", {"target", "sampleName"}), [this](const Session & session) {
		auto target = pathPar(session, "target");
		auto sample = pathPar(session, "sampleName");
		if (!project_->hasSample(target, sample)) {
			respondNotFound(session, "no such sample as " + sample + " for " + target);
			return;
		}
		respondPage(session, "sampleMainPage.js");
	}));
	ret.emplace_back(postRoute(seg("sampSeqData", {"target"}), [this, targetOr404](const Session & session, const Json::Value & post) {
		targetOr404(session, [this, &session, &post](const std::string & target) {
			auto sample = post["sampleName"].asString();
			if (!project_->hasSample(target, sample)) {
				respondNotFound(session, "no such sample as " + sample + " for " + target);
				return;
			}
			auto uid = target + "_" + sample;
			registerSeqs(uid, project_->sampleSeqsFasta(target, sample));
			auto sesUid = startSeqCacheSession();
			auto out = seqsBySession_[sesUid]->getJson(uid);
			out["sessionUID"] = njh::json::toJson(sesUid);
			respondJson(session, out);
		});
	}));

	//group
	ret.emplace_back(getRoute(seg("groupInfoPage", {"target", "groupName"}), [this](const Session & session) {
		auto target = pathPar(session, "target");
		auto group = pathPar(session, "groupName");
		if (!project_->hasGroup(target, group)) {
			respondNotFound(session, "no such group as " + group + " for " + target);
			return;
		}
		respondPage(session, "groupInfoPage.js");
	}));
	ret.emplace_back(getRoute(seg("groupsPopInfos", {"target", "groupName"}), [this](const Session & session) {
		auto target = pathPar(session, "target");
		auto group = pathPar(session, "groupName");
		if (!project_->hasGroup(target, group)) {
			respondNotFound(session, "no such group as " + group + " for " + target);
			return;
		}
		Json::Value out;
		out["groupNames"] = njh::json::toJson(project_->subGroupNames(target, group));
		out["popInfo"] = tableToJsonByRow(project_->groupInfo(target, group), "g_GroupName", VecStr{});
		respondJson(session, out);
	}));

	//sub group
	const VecStr subGroupPars{"target", "groupName", "subGroupName"};
	ret.emplace_back(getRoute(seg("groupMainPage", subGroupPars), [this, subGroupOr404](const Session & session) {
		subGroupOr404(session, [this, &session](const std::string &, const std::string &, const std::string &) {
			respondPage(session, "groupMainPage.js");
		});
	}));
	ret.emplace_back(getRoute(seg("groupSampleNames", subGroupPars), [this, subGroupOr404](const Session & session) {
		subGroupOr404(session, [this, &session](const std::string & target, const std::string & group, const std::string & subGroup) {
			Json::Value out;
			out["groupSamples"] = njh::json::toJson(project_->subGroupSampleNames(target, group, subGroup));
			respondJson(session, out);
		});
	}));
	ret.emplace_back(postRoute(seg("groupSampInfo", subGroupPars), [this, subGroupOr404, toStrVec](const Session & session, const Json::Value & post) {
		subGroupOr404(session, [this, &session, &post, &toStrVec](const std::string & target, const std::string & group, const std::string & subGroup) {
			respondJson(session, sampInfoJson(project_->subGroupTables(target, group, subGroup).sampInfo_, toStrVec(post["sampNames"]), true));
		});
	}));
	ret.emplace_back(postRoute(seg("groupPopInfo", subGroupPars), [this, subGroupOr404, toStrVec, groupPopInfoExclude](const Session & session, const Json::Value & post) {
		subGroupOr404(session, [this, &session, &post, &toStrVec, &groupPopInfoExclude](const std::string & target, const std::string & group, const std::string & subGroup) {
			respondJson(session, popInfoJson(project_->subGroupTables(target, group, subGroup).popInfo_, toStrVec(post["popUIDs"]), groupPopInfoExclude));
		});
	}));
	ret.emplace_back(postRoute(seg("groupHapIdTable", subGroupPars), [this, subGroupOr404, toStrVec](const Session & session, const Json::Value & post) {
		subGroupOr404(session, [this, &session, &post, &toStrVec](const std::string & target, const std::string & group, const std::string & subGroup) {
			respondJson(session, hapIdTabJson(project_->subGroupTables(target, group, subGroup).hapIdTab_, toStrVec(post["popUIDs"]), toStrVec(post["samples"])));
		});
	}));
	ret.emplace_back(postRoute(seg("groupPopSeqData", subGroupPars), [this, subGroupOr404](const Session & session, const Json::Value & post) {
		subGroupOr404(session, [this, &session, &post](const std::string & target, const std::string &, const std::string &) {
			respondJson(session, popSeqsJson(target, post));
		});
	}));
	return ret;
}

}  // namespace njhseq
