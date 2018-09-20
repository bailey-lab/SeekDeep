/*
 * PrimersAndMids.cpp
 *
 *  Created on: Nov 27, 2016
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
#include "PrimersAndMids.hpp"

namespace bibseq {

PrimersAndMids::Target::lenCutOffs::lenCutOffs(uint32_t minLen, uint32_t maxLen,
		bool mark) :
		minLenChecker_(ReadCheckerLenAbove(minLen, mark)), maxLenChecker_(
				ReadCheckerLenBelow(maxLen, mark)) {
}


PrimersAndMids::Target::Target(const std::string & name,
		const std::string & forPrimer, const std::string & revPrimer) :
		info_(name, forPrimer, revPrimer),
		overlapStatus_(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP) {
}

void PrimersAndMids::Target::addLenCutOff(uint32_t minLen, uint32_t maxLen,
		bool mark) {
	lenCuts_ = std::make_shared<lenCutOffs>(minLen, maxLen, mark);
}

void PrimersAndMids::Target::setRefKInfos(uint32_t klen, bool setRevComp) {
	refKInfos_.clear();
	for (const auto & ref : refs_) {
		refKInfos_.emplace_back(ref.seq_, klen, setRevComp);
	}
}

void PrimersAndMids::Target::addSingleRef(const seqInfo & ref) {
	if (info_.primerPairName_ == ref.name_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": error, adding ref that doesn't the target name" << "\n";
		ss << "TargetName: " << info_.primerPairName_ << "\n";
		ss << "AddingRefName: " << ref.name_ << "\n";
	}
	refs_.emplace_back(ref);
}
void PrimersAndMids::Target::addMultileRef(const std::vector<seqInfo> & refs) {
	//mutliple refs so no name matchig
	addOtherVec(refs_, refs);
}


PrimersAndMids::PrimersAndMids(
		const std::unordered_map<std::string, Target> & targets) :
		targets_(targets) {

}





PrimersAndMids::PrimersAndMids(const bfs::path & idFileFnp) : idFile_(idFileFnp) {
	if (!bfs::exists(idFile_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, "
				<< bib::bashCT::boldRed(idFile_.string()) << " doesn't exist\n";
		throw std::runtime_error { ss.str() };
	}
	auto firstLine = bib::files::getFirstLine(idFileFnp);
	bib::strToLower(firstLine);
	if (!bib::beginsWith(firstLine, "gene")
			&& !bib::beginsWith(firstLine, "target")
				&& !bib::beginsWith(firstLine, "id")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error the id file "
				<< bib::bashCT::boldRed(idFile_.string())
				<< " should start with either target, gene, or id (case insensitive)\n";
		ss << "line: " << firstLine << "\n";
		throw std::runtime_error { ss.str() };
	}
	bool readingPrimers = false;
	bool readingMids = false;
	InputStream idFile{InOptions{idFileFnp}};
	std::string line = "";
	while (bib::files::crossPlatGetline(idFile, line)) {
		auto lowerLine = bib::strToLowerRet(line);
		auto lowerLineToks = bib::tokenizeString(lowerLine, "whitespace");

		if (lowerLineToks.size() > 1 &&
				("gene" == lowerLineToks[0]
				|| "target" == lowerLineToks[0])) {
			readingPrimers = true;
			readingMids = false;
		} else if (lowerLineToks.size() > 1 && "id" == lowerLineToks[0]) {
			readingPrimers = false;
			readingMids = true;
		} else if (readingPrimers) {
			auto toks = bib::tokenizeString(line, "whitespace");
			if (3 != toks.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error in id file "
						<< bib::bashCT::boldRed(idFile_.string())
						<< "primer line should contain 3 items not " << toks.size() << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error { ss.str() };
			}
			addTarget(toks[0], toks[1], toks[2]);
		} else if (readingMids) {
			auto toks = bib::tokenizeString(line, "whitespace");
			if (2 != toks.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error in id file "
						<< bib::bashCT::boldRed(idFile_.string())
						<< "barcode line should contain 2 items not " << toks.size()
						<< "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error { ss.str() };
			}
			addMid(toks[0], toks[1]);
		}
	}
}

bool PrimersAndMids::hasTarget(const std::string & target) const {
	return targets_.end() != targets_.find(target);
}

VecStr PrimersAndMids::getTargets() const {
	auto ret = bib::getVecOfMapKeys(targets_);
	bib::sort(ret);
	return ret;
}

VecStr PrimersAndMids::getMids() const {
	auto ret = bib::getVecOfMapKeys(mids_);
	bib::sort(ret);
	return ret;
}

bool PrimersAndMids::hasMid(const std::string & mid) const {
	return mids_.end() != mids_.find(mid);
}

void PrimersAndMids::addTarget(const std::string & primerName,
		const std::string & forPrimer, const std::string & revPrimer) {
	if (hasTarget(primerName)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, already have "
				<< bib::bashCT::boldRed(primerName) << "\n";
		throw std::runtime_error { ss.str() };
	}
	targets_.emplace(primerName, Target(primerName, forPrimer, revPrimer));
}

void PrimersAndMids::addMid(const std::string & midNmae,
		const std::string & barcode) {
	if (hasMid(midNmae)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, already have "
				<< bib::bashCT::boldRed(midNmae) << "\n";
		throw std::runtime_error { ss.str() };
	}
	mids_.emplace(midNmae, MidDeterminator::MidInfo(midNmae, barcode));
}

bool PrimersAndMids::hasMultipleTargets() const {
	return targets_.size() > 1;
}

bool PrimersAndMids::containsMids() const {
	return !mids_.empty();
}

bool PrimersAndMids::containsTargets() const {
	return !targets_.empty();
}

bool PrimersAndMids::screeningForPossibleContamination() const {
	bool hasRefs = false;
	for(const auto & tar : targets_){
		if(!tar.second.refs_.empty()){
			hasRefs = true;
			break;
		}
	}
	return hasRefs;
}





void PrimersAndMids::writeIdFile(const OutOptions & outOpts) const {
	std::ofstream idFileOut;
	openTextFile(idFileOut, outOpts);
	idFileOut << "target\tforward\treverse\n";
	auto pKeys = getTargets();
	bib::sort(pKeys);
	for (const auto & pKey : pKeys) {
		idFileOut << pKey << "\t" << targets_.at(pKey).info_.forwardPrimer_ << "\t"
				<< targets_.at(pKey).info_.reversePrimer_ << "\n";
	}
	if (containsMids()) {
		idFileOut << "id\tbarcode\n";
		auto mKeys = getMids();
		bib::sort(mKeys);
		for (const auto & mKey : mKeys) {
			idFileOut << mKey << "\t" << mids_.at(mKey).bar_->motifOriginal_ << "\n";
		}
	}
}

void PrimersAndMids::writeIdFile(const OutOptions & outOpts,
		const VecStr & targets) const {
	std::ofstream idFileOut;
	openTextFile(idFileOut, outOpts);
	idFileOut << "target\tforward\treverse\n";
	auto pKeys = getTargets();
	bib::sort(pKeys);
	for (const auto & pKey : pKeys) {
		if (bib::in(pKey, targets)) {
			idFileOut << pKey << "\t" << targets_.at(pKey).info_.forwardPrimer_
					<< "\t" << targets_.at(pKey).info_.reversePrimer_ << "\n";
		}
	}
	if (containsMids()) {
		idFileOut << "id\tbarcode\n";
		auto mKeys = getMids();
		bib::sort(mKeys);
		for (const auto & mKey : mKeys) {
			idFileOut << mKey << "\t" << mids_.at(mKey).bar_->motifOriginal_ << "\n";
		}
	}
}

table PrimersAndMids::genLenCutOffs(const VecStr & targets) const {
	table ret(VecStr { "target", "minlen", "maxlen" });
	for (const auto & tar : targets) {
		if (!hasTarget(tar)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, doesn't contain " << tar
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		if (nullptr != targets_.at(tar).lenCuts_) {
			ret.addRow(tar, targets_.at(tar).lenCuts_->minLenChecker_.minLen_,
					targets_.at(tar).lenCuts_->maxLenChecker_.maxLen_);
		}
	}
	return ret;
}

table PrimersAndMids::genOverlapStatuses(const VecStr & targets) const {
	table ret(VecStr { "target", "status"});
	for (const auto & tar : targets) {
		if (!hasTarget(tar)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, doesn't contain " << tar
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		if(PairedReadProcessor::ReadPairOverLapStatus::NONE != targets_.at(tar).overlapStatus_){
			ret.addRow(tar, PairedReadProcessor::getOverlapStatusStr(targets_.at(tar).overlapStatus_));
		}
	}
	return ret;
}


std::vector<seqInfo> PrimersAndMids::getRefSeqs(const VecStr & targets) const {
	std::vector<seqInfo> ret;
	for (const auto & tar : targets) {
		if (!hasTarget(tar)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error, doesn't contain " << tar
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		if (!targets_.at(tar).refs_.empty()) {
			addOtherVec(ret, targets_.at(tar).refs_);
		}
	}
	return ret;
}

void PrimersAndMids::checkMidNamesThrow() const {
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ": error\n";
	bool failed = false;
	std::regex pat{"(^MID)([0-9]+$)"};
	std::smatch match;
	for (const auto & mid : mids_) {
		if(!std::regex_match(mid.second.midName_, match, pat)){
			ss << "MID names need to begin with MID and end with a number: failed " << bib::bashCT::red
				<< mid.second.midName_ << bib::bashCT::reset << "\n";
			failed = true;
		}
	}
	if (failed) {
		throw std::runtime_error { ss.str() };
	}
}


std::map<std::string, PrimersAndMids::Target::lenCutOffs> PrimersAndMids::readInLenCutOffs(const bfs::path & lenCutOffsFnp){
	std::map<std::string, PrimersAndMids::Target::lenCutOffs> ret;
	table lenCutTab = table(lenCutOffsFnp.string(), "whitespace", true);
	bib::for_each(lenCutTab.columnNames_,
			[](std::string & str) {stringToLower(str);});
	lenCutTab.setColNamePositions();
	lenCutTab.checkForColumnsThrow(VecStr{"target", "minlen", "maxlen"}, __PRETTY_FUNCTION__);
	for (const auto & row : lenCutTab.content_) {
		ret.emplace(row[lenCutTab.getColPos("target")],
				PrimersAndMids::Target::lenCutOffs { estd::stou(
						row[lenCutTab.getColPos("minlen")]), estd::stou(
						row[lenCutTab.getColPos("maxlen")]) });
	}
	return ret;
}

void PrimersAndMids::checkIfMIdsOrPrimersReadInThrow(
		const std::string & funcName) const {
	if (0 == getMids().size() && 0 == getTargets().size()) {
		std::stringstream ss;
		ss << funcName << ", failed to read either targets or primers from "
				<< idFile_ << "\n";
		throw std::runtime_error { ss.str() };
	}
}

void PrimersAndMids::initAllAddLenCutsRefs(const InitPars & pars){
	//init mids
	if(containsMids()){
		initMidDeterminator();
		mDeterminator_->setAllowableMismatches(pars.barcodeErrors_);
		mDeterminator_->setMidEndsRevComp(pars.midEndsRevComp_);
	}

	//init primers
	if(containsTargets()){
		initPrimerDeterminator();
	}

	//add in any length cuts if any
	if("" != pars.lenCutOffFilename_){
		addLenCutOffs(pars.lenCutOffFilename_);
	}

	//add in ref sequences if any
	if("" != pars.comparisonSeqFnp_){
		addRefSeqs(pars.comparisonSeqFnp_);
		setRefSeqsKInfos(pars.compKmerLen_, true);
	}
}

void PrimersAndMids::initMidDeterminator(){
	if(mids_.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error mids_ are empty, can't init mid determinator" << "\n";
		throw std::runtime_error{ss.str()};
	}
	mDeterminator_ = std::make_unique<MidDeterminator>(mids_);
}

void PrimersAndMids::initPrimerDeterminator(){
	if(targets_.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error targets_ are empty, can't init primer determinator" << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, PrimerDeterminator::primerInfo> pInfos;
	for(const auto & tar : targets_){
		pInfos.emplace(tar.second.info_.primerPairName_, tar.second.info_);
	}
	pDeterminator_ = std::make_unique<PrimerDeterminator>(pInfos);
}

void PrimersAndMids::addLenCutOffs(const bfs::path & lenCutOffsFnp){
	bool failedOverlapStatusProcessing = false;
	std::stringstream errorStream;

	std::map<std::string, PrimersAndMids::Target::lenCutOffs> ret;
	table lenCutTab = table(lenCutOffsFnp, "whitespace", true);
	bib::for_each(lenCutTab.columnNames_,
			[](std::string & str) {stringToLower(str);});
	lenCutTab.setColNamePositions();
	lenCutTab.checkForColumnsThrow(VecStr{"target", "minlen", "maxlen"}, __PRETTY_FUNCTION__);
	for (const auto & row : lenCutTab.content_) {
		if (!bib::in(row[lenCutTab.getColPos("target")], targets_)) {
			failedOverlapStatusProcessing = true;
			errorStream << __PRETTY_FUNCTION__ << ", error found "
					<< row[lenCutTab.getColPos("target")]
					<< " but it isn't a listed target in " << idFile_ << "\n";
			errorStream << "options are: "
					<< bib::conToStr(getVectorOfMapKeys(targets_), ", ") << "\n";
		}
		if (bib::in(row[lenCutTab.getColPos("target")], ret)) {
			failedOverlapStatusProcessing = true;
			errorStream << __PRETTY_FUNCTION__ << ", error already have "
					<< row[lenCutTab.getColPos("target")] << " in table" << "\n";
		}
		ret.emplace(row[lenCutTab.getColPos("target")],
				PrimersAndMids::Target::lenCutOffs { estd::stou(
						row[lenCutTab.getColPos("minlen")]), estd::stou(
						row[lenCutTab.getColPos("maxlen")]) });
	}

	for (const auto & tar : getTargets()) {
		if (!bib::in(tar, ret)) {
			std::cerr << __PRETTY_FUNCTION__
					<< ", warning, didn't indicate len cuts for " << tar
					<< " cut offs will be automatically determined based on median read length for reads found for target"
					<< "\n";
		}
	}

	if (failedOverlapStatusProcessing) {
		throw std::runtime_error { errorStream.str() };
	}
	for(const auto & tlenCut : ret){
		targets_.at(tlenCut.first).addLenCutOff(
				tlenCut.second.minLenChecker_.minLen_,
				tlenCut.second.maxLenChecker_.maxLen_);
	}
}

void PrimersAndMids::addOverLapStatuses(const bfs::path & overlapStatuses){
	table overlapStatusTab(overlapStatuses, "whitespace", true);
	std::unordered_map<std::string, PairedReadProcessor::ReadPairOverLapStatus> targetStatus;
	bib::for_each(overlapStatusTab.columnNames_,
			[](std::string & str) {stringToLower(str);});
	overlapStatusTab.setColNamePositions();
	overlapStatusTab.checkForColumnsThrow(VecStr{"target", "status"}, __PRETTY_FUNCTION__);
	bool failedOverlapStatusProcessing = false;
	std::stringstream errorStream;
	for (const auto & row : overlapStatusTab.content_) {
		if (!bib::in(row[overlapStatusTab.getColPos("target")], targets_)) {
			failedOverlapStatusProcessing = true;
			errorStream << __PRETTY_FUNCTION__ << ", error found "
					<< row[overlapStatusTab.getColPos("target")]
					<< " but it isn't a listed target in " << idFile_ << "\n";
			errorStream << "options are: "
					<< bib::conToStr(getVectorOfMapKeys(targets_), ", ") << "\n";
		}else 	if (bib::in(row[overlapStatusTab.getColPos("target")], targetStatus)) {
			failedOverlapStatusProcessing = true;
			errorStream << __PRETTY_FUNCTION__ << ", error already have "
					<< row[overlapStatusTab.getColPos("target")] << " in table" << "\n";
		} else {
			if ("NOOVERLAP"
					== bib::strToUpperRet(row[overlapStatusTab.getColPos("status")])) {
				targetStatus[row[overlapStatusTab.getColPos("target")]] =
						PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP;
			} else if ("R1BEGINSINR2"
					== bib::strToUpperRet(row[overlapStatusTab.getColPos("status")])) {
				targetStatus[row[overlapStatusTab.getColPos("target")]] =
						PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2;
			} else if ("R1ENDSINR2"
					== bib::strToUpperRet(row[overlapStatusTab.getColPos("status")])) {
				targetStatus[row[overlapStatusTab.getColPos("target")]] =
						PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2;
			} else {
				failedOverlapStatusProcessing = true;
				errorStream << __PRETTY_FUNCTION__
						<< ", error status should be NOOVERLAP, R1StartsInR2, or R1EndsInR2, not "
						<< bib::strToUpperRet(row[overlapStatusTab.getColPos("status")])
						<< "\n";
			}
		}
	}
	for(const auto & tar : getTargets()){
		if(!bib::in(tar, targetStatus)){
			failedOverlapStatusProcessing = true;
			errorStream << __PRETTY_FUNCTION__  << ", error, didn't indicate overlap status for " << tar << "\n";
		}
	}

	if (failedOverlapStatusProcessing) {
		throw std::runtime_error { errorStream.str() };
	}

	for(const auto & ts : targetStatus){
		targets_.at(ts.first).overlapStatus_ = ts.second;
	}
}
void PrimersAndMids::setRefSeqsKInfos(uint32_t klen, bool setRevComp){
	for(auto & tar : targets_){
		tar.second.setRefKInfos(klen, setRevComp);
	}
}

void PrimersAndMids::addRefSeqs(const bfs::path & refSeqsDir){
	auto fastaFiles = bib::files::listAllFiles(refSeqsDir.string(), false, {
			std::regex { R"(.*\.fasta$)" } });
	VecStr found;
	auto tars = getTargets();
	bool foundOtherRefSeqs = false;
	std::stringstream errorStream;
	std::set<std::string> alreadyAdded;
	for (const auto & ff : fastaFiles) {
		auto tarName = bfs::basename(ff.first);
		found.emplace_back(tarName);
		if (!bib::in(tarName, targets_)) {
			foundOtherRefSeqs = true;
			errorStream << __PRETTY_FUNCTION__ << ", error found " << tarName << " but it isn't a listed target in " << idFile_ << "\n";
			errorStream << "options are: " << bib::conToStr(getVectorOfMapKeys(targets_), ", ") << "\n";
		} else if (bib::in(tarName, alreadyAdded)) {
			foundOtherRefSeqs = true;
			errorStream << __PRETTY_FUNCTION__ << ", error already have " << tarName << " in table" << "\n";
		} else {
			auto inOpts = SeqIOOptions::genFastaIn(ff.first.string(), false);
			SeqInput reader(inOpts);
			auto seqs = reader.readAllReads<seqInfo>();
			if (seqs.size() == 1) {
				targets_.at(tarName).addSingleRef(seqs.front());
			} else {
				targets_.at(tarName).addMultileRef(seqs);
			}
		}
	}
	for (const auto & tar : tars) {
		if (!bib::in(tar, found)) {
			std::cerr << __PRETTY_FUNCTION__
					<< ", warning, didn't find any ref sequences for " << tar
					<< " when loading sequences for others though" << "\n";
		}
	}
	if (foundOtherRefSeqs) {
		//for now leaving this out of reporting as projects with mixed targets will generate what is probably unnecessary warnings
		//std::cerr << errorStream.str() << std::endl;
	}
}


void PrimersAndMids::addDefaultLengthCutOffs(uint32_t minLength, uint32_t maxLength){
	for(auto & tar : targets_ ){
		if(nullptr == tar.second.lenCuts_){
			tar.second.addLenCutOff(minLength, maxLength);
		}
	}
}


}  // namespace bibseq
