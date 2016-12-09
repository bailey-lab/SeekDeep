/*
 * PrimersAndMids.cpp
 *
 *  Created on: Nov 27, 2016
 *      Author: nick
 */

#include "PrimersAndMids.hpp"

namespace bibseq {

PrimersAndMids::Target::lenCutOffs::lenCutOffs(uint32_t minLen, uint32_t maxLen,
		bool mark) :
		minLenChecker_(ReadCheckerLenAbove(minLen, mark)), maxLenChecker_(
				ReadCheckerLenBelow(maxLen, mark)) {
}


PrimersAndMids::Target::Target(const std::string & name,
		const std::string & forPrimer, const std::string & revPrimer) :
		info_(name, forPrimer, revPrimer) {
}

void PrimersAndMids::Target::addLenCutOff(uint32_t minLen, uint32_t maxLen,
		bool mark) {
	lenCuts_ = std::make_unique<lenCutOffs>(minLen, maxLen, mark);
}

void PrimersAndMids::Target::addSingleRef(const seqInfo & ref) {
	if (info_.primerPairName_ == ref.name_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": error, adding ref that doesn't the target name" << "\n";
		ss << "TargetName: " << info_.primerPairName_ << "\n";
		ss << "AddingRefName: " << ref.name_ << "\n";
	}
}
void PrimersAndMids::Target::addMultileRef(const std::vector<seqInfo> & ref) {
	//mutliple refs so no name matchig
	addOtherVec(refs_, ref);
}

PrimersAndMids::PrimersAndMids(const bfs::path & idFileFnp) :
		idFile_(idFileFnp) {
	if (!bfs::exists(idFile_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, "
				<< bib::bashCT::boldRed(idFile_.string()) << " doesn't exist\n";
		throw std::runtime_error { ss.str() };
	}
	auto firstLine = bib::files::getFirstLine(idFileFnp.string());
	bib::strToLower(firstLine);
	if (!bib::beginsWith(firstLine, "gene")
			&& !bib::beginsWith(firstLine, "target")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error the id file "
				<< bib::bashCT::boldRed(idFile_.string())
				<< " should start with either target or gene (case insensitive)\n";
		ss << "line: " << firstLine << "\n";
		throw std::runtime_error { ss.str() };
	}
	bool readingPrimers = false;
	bool readingMids = false;
	std::ifstream idFile(idFileFnp.string());
	std::string line = "";
	while (bib::files::crossPlatGetline(idFile, line)) {
		auto lowerLine = bib::strToLowerRet(line);
		if (bib::beginsWith(lowerLine, "gene")
				|| bib::beginsWith(lowerLine, "target")) {
			readingPrimers = true;
			readingMids = false;
		} else if (bib::beginsWith(lowerLine, "id")) {
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
	return bib::getVecOfMapKeys(targets_);
}

VecStr PrimersAndMids::getMids() const {
	return bib::getVecOfMapKeys(mids_);
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

}  // namespace bibseq
