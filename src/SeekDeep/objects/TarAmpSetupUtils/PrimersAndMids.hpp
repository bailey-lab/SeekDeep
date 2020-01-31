#pragma once

/*
 * PrimersAndMids.hpp
 *
 *  Created on: Nov 27, 2016
 *      Author: nick
 */
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
#include <njhseq.h>

#include "SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp"

namespace njhseq {

class PrimersAndMids {
public:

	struct InitPars{
		InitPars(){
			mPars_.allowableErrors_ = 1;
		}
		bfs::path idFile_;

		bfs::path lenCutOffFilename_ = "";

	  bfs::path comparisonSeqFnp_ = "";
	  uint32_t compKmerLen_ = 5;
	  double compKmerSimCutOff_ = 0.50;

	  MidDeterminator::MidDeterminePars mPars_;

	  bfs::path overlapStatusFnp_ = "";
	  bool noOverlapProcessForNoOverlapStatusTargets_ = false;

	};

	class Target {
	public:
		class lenCutOffs {
		public:
			lenCutOffs(uint32_t minLen, uint32_t maxLen, bool mark = true);
			ReadCheckerLenAbove minLenChecker_;
			ReadCheckerLenBelow maxLenChecker_;
		};
		Target(const std::string & name, const std::string & forPrimer, const std::string & revPrimer);

		PrimerDeterminator::primerInfo info_;
		std::vector<seqInfo> refs_;
		std::vector<kmerInfo> refKInfos_;

		std::shared_ptr<lenCutOffs> lenCuts_;

		void addLenCutOff(uint32_t minLen, uint32_t maxLen, bool mark = true);

		void addSingleRef(const seqInfo & ref);
		void addMultileRef(const std::vector<seqInfo> & refs);

		void setRefKInfos(uint32_t klen, bool setRevComp);

		std::vector<PairedReadProcessor::ReadPairOverLapStatus> overlapStatuses_ {PairedReadProcessor::ReadPairOverLapStatus::NONE};

	};

	PrimersAndMids(const bfs::path & idFileFnp);
	PrimersAndMids(const std::unordered_map<std::string, Target> & targets);


	void checkIfMIdsOrPrimersReadInThrow(const std::string & funcName) const;

	const bfs::path idFile_;

	std::unordered_map<std::string, Target> targets_;
	std::unordered_map<std::string, MidDeterminator::MID> mids_;

	std::unique_ptr<MidDeterminator> mDeterminator_;
	std::unique_ptr<PrimerDeterminator> pDeterminator_;

	void initAllAddLenCutsRefs(const InitPars & pars);

	void initMidDeterminator(const MidDeterminator::MidDeterminePars & midSearchPars);
	void initPrimerDeterminator();

	bool hasTarget(const std::string & target) const;

	VecStr getTargets() const;

	VecStr getMids() const;

	bool hasMid(const std::string & mid) const;
	void addTarget(const std::string & primerName, const std::string & forPrimer,
			const std::string & revPrimer);

	void addMid(const std::string & midNmae, const std::string & barcode);
	void addMid(const std::string & midNmae, const std::string & forBarcode, const std::string & revBarcode);

	bool hasMultipleTargets() const;

	bool containsMids() const;
	bool containsTargets() const;
	bool screeningForPossibleContamination() const;

	void writeIdFile(const OutOptions & outOpts) const;
	void writeIdFile(const OutOptions & outOpts, const VecStr & targets) const;

	table genLenCutOffs(const VecStr & targets) const;
	table genOverlapStatuses(const VecStr & targets) const;



	std::vector<seqInfo> getRefSeqs(const VecStr & targets) const;
	void checkMidNamesThrow() const;

	static std::map<std::string, PrimersAndMids::Target::lenCutOffs> readInLenCutOffs(
			const bfs::path & lenCutOffsFnp);

	void addLenCutOffs(const bfs::path & lenCutOffsFnp);
	void addOverLapStatuses(const bfs::path & overlapStatuses);
	void addOverLapStatuses(const std::vector<PairedReadProcessor::ReadPairOverLapStatus> & allStatus);

	void addRefSeqs(const bfs::path & refDirectory);
	void setRefSeqsKInfos(uint32_t klen, bool setRevComp);

	void addDefaultLengthCutOffs(uint32_t minLength, uint32_t maxLength);

};

}  // namespace njhseq

