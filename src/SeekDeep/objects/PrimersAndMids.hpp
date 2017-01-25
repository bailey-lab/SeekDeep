#pragma once

/*
 * PrimersAndMids.hpp
 *
 *  Created on: Nov 27, 2016
 *      Author: nick
 */

#include <bibseq.h>

namespace bibseq {

class PrimersAndMids {
public:

	class Target {
	public:
		class lenCutOffs {
		public:
			lenCutOffs(uint32_t minLen, uint32_t maxLen, bool mark = true);
			ReadCheckerLenAbove minLenChecker_;
			ReadCheckerLenBelow maxLenChecker_;
		};
		Target(const std::string & name, const std::string & forPrimer,
				const std::string & revPrimer);

		PrimerDeterminator::primerInfo info_;
		std::vector<seqInfo> refs_;

		std::unique_ptr<lenCutOffs> lenCuts_;

		void addLenCutOff(uint32_t minLen, uint32_t maxLen, bool mark = true);

		void addSingleRef(const seqInfo & ref);
		void addMultileRef(const std::vector<seqInfo> & refs);
	};

	PrimersAndMids(const bfs::path & idFileFnp);

	bfs::path idFile_;

	std::unordered_map<std::string, Target> targets_;
	std::unordered_map<std::string, MidDeterminator::MidInfo> mids_;

	bool hasTarget(const std::string & target) const;

	VecStr getTargets() const;

	VecStr getMids() const;

	bool hasMid(const std::string & mid) const;
	void addTarget(const std::string & primerName, const std::string & forPrimer,
			const std::string & revPrimer);

	void addMid(const std::string & midNmae, const std::string & barcode);

	bool hasMultipleTargets() const;

	bool containsMids() const;

	void writeIdFile(const OutOptions & outOpts) const;

	void writeIdFile(const OutOptions & outOpts, const VecStr & targets) const;

	table genLenCutOffs(const VecStr & targets) const;

	std::vector<seqInfo> getRefSeqs(const VecStr & targets) const;

	void checkMidNamesThrow() const;


	static std::map<std::string, PrimersAndMids::Target::lenCutOffs> readInLenCutOffs(
			const bfs::path & lenCutOffsFnp);

};

}  // namespace bibseq

