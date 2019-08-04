#pragma once

/*
 * ControlBencher.hpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */


#include <njhseq.h>
#include "SeekDeep/objects/ControlBenchmarking/ControlMixSetUp.hpp"

namespace njhseq {

class ControlBencher {
public:
	struct ControlBencherPars{
		bfs::path mixSetUpFnp_;
		bfs::path samplesToMixFnp_;
	};

	ControlBencher(const ControlBencherPars & pars);

	const ControlBencherPars pars_;

	std::unordered_map<std::string, ControlMixSetUp> mixSetups_;
	std::unordered_map<std::string, std::string> samplesToMix_;


	std::set<std::string> getAllStrains() const;
	VecStr getSamples() const;


	void checkForStrainsThrow(const std::set<std::string> & names,
			const std::string & funcName) const;

	struct benchResults{
		std::unordered_map<std::string, std::string> resSeqToExpSeq_;
		uint32_t recoveredHaps_ = 0;
		uint32_t falseHaps_ = 0;
		uint32_t expectedHapCnt_ = 0;
		double sumOfSquares_ = 0;

		double RMSE() const{
			return 0 == sumOfSquares_ ? 0 : std::sqrt(sumOfSquares_/static_cast<double>(recoveredHaps_));
		}
		uint32_t totalHaps() const {
			return falseHaps_ + recoveredHaps_;
		}
		double falseHapRate() const{
			return falseHaps_/static_cast<double>(totalHaps());
		}
		double hapRecoveryRate() const{
			return recoveredHaps_/static_cast<double>(expectedHapCnt_);
		}

		std::map<std::string, std::map<std::string, comparison>> falseHapsCompsToExpected;
		std::map<std::string, std::map<std::string, comparison>> falseHapsCompsToOthers;

	};


	template<typename RESULTSEQ, typename REFSEQ>
	static benchResults benchmark(const std::vector<RESULTSEQ> & resultSeqs, const std::vector<REFSEQ> & expectedSeqs,
	const std::unordered_map<std::string, double> & expectedSeqFracs,
	const std::unordered_map<std::string, std::string> & expectedSeqNameKey){
		benchResults ret;
		ret.expectedHapCnt_ = expectedSeqs.size();
		for(const auto & resultSeq : resultSeqs){
			const auto & resSeq = getSeqBase(resultSeq);
			std::string matchingRef = "";
			double expectedFrac = 0;
			for(const auto & expectedSeq : expectedSeqs){
				const auto & expSeq = getSeqBase(expectedSeq);
				if(resSeq.seq_ == expSeq.seq_){
					matchingRef =  njh::mapAt(expectedSeqNameKey,expSeq.name_);
					expectedFrac = njh::mapAt(expectedSeqFracs,expSeq.name_);
					break;
				}
			}
			if("" != matchingRef){
				++ret.recoveredHaps_;
				ret.sumOfSquares_+= std::pow(resultSeq.frac_ - expectedFrac, 2.0);
			}else{
				++ret.falseHaps_;
			}
			ret.resSeqToExpSeq_[resSeq.name_] = matchingRef;
		}
		return ret;
	}

	template<typename RESULTSEQ, typename REFSEQ>
	static benchResults benchmark(const std::vector<RESULTSEQ> & resultSeqs, const std::vector<REFSEQ> & expectedSeqs,
	const std::unordered_map<std::string, double> & expectedSeqFracs,
	const std::unordered_map<std::string, std::string> & expectedSeqNameKey,
	aligner & alignerObj){
		benchResults ret;
		ret.expectedHapCnt_ = expectedSeqs.size();
		for(const auto & resultSeq : resultSeqs){
			const auto & resSeq = getSeqBase(resultSeq);
			std::string matchingRef = "";
			double expectedFrac = 0;
			for(const auto & expectedSeq : expectedSeqs){
				const auto & expSeq = getSeqBase(expectedSeq);
				if(resSeq.seq_ == expSeq.seq_){
					matchingRef =  njh::mapAt(expectedSeqNameKey,expSeq.name_);
					expectedFrac = njh::mapAt(expectedSeqFracs,expSeq.name_);
					break;
				}
			}
			if("" != matchingRef){
				++ret.recoveredHaps_;
				ret.sumOfSquares_+= std::pow(resultSeq.frac_ - expectedFrac, 2.0);
			}else{
				++ret.falseHaps_;
				for(const auto & exp : expectedSeqs){
					alignerObj.alignCacheGlobal(exp, resultSeq);
					alignerObj.profileAlignment(exp, resultSeq, false, false, false);
					ret.falseHapsCompsToExpected[resSeq.name_][getSeqBase(exp).name_] = alignerObj.comp_;
				}
				for(const auto & otherRes : resultSeqs){
					if(getSeqBase(otherRes).name_ == getSeqBase(resultSeq).name_){
						continue;
					}
					alignerObj.alignCacheGlobal(otherRes, resultSeq);
					alignerObj.profileAlignment(otherRes, resultSeq, false, false, false);
					ret.falseHapsCompsToOthers[resSeq.name_][getSeqBase(otherRes).name_] = alignerObj.comp_;
				}
			}
			ret.resSeqToExpSeq_[resSeq.name_] = matchingRef;
		}
		return ret;
	}

};


}  // namespace njhseq



