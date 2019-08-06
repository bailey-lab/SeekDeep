/*
 * ControlBencher.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */



#include "ControlBencher.hpp"
namespace njhseq {
ControlBencher::ControlBencher(const ControlBencherPars & pars):pars_(pars){
	//read in mixture set ups
	mixSetups_  = ControlMixSetUp::readInSetUps(pars_.mixSetUpFnp_);

	//read in samples to mixture
	table sampleToMixtureTab(pars_.samplesToMixFnp_, "\t", true);
	sampleToMixtureTab.checkForColumnsThrow(VecStr{"sample", "MixName"}, __PRETTY_FUNCTION__);
	if(sampleToMixtureTab.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << pars_.samplesToMixFnp_ << " is empty " << "\n";
		throw std::runtime_error{ss.str()};
	}
	for(const auto & row : sampleToMixtureTab){
		samplesToMix_[row[sampleToMixtureTab.getColPos("sample")]] = row[sampleToMixtureTab.getColPos("MixName")];
	}
	VecStr missingMixs;
	for(const auto & sampToMix : samplesToMix_){
		if(!njh::in(sampToMix.second, mixSetups_)){
			missingMixs.emplace_back(sampToMix.second);
		}
	}
	if(!missingMixs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " error missing the following mixture information " << njh::conToStr(missingMixs, ",") << " from " << pars_.mixSetUpFnp_ << "\n";
		throw std::runtime_error{ss.str()};
	}
}

VecStr ControlBencher::getSamples() const {
	auto ret = getVectorOfMapKeys(samplesToMix_);
	njh::sort(ret);
	return ret;
}


void ControlBencher::removeStrain(const std::string & name){
	for(auto & mix : mixSetups_){
		mix.second.removeStrain(name);
	}
}

void ControlBencher::removeStrains(const VecStr & names){
	for(const auto & name : names){
		removeStrain(name);
	}
}

void ControlBencher::checkForStrainsThrow(const std::set<std::string> & names,
		const std::string & funcName) const {
	VecStr missingStrains;
	auto allStrains = getAllStrains();
	for (const auto & strain : allStrains) {
		if (!njh::in(strain, names)) {
			missingStrains.emplace_back(strain);
		}
	}
	if (!missingStrains.empty()) {
		std::stringstream ss;
		ss << funcName << ", error "
				<< "the following strains are missing from expected seqs: "
				<< njh::conToStrEndSpecial(missingStrains, ", ", " and ") << "\n";
		throw std::runtime_error { ss.str() };
	}
}

std::set<std::string> ControlBencher::getAllStrains() const {
	std::set < std::string > allStrains;
	for (const auto & mixSetup : mixSetups_) {
		auto strains = mixSetup.second.getStrains();
		allStrains.insert(strains.begin(), strains.end());
	}
	return allStrains;
}


}  // namespace njhseq

