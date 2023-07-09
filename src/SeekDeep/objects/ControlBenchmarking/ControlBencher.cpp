/*
 * ControlBencher.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */



#include "ControlBencher.hpp"

#include <njhseq/objects/dataContainers/tables.h>

#include <utility>

namespace njhseq {
ControlBencher::ControlBencher(ControlBencherPars pars):pars_(std::move(pars)){
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
		if(njh::in(name, mix.second.getStrains())){
			mix.second.removeStrain(name);
		}
	}
}

void ControlBencher::removeStrains(const VecStr & names){
	for(const auto & name : names){
		removeStrain(name);
	}
}

void ControlBencher::writeMixSetUpsInSamples(const OutOptions &outOptions) {
	OutputStream out(outOptions);
	std::set<std::string> mixsPresent;
	for (const auto &sampToMix: samplesToMix_) {
		mixsPresent.emplace(sampToMix.second);
	}
	out << "MixName\tstrain\trelative_abundance" << std::endl;
	for (const auto &mix: mixsPresent) {
		for (const auto &strain: njh::mapAt(mixSetups_, mix).relativeAbundances_) {
			out << mix << "\t" << strain.first << "\t" << strain.second << std::endl;
		}
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


double ControlBencher::benchResults::RMSE() const{
	return 0 == sumOfSquares_ ? 0 : std::sqrt(sumOfSquares_/static_cast<double>(recoveredHaps_));
}
uint32_t ControlBencher::benchResults::totalHaps() const {
	return falseHaps_ + recoveredHaps_;
}
double ControlBencher::benchResults::falseHapRate() const{
	return totalHaps() > 0 ? falseHaps_/static_cast<double>(totalHaps()) : 0;
}
double ControlBencher::benchResults::hapRecoveryRate() const{
	return expectedHapCnt_ > 0 ? recoveredHaps_/static_cast<double>(expectedHapCnt_) : 0;
}


}  // namespace njhseq

