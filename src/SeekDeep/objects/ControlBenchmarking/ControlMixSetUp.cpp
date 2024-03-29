/*
 * ControlMixSetUp.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */


#include "ControlMixSetUp.hpp"
#include <njhseq/objects/dataContainers/tables.h>
namespace njhseq {

ControlMixSetUp::ControlMixSetUp(const std::string & name,
		const std::unordered_map<std::string, double> & relativeAbundances) :
		name_(name), rawRelativeAbundances_(relativeAbundances) {
	if (relativeAbundances.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< " relativeAbundances can't be empty" << "\n";
		throw std::runtime_error { ss.str() };
	}

	double total = 0;
	for (const auto & relAbund : rawRelativeAbundances_) {
		if (relAbund.second <= 0) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error "
					<< " abundances need to be greater than zero, " << relAbund.first
					<< ": " << relAbund.second << "\n";
			throw std::runtime_error { ss.str() };
		}
		total += relAbund.second;
	}
	for (const auto & relAbund : rawRelativeAbundances_) {
		relativeAbundances_[relAbund.first] = relAbund.second/total;
	}

}

void ControlMixSetUp::removeStrain(const std::string & name){
	if(!njh::in(name, rawRelativeAbundances_) || !njh::in(name, relativeAbundances_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no strain named: " << name << " options are: " << njh::conToStr(njh::getVecOfMapKeys(rawRelativeAbundances_), ",")<< "\n";
		throw std::runtime_error{ss.str()};
	}
	rawRelativeAbundances_.erase(name);
	relativeAbundances_.erase(name);

	double total = 0;
	for (const auto & relAbund : rawRelativeAbundances_) {
		if (relAbund.second <= 0) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error "
					<< " abundances need to be greater than zero, " << relAbund.first
					<< ": " << relAbund.second << "\n";
			throw std::runtime_error { ss.str() };
		}
		total += relAbund.second;
	}
	for (const auto & relAbund : rawRelativeAbundances_) {
		relativeAbundances_[relAbund.first] = relAbund.second/total;
	}
}

VecStr ControlMixSetUp::getStrains() const {
	return njh::getVecOfMapKeys(rawRelativeAbundances_);
}

std::unordered_map<std::string, ControlMixSetUp> ControlMixSetUp::readInSetUps(
		const bfs::path & mixtureSetUpFnp) {
	//read in mixture setup
	table mixtureSetupTab(mixtureSetUpFnp, "\t", true);
	mixtureSetupTab.checkForColumnsThrow(VecStr { "MixName", "strain", "relative_abundance" }, __PRETTY_FUNCTION__);

	if (mixtureSetupTab.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << mixtureSetUpFnp << " is empty "
				<< "\n";
		throw std::runtime_error { ss.str() };
	}

	std::unordered_map<std::string, ControlMixSetUp> mixSetups;
	std::unordered_map<std::string, std::unordered_map<std::string, double>> mixInfos;
	std::unordered_map<std::string, std::unordered_map<std::string, double>> mixInfosWithZeros;

	for (const auto & row : mixtureSetupTab) {
		std::string mixname = row[mixtureSetupTab.getColPos("MixName")];
		std::string strain = row[mixtureSetupTab.getColPos("strain")];
		double relative_abundance = njh::StrToNumConverter::stoToNum<double>(
				row[mixtureSetupTab.getColPos("relative_abundance")]);
		if (njh::in(strain, mixInfos[mixname])) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have " << strain
					<< " for mixture " << mixname << "\n";
			throw std::runtime_error { ss.str() };
		}
		if (0 != relative_abundance) {
			mixInfos[mixname][strain] = relative_abundance;
		} else {
			mixInfosWithZeros[mixname][strain] = relative_abundance;
		}
	}
	std::set<std::string> mixturesNamesWithAllZeros;
	for (const auto & mixInfo : mixInfosWithZeros) {
		if (!njh::in(mixInfo.first, mixInfos)) {
			mixturesNamesWithAllZeros.emplace(mixInfo.first);
		}
	}
	if (!mixturesNamesWithAllZeros.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "the following mixtures had strains all with relative abundances of 0: "
				<< njh::conToStrEndSpecial(mixturesNamesWithAllZeros, ", ", " and ")
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	for (const auto & mixInfo : mixInfos) {
		mixSetups.emplace(mixInfo.first,
				ControlMixSetUp(mixInfo.first, mixInfo.second));
	}
	return mixSetups;
}

}  // namespace njhseq





