/*
 * TarAmpPEAnalysisSetup.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */


#include "TarAmpPEAnalysisSetup.hpp"


namespace bibseq {

TarAmpPEAnalysisSetup::TarAmpPEAnalysisSetup(const bfs::path & dir) :
		dir_(dir) {
	bib::files::makeDir(bib::files::MkdirPar(dir_.string()));
	infoDir_ = bib::files::makeDir(dir_.string(), bib::files::MkdirPar("info"));
	logsDir_ = bib::files::makeDir(dir_.string(), bib::files::MkdirPar("logs"));
	idsDir_ = bib::files::makeDir(infoDir_.string(), bib::files::MkdirPar("ids"));
	reportsDir_ = bib::files::makeDir(dir_.string(),
			bib::files::MkdirPar("reports"));
	serverConfigsDir_ = bib::files::makeDir(dir_.string(), bib::files::MkdirPar("serverConfigs"));
}


TarAmpPEAnalysisSetup::Sample::Sample(const std::string & name) :
		name_(name) {
}

void TarAmpPEAnalysisSetup::Sample::addRep(const std::string & rep) {
	if (bib::in(rep, reps_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, sample " << name_
				<< " already has rep name " << rep << "\n";
		throw std::runtime_error { ss.str() };
	}
	reps_.emplace_back(rep);
}

void TarAmpPEAnalysisSetup::Sample::addReps(const VecStr & reps) {
	for (const auto & rep : reps) {
		addRep(rep);
	}
}

VecStr TarAmpPEAnalysisSetup::Sample::getReps() const {
	return reps_;
}

TarAmpPEAnalysisSetup::Samples::Samples(const std::string & target): target_(target){};


bool TarAmpPEAnalysisSetup::Samples::hasSample(const std::string & sample) {
	return samples_.end() != samples_.find(sample);
}

void TarAmpPEAnalysisSetup::Samples::addSample(const std::string & sample) {
	if (hasSample(sample)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, already have " << sample << " for " << target_ << "\n";
		ss << "With input: " << bib::conToStr(samples_.at(sample).reps_ , ", ")<< "\n";
		throw std::runtime_error { ss.str() };
	}
	samples_.emplace(sample, Sample{sample});
}

void TarAmpPEAnalysisSetup::Samples::addSample(const std::string & sample, const VecStr & reps) {
	addSample(sample);
	samples_.find(sample)->second.addReps(reps);
}

VecStr TarAmpPEAnalysisSetup::Samples::getSamples() const {
	return bib::getVecOfMapKeys(samples_);
}

std::vector<std::string> TarAmpPEAnalysisSetup::Samples::getReps() const {
	std::vector<std::string> ret;
	for (const auto & samp : samples_) {
		ret.insert(ret.end(), samp.second.reps_.begin(), samp.second.reps_.end());
	}
	return ret;
}

std::set<std::string> TarAmpPEAnalysisSetup::getSamples() const {
	std::set<std::string> ret;
	for (const auto & tar : samplesForTargets_) {
		auto samps = tar.second.getSamples();
		ret.insert(samps.begin(), samps.end());
	}
	return ret;
}

std::vector<std::string> TarAmpPEAnalysisSetup::getReps() const {
	std::set<std::string> ret;
	for (const auto & tar : samplesForTargets_) {
		auto reps = tar.second.getReps();
		ret.insert(reps.begin(), reps.end());
	}
	return VecStr{ret.begin(), ret.end()};
}

void TarAmpPEAnalysisSetup::addGroupingMetaData(const bfs::path & groupingsFileFnp) {
	groupMetaData_ = std::make_unique<MultipleGroupMetaData>(
			bib::files::normalize(groupingsFileFnp), getSamples());
}

void TarAmpPEAnalysisSetup::addGroupingsFile()const{
	if(nullptr != groupMetaData_){
		bfs::copy(groupMetaData_->groupingsFile_, bib::files::make_path(infoDir_, "groupingMetaData.tab.txt"));
	}
}

void TarAmpPEAnalysisSetup::writeSampleNamesFile() const {
	std::ofstream sampleNamesFile;
	openTextFile(sampleNamesFile,
			OutOptions(bib::files::make_path(infoDir_, "sampNames").string(), ".tab.txt"));
	auto tarKeys = getVectorOfMapKeys(samplesForTargets_);
	bib::sort(tarKeys);
	for (const auto & tarKey : tarKeys) {
		auto sampKeys = bib::getVecOfMapKeys(samplesForTargets_.at(tarKey).samples_);
		bib::sort(sampKeys);
		for (const auto & sampKey : sampKeys) {
			sampleNamesFile << tarKey << "\t" << sampKey;
			for (auto rep : samplesForTargets_.at(tarKey).samples_.at(sampKey).reps_) {
				if(!bib::beginsWith(rep, "MID")){
					rep = "MID" + rep;
				}
				sampleNamesFile << "\t" << rep;
			}
			sampleNamesFile << std::endl;
		}
	}
}

VecStr TarAmpPEAnalysisSetup::getTargets() const {
	return bib::getVecOfMapKeys(samplesForTargets_);
}

void TarAmpPEAnalysisSetup::addSamplesNames(const bfs::path & samplesNamesFnp){
	VecStr warnings;
	table samplesNamesInputTab(samplesNamesFnp.string());
	std::vector<VecStr> finalSamps;
	uint32_t rowCount = 0;
	for(const auto & row : samplesNamesInputTab.content_){
		++rowCount;
		if(row.empty()){
			continue;
		}
		if(0 != row.front().size() && '#' == row.front().front()){
			continue;
		}
		if(row.size() < 3){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": error each row should have at least 3 columns, not " << row.size() << " for row \n";
			ss << "Row Number(1 based pos): " << rowCount << " - " << bib::conToStr(row, "\t") << "\n";
			warnings.emplace_back(ss.str());
			continue;
		}
		std::vector<std::string> nonEmptyRow;
		for(const auto & colPos : iter::range(row.size())){
			if(colPos < 2 && row[colPos].empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error column " << colPos << " can't be an empty string" << "\n";
				ss << "Row Number(1 based pos): " << rowCount << " - " << bib::conToStr(row, "\t") << "\n";
				warnings.emplace_back(ss.str());
				continue;
			}
			if(!row[colPos].empty()){
				nonEmptyRow.emplace_back(row[colPos]);
			}
		}
		finalSamps.emplace_back(nonEmptyRow);
	}
	if(!warnings.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": There are " << warnings.size() << " errors." << "\n";
		for(const auto & warn : warnings){
			ss << warn << std::endl;
		}
		throw std::runtime_error{ss.str()};
	}

	for (const auto & row : finalSamps) {
		std::string target = row[0];
		std::string sample = row[1];
		VecStr reps = getSubVector(row, 2);
		for (const auto & rep : reps) {
			for (const auto & tar : samplesForTargets_) {
				if(tar.first != target){
					continue;
				}
				for (const auto & samp : tar.second.samples_) {
					if (bib::in(rep, samp.second.reps_)) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << " : Error, already contains " << rep
								<< " in " << samp.first << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			}
		}
		auto search = samplesForTargets_.find(target);
		if(samplesForTargets_.end() == search){
			samplesForTargets_.emplace(target, Samples(target));
		}
		samplesForTargets_.at(target).addSample(sample, reps);
	}
}


}  // namespace bibseq

