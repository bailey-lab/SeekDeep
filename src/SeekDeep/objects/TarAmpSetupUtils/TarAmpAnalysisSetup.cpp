/*
 * TarAmpPEAnalysisSetup.cpp
 *
 *  Created on: Nov 25, 2016
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
#include "TarAmpAnalysisSetup.hpp"
#include "SeekDeep/objects/IlluminaUtils/ReadPairsOrganizer.hpp"
#include "SeekDeep/objects/TarAmpSetupUtils/SampleFileNameGenerator.hpp"

namespace njhseq {





bool TarAmpAnalysisSetup::TarAmpPars::checkForOutDir(VecStr & warnings) const {
	if (bfs::exists(outDir)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, out directory " << "\""
				<< njh::bashCT::boldRed(outDir.string()) << "\"" << " already exists"
				<< ", must remove by hand will not overwrite\n";
		warnings.emplace_back(ss.str());
		return false;
	}
	return true;
}

bool TarAmpAnalysisSetup::TarAmpPars::checkIfFnpExists(const bfs::path & fnp,
		VecStr & warnings) {
	if (!bfs::exists(fnp)) {
		std::stringstream ss;
		ss << njh::bashCT::boldRed(fnp.string()) << " doesn't exist" << "\n";
		warnings.emplace_back(ss.str());
		return false;
	}
	return true;
}


bool TarAmpAnalysisSetup::TarAmpPars::techIs454() const {
	return "454" == technology;
}

bool TarAmpAnalysisSetup::TarAmpPars::techIsIllumina() const {
	return "illumina" == technology;
}

bool TarAmpAnalysisSetup::TarAmpPars::techIsIlluminaSingleEnd() const {
	return "illumina-singleend" == technology;
}

bool TarAmpAnalysisSetup::TarAmpPars::techIsIonTorrent() const {
	return "iontorrent" == technology;
}

bool TarAmpAnalysisSetup::TarAmpPars::teschIsNanopore() const {
	return "nanopore" == technology;
}


bool TarAmpAnalysisSetup::TarAmpPars::checkForRequiredFnpPars(
		VecStr & warnings) const {
	bool status = true;
	status = status && checkForOutDir(warnings);
	if("" != samplesNamesFnp){
		status = status && checkIfFnpExists(samplesNamesFnp, warnings);
	}
	status = status && checkIfFnpExists(inputDir, warnings);
	status = status && checkIfFnpExists(idFile, warnings);
	return status;
}

bool TarAmpAnalysisSetup::TarAmpPars::allChecks(VecStr & warnings) const {
	bool status = true;
	status = status && checkForRequiredFnpPars(warnings);
	status = status && checkForOptionalFnpPars(warnings);
	return status;
}



bool TarAmpAnalysisSetup::TarAmpPars::checkForOptionalFnpPars(VecStr & warnings) const{
	bool status = true;
	if(refSeqsDir != ""){
		status = status && checkIfFnpExists(refSeqsDir, warnings);
		if(bfs::exists(refSeqsDir) && !bfs::is_directory(refSeqsDir)){
			std::stringstream ss;
			ss << njh::bashCT::boldRed(refSeqsDir.string()) << " is not a directory, the ref seqs directory should be, well, a directory" << "\n";
			warnings.emplace_back(ss.str());
			status = false;
		}
	}
	if(lenCutOffsFnp != ""){
		status = status && checkIfFnpExists(lenCutOffsFnp, warnings);
	}
	if(groupMeta != ""){
		status = status && checkIfFnpExists(groupMeta, warnings);
	}
	return status;
}




TarAmpAnalysisSetup::TarAmpAnalysisSetup(const TarAmpPars & inputPars) :
		pars_(inputPars) {
	dir_ = pars_.outDir;
	njh::files::makeDir(njh::files::MkdirPar(dir_.string()));
	VecStr warnings;



	infoDir_ = njh::files::makeDir(dir_.string(), njh::files::MkdirPar("info"));
	logsDir_ = njh::files::makeDir(dir_.string(), njh::files::MkdirPar("logs"));
	idsDir_ = njh::files::makeDir(infoDir_.string(), njh::files::MkdirPar("ids"));
	refsDir_ = njh::files::makeDir(infoDir_.string(), njh::files::MkdirPar("refs"));
	reportsDir_ = njh::files::makeDir(dir_.string(),njh::files::MkdirPar("reports"));
	serverConfigsDir_ = njh::files::makeDir(dir_.string(), njh::files::MkdirPar("serverConfigs"));

	//add ids and primers
	idsMids_ = std::make_unique<PrimersAndMids>(pars_.idFile);

	if("" == pars_.samplesNamesFnp){
		if("" != pars_.samplesNamesWithBarcodeInfoFnp){
			//
			SampleFileNameGenerator fileGen(pars_.idFile, pars_.samplesNamesWithBarcodeInfoFnp);
			//sample name file
			OutOptions sampleNamesOutOpts(njh::files::make_path(infoDir_, "input_sampleNames.tab.txt"));
			fileGen.writeSampleNameFile(sampleNamesOutOpts);
			pars_.samplesNamesFnp = sampleNamesOutOpts.outFilename_;
			//barcode/primer file
			OutOptions idFileOutOpts(njh::files::make_path(infoDir_, "input_ids.tab.txt"));
			fileGen.writeBarcodePrimerFile(idFileOutOpts);
			pars_.idFile = idFileOutOpts.outFilename_;
		} else if("" != pars_.samplesNamesByLibraryNameFnp){
			SampleFileNameGenerator fileGen(pars_.idFile, pars_.samplesNamesByLibraryNameFnp, true);
			//sample name file
			OutOptions sampleNamesOutOpts(njh::files::make_path(infoDir_, "input_sampleNames.tab.txt"));
			fileGen.writeSampleNameFile(sampleNamesOutOpts);
			pars_.samplesNamesFnp = sampleNamesOutOpts.outFilename_;
		} else {
			//if samples fnp not set, try to guess samples
			auto samps = GuessPossibleSamps(pars_);
			auto guessSamplesFnp = njh::files::make_path(infoDir_, "guessedSamples.tab.txt");

			if(idsMids_->containsMids()){
				table guessedSamplesWithMids(VecStr{"#index", "sample", "run1"});
				auto sampleCol = samps.getColPos("sample");
				for(const auto & row : samps){
					for(const auto & mid : idsMids_->getMids()){
						guessedSamplesWithMids.addRow(row[sampleCol], njh::pasteAsStr(row[sampleCol], mid), mid);
					}
				}
				auto guessSamplesWithMidsFnp = njh::files::make_path(infoDir_, "guessedSamplesWithMids.tab.txt");
				OutputStream guessedSamplesWithMidsOut(guessSamplesWithMidsFnp);
				guessedSamplesWithMids.outPutContents(guessedSamplesWithMidsOut,"\t");
				pars_.samplesNamesFnp = guessSamplesWithMidsFnp;
			}else{
				OutputStream guessedSamplesOut(guessSamplesFnp);
				samps.outPutContents(guessedSamplesOut,"\t");
				pars_.samplesNamesFnp = guessSamplesFnp;
			}
		}
	}



	if(idsMids_->containsMids()){
		// set this automatically
		pars_.byIndex = true;
	}
	if(idsMids_->containsMids()){
		idsMids_->checkMidNamesThrow();
	}
	//add overlap status
	if("" != pars_.overlapStatusFnp){
		addOverlapStatus(pars_.overlapStatusFnp);
	}else if(!pars_.defaultStatuses_.empty()){
		idsMids_->addOverLapStatuses(pars_.defaultStatuses_);
	}


	//add sample names
	addSamplesNames(pars_.samplesNamesFnp);
	writeSampleNamesFile();
	if(pars_.byIndex){
		addIndexToTargetsNames(pars_.targetsToIndexFnp);
		if("" == pars_.targetsToIndexFnp){
			table indexToTargetsTab(VecStr{"index", "	targets"});
			for(const auto & indexToTar : indexToTars_){
				indexToTargetsTab.addRow(indexToTar.first, njh::conToStr(indexToTar.second, ","));
			}
			OutputStream indexToTargetsOut(OutOptions(njh::files::make_path(infoDir_, "indexToTargets.tab.txt")));
			indexToTargetsTab.outPutContents(indexToTargetsOut, "\t");
		}else{
			bfs::copy_file(pars_.targetsToIndexFnp, njh::files::make_path(infoDir_, "indexToTargets.tab.txt"));
		}
	}
	//add ref seqs if provided
	if("" != pars_.refSeqsDir){
		addRefSeqs(pars_.refSeqsDir);
	}
	//add len cut offs if provided
	if("" != pars_.lenCutOffsFnp){
		addLenCutOffs(pars_.lenCutOffsFnp);
	}

	//check for incompatible set ups
//	if (!idsMids_->containsMids() && pars_.byIndex) {
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__
//				<< ": error, if input files don't have MID still, then the sample file"
//				<< pars_.samplesNamesFnp << " should be by target and not index" << "\n";
//		throw std::runtime_error { ss.str() };
//	}
//
//	if (idsMids_->containsMids() && !pars_.byIndex) {
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__
//				<< ": error, if input files have MIDs still, then the sample file"
//				<< pars_.samplesNamesFnp << " should be by index and not by target" << "\n";
//		throw std::runtime_error { ss.str() };
//	}

	//add meta data if available
	if("" != pars_.groupMeta){
		addGroupingMetaData(pars_.groupMeta);
		bfs::copy_file(groupMetaData_->groupingsFile_,
				njh::files::make_path(infoDir_, "groupMeta.tab.txt"));
	}


	auto targets = getTargets();
	VecStr failedTargets;
	for(const auto & tar : targets){
		if(std::string::npos != tar.find("_")){
			failedTargets.emplace_back(tar);
		}
	}
	if(!failedTargets.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": error, target names shouldn't have underscores in them (this makes downstream patern matching and such easier) replace with -" << "\n";
		ss << "The following targets should be renamed" << "\n";
		ss << njh::conToStr(failedTargets, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}


}

void TarAmpAnalysisSetup::addIndexToTargetsNames(
		const bfs::path & targetsToIndexFnp) {
	if("" == targetsToIndexFnp){
		auto tars = idsMids_->getTargets();
		auto indxs = getIndexes();
		for(const auto & indx : indxs){
			indexToTars_[indx] = tars;
		}
	}else{
		table targetsToIndexTab(targetsToIndexFnp.string(), "\t", true);
		std::set<std::string> tarsSet;
		for (const auto & row : targetsToIndexTab.content_) {
			auto tars = tokenizeString(row[targetsToIndexTab.getColPos("targets")],
					",");
			auto index = row[targetsToIndexTab.getColPos("index")];
			if (njh::in(index, indexToTars_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error already have targets for index: "
						<< index << "\n";
				ss << "Check to see if there are repeat values in " << targetsToIndexFnp
						<< "\n";
				throw std::runtime_error { ss.str() };
			}
			njh::sort(tars);
			indexToTars_[index] = tars;
			tarsSet.insert(tars.begin(), tars.end());
		}
		auto indexes = getIndexes();
		VecStr missingIndex;
		for (const auto & index : indexes) {
			if (!njh::in(index, indexToTars_)) {
				missingIndex.emplace_back(index);
			}
		}
		if (!missingIndex.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ": error, missing the following index from the index to targets file "
					<< targetsToIndexFnp << "\n";
			ss << "Targets: " << njh::conToStr(missingIndex, ", ") << "\n";
			ss << "Indexes Read In: "
					<< njh::conToStr(njh::getVecOfMapKeys(indexToTars_), ", ") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}


}



TarAmpAnalysisSetup::Sample::Sample(const std::string & name) :
		name_(name) {
}

void TarAmpAnalysisSetup::Sample::addRep(const std::string & rep) {
	if (njh::in(rep, reps_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, sample " << name_
				<< " already has rep name " << rep << "\n";
		throw std::runtime_error { ss.str() };
	}
	reps_.emplace_back(rep);
}

void TarAmpAnalysisSetup::Sample::addReps(const VecStr & reps) {
	for (const auto & rep : reps) {
		addRep(rep);
	}
}

VecStr TarAmpAnalysisSetup::Sample::getReps() const {
	return reps_;
}

TarAmpAnalysisSetup::Samples::Samples(const std::string & target): target_(target){};


bool TarAmpAnalysisSetup::Samples::hasSample(const std::string & sample) {
	return samples_.end() != samples_.find(sample);
}

void TarAmpAnalysisSetup::Samples::addSample(const std::string & sample) {
	if (hasSample(sample)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, already have " << sample << " for " << target_ << "\n";
		ss << "With input: " << njh::conToStr(samples_.at(sample).reps_ , ", ")<< "\n";
		throw std::runtime_error { ss.str() };
	}
	samples_.emplace(sample, Sample{sample});
}

void TarAmpAnalysisSetup::Samples::addSample(const std::string & sample, const VecStr & reps) {
	addSample(sample);
	samples_.find(sample)->second.addReps(reps);
}

VecStr TarAmpAnalysisSetup::Samples::getSamples() const {
	return njh::getVecOfMapKeys(samples_);
}

std::vector<std::string> TarAmpAnalysisSetup::Samples::getReps() const {
	std::vector<std::string> ret;
	for (const auto & samp : samples_) {
		ret.insert(ret.end(), samp.second.reps_.begin(), samp.second.reps_.end());
	}
	return ret;
}

std::set<std::string> TarAmpAnalysisSetup::getSamples() const {
	std::set<std::string> ret;
	for (const auto & tar : samples_) {
		auto samps = tar.second.getSamples();
		ret.insert(samps.begin(), samps.end());
	}
	return ret;
}

std::vector<std::string> TarAmpAnalysisSetup::getReps() const {
	std::set<std::string> ret;
	for (const auto & tar : samples_) {
		auto reps = tar.second.getReps();
		ret.insert(reps.begin(), reps.end());
	}
	return VecStr{ret.begin(), ret.end()};
}

void TarAmpAnalysisSetup::addGroupingMetaData(const bfs::path & groupingsFileFnp) {
	groupMetaData_ = std::make_unique<MultipleGroupMetaData>(
			njh::files::normalize(groupingsFileFnp), getSamples());
}

void TarAmpAnalysisSetup::addGroupingsFile()const{
	if(nullptr != groupMetaData_){
		bfs::copy_file(groupMetaData_->groupingsFile_, njh::files::make_path(infoDir_, "groupingMetaData.tab.txt"));
	}
}

void TarAmpAnalysisSetup::writeSampleNamesFile() const {
	std::ofstream sampleNamesFile;
	openTextFile(sampleNamesFile,
			OutOptions(njh::files::make_path(infoDir_, "sampNames").string(),
					".tab.txt"));
	auto keys = getVectorOfMapKeys(samples_);
	njh::sort(keys);
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for (const auto & tarKey : keys) {
		if(!pars_.byIndex && 1 == idsMids_->targets_.at(tarKey).overlapStatuses_.size()
				&& PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP
						== idsMids_->targets_.at(tarKey).overlapStatuses_.front() &&
						pars_.techIsIllumina())  {

//			std::cout << tarKey << std::endl;
//			std::cout << "\t" << idsMids_->targets_.at(tarKey).overlapStatuses_.size() << std::endl;
//			std::cout << "\t" << static_cast<uint32_t>(idsMids_->targets_.at(tarKey).overlapStatuses_.front()) << std::endl;

			auto sampKeys = njh::getVecOfMapKeys(samples_.at(tarKey).samples_);
			njh::sort(sampKeys);

			for (const auto & sampKey : sampKeys) {
				sampleNamesFile << tarKey << "-R1" << "\t" << sampKey;
				for (auto rep : samples_.at(tarKey).samples_.at(sampKey).reps_) {
					if (!njh::beginsWith(rep, "MID")) {
						rep = "MID" + rep;
					}
					sampleNamesFile << "\t" << rep;
				}
				sampleNamesFile << std::endl;
			}
			for (const auto & sampKey : sampKeys) {
				sampleNamesFile << tarKey<< "-R2" << "\t" << sampKey;
				for (auto rep : samples_.at(tarKey).samples_.at(sampKey).reps_) {
					if (!njh::beginsWith(rep, "MID")) {
						rep = "MID" + rep;
					}
					sampleNamesFile << "\t" << rep;
				}
				sampleNamesFile << std::endl;
			}
		}else{
			auto sampKeys = njh::getVecOfMapKeys(samples_.at(tarKey).samples_);
			njh::sort(sampKeys);
			for (const auto & sampKey : sampKeys) {
				sampleNamesFile << tarKey << "\t" << sampKey;
				for (auto rep : samples_.at(tarKey).samples_.at(sampKey).reps_) {
					if (!njh::beginsWith(rep, "MID")) {
						rep = "MID" + rep;
					}
					sampleNamesFile << "\t" << rep;
				}
				sampleNamesFile << std::endl;
			}
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

VecStr TarAmpAnalysisSetup::getTargets() const {
	if(pars_.byIndex){
		std::set<std::string> tarsSet;
		for(const auto & indToTar : indexToTars_){
			tarsSet.insert(indToTar.second.begin(), indToTar.second.end());
		}
		return VecStr{tarsSet.begin(), tarsSet.end()};
	}else{
		return njh::getVecOfMapKeys(samples_);
	}
}

VecStr TarAmpAnalysisSetup::getIndexes() const {
	return njh::getVecOfMapKeys(samples_);
}

void TarAmpAnalysisSetup::addSamplesNames(const table & samplesNamesInputTab){
	VecStr warnings;
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
			ss << "Row Number(1 based pos): " << rowCount << " - " << njh::conToStr(row, "\t") << "\n";
			warnings.emplace_back(ss.str());
			continue;
		}
		std::vector<std::string> nonEmptyRow;
		for(const auto colPos : iter::range(row.size())){
			if(colPos < 2 && row[colPos].empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error column " << colPos << " can't be an empty string" << "\n";
				ss << "Row Number(1 based pos): " << rowCount << " - " << njh::conToStr(row, "\t") << "\n";
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
	std::stringstream errorStream;
	errorStream << __PRETTY_FUNCTION__ << ": error\n";
	bool failed = false;
	std::regex pat{"(^MID)([0-9]+$)"};
	std::smatch match;
	for (const auto & row : finalSamps) {
		std::string targetRaw = row[0];
		VecStr targets{targetRaw};
		if("all" == njh::strToLowerRet(targetRaw)){
			if(nullptr == idsMids_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "error, idsMids_ cannot be nullptr"<< "\n";
				throw std::runtime_error{ss.str()};
			}
			targets = idsMids_->getTargets();
		}
		std::string sample = row[1];
		VecStr reps = getSubVector(row, 2);
		for(const auto & target : targets){
			if(pars_.byIndex){
				for(const auto & rep : reps){
					if(!std::regex_match(rep, match, pat)){
						errorStream << "Error for " << sample << " in " << target << " rep name needs begins with MID and end with a number" << "\n";
						errorStream << "Failed: " << njh::bashCT::red << rep << njh::bashCT::reset << "\n";
						failed = true;
					}
				}
			}
			for (const auto & rep : reps) {
				for (const auto & tar : samples_) {
					if(tar.first != target){
						continue;
					}
					for (const auto & samp : tar.second.samples_) {
						if (njh::in(rep, samp.second.reps_)) {
							errorStream << __PRETTY_FUNCTION__ << " : Error for " << sample << " in " << target << ", already have " << rep
									<< " in " << samp.first << "\n";
						}
					}
				}
			}
			auto search = samples_.find(target);
			if(samples_.end() == search){
				samples_.emplace(target, Samples(target));
			}
			samples_.at(target).addSample(sample, reps);
			//std::cout << "target: " << target << ", sample: " << sample << ", reps: " << njh::conToStr(reps, ",") << std::endl;
		}
	}
	if(failed){
		throw std::runtime_error{errorStream.str()};
	}
}

void TarAmpAnalysisSetup::addSamplesNames(const bfs::path & samplesNamesFnp){
	table samplesNamesInputTab(samplesNamesFnp.string());
	addSamplesNames(samplesNamesInputTab);
}


void TarAmpAnalysisSetup::addRefSeqs(const bfs::path & refSeqsDir){
	/**@todo add ref seq functions and checks*/
	auto fastaFiles = njh::files::listAllFiles(refSeqsDir.string(), false, {
			std::regex { R"(.*\.fasta$)" } });
	VecStr found;
	auto tars = getTargets();
	for (const auto & ff : fastaFiles) {
		auto tarName = bfs::basename(ff.first);
		found.emplace_back(tarName);
		if (!njh::in(tarName, tars)) {
			forRefSeqs_.notMatching_.emplace_back(tarName);
		} else {
			auto inOpts = SeqIOOptions::genFastaIn(ff.first.string(), false);
			SeqInput reader(inOpts);
			auto seqs = reader.readAllReads<seqInfo>();
			if (seqs.size() == 1) {
				idsMids_->targets_.at(tarName).addSingleRef(seqs.front());
			} else {
				idsMids_->targets_.at(tarName).addMultileRef(seqs);
			}
		}
	}
	for (const auto & tar : tars) {
		if (!njh::in(tar, found)) {
			forRefSeqs_.missing_.emplace_back(tar);
		}
	}
}


void TarAmpAnalysisSetup::addLenCutOffs(const bfs::path & lenCutOffsFnp){
	auto multipleLenCutOffs = PrimersAndMids::readInLenCutOffs(lenCutOffsFnp);
	auto targets = getTargets();
	for(const auto & cutOff : multipleLenCutOffs){
		if(!njh::in(cutOff.first, targets)){
			forLenCutOffs_.notMatching_.emplace_back(cutOff.first);
		}
	}
	for (const auto & tar : idsMids_->getTargets()) {
		if (!njh::in(tar, multipleLenCutOffs)) {
			forLenCutOffs_.missing_.emplace_back(tar);
		} else {
			idsMids_->targets_.at(tar).addLenCutOff(
					multipleLenCutOffs.at(tar).minLenChecker_.minLen_,
					multipleLenCutOffs.at(tar).maxLenChecker_.maxLen_);
		}
	}
}

void TarAmpAnalysisSetup::addOverlapStatus(const bfs::path & overlapStatusFnp){
	idsMids_->addOverLapStatuses(overlapStatusFnp);
}




std::vector<VecStr> TarAmpAnalysisSetup::getTarCombos() const{
	std::unordered_map<std::string, VecStr> targetsForReps;
	if (pars_.byIndex) {
		targetsForReps = indexToTars_;
	} else {
		for (const auto & tars : samples_) {
			for (const auto & rep : tars.second.getReps()) {
				targetsForReps[rep].emplace_back(tars.first);
			}
		}
	}

	std::vector<VecStr> tarCombos;
	for (auto & sampTars : targetsForReps) {
		njh::sort(sampTars.second);
		if (!njh::in(sampTars.second, tarCombos)) {
			tarCombos.emplace_back(sampTars.second);
		}
	}
	return tarCombos;
}

void TarAmpAnalysisSetup::writeOutIdFiles() {
	//now write id files
	std::vector<VecStr> tarCombos = getTarCombos();
	tarsToTargetSubSets_.clear();

	for (const auto & tarCombo : tarCombos) {

		auto collapseStr = njh::conToStr(tarCombo, "_");
		auto collapseIdx = estd::to_string(tarsToTargetSubSets_.size());
		tarsToTargetSubSets_[collapseStr] = collapseIdx;
		auto refs = idsMids_->getRefSeqs(tarCombo);
		auto lens = idsMids_->genLenCutOffs(tarCombo);
		auto overlapStatuses = idsMids_->genOverlapStatuses(tarCombo);

		if (!lens.empty()) {
			auto lensOutOpts = TableIOOpts::genTabFileOut(
					njh::files::make_path(idsDir_,
							collapseIdx + "_lenCutOffs.tab.txt"));
			lens.outPutContents(lensOutOpts);
		}

		if (!overlapStatuses.empty()) {
			auto overlapStatusesOpts = TableIOOpts::genTabFileOut(
					njh::files::make_path(idsDir_,
							collapseIdx + "_overlapStatus.tab.txt"));
			overlapStatuses.outPutContents(overlapStatusesOpts);
		}

		idsMids_->writeIdFile(
				OutOptions(
						njh::files::make_path(idsDir_, collapseIdx + ".id.txt")),
				tarCombo);
	}
	for(const auto & tar : idsMids_->targets_){
		if(!tar.second.refs_.empty()){
			SeqOutput::write(tar.second.refs_,
					SeqIOOptions::genFastaOut(
							njh::files::make_path(refsDir_, tar.first)));
		}
	}
}

VecStr TarAmpAnalysisSetup::getExpectantInputNames()const{
	VecStr ret;
	if(pars_.byIndex){
		ret = njh::getVecOfMapKeys(samples_);
	}else{
		ret = getReps();
	}
	return ret;
}


void TarAmpAnalysisSetup::setUpPopClusteringDirs(bool verbose) const {
	auto tars = getTargets();
	if (pars_.byIndex) {
		auto topPopDir = njh::files::makeDir(dir_.string(), njh::files::MkdirPar("popClustering"));
		for (const auto & tar : tars) {
			if (1 == idsMids_->targets_.at(tar).overlapStatuses_.size() &&
					PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP == idsMids_->targets_.at(tar).overlapStatuses_.front() &&
					pars_.techIsIllumina() ) {
				setUpSampleDirs(
						njh::files::make_path(infoDir_, "sampNames.tab.txt").string(),
						njh::files::make_path(topPopDir, tar + "-R1").string(), false, verbose);
				setUpSampleDirs(
						njh::files::make_path(infoDir_, "sampNames.tab.txt").string(),
						njh::files::make_path(topPopDir, tar + "-R2").string(), false, verbose);
			} else {
				setUpSampleDirs(
						njh::files::make_path(infoDir_, "sampNames.tab.txt").string(),
						njh::files::make_path(topPopDir, tar).string(), false, verbose);
			}
		}
	} else {
		setUpSampleDirs(
							njh::files::make_path(infoDir_, "sampNames.tab.txt").string(),
							njh::files::make_path(dir_, "popClustering").string(), true, verbose);
	}
}



table GuessPossibleSamps(const TarAmpAnalysisSetup::TarAmpPars & pars){
	table sampleNames(VecStr{"#target", "sample", "run1"});
	std::regex inputfilePatReg{ pars.inputFilePat };
	auto files = njh::files::listAllFiles(pars.inputDir.string(), false, {inputfilePatReg});

	if("" != pars.replicatePattern){
		std::regex replicatePatternReg{ pars.replicatePattern };
		std::cout << "replicatePatternReg.mark_count(): " << replicatePatternReg.mark_count() << std::endl;
		std::map<std::string, VecStr> replicates;
		if (pars.techIsIllumina()) {
			ReadPairsOrganizer rpOrganizer(VecStr{});
			rpOrganizer.processFiles(files);
			for(const auto & samp : rpOrganizer.readPairsUnrecognized_){
				if(njh::in(samp.first, pars.ignoreSamples)){
					continue;
				}
				std::smatch match;
				if(std::regex_match(samp.first, match, replicatePatternReg)){
					std::cout << match[1] << ": " << match[2] << std::endl;
					replicates[match[1]].emplace_back(samp.first);
				}else{
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "failed to find replicate pattern: " << pars.replicatePattern << " in " << samp.first<< "\n";
					throw std::runtime_error{ss.str()};
				}
			}
		} else {
			VecStr samplesFiles;
			for(const auto & f : files){
				samplesFiles.emplace_back(f.first.filename().string());
			}
			for(const auto & sf : samplesFiles){
				auto sampName = sf.substr(0, sf.rfind(".fast"));
				if(njh::in(sampName, pars.ignoreSamples)){
					continue;
				}
				std::smatch match;
				if(std::regex_match(sampName, match, replicatePatternReg)){
					std::cout << match[1] << ": " << match[2] << std::endl;
					replicates[match[1]].emplace_back(sampName);
				}else{
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "failed to find replicate pattern: " << pars.replicatePattern << " in " << sampName<< "\n";
					throw std::runtime_error{ss.str()};
				}
			}
		}
		uint32_t maxReplicate = 1;
		for(auto & rep : replicates){
			njh::sort(rep.second);
			maxReplicate = std::max<uint32_t>(maxReplicate, rep.second.size());
		}
		VecStr colNames{"#target", "sample"};
		for(uint32_t run = 1; run <= maxReplicate; ++run){
			colNames.emplace_back(njh::pasteAsStr("run", run));
		}
		sampleNames = table(colNames);
		for(const auto & rep : replicates){
			VecStr rowToAdd{"all", rep.first};
			for(const auto & repName : rep.second){
				rowToAdd.emplace_back(repName);
			}
			for(uint32_t rest = 0; rest < maxReplicate - rep.second.size(); ++rest){
				rowToAdd.emplace_back("");
			}
			sampleNames.addRow(rowToAdd);
		}

	}else{
		if (pars.techIsIllumina()) {
			ReadPairsOrganizer rpOrganizer(VecStr{});
			rpOrganizer.processFiles(files);
			for(const auto & samp : rpOrganizer.readPairsUnrecognized_){
				if(njh::in(samp.first, pars.ignoreSamples)){
					continue;
				}
				sampleNames.addRow("all", samp.first, samp.first);
			}
		} else {
			VecStr samplesFiles;
			for(const auto & f : files){
				samplesFiles.emplace_back(f.first.filename().string());
			}
			for(const auto & sf : samplesFiles){
				auto sampName = sf.substr(0, sf.rfind(".fast"));
				if(njh::in(sampName, pars.ignoreSamples)){
					continue;
				}
				sampleNames.addRow("all", sampName, sampName);
			}
		}
	}

	//sampleNames.outPutContents(std::cout, "\t");
	return sampleNames;
}




}  // namespace njhseq

