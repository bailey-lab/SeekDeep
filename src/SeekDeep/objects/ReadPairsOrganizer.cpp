/*
 * ReadPairsOrganizer.cpp
 *
 *  Created on: Jan 23, 2017
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

#include "ReadPairsOrganizer.hpp"


namespace njhseq {

std::regex ReadPairsOrganizer::illuminaPat_{"(.*?)((_S[0-9]+)?(_L[0-9]+)?_(R[12])(_[0-9]+)?\\.fastq(\\.gz)?)"};


ReadPairsOrganizer::ReadPairsOrganizer(const VecStr & expectedSamples) :
		expectedSamples_(expectedSamples) {

}


void ReadPairsOrganizer::processFiles(const std::map<bfs::path, bool> & files) {

	for (const auto & f : files) {
		auto filename = f.first.filename().string();
		auto underPos = filename.find("_");
		if (0 == underPos) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
					<< f.first
					<< ", shouldn't start with an _, can't determine sample name\n";
			throw std::runtime_error { ss.str() };
		}
		if (std::string::npos == underPos) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
					<< f.first
					<< ", should contain an _ to separate sample name from mate pair designations\n";
			ss << "Examples: Samp1_R1.fastq Samp1_R2.fastq Samp2_R1.fastq ..."
					<< "\n";
			throw std::runtime_error { ss.str() };
		}



		if (doNotGuessSampleNames_) {
			VecStr matchingSamples;
			for(const auto & sampName : expectedSamples_){
				if(njh::beginsWith(filename, sampName)){
					matchingSamples.emplace_back(sampName);
				}
			}
			if(matchingSamples.empty()){
				std::string sampName = filename.substr(0, underPos);
				std::smatch matchRes;
				if(std::regex_match(filename, matchRes, illuminaPat_)){
					sampName = matchRes[1];
				}
				readPairsUnrecognized_[sampName].emplace_back(f.first.string());
			} else if(matchingSamples.size() == 1){
				readPairs_[matchingSamples.front()].emplace_back(f.first.string());
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << filename << " matched multiple input names: " << njh::conToStrEndSpecial(matchingSamples, ", ", " and ") << "\n";
				throw std::runtime_error{ss.str()};
			}
		} else {
			std::smatch matchRes;
			std::string sampName = filename.substr(0, underPos);
			if(std::regex_match(filename, matchRes, illuminaPat_)){
				sampName = matchRes[1];
			}
			if (njh::in(sampName, expectedSamples_)
					|| njh::in("MID" + sampName, expectedSamples_)) {
				readPairs_[sampName].emplace_back(f.first.string());
			} else {
				readPairsUnrecognized_[sampName].emplace_back(f.first.string());
			}
		}
	}
}

std::unordered_map<std::string, std::pair<VecStr, VecStr>> ReadPairsOrganizer::processReadPairs() {
	std::unordered_map<std::string, std::pair<VecStr, VecStr>> readsByPairs;
	for (const auto & reads : readPairs_) {
		for (const auto & read : reads.second) {
			auto filename = bfs::path(read).filename().string();
			auto lastUnderPos = filename.rfind("_");
			auto periodPos = filename.find(".", lastUnderPos);
			if(std::string::npos == lastUnderPos){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
						<< read <<  ", should contain an _ before the read mate designation\n";
				throw std::runtime_error{ss.str()};
			}
			if(std::string::npos == periodPos){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
						<< read <<  ", should contain an . after the read designation, normally for the file extension\n";
				throw std::runtime_error{ss.str()};
			}
			std::string pairNum = filename.substr(lastUnderPos + 1, periodPos - lastUnderPos - 1);
			while("R1" != pairNum
					&& "1" != pairNum
					&& "R2" != pairNum
					&& "2" != pairNum){
				//find next _ to try to determine read pairs designation
				periodPos = lastUnderPos;
				lastUnderPos = filename.rfind("_", lastUnderPos - 1);
				if(std::string::npos == lastUnderPos){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ": Error, in processing file name "
							<< read <<  ", couldn't find the mate designation \n";
					ss << "Examples: Samp1_R1.fastq Samp1_R2.fastq Samp2_R1.fastq ..." << "\n";
					ss << "Examples: Samp1_R1_001.fastq Samp1_R2_001.fastq Samp2_R1_001.fastq ..." << "\n";
					throw std::runtime_error{ss.str()};
				}
				pairNum = filename.substr(lastUnderPos + 1, periodPos - lastUnderPos - 1);
			}
			if("R1" == pairNum || "1" == pairNum){
				readsByPairs[reads.first].first.emplace_back(read);
			}else if("R2" == pairNum || "2" == pairNum){
				readsByPairs[reads.first].second.emplace_back(read);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, cannot determine pair number for "
						<< filename << ", found designation of " << pairNum << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
	}
	//file checks
	for(auto & reads : readsByPairs){
		//check size
		if(reads.second.first.size() != reads.second.second.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, read paired file numbers for " << reads.first << " don't match up \n";
			ss << "R1 file number: " << reads.second.first.size() << ", files: " << njh::conToStr(reads.second.first, ", ") << "\n";
			ss << "R2 file number: " << reads.second.second.size() << ", files: " << njh::conToStr(reads.second.second, ", ") << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::vector<size_t> needToErrase;
		//sort
		njh::sort(reads.second.first);
		njh::sort(reads.second.second);
		for(const auto pos : iter::range(reads.second.first.size())){

			auto lastUnderPos1 = reads.second.first[pos].rfind("_");
			{
				auto periodPos = reads.second.first[pos].find(".", lastUnderPos1);
				std::string pairNum = reads.second.first[pos].substr(lastUnderPos1 + 1,
						periodPos - lastUnderPos1 - 1);
				while ("R1" != pairNum && "1" != pairNum) {
					//find next _ to try to determine read pairs designation
					periodPos = lastUnderPos1;
					lastUnderPos1 = reads.second.first[pos].rfind("_", lastUnderPos1 - 1);
					pairNum = reads.second.first[pos].substr(lastUnderPos1 + 1,
							periodPos - lastUnderPos1 - 1);
				}
			}
			auto lastUnderPos2 = reads.second.second[pos].rfind("_");
			{
				auto periodPos = reads.second.second[pos].find(".", lastUnderPos2);
				std::string pairNum = reads.second.second[pos].substr(lastUnderPos2 + 1,
						periodPos - lastUnderPos2 - 1);
				while ("R2" != pairNum && "2" != pairNum) {
					//find next _ to try to determine read pairs designation
					periodPos = lastUnderPos2;
					lastUnderPos2 = reads.second.second[pos].rfind("_",
							lastUnderPos2 - 1);
					pairNum = reads.second.second[pos].substr(lastUnderPos2 + 1,
							periodPos - lastUnderPos2 - 1);
				}
			}
			auto nameStub1 = reads.second.first[pos].substr(0, lastUnderPos1);
			auto nameStub2 = reads.second.second[pos].substr(0, lastUnderPos2);
			//check name
			if(nameStub1 != nameStub2){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, read paired file names don't match up for " << reads.first << "\n";
				ss << "for " << reads.second.first[pos] << " and " << reads.second.second[pos] << "\n";
				ss << "stub1: " << nameStub1 << "\n";
				ss << "stub2: " << nameStub2 << "\n";
				throw std::runtime_error{ss.str()};
			}
			//check to see if they are empty
			if(0 == bfs::file_size(reads.second.first[pos]) ||
					0 == bfs::file_size(reads.second.second[pos])){
				needToErrase.emplace_back(pos);
			}
		}
		if(!needToErrase.empty()){
			//sort the position so the last positions come first so positions don't get invalidated after erasing (back -> front errasing)
			std::sort(needToErrase.rbegin(), needToErrase.rend());
			for(const auto pos : needToErrase){
				reads.second.first.erase(reads.second.first.begin() + pos);
				reads.second.second.erase(reads.second.second.begin() + pos);
			}
		}
	}
	return readsByPairs;
}


}  // namespace njhseq
