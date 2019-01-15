#pragma once
/*
 * ReadPairsOrganizer.hpp
 *
 *  Created on: Jan 23, 2017
 *      Author: nick
 */
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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


namespace njhseq {



class ReadPairsOrganizer {
public:
	ReadPairsOrganizer(const VecStr & expectedSamples);

	VecStr expectedSamples_;
	bool doNotGuessSampleNames_{false};
	std::regex illuminaPat_{"(.*?)((_S[0-9]+)?_(R[12])(_[0-9]+)?\\.fastq(\\.gz)?)"};

	std::unordered_map<std::string, VecStr> readPairs_;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized_;

	void processFiles(const std::map<bfs::path, bool> & files);
	std::unordered_map<std::string, std::pair<VecStr, VecStr>> processReadPairs();
};


}  // namespace njhseq
