#pragma once
/*
 * ReadPairsOrganizer.hpp
 *
 *  Created on: Jan 23, 2017
 *      Author: nick
 */



#include <bibseq.h>


namespace bibseq {



class ReadPairsOrganizer {
public:
	ReadPairsOrganizer(const VecStr & expectedSamples);

	VecStr expectedSamples_;
	std::unordered_map<std::string, VecStr> readPairs_;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized_;

	void processFiles(const std::map<bfs::path, bool> & files);
	std::unordered_map<std::string, std::pair<VecStr, VecStr>> processReadPairs() ;
};


}  // namespace bibseq
