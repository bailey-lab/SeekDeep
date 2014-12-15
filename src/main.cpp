//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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


#include "programs/SeekDeep.h"

int main(int argc, char* argv[]) {
	//std::cout << "Start " << std::endl;
	try{
  bibseq::SeekDeepRunner seqRunner;
  if (argc > 1) {
    return seqRunner.run(argc, argv);
  }

  seqRunner.listPrograms(std::cout);
	}catch (std::exception & e) {
		std::cout << e.what() << std::endl;
	}
  return 0;
}
