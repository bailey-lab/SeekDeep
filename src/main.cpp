//
//  main.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/11/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
#include "SeekDeep.h"




int main(int argc, char* argv[]) {
	try {
		bibseq::SeekDeepRunner seqRunner;
		return seqRunner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
