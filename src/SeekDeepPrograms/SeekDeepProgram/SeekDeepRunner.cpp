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
//
//  main.cpp
//  SeekDeep
//
//  Created by Nicholas Hathaway on 8/11/13.
//
#include "SeekDeepRunner.hpp"
#include <njhcpp.h>
#include <SeekDeepPrograms/SeekDeepServerRunner.h>
#include "SeekDeepPrograms/SeekDeepUtils.h"

namespace njhseq {

SeekDeepRunner::SeekDeepRunner() :
		njh::progutils::OneRing(
				{ addRing<SeekDeepUtilsRunner>(),
					addRing<SeekDeepServerRunner>()
				},
				{
						addFunc("extractor", extractor, false),
						addFunc("extractorPairedEnd", extractorPairedEnd, false),
						addFunc("processClusters", processClusters,false),
						addFunc("qluster", clusterDown, false),
						addFunc("clusterDown",clusterDown, true),
						addFunc("makeSampleDirectories", makeSampleDirectories, false)
				}, "SeekDeep", "2", "6", "3") {
}

//

}  // namespace njhseq
