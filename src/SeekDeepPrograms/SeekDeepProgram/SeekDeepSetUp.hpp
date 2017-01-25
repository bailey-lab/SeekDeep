#pragma once
//
//  SeekDeepSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/13.
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
#include <bibseq.h>
#include <bibseq/programUtils/seqSetUp.hpp>

#include "SeekDeep.h"

namespace bibseq {



class SeekDeepSetUp : public seqSetUp {

 public:
	using seqSetUp::seqSetUp;

	void setUpExtractor(extractorPars & pars);
	void setUpClusterDown(clusterDownPars & pars);
	void setUpMultipleSampleCluster(processClustersPars & pars);
	void setUpMakeSampleDirectories(makeSampleDirectoriesPars & pars);

};
}  // namespace bibseq


