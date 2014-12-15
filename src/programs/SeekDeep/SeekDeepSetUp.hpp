#pragma once
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
#include <bibseq.h>

#include <bibseq/programUtils/seqSetUp.hpp>

namespace bibseq {

class SeekDeepSetUp : public seqSetUp {

 public:
  // constructors
  SeekDeepSetUp(int argc, char* argv[]) : seqSetUp(argc, argv) {}
  SeekDeepSetUp(const bib::progutils::commandLineArguments& inputCommands)
      : seqSetUp(inputCommands) {}
  SeekDeepSetUp(const MapStrStr& inputCommands)
      : seqSetUp(inputCommands) {}

  void setUpExtractor(std::string& idFileName, bool& multiplex, bool& condensed,
                      int& minLen, int& maxLength, int& within,
                      int& qualityWindowLength, int& qualityWindowStep,
                      int& qualityWindowThres, bool& unknown,
                      int& unknownPrimerSize, int& unknownMinSize,
                      bool& findReversePrimer, double& queryCoverageCutoff,
                      double& percentIdentityCutoff, int& numberOfNs,
                      bool& pyroExtract, bool& reversePrimerToLowerCase,
                      bool& flowFiltering, int& maxFlows, bool& checkComplement,
                      bool& screenForPossibleContaimination,
                      std::string& compareSeq, std::string& idFileDelim,
                      int& smallFragmentCutOff, bool& qualWindowTrim);
  void setUpClusterDown(std::string& qualRep, std::string& parameters,
                        bool& extra,
                        std::map<int, std::vector<double>>& iteratorMap,
                        bool& smallestFirst, bool& markChimeras, int& parFreqs,
                        bool& bestMatch, int& bestMatchCheck, bool& snapShots,
                        std::string& sortBy, bool& additionalOut,
                        std::string& additonalOutLocationFile,
                        bool& collapsingTandems,
                        bool& kmerCheckingOnAlignProfile,
                        bool& condesnedCollpase, bool& removeLowQualBases,
                        int& lowQualityCutOff, bool& adjustHomopolyerRuns);
  void setUpMultipleSampleCluster(
      std::string& parameters, bool& extra, int& cutOff,
      std::map<int, std::vector<double>>& iteratorMap, bool& population,
      double& fracCutoff, bool& smallestFirst, bool& bestMatch,
      int& bestMatchCheck, bool& checkChimeras, int& parFreqs,
      std::string& parametersPopulation, bool& differentPar,
      std::map<int, std::vector<double>>& popIteratorMap, bool& popBoth,
      bool& keepChimeras, std::string& experimentName, uint32_t& runsRequired);
  void setUpCollapseTandems(double& freqCutoff, bool& extra,
                            bool& additionalOut,
                            std::string& additionalOutLocationFile);
  void setUpMarkChimeras(double& parentFreqs);
  void setUpMakeSampleDirectories(std::string& sampleNameFilename);
  void setUpPopulationClustering(
      std::string& directory, VecStr& contains, bool& specific, bool& recursive,
      std::string& parameters, bool& extra, int& cutOff,
      std::map<int, std::vector<double>>& popIteratorMap, double& fracCutoff,
      bool& smallestFirst, bool& bestMatch, int& bestMatchCheck,
      bool& checkChimeras, int& parFreqs);
  void setUpCompareTwoReplicates(
      std::string& filename1, std::string& qualFilename1, std::string& format1,
      std::string& filename2, std::string& qualFilename2, std::string& format2,
      std::string& parameters, bool& extra, int& cutOff,
      std::map<int, std::vector<double>>& iteratorMap, bool& smallestFirst,
      bool& markChimeras, int& parFreq, bool& bestMatch, int& bestMatchCheck);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "SeekDeepSetUp.cpp"
#endif
