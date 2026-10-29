/*
 * bamSampleNames.cpp
 *
 *  Created on: Oct 28, 2026
 *      Author: nick
 */

#include "bamSampleNames.hpp"

#include <api/BamReader.h>

namespace njhseq {

std::set<std::string> getSampleNamesFromBamHeader(const bfs::path & bamFnp){
  std::set<std::string> samples;
  BamTools::BamReader headerReader;
  if(!headerReader.Open(bamFnp.string())){
    return samples;
  }
  for(const auto & rg : headerReader.GetHeader().ReadGroups){
    if(rg.HasSample() && !rg.Sample.empty()){
      samples.emplace(rg.Sample);
    }
  }
  headerReader.Close();
  return samples;
}

bool setSampleNameFromBamHeader(const bfs::path & bamFnp, std::string & sampleName, const std::string & sampleFlag){
  auto samples = getSampleNamesFromBamHeader(bamFnp);
  if(1 == samples.size()){
    sampleName = *samples.begin();
    return true;
  }
  if(samples.size() > 1){
    std::cerr << "Warning, found multiple SM values in the header of " << bamFnp
              << " (" << njh::conToStr(samples, ",") << "), not defaulting sample name to SM, use " << sampleFlag << " to set it" << std::endl;
  }
  return false;
}

}  // namespace njhseq
