/*
 * setUpPars.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */




#include "setUpPars.hpp"

namespace bibseq {


extractorPars::extractorPars(){
  rPrimerErrors.hqMismatches_ = 4;
  rPrimerErrors.distances_.query_.coverage_ = .50;
  rPrimerErrors.largeBaseIndel_ = .99;
  rPrimerErrors.oneBaseIndel_ = 2;
  rPrimerErrors.twoBaseIndel_ = 1;

  fPrimerErrors.hqMismatches_ = 2;
  fPrimerErrors.distances_.query_.coverage_ = 1;
  fPrimerErrors.largeBaseIndel_ = .99;
  fPrimerErrors.oneBaseIndel_ = 2;
  fPrimerErrors.twoBaseIndel_ = 1;
}


ExtractorPairedEndPars::ExtractorPairedEndPars(){


  primerErrors.hqMismatches_ = 2;
  primerErrors.lqMismatches_ = 5; /**@todo incorporate this*/
  primerErrors.distances_.query_.coverage_ = 1;
  primerErrors.largeBaseIndel_ = 0.99;
  primerErrors.oneBaseIndel_ = 0.5;
  primerErrors.twoBaseIndel_ = 0.5;
}

}  // namespace bibseq
