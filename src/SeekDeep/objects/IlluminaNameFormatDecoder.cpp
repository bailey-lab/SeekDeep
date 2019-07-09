/*
 * IlluminaNameFormatDecoder.cpp
 *
 *  Created on: Jul 8, 2019
 *      Author: nicholashathaway
 */




#include "IlluminaNameFormatDecoder.hpp"

namespace njhseq {
std::regex IlluminaNameFormatDecoder::NameRegPatStr_{"([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([AGCTN+]+))"};

}  // namespace njhseq
