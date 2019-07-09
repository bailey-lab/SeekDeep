#pragma once

/*
 * IlluminaNameFormatDecoder.hpp
 *
 *  Created on: Jul 8, 2019
 *      Author: nicholashathaway
 */


#include <njhseq.h>
namespace njhseq {

class IlluminaNameFormatDecoder{
public:

	IlluminaNameFormatDecoder(const std::string & name){
		std::regex_match(name, match_, NameRegPatStr_);
	}
	static std::regex NameRegPatStr_;
	std::smatch match_;

	/** index in match to what value in the illumina
	 * 0 full name
	 * 1 <instrument>
	 * 2 <run number>
	 * 3 <flowcell ID>
	 * 4 <lane>
	 * 5 <tile>
	 * 6 <x_pos>
	 * 7 <y_pos>
	 * space
	 * 8 <read> (1 or 2 for first mate or second match)
	 * 9 <is filtered> Y or N for yes or no
	 *10 <control number>
	 *11 <sample number> barcode normally eg AGGCGT
	 */

	std::string getIndexValue(uint32_t idx) const{
		if(idx <match_.size()){
			return match_[idx];
		}
		return "";
	}

	std::string getSampleNumber() const{
		return getIndexValue(11);
	}
};



} //namespace njhseq

