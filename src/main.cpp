//
//  main.cpp
//
//  Created by Nicholas Hathaway on 8/11/13.
//

//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2019 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include <utility>

#include "SeekDeepPrograms.h"

//namespace njhseq {
//
//
//int test(){
//	{
//		VecStr targets{"t1", "t2", "t100", "t10", "t22", "t99", "t9"};
//		std::cout << njh::conToStr(targets) << std::endl;
//		njh::naturalSortNameSet(targets);
//		std::cout << njh::conToStr(targets) << std::endl;
//	}
//
//	return 0;
//}
//
//}  // namespace njhseq
//





int main(int argc, char* argv[]) {
//	return njhseq::test();
//	{
//		std::string name = "M01380:70:000000000-B9CJY:1:2118:15005:24494 1:N:0:hu136 AGGAGTCC|0|TAGATCGC|0";
//		//std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-z0-9_-+]+):([A-z0-9_-+]+) ([A-z0-9_|+-]+)( .*)?";
//		uint32_t BackUpIlluminaSampleNumberPos_ = 12;
//	//	std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-Za-z0-9_-+]+):([A-Za-z0-9_-+]+) ([A-Za-z0-9_|+-]+)";
//		//std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-Za-z0-9_-+:]+) ([A-z0-9_|+-]+)( .*)?";
//		std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-Za-z0-9_+:-]+) ([A-z0-9_|+-]+)( .*)?";
//
//		std::cout << "name: " << name << std::endl;
//
//		std::smatch mat;
//		std::cout << njh::colorBool(std::regex_match(name,mat,std::regex{BackUpIlluminaSampleRegPatStr_} )) << std::endl;
//		std::cout << "mat.size(): " << mat.size() << std::endl;
//		if(std::regex_match(name,mat,std::regex{BackUpIlluminaSampleRegPatStr_} )){
//			std::cout << mat[BackUpIlluminaSampleNumberPos_] << std::endl;
//			for( uint32_t pos = 0; pos < mat.size(); ++pos){
//				std::cout << "pos: " << pos << " " << mat[pos] << std::endl;
//			}
//		}
//	}
//	{
//		std::string name = "M01380:109:000000000-BJD3T:1:1102:23389:3505 1:N:0:hu99b:5_target_multiplex TTCTGGGT|1|TGTTCTCT|0 Pf_ama1d1_fwd|0|24|";
//		//std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-z0-9_-+]+):([A-z0-9_-+]+) ([A-z0-9_|+-]+)( .*)?";
//	//	std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-Za-z0-9_-+]+):([A-Za-z0-9_-+]+) ([A-Za-z0-9_|+-]+)";
//		std::string BackUpIlluminaSampleRegPatStr_ = "([A-Za-z0-9_]+):([0-9]+):([A-Za-z0-9-]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+) ([12]):([NY]):([0-9]):([A-Za-z0-9_+:-]+) ([A-z0-9_|+-]+)( .*)?";
//		uint32_t BackUpIlluminaSampleNumberPos_ = 12;
//
//		std::cout << "name: " << name << std::endl;
//		std::smatch mat;
//		std::cout << njh::colorBool(std::regex_match(name,mat,std::regex{BackUpIlluminaSampleRegPatStr_} )) << std::endl;
//		std::cout << "mat.size(): " << mat.size() << std::endl;
//		if(std::regex_match(name,mat,std::regex{BackUpIlluminaSampleRegPatStr_} )){
//			std::cout << mat[BackUpIlluminaSampleNumberPos_] << std::endl;
//			for(uint32_t pos = 0; pos < mat.size(); ++pos){
//				std::cout << "pos: " << pos << " " << mat[pos] << std::endl;
//			}
//		}
//	}
//
//
//	return 0;
	try {
		njhseq::SeekDeepRunner seqRunner;
		return seqRunner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}


