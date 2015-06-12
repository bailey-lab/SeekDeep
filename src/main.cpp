//
//  main.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/11/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//


#include "SeekDeep.h"

int main(int argc, char* argv[]) {
	try {
		bibseq::SeekDeepRunner seqRunner;
		return seqRunner.run(argc, argv);
	} catch (std::exception & e) {
		std::cout << e.what() << std::endl;
	}
	return 0;
}
