#pragma once

/*
 * ControlMixSetUp.hpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */

#include <njhseq/common.h>

#include <njhseq/objects/Meta/MetaDataInName.hpp>
#include <njhseq/objects/dataContainers/tables/table.hpp>

namespace njhseq {

class ControlMixSetUp {
public:
	ControlMixSetUp(const std::string & name,
			const std::unordered_map<std::string, double> & relativeAbundances);


	std::string name_;
	std::unordered_map<std::string, double> rawRelativeAbundances_;
	std::unordered_map<std::string, double> relativeAbundances_;

	MetaDataInName meta_;

	VecStr getStrains() const;

	void removeStrain(const std::string & name);


	static std::unordered_map<std::string, ControlMixSetUp> readInSetUps(const bfs::path & mixtureSetUpFnp);

	static std::unordered_map<std::string, ControlMixSetUp> processSetUps(const table & mixtureSetUp);


};

}  // namespace njhseq




