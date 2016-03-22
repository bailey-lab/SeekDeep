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
function addPanelWithDiv(addToSelector, divId, panelTitle, panelType){
	var top = d3.select(addToSelector)
				.append("div")
					.attr("class", "panel " + (panelType || "panel-primary"));
	top.append("div")
		.attr("class", "panel-heading")
			.append("h1")
				.text(panelTitle);
	top.append("div")
		.attr("class", "panel-body")
			.append("div")
				.attr("id", divId);
	return top;
}

function addPanelOnlyHead(addToSelector, panelTitle, panelType){
	var top = d3.select(addToSelector)
				.append("div")
					.attr("class", "panel " + (panelType || "panel-primary"));
	top.append("div")
		.attr("class", "panel-heading")
			.append("h1")
				.text(panelTitle);
	return top;
}



function addMainDiv(addToSelector, addJumbo){
	var main = d3.select(addToSelector).append("div")
		.attr("class", "container theme-showcase")
		.attr("id", "mainContent")
		.attr("role", "main");
	if(addJumbo){
		main.append("div")
			.attr("class", "jumbotron")
			.append("h1")
				.attr("id", "jumboTitle");
	}
	return main;
}

