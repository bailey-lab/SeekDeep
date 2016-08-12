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
function createProjectNavBarSkeleton(wrappingNavSelector, rName, projectName){
	addFixedTopNavSkeleton(wrappingNavSelector, "SeekDeep", "mainNav", "projectNav");
 	addNavLink("#mainNav", "Home", "/" + rName, "#siteHomeLink");
	addNavDrop("#mainNav", "Projects", "projectDrop");
	addNavLink("#mainNav", projectName, "/" + rName + "/mainProjectPage/"+ projectName, "projectHomeLink", true);
	addNavLink("#projectNav", "Extraction Info", "/" + rName + "/showExtractionInfo/" + projectName, "extractionLink");
	addNavDrop("#projectNav", "Groups", "groupsDrop");
	addNavDrop("#projectNav", "Samples", "samplesDrop");
}

function populateProjectNavBar(wrappingNavSelector, rName, projectName, groupName, sampName){
   	//get project names
	var projectsNames;
	ajax("/" + rName + "/projectsNames", function(mn){ projectsNames = mn; });
	var projLinkPre = "/" + rName + "/mainProjectPage/";
	
	//get group names
	var groupNames;
	ajax("/" + rName + "/getGroupNames/" + projectName, function(gn){ groupNames = gn; });
    var linkPreGroup = "/" + rName + "/showSubGroupsPage/" + projectName + "/";
	//get sample names
	var sampNames;
	ajax("/" + rName + "/sampleNames/" + projectName, function(mn){ sampNames = mn; });

    var linkPreSample = "/" + rName + "/individualSamplePage/" + projectName + "/";
	
    d3.select("#projectDrop")
		.selectAll("li")
		.data(projectsNames)
		.enter()
			.append("li")
				.attr("class", function(d){
					if(d == projectName){
						return "active";
					}else{
						return "";
					}
				})
			.append("a")
				.attr("href", function(d){ return projLinkPre + d;})
				.text(function(d){return d;});
    
	d3.select("#groupsDrop")
		.selectAll("li")
		.data(groupNames)
		.enter()
			.append("li")
				.attr("class", function(d){
					if(d == groupName){
						return "active";
					}else{
						return "";
					}
				})
			.append("a")
				.attr("href",function(d){return linkPreGroup + d;} )
				.text(function(d){return d;});
    
	d3.select("#samplesDrop")
		.selectAll("li")
		.data(sampNames)
		.enter()
			.append("li")
				.attr("class", function(d){
					if(d == sampName){
						return "active";
					}else{
						return "";
					}
				})
			.append("a")
				.attr("href",function(d){return linkPreSample + d;} )
				.text(function(d, i){return sampNames[i];});
}


function createProjectNavBar(wrappingNavSelector, rName, projectName, groupName, sampName){
	createProjectNavBarSkeleton(wrappingNavSelector, rName, projectName);
	populateProjectNavBar(wrappingNavSelector, rName, projectName, groupName, sampName);
}


function createMainPageNavBarSkeleton(wrappingNavSelector, rName){
	addFixedTopNavSkeleton(wrappingNavSelector, "SeekDeep", "mainNav", "projectNav");
 	addNavLink("#mainNav", "Home", "/" + rName, "#siteHomeLink");
	addNavDrop("#mainNav", "Projects", "projectDrop");
}

function populateMainPageNavBar(wrappingNavSelector, rName){
   	//get project names
	var projectsNames;
	ajax("/" + rName + "/projectsNames", function(mn){ projectsNames = mn; });
	var projLinkPre = "/" + rName + "/mainProjectPage/";

    d3.select("#projectDrop")
		.selectAll("li")
		.data(projectsNames)
		.enter()
			.append("li")
			.append("a")
				.attr("href", function(d){ return projLinkPre + d;})
				.text(function(d){return d;});
}


function createMainPageNavBar(wrappingNavSelector, rName){
	createMainPageNavBarSkeleton(wrappingNavSelector, rName);
	populateMainPageNavBar(wrappingNavSelector, rName);
}

