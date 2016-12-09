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
function createProjectNavBarSkeleton(wrappingNavSelector, names){
	//get rootName
	var rName = getRootName();
	var projectName = names["projectName"];
	addFixedTopNavSkeleton(wrappingNavSelector, "SeekDeep", "mainNav", "projectNav");
 	addNavLink("#mainNav", "Home", "/" + rName, "#siteHomeLink");
	addNavDrop("#mainNav", "Projects", "projectDrop");
	addNavLink("#mainNav", projectName, "/" + rName + "/mainProjectPage/"+ projectName, "projectHomeLink", true);
	addNavLink("#projectNav", "Extraction Info", "/" + rName + "/extractionInfo/" + projectName, "extractionLink");
	addNavDrop("#projectNav", "Groups", "groupsDrop");
	addNavDrop("#projectNav", "Samples", "samplesDrop");
}

//names should have fields should have at least projectName and possibly groupName,sampName
//names should also have array fields projectNames, groupNames, sampNames
function populateProjectNavBar(wrappingNavSelector, names){
	//get rootName
	var rName = getRootName();
   	//get project names
	var projectName = names["shortName"];
	var projectsNames = names["projects"];
	var projLinkPre = "/" + rName + "/mainProjectPage/";
	//get group names
	var groupName = (names.hasOwnProperty("groupName")) ?   names["groupName"] : "";
	var groupNames = names["groups"];
    var linkPreGroup = "/" + rName + "/groupInfoPage/" + projectName + "/";
	//get sample names
    var sampName = (names.hasOwnProperty("sampleName")) ?   names["sampleName"] : "";
	var sampNames = names["samples"];
    var linkPreSample = "/" + rName + "/sampleMainPage/" + projectName + "/";
	
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


function createProjectNavBar(wrappingNavSelector, names){
	createProjectNavBarSkeleton(wrappingNavSelector, names);
	populateProjectNavBar(wrappingNavSelector, names);
}



function createMainPageNavBarSkeleton(wrappingNavSelector){
	var rName = getRootName();
	addFixedTopNavSkeleton(wrappingNavSelector, "SeekDeep", "mainNav", "projectNav");
 	addNavLink("#mainNav", "Home", "/" + rName, "#siteHomeLink");
	addNavDrop("#mainNav", "Projects", "projectDrop");
}

function populateMainPageNavBar(wrappingNavSelector, projectsNames){
	var rName = getRootName();
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


function createMainPageNavBar(wrappingNavSelector, projectsNames){
	createMainPageNavBarSkeleton(wrappingNavSelector);
	populateMainPageNavBar(wrappingNavSelector, projectsNames);
}

function getNameUrls(projectName){
	var rName = getRootName();
	var infoUrls = ["/" + rName + "/projectsNames/"];
	infoUrls.push("/" + rName + "/projectName/" + projectName);
	infoUrls.push("/" + rName + "/sampleNames/" + projectName);
	infoUrls.push("/" + rName + "/getGroupNames/" + projectName)
	return infoUrls;
}


