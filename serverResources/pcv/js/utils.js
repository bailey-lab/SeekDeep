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
	var sampNameEncoding;
	ajax('/' + rName + '/sampleNamesEncoding/' + projectName, function(sne){ sampNameEncoding = sne; });
	ajax("/" + rName + "/sampleNames/" + projectName, function(mn){ sampNames = mn; });
	var encodedSampNames =[];
	for (sampN in sampNames){
		encodedSampNames.push(sampNameEncoding[sampNames[sampN]]);
	}
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
		.data(encodedSampNames)
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

