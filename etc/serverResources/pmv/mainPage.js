$(document).ready(function(){
	//get current name from window location
	var rName = getRootName();
	addDiv("body", "topNav");
	addMainDiv("body", false);
	setHeadTitle(rName);
	addPanelWithDiv("#mainContent", "pmoHeader", "PMO");
	addPanelWithDiv("#mainContent", "projectsLinks", "All Targets");
	addPanelWithDiv("#mainContent", "rainbowLinks", "Haplotype Rainbows");
	addPanelWithDiv("#mainContent", "sampleLinks", "Samples Across All Targets");
	addPanelOnlyHead("#mainContent", "Targets");
	addDiv("#mainContent", "targetsTable");
	addPanelOnlyHead("#mainContent", "Library Samples");
	addDiv("#mainContent", "librarySamplesTable");
	getJSON('/' + rName + '/pmoHeader').then(function (header) {
		var method = header["generation_method"];
		var lines = ["PMO version " + header["pmo_version"]];
		if(header["creation_date"]){
			lines.push("Created " + header["creation_date"]);
		}
		if(method){
			lines.push("By " + method["program_name"] + " " + method["program_version"]);
		}
		d3.select("#pmoHeader").append("div").text(lines.join(", "));
	}).catch(logRequestError);
	getJSON('/' + rName + '/projectsNames').then(function (projectsNames) {
		createMainPageNavBar("#topNav", projectsNames["projects"]);
		var linkPre = "/" + rName + "/mainProjectPage/";
		createLinksTable("#projectsLinks", linkPre, projectsNames["projects"], 8, "#999", "#FFF");
	}).catch(logRequestError);
	getJSON('/' + rName + '/metaGroups').then(function (metaGroups) {
		var names = ["All samples"];
		var urls = [""];
		Object.keys(metaGroups["groups"]).sort().forEach(function(group){
			metaGroups["groups"][group].forEach(function(subGroup){
				names.push(group + ": " + subGroup);
				urls.push("/" + encodeURIComponent(group) + "/" + encodeURIComponent(subGroup));
			});
		});
		createLinksTableDifferentLinks("#rainbowLinks", "/" + rName + "/haplotypeRainbowPage", names, urls, 6, "#999", "#FFF");
	}).catch(logRequestError);
	getJSON('/' + rName + '/librarySampleNames').then(function (samples) {
		createLinksTable("#sampleLinks", "/" + rName + "/sampleAcrossTargetsPage/", samples["samples"], 8, "#999", "#FFF");
	}).catch(logRequestError);
	getJSON('/' + rName + '/targetsInfo').then(function (targetsInfo) {
		new njhTable("#targetsTable", targetsInfo, "pmoTargets", false);
	}).catch(logRequestError);
	getJSON('/' + rName + '/librarySamplesInfo').then(function (samplesInfo) {
		new njhTable("#librarySamplesTable", samplesInfo, "pmoLibrarySamples", false);
	}).catch(logRequestError);
});
