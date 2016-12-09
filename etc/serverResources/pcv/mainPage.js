$(document).ready(function(){
	//get current name from window location
	var rName = getRootName();
	addDiv("body", "topNav");
	addMainDiv("body", false);
	setHeadTitle(rName);
	addPanelWithDiv("#mainContent","projectsLinks", "All Projects");	
	getJSON('/' + rName + '/projectsNames').then(function (projectsNames) {
		createMainPageNavBar("#topNav", projectsNames["projects"]);
		var cols = 10;
		var linkPre = "/" + rName + "/mainProjectPage/";
		var mouseOverC = "#999";
		var mouseLeaveC = "#FFF";
		var addTo = "#projectsLinks";
		createLinksTable(addTo, linkPre, projectsNames["projects"],8, mouseOverC, mouseLeaveC);
	}).catch(logRequestError);
	
});
    	
    	