$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var projectName = locSplit[4];
	var groupName = locSplit.pop();
	var infoUrls = getNameUrls(projectName);
	var gifLoading = prsentDivGifLoading();
	Promise.all(infoUrls.map(getJSON)).then(function(projNames) {
		addDiv("body", "topNav");
		var names = {shortName:projectName};
		projNames.map(function(name){$.extend(names,name);});
		createProjectNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle(names["shortName"] + " " + groupName);
		$("#jumboTitle").html(names["projectName"] + " " + groupName);
		addPanelWithDiv("#mainContent","sampLinks", "See Group");
		addPanelOnlyHead("#mainContent", "Population Infomation Table");
		addDiv("#mainContent", "popTable");
		//get group names 
		
		getJSON("/" + rName + "/groupsPopInfos/" + projectName + "/" + groupName).then(function(groupInfo){
			var cols = 10;
			var linkPre = "/" + rName + "/groupMainPage/" + projectName + "/" + groupName + "/";
			var mouseOverC = "#999";
			var mouseLeaveC = "#FFF";
			var addTo = "#sampLinks";
			createLinksTable(addTo, linkPre, groupInfo["groupNames"],cols, mouseOverC, mouseLeaveC);
			//create the population table and populate it 
			var popTable = new njhTable("#popTable", groupInfo["popInfo"], names["projectName"] + "_" + groupName  + "_popInfo", true);
		}).catch(logRequestError).then(function() {
		  //stop loading info
		});
	}).catch(logRequestError).then(function() {
		removeAllDivGifLoading();
	  //stop loading page
	});
});