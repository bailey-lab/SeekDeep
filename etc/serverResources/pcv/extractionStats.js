
$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var projectName = locSplit[4];
	
	var infoUrls = getNameUrls(projectName);
	var gifLoading = prsentDivGifLoading();
	Promise.all(infoUrls.map(getJSON)).then(function(projNames) {
		addDiv("body", "topNav");
		var names = {shortName:projectName};
		projNames.map(function(name){$.extend(names,name);});
		createProjectNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle(names["projectName"]);
		$("#jumboTitle").html(names["projectName"] + ": Extraction");
		//change title to current name
		addPanelOnlyHead("#mainContent", "Sample Extraction Information");	
		addDiv("#mainContent", "profileTable");
		//addDiv("#mainContent", "sampNameMenu");
		addPanelOnlyHead("#mainContent", "Extraction By File");
		addDiv("#mainContent", "statsTable");
		var extractionTabsUrl = ["/" + rName + "/getExtractionProfileInfo/" + projectName, 
		                      "/" + rName + "/getExtractionStatsInfo/" + projectName];
		Promise.all(extractionTabsUrl.map(getJSON)).then(function(extractionTabs) {
			var mainProfileInfoTab = extractionTabs[0];
			var mainStatsInfoTab = extractionTabs[1];
			if(mainProfileInfoTab !== null){
				var sampleTable =  new njhTable("#profileTable", mainProfileInfoTab, names["projectName"] + "_extractionProfileInfo", true);	
				var popTable = new njhTable("#statsTable", mainStatsInfoTab, names["projectName"] + "_extractionStatsInfo", true);
			}
		}).catch(logRequestError).then(function() {
		  //stop loading 
		});
	}).catch(logRequestError).then(function() {
		removeAllDivGifLoading();
	  //stop loading page
	});
});
