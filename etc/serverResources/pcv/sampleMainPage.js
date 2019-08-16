$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var projectName = locSplit[4];
	var sampName = locSplit.pop();

	var infoUrls = getNameUrls(projectName);
	var gifLoading = prsentDivGifLoading();
	Promise.all(infoUrls.map(getJSON)).then(function(projNames) {
		addDiv("body", "topNav");
		var names = {shortName:projectName, sampleName:sampName};
		projNames.map(function(name){$.extend(names,name);});
		createProjectNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle(sampName);
		$("#jumboTitle").html(sampName);
		addPanelOnlyHead("#mainContent", "Sample Info Table");
		addDiv("#mainContent", "sampTable");
		addPanelOnlyHead("#mainContent", "Haplotypes Relative Abundances");
		addDiv("#mainContent", "sampNameMenu");
		addDiv("#mainContent", "sampleChartMaster");
		addPanelOnlyHead("#mainContent", "Haplotypes Sequences");
		addDiv("#mainContent", "dnaViewer");
		var sampNames = [sampName];
		var mainPopInfoTab;
		postJSON("/" + rName + "/sampInfo/" + projectName, {"sampNames":sampNames}).then(function(sampInfo){
			//create sample table
			var sampleTable =  new njhTable("#sampTable", sampInfo, names["projectName"] + "_" + sampName + "_sampInfo", false);
			var sampleChart = new njhSampleChart("#sampleChartMaster", sampInfo, names["projectName"] + "_" + sampName + "_sampChart","s_Sample", "c_AveragedFrac","h_popUID", ["s_Sample", "h_popUID", "c_clusterID", "c_AveragedFrac", "c_ReadCnt"]);
			//get the seq and color data for the sequence view of the population sequences
			postJSON("/" + rName + "/sampSeqData/" + projectName, {"sampleName":sampName}).then(function(sampSeqData){
				//create SeqViewer for the population final sequences
				var sesUid = sampSeqData["sessionUID"];
				var SeqViewer = new njhSeqView("#dnaViewer", sampSeqData);
				setUpCloseSession(sesUid);
			}).catch(logRequestError).then(function(){
				//stop loading page
			});
		}).catch(logRequestError).then(function(){
			//stop loading page
		});
	}).catch(logRequestError).then(function() {
		removeAllDivGifLoading();
	  //stop loading page
	});
});
