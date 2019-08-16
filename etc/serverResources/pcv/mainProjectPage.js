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
		$("#jumboTitle").html(names["projectName"]);
		addPanelWithDiv("#mainContent","extractionLinkDiv", "See Extraction Stats");
		d3.select("#extractionLinkDiv").append("a")
			.attr("id", "extractionLink")
			.text("See Extraction Info")
			.attr("href", function(){ return "/" + rName + "/extractionInfo/" + projectName;});;
		addPanelWithDiv("#mainContent", "groupLinks", "See Group Infos");
		addPanelWithDiv("#mainContent", "sampLinks", "See Sample Infos");
		addPanelOnlyHead("#mainContent", "Sample Info Table");
		addDiv("#mainContent", "sampTable");
		var sampHeader = addPanelOnlyHead("#mainContent", "Haplotypes Across Samples");
		sampHeader.attr("id", "sampleHeader");
		addDiv("#mainContent", "sampNameMenu");
		addDiv("#mainContent", "sampleChartMaster");
		addPanelOnlyHead("#mainContent", "Haplotypes Sequences");
		addDiv("#mainContent", "dnaViewer");
		addPanelOnlyHead("#mainContent", "Haplotype Information Table");
		addDiv("#mainContent", "popTable");
		addPanelOnlyHead("#mainContent", "Haplotype ID Table");
		addDiv("#mainContent", "hapIdTab");
		//get sample names and the table with the sample names
		//ajaxPost('/' + rName + '/sampInfo' + "/" + projectName, {"sampNames":sampNames}, function(tab){ mainPopInfoTab = tab; });
		var cols = 10;
		var linkPre = "/" + rName + "/sampleMainPage/" + projectName + "/";
		var mouseOverC = "#999";
		var mouseLeaveC = "#FFF";
		var addTo = "#sampLinks";
		createLinksTable(addTo, linkPre, names["samples"],8, mouseOverC, mouseLeaveC);
		var linkPreGroup = "/" + rName + "/groupInfoPage/" + projectName + "/";
		var addToGroup = "#groupLinks";
		createLinksTable(addToGroup, linkPreGroup, names["groups"],cols, mouseOverC, mouseLeaveC);
		postJSON("/" + rName + "/sampInfo/" + projectName, {"sampNames":names["samples"]}).then(function(sampInfo){
			//create sample table
			var sampleTable = new njhTable("#sampTable", sampInfo, names["projectName"] + "_sampInfo", false);
			var sampleChart = new njhSampleChart("#sampleChartMaster", sampInfo, names["projectName"] +  "_sampChart","s_Sample", "c_AveragedFrac","h_popUID", ["s_Sample", "h_popUID", "c_clusterID", "c_AveragedFrac", "c_ReadCnt"]);
			var popUrls = ["/" + rName + "/popInfo/" + projectName]
			popUrls.push("/" + rName + "/popSeqData/" + projectName);
			popUrls.push("/" + rName + "/hapIdTable/" + projectName);
			Promise.all(popUrls.map(function(popUrl){
				return postJSON(popUrl, {popUIDs:sampInfo["popUIDs"], samples:names["samples"]});
			})).then(function(popData){
				//0 is popIno, 1 is popSeqs
				var popInfoTab = popData[0];
				var popSeqData = popData[1];
				var hapIdTabInfo = popData[2];
				//set up seqs
				var sesUid = popSeqData["sessionUID"];
				var SeqViewer = new njhSeqView("#dnaViewer", popSeqData);
				setUpCloseSession(sesUid);
				var popTable = new njhTable("#popTable", popInfoTab, names["projectName"] + "_popInfo", true);
				var hapIdTab = new njhTable("#hapIdTab", hapIdTabInfo, names["projectName"] + "_hapIdTable", false);
				hapIdTab.enactPoorMansHeatMap();
				function updateChartOnClick() {
					//get all currently checked sample boxes and then update the current samples
				    var allVals = [];
				    $('#sampNameMenu :checked').each(function() {
				      allVals.push($(this).val());
				    });
				    var currentSampNames = _.intersection(names["samples"], allVals);
				    var currentPopInfoTab = {};
				    var currentPopSeqData = {};
				    gifLoading = prsentDivGifLoading();
				    postJSON("/" + rName + "/sampInfo/" + projectName,
				    		{"sampNames":currentSampNames, sessionUID:sesUid}).then(function(sampInfo){
				    	Promise.all(popUrls.map(function(popUrl){
							return postJSON(popUrl, {popUIDs:sampInfo["popUIDs"], sessionUID:sesUid, samples:currentSampNames});
						})).then(function(popData){
							//for popData 0 is popIno, 1 is popSeqs
							sampleTable.updateWithData(sampInfo);
						 	sampleChart.updateWithData(sampInfo);
						 	popTable.updateWithData(popData[0]);
						 	SeqViewer.updateData(popData[1]);
						 	d34.select("#hapIdTab").selectAll("*").remove();
						 	hapIdTab = new njhTable("#hapIdTab", popData[2], names["projectName"] + "_hapIdTable", false);
						 	hapIdTab.enactPoorMansHeatMap();

						 	$("#sampleHeader").scrollView(60, 0);
						}).catch(logRequestError).then(function(){
							//stop loading population info
						});
				    }).catch(logRequestError).then(function(){
						//stop loading sample table info
				    	removeAllDivGifLoading();
					});
				};
				//create samp menu
				var sampMenu = new njhCheckboxMenu("#sampNameMenu", names["samples"], updateChartOnClick);

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
