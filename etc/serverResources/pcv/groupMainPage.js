$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var projectName = locSplit[4];
	var subGroupName = locSplit.pop();
	var groupName = locSplit.pop();

	var infoUrls = getNameUrls(projectName);
	infoUrls.push("/" + rName + "/groupSampleNames/" + projectName + "/" + groupName + "/" + subGroupName);
	var gifLoading = prsentDivGifLoading();
	Promise.all(infoUrls.map(getJSON)).then(function(projNames) {
		addDiv("body", "topNav");
		var names = {shortName:projectName, groupName:groupName};
		projNames.map(function(name){$.extend(names,name);});
		createProjectNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle(names["shortName"] + "_" + groupName + "_" + subGroupName);
		$("#jumboTitle").html(names["projectName"] + " " + groupName + ":" + subGroupName);
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
		createLinksTable(addTo, linkPre, names["groupSamples"],8, mouseOverC, mouseLeaveC);
		postJSON("/" + rName + "/groupSampInfo/" + projectName + "/" + groupName + "/" + subGroupName, {"sampNames":names["groupSamples"]}).then(function(sampInfo){
			//create sample table
			var sampleTable = new njhTable("#sampTable", sampInfo, names["projectName"] + "_sampInfo", false);
			var sampleChart = new njhSampleChart("#sampleChartMaster", sampInfo, names["projectName"] + "_" + groupName + "_" + subGroupName +  "_sampChart","s_Sample", "c_AveragedFrac","h_popUID", ["s_Sample", "h_popUID", "c_clusterID", "c_AveragedFrac", "c_ReadCnt"]);
			var popUrls = ["/" + rName + "/groupPopInfo/" + projectName + "/" + groupName + "/" + subGroupName]
			popUrls.push("/" + rName + "/groupPopSeqData/" + projectName + "/" + groupName + "/" + subGroupName);
			popUrls.push("/" + rName + "/groupHapIdTable/" + projectName + "/" + groupName + "/" + subGroupName);
			Promise.all(popUrls.map(function(popUrl){
				return postJSON(popUrl, {popUIDs:sampInfo["popUIDs"],"samples":names["groupSamples"]});
			})).then(function(popData){
				//0 is popIno, 1 is popSeqs
				var popInfoTab = popData[0];
				var popSeqData = popData[1];
				var hapIdTabInfo = popData[2];
				//set up seqs
				var sesUid = popSeqData["sessionUID"];
				var SeqViewer = new njhSeqView("#dnaViewer", popSeqData);
				setUpCloseSession(sesUid);
				var popTable = new njhTable("#popTable", popInfoTab, names["projectName"] + "_" + groupName + "_" + subGroupName  + "_popInfo", true);
				var hapIdTab = new njhTable("#hapIdTab", hapIdTabInfo, names["projectName"] + "_" + groupName + "_" + subGroupName + "_hapIdTable", false);
				hapIdTab.enactPoorMansHeatMap();
				function updateChartOnClick() {
					//get all currently checked sample boxes and then update the current samples
				    var allVals = [];
				    $('#sampNameMenu :checked').each(function() {
				      allVals.push($(this).val());
				    });
				    var currentSampNames = _.intersection(names["groupSamples"], allVals);
				    var currentPopInfoTab = {};
				    var currentPopSeqData = {};
				    gifLoading = prsentDivGifLoading();
				    postJSON("/" + rName + "/groupSampInfo/" + projectName + "/" + groupName + "/" + subGroupName,
				    		{"sampNames":currentSampNames, sessionUID:sesUid}).then(function(sampInfo){
				    	Promise.all(popUrls.map(function(popUrl){
							return postJSON(popUrl, {popUIDs:sampInfo["popUIDs"], sessionUID:sesUid,"samples":currentSampNames });
						})).then(function(popData){
							//for popData 0 is popIno, 1 is popSeqs
							sampleTable.updateWithData(sampInfo);
						 	sampleChart.updateWithData(sampInfo);
						 	popTable.updateWithData(popData[0]);
						 	SeqViewer.updateData(popData[1]);
						 	d34.select("#hapIdTab").selectAll("*").remove();
						 	hapIdTab = new njhTable("#hapIdTab", popData[2], names["projectName"] + "_" + groupName + "_" + subGroupName + "_hapIdTable", false);
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
				var sampMenu = new njhCheckboxMenu("#sampNameMenu", names["groupSamples"], updateChartOnClick);

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
