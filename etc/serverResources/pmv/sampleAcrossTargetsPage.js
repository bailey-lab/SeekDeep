$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var sampleUrlName = locSplit.pop();
	var sampName = decodeURIComponent(sampleUrlName);
	var gifLoading = prsentDivGifLoading();
	addDiv("body", "topNav");
	addMainDiv("body", true);
	setHeadTitle(sampName);
	$("#jumboTitle").html(sampName);
	addPanelWithDiv("#mainContent", "sampleSummary", "Summary");
	addPanelOnlyHead("#mainContent", "Haplotypes in Each Target");
	addDiv("#mainContent", "sampleChartMaster");
	addPanelOnlyHead("#mainContent", "Haplotypes Table");
	addDiv("#mainContent", "sampTable");
	addPanelWithDiv("#mainContent", "notDetectedLinks", "Targets Without Haplotypes");
	addPanelOnlyHead("#mainContent", "Haplotypes Sequences");
	addDiv("#mainContent", "dnaViewer");
	Promise.all([getJSON('/' + rName + '/projectsNames'), getJSON('/' + rName + '/librarySampleNames')]).then(function (navNames) {
		createMainPageNavBar("#topNav", navNames[0]["projects"]);
		addNavDrop("#projectNav", "Samples", "samplesDrop");
		var sampleLinkPre = "/" + rName + "/sampleAcrossTargetsPage/";
		var samplesDrop = d3.select("#samplesDrop");
		//filter box kept at the top while the list scrolls
		var filter = samplesDrop.append("li")
			.style("position", "sticky")
			.style("top", "0")
			.style("z-index", "1")
			.style("background", "#fff")
			.style("padding", "4px 8px")
			.append("input")
				.attr("type", "search")
				.attr("class", "form-control input-sm")
				.attr("placeholder", "Filter samples");
		var noMatches = samplesDrop.append("li")
			.attr("class", "dropdown-header")
			.style("display", "none")
			.text("No matching samples");
		var sampleItems = samplesDrop.selectAll("li.sampleItem")
			.data(navNames[1]["samples"])
			.enter()
				.append("li")
					.attr("class", function(d){ return d == sampName ? "sampleItem active" : "sampleItem"; });
		sampleItems.append("a")
			.attr("href", function(d){ return sampleLinkPre + encodeURIComponent(d); })
			.text(function(d){ return d; });
		//case insensitive substring match, Enter opens the first match
		filter.on("input", function(){
			var term = this.value.trim().toLowerCase();
			var shown = 0;
			sampleItems.style("display", function(d){
				var match = "" === term || d.toLowerCase().indexOf(term) >= 0;
				if(match){
					++shown;
				}
				return match ? null : "none";
			});
			noMatches.style("display", 0 === shown ? null : "none");
		});
		filter.node().addEventListener("keydown", function(e){
			if("Enter" === e.key){
				e.preventDefault();
				var first = sampleItems.filter(function(){ return "none" !== this.style.display; }).select("a").node();
				if(first){
					window.location.href = first.getAttribute("href");
				}
			}
		});
		$("#samplesDrop").parent().on("shown.bs.dropdown", function(){
			filter.node().focus();
		});
	}).catch(logRequestError);
	getJSON('/' + rName + '/sampleAcrossTargets/' + sampleUrlName).then(function(info){
		var rows = info["tab"];
		var targets = _.uniq(rows.map(function(row){ return row["p_name"]; }));
		var coiByTarget = {};
		var totalReads = 0;
		rows.forEach(function(row){
			coiByTarget[row["p_name"]] = row["s_COI"];
			totalReads += row["c_ReadCnt"];
		});
		var cois = Object.keys(coiByTarget).map(function(target){ return coiByTarget[target]; });
		var summary = d3.select("#sampleSummary");
		summary.append("div").text("Haplotypes detected in " + targets.length + " of " + info["targetCount"] + " targets, " + totalReads + " reads");
		if(cois.length > 0){
			summary.append("div").text("Haplotypes per target: mean " + (cois.reduce(function(a, b){ return a + b; }, 0) / cois.length).toFixed(2)
				+ ", max " + Math.max.apply(null, cois)
				+ ", targets with more than one haplotype " + cois.filter(function(coi){ return coi > 1; }).length);
		}
		Object.keys(info["meta"]).sort().forEach(function(field){
			summary.append("div").text(field + ": " + info["meta"][field]);
		});
		var targetLinkPre = "/" + rName + "/mainProjectPage/";
		createLinksTable("#notDetectedLinks", targetLinkPre, info["targetsNotDetected"], 8, "#999", "#FFF");
		if(rows.length > 0){
			new njhSampleChart("#sampleChartMaster", info, sampName + "_acrossTargets", "p_name", "c_AveragedFrac", "c_clusterID",
				["p_name", "h_popUID", "c_clusterID", "c_AveragedFrac", "c_ReadCnt", "h_SampCnt", "h_SampFrac"]);
			new njhTable("#sampTable", info, sampName + "_acrossTargets", false);
			getJSON('/' + rName + '/sampleAcrossTargetsSeqs/' + sampleUrlName).then(function(seqData){
				new njhSeqView("#dnaViewer", seqData);
				setUpCloseSession(seqData["sessionUID"]);
			}).catch(logRequestError);
		}
	}).catch(logRequestError).then(function(){
		removeAllDivGifLoading();
	});
});
