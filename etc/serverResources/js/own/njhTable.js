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



function njhTable(masterDivId, tableMasterData, tableDownloadStubName, addChart){
	this.masterDivId = masterDivId;
	this.tableMasterData = tableMasterData;
	this.tableDownloadStubName = tableDownloadStubName;
	d3.select(masterDivId).attr("style", "margin-top:10px; margin-bottom:10px;")
		.attr("class", "njhTable");
	//create internal table divs
	//d3.select(masterDivId).attr("class", "njhTableMenu");
	d3.select(masterDivId).append("div").attr("class", "njhTableMenuOrganized");
	this.tabDiv = createTable(this.masterDivId);
	d3.select(masterDivId).append("div").attr("class", "njhTableChart");
	//populate table 
	updateTable(this.tabDiv,this.tableMasterData["tab"],this.tableMasterData["columnNames"]);
	//var menu = d3.select(this.masterDivId + " .njhTableMenu");
	//this.menu = new njhCheckboxMenu(this.masterDivId + " .njhTableMenu", this.tableMasterData["columnNames"],this.updateTableOnClick.bind(this) );
	
	var menuOrganized = d3.select(this.masterDivId + " .njhTableMenuOrganized");
	this.menuOrganized = new njhCheckboxMenuOrganized(this.masterDivId + " .njhTableMenuOrganized", this.tableMasterData["columnNames"],this.updateTableOnClickOrganized.bind(this) );
	
	var self = this;
	//add download button for table 
	menuOrganized.append("br");
	menuOrganized.append("button").text("Download Table").attr("class", "njhTabSaveButton btn btn-success");
	menuOrganized.select(".njhTabSaveButton").append("a").attr("class", "njhTabDownLink");
	menuOrganized.select(".njhTabSaveButton").on("click", function() {
		var allVals = [];
		menuOrganized.selectAll("input:checked").each(function() {
			allVals.push($(this).val());
		});
		var currentColumnNames = _.intersection(self.tableMasterData["columnNames"], allVals);
		var mainTable = [];
		mainTable.push(currentColumnNames);
		//
		for ( i = 0; i < self.tableMasterData["tab"].length; i++) {
			var currentRow = [];
			for (colNum in currentColumnNames) {
				currentRow.push(self.tableMasterData["tab"][i][currentColumnNames[colNum]]);
			}
			mainTable.push(currentRow);
		}
		var dataSrc = 'data:text/csv;base64,' + btoa(d3.tsv.format(mainTable));
		var downLink = menuOrganized.select(".njhTabDownLink");
		downLink.attr("download", self.tableDownloadStubName + ".tab.csv");
		downLink.attr("href", dataSrc);
		downLink.node().click();
	});
	//create chart of numeric columns if needed
	this.chart;
	if (addChart) {
		this.addChart();
	}
	if(this.tableMasterData["hideOnStartColNames"].length > 0){
		this.tableMasterData["hideOnStartColNames"].forEach(function(d){
			menuOrganized.select("#" + String(d).replaceAll(".", "\\.")).property("checked", false);
		});
		this.updateTableOnClickOrganized();
	}
}



njhTable.prototype.updateTableOnClickOrganized = function() {
	var allVals = [];
	d3.selectAll(this.masterDivId + " .njhTableMenuOrganized input:checked").each(function() {
		allVals.push($(this).val());
	});
	var currentColumnNames = _.intersection(this.tableMasterData["columnNames"], allVals);
	
	updateTable(this.tabDiv, this.tableMasterData["tab"], currentColumnNames);
	if(this.chart){
		var showCols = [];
		var hidCols = [];
		for(col in this.tableMasterData["numericColNames"]){
			//console.log(this.tableMasterData["numericColNames"][col]);
			if(arrayContains(currentColumnNames,this.tableMasterData["numericColNames"][col])){
				showCols.push(this.tableMasterData["numericColNames"][col]);
				//this.chart.show([this.tableMasterData["numericColNames"][col]]);
			}else{
				hidCols.push(this.tableMasterData["numericColNames"][col]);
				//this.chart.hide([this.tableMasterData["numericColNames"][col]]);
			}
		}
		this.chart.show(showCols);
		this.chart.hide(hidCols);
	}
}; 

njhTable.prototype.addChart = function(){
	this.addedChart = true;
	var xLabs = [];
	for (obj in this.tableMasterData["tab"]) {
		xLabs.push(this.tableMasterData["tab"][obj][this.tableMasterData["mainColName"]]);
	}
	this.chart = c3.generate({
		data : {
			json : this.tableMasterData["tab"],
			keys : {
				value : this.tableMasterData["numericColNames"],
			},
			type : 'bar'
		},
		axis : {
			x : {
				type : 'category',
				categories : xLabs,
				tick : {
					rotate : 75
				},
				height : 130
			}
		},
		bindto : this.masterDivId + " .njhTableChart"
	});
	if(this.tableMasterData["hideOnStartColNames"].length > 0){
		this.chart.hide(this.tableMasterData["hideOnStartColNames"]);
	}
};

njhTable.prototype.updateWithData = function(updatedDataTab){
	//as long as there isn't different column names this should work, other wise the table menu has to change
	/**@todo might want to add in the ability to update with new column names*/
	this.tableMasterData = updatedDataTab;
	//this.updateTableOnClick();
	this.updateTableOnClickOrganized();
	if(this.addedChart){
		this.addChart();
	}
};




