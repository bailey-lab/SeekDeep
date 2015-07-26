

function njhTable(masterDivId, tableMasterData, tableDownloadStubName, addChart){
	this.masterDivId = masterDivId;
	this.tableMasterData = tableMasterData;
	this.tableDownloadStubName = tableDownloadStubName;
	//create internal table divs
	//d3.select(masterDivId).append("div").attr("class", "njhTableMenu");
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
	menuOrganized.append("button").text("Download Table").attr("class", "njhTabSaveButton");
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
	/*menu.append("br");
	menu.append("button").text("Download Table").attr("class", "njhTabSaveButton");
	menu.select(".njhTabSaveButton").append("a").attr("class", "njhTabDownLink");
	menu.select(".njhTabSaveButton").on("click", function() {
		var allVals = [];
		menu.selectAll("input:checked").each(function() {
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
		var downLink = menu.select(".njhTabDownLink");
		downLink.attr("download", self.tableDownloadStubName + ".tab.csv");
		downLink.attr("href", dataSrc);
		downLink.node().click();
	}); */
	//create chart of numeric columns if needed
	this.chart;
	if (addChart) {
		this.addChart();
	}
}

/*njhTable.prototype.updateTableOnClick = function() {
	var allVals = [];
	d3.selectAll(this.masterDivId + " .njhTableMenu input:checked").each(function() {
		allVals.push($(this).val());
	});
	var currentColumnNames = _.intersection(this.tableMasterData["columnNames"], allVals);
	updateTable(this.tabDiv, this.tableMasterData["tab"], currentColumnNames);
}; */

njhTable.prototype.updateTableOnClickOrganized = function() {
	var allVals = [];
	d3.selectAll(this.masterDivId + " .njhTableMenuOrganized input:checked").each(function() {
		allVals.push($(this).val());
	});
	var currentColumnNames = _.intersection(this.tableMasterData["columnNames"], allVals);
	updateTable(this.tabDiv, this.tableMasterData["tab"], currentColumnNames);
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
	this.chart.hide(this.tableMasterData["hideOnStartColNames"]);
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




