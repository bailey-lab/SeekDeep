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


function njhMenuItem(idName, displayName, func){
	this.idName = idName;
	this.displayName = displayName;
	this.func = func;
};

function createSeqMenu(idNameOfParentDiv, menuContent){
	//console.log(idNameOfParentDiv);
	var nav = d3.select(idNameOfParentDiv)
				.append("nav")
					.attr("id", "myNavBar")
					.attr("class", "navbar navbar-default")
					.attr("role", "navigation")
						.append("div")
							.attr("class", "container")
							.style("margin-left", 0)
								.append("div")
									.attr("class", "navbar-header")
										.append("ul")
										.attr("class", "nav navbar-nav");
	var menuKeys = Object.keys(menuContent);
	for (mk in menuKeys){
		var currentMenuItem = nav.append("li").attr("class", "dropdown");
		currentMenuItem.append("a")
			.attr("href", "#")
			.attr("data-toggle", "dropdown")
			.attr("class", "dropdown-toggle")
			.text(menuKeys[mk])
				.append("b").attr("class", "caret");
		var currentMenuItemOptions = currentMenuItem.append("ul").attr("class", "dropdown-menu").attr("id", menuKeys[mk] + "Drops");
		for (it in menuContent[menuKeys[mk]]){
			currentMenuItemOptions.append("li")
				.append("a")
					.attr("href", "javascript:void(0)")
					.attr("id",menuContent[menuKeys[mk]][it].idName)
					.text(menuContent[menuKeys[mk]][it].displayName)
					.on("click", menuContent[menuKeys[mk]][it].func);
		}
		//console.log(menuContent[menuKeys[mk]]);
	}
}

    
	function Canvas(i){
        this.canvas = $(i)[0];
        this.context = this.canvas.getContext("2d");
    }
	
	function SeqPainter(cellWidth, cellHeight, numOfSeqs, numOfBases, nameOffSet, baseColors){
		this.needToPaint = true;
		this.cw = cellWidth;
		this.ch = cellHeight;
		this.nSeqs = numOfSeqs;
		this.nBases = numOfBases;
		this.bColors = baseColors;
		this.nameOffSet = nameOffSet;
	}
	
	//draw given seq
	SeqPainter.prototype.paintSeq = function(seqContext, index, seq, start){
    	seqContext.textAlign = "center";
    	seqContext.textBaseline = "middle";
    	seqContext.font = "bold 15px Arial, sans-serif";
    	for(var wi = start; wi - start < this.nBases; wi++){
            if(wi < seq.length){
            	seqContext.fillStyle = this.bColors[seq[wi]];
            	seqContext.fillRect((wi - start) * this.cw + this.nameOffSet, index*this.ch, this.cw, this.ch);
            	seqContext.fillStyle = "#000000";
	        	seqContext.fillText(seq[wi], (wi -start)*this.cw + this.cw/2.0 + this.nameOffSet, index*this.ch + this.ch/2.0);
            }else{
            	seqContext.strokeStyle = "#000000";
    			seqContext.lineWidth   = 1;
            	seqContext.fillStyle = "#EEEEEE";
            	seqContext.fillRect((wi -start) * this.cw + this.nameOffSet, index*this.ch, this.cw, this.ch);
            	seqContext.strokeRect((wi -start) * this.cw + this.nameOffSet, index*this.ch, this.cw, this.ch);
            }
        }
	};
	
	//draw all necessary seqs
	SeqPainter.prototype.paintSeqs = function(seqContext, seqData, sStart, bStart){
    	if(this.needToPaint){
    		//console.log("painting");
    		//box around seqs
    		seqContext.strokeStyle = "#000000";
			seqContext.lineWidth   = 1;
	    	seqContext.strokeRect(this.nameOffSet, 0,this.cw * this.nBases,
	    		 Math.min(this.ch * this.nSeqs,this.ch * seqData.length ));
	    	//write names
	    	seqContext.textAlign = "left";
	    	seqContext.textBaseline = "middle";
	    	seqContext.font = "bold 15px Arial, sans-serif";
	        seqContext.strokeStyle = "#000000";
			seqContext.lineWidth   = 1;
	    	for(var hi = sStart;( hi -sStart < this.nSeqs ) && (hi < seqData.length); hi++){
	    		//console.log(hi);
	    		//console.log(this.nSeqs);
	    		seqContext.fillStyle = "#EEEEEE";
		    	seqContext.fillRect(0, (hi - sStart)*this.ch, this.nameOffSet, this.ch);
		    	seqContext.strokeRect(0, (hi - sStart)*this.ch, this.nameOffSet, this.ch);
		    	seqContext.fillStyle = "#000000";
				seqContext.fillText(seqData[hi]["name"], 2, (hi - sStart)*this.ch + this.ch/2.0);
	    	}	    	
	    	for(var hi = sStart;( hi -sStart < this.nSeqs ) && (hi < seqData.length) ; hi++){
	    		this.paintSeq(seqContext, hi - sStart,
	    		 seqData[hi]["seq"], bStart);
	    	}
	    	this.needToPaint = false;
	    	//console.log("end-painting");
    	}
	};
	//paint seq position
	SeqPainter.prototype.placeBasePos = function(seqContext, sStart, bStart){
    	seqContext.textAlign = "center";
    	seqContext.textBaseline = "middle";
    	seqContext.fillStyle = "#EEEEEE";
	    seqContext.fillRect(this.nameOffSet - this.cw, (this.nSeqs)*this.ch + 2 , this.cw * 2, this.ch);
	    seqContext.fillRect(this.nameOffSet - this.cw + this.nBases * this.cw -this.cw, (this.nSeqs)*this.ch + 2 , this.cw * 2, this.ch);
	    seqContext.fillRect(this.nameOffSet + this.nBases * this.cw + 2, 0 , this.cw * 2, this.ch);
	    seqContext.fillRect(this.nameOffSet + this.nBases * this.cw + 2, (this.nSeqs)*this.ch - this.ch , this.cw * 2, this.ch);
	    seqContext.fillStyle = "#000000";
    	seqContext.font = "bold 15px Arial, sans-serif";
    	//number of bases
      	seqContext.fillText(bStart, this.nameOffSet, (this.nSeqs)*this.ch +2 + this.ch/2.0);
      	seqContext.fillText(bStart + this.nBases -1, this.nameOffSet + this.nBases * this.cw -this.cw, (this.nSeqs)*this.ch +2 + this.ch/2.0);
   		//number of seqs
      	seqContext.fillText(sStart, this.nameOffSet + this.nBases * this.cw + (2.5/4) * this.cw , this.ch/2.0  );
   		seqContext.fillText(sStart + this.nSeqs - 1, this.nameOffSet + this.nBases * this.cw + (2.5/4) * this.cw, (this.nSeqs)*this.ch - this.ch+ this.ch/2.0  );
   };
   
   SeqPainter.prototype.paintSelectedSeq = function(seqContext, seq, currentBase){
   		seqContext.font = "bold 15px Arial, sans-serif";
   		seqContext.textAlign = "left";
    	seqContext.textBaseline = "middle";
   		var logInfo ="";
   		if(currentBase < 0){
   	   		logInfo = "name: " + 
        	seq["name"];
   		}else{
   	   		logInfo = "name: " + 
        	seq["name"]
        	+ " base: "  + seq["seq"][currentBase] 
        	+ " qual: " +  seq["qual"][currentBase]
        	+ " pos: " + currentBase;
   		}
        var tWidth = seqContext.measureText(logInfo).width;
        seqContext.fillStyle = "#FFFFFF";
	    seqContext.fillRect(this.nameOffSet + this.cw + 2, (this.nSeqs)*this.ch + 2 , this.nBases * this.cw - (3 *this.cw), this.ch);
        seqContext.fillStyle = "#EEEEEE";
	    seqContext.fillRect(this.nameOffSet + this.cw + 2, (this.nSeqs)*this.ch + 2 , tWidth, this.ch);
	    seqContext.fillStyle = "#000000";
	    //console.log(tWidth);
	    seqContext.fillText(logInfo,this.nameOffSet + this.cw + 2, (this.nSeqs)*this.ch + 2  + this.ch/2 );
   };

	function njhSeqView(viewName, seqData, cellWidth, cellHeight, baseColors,addQualChart, minTreeLink, protein){
		//need to add style and html, currently just there
		//retrieve html elements 
		this.topDivName = viewName;
		this.topDiv = document.getElementById(viewName);
		$(this.topDiv).addClass("SeqView");
		this.uid = seqData["uid"];
		this.selected = new Set(seqData["selected"]);
		//this.selected = seqData["selected"];
		//this.selected = [1,4,5];
		//this.selected.add(1);
		//this.selected.add(4);
		//console.log(this.selected);
		this.menuDiv = d3.select(viewName).append("div")
				.attr("class", "njhSeqViewMenu");
		this.setUpDefaultMenu(protein);
		if(minTreeLink != undefined && !protein){
			d3.select(viewName).append("div")
				.attr("id", "minTreeDiv")
				.style("border", "2px solid black")
				.style("padding", "2px")
				.style("margin", "5px")
				.style("width", "100px")
				.style("float", "left")
				.append("a")
					.attr("href", minTreeLink)
					.text("Show Minimum Spanning Tree");
		}
		d3.select(viewName).append("div")
			.attr("class", "downFastaDiv")
			.style("margin", "5px")
			.style("float", "left");
		d3.select(viewName).append("div")
		.attr("class", "deselectDiv")
		.style("margin", "5px")
		.style("float", "left");
		this.masterDivD3 = d3.select(viewName).append("div").attr("class", "SeqViewCanvasDiv");
		this.masterDivD3.append("div").attr("class", "rightSlider");
		this.masterDivD3.append("canvas").attr("class", "canvas");
		this.masterDivD3.append("div").attr("class", "bottomSlider");
		this.masterDivD3.append("div").attr("class", "pop-up").append("p").attr("id", "info");
		this.masterDivD3.append("div").attr("class", "select");
		d3.select(viewName).append("div").attr("class", "qualChart");
		d3.select(viewName).append("div")
			.attr("id", "minTreeChartTop");
		var self = this;
		var linkButton = d3.select(this.topDivName + " .downFastaDiv")
			.append("button")
			.text("Download Fasta")
			.attr("class", "fastaSaveButton btn btn-success");
		var deselectButton = d3.select(this.topDivName + " .deselectDiv")
			.append("button")
			.text("Un-select All")
			.attr("class", "deselectAllBut btn btn-info");
		deselectButton.on("click", function(){
			self.selected.clear();
			self.updateSelectors();
		});
		linkButton.append("a")
				.attr("class", "fastaDownLink");
		
		linkButton.on("click", function(){
				    var mainTable = [];
				    //
				    if (self.selected.size > 0){
				    	var sels = setToArray(self.selected);
					    for (i in sels) {
					    	//console.log(sels[i]);
							mainTable.push([">" + self.seqData["seqs"][sels[i]]["name"]]);
							mainTable.push([self.seqData["seqs"][sels[i]]["seq"]]);
						}
				    }else{
					    for (i = 0; i <self.seqData["seqs"].length ; i++) { 
							mainTable.push([">" + self.seqData["seqs"][i]["name"]]);
							mainTable.push([self.seqData["seqs"][i]["seq"]]);
						}
				    }

				  	var fastaData = 'data:text/plain;base64,'+ btoa(d3.tsv.format(mainTable));
				  	linkButton.select(".fastaDownLink").attr("download", self.seqData["uid"] + ".fasta");
				  	linkButton.select(".fastaDownLink").attr("href", fastaData).node().click();
				});
		//this.masterDiv = document.getElementById(viewName);
		this.masterDiv = this.masterDivD3.node();
		
		this.canvas = $(".canvas", this.masterDiv)[0];
		this.context = this.canvas.getContext('2d');
		
		this.rSlider = $(".rightSlider", this.masterDiv)[0];
		this.bSlider = $(".bottomSlider", this.masterDiv)[0];
		this.popUp = $(".pop-up", this.masterDiv)[0];
		this.sel = $(".select", this.masterDiv)[0];
		
		this.seqData = seqData;
		this.seqStart = 0;
		this.baseStart = 0;
		this.currentSeq = 0;
		this.currentBase = 0;
		// set up of sizes
		$(this.masterDiv).width((window.innerWidth - 10) * 0.98);
		$(this.masterDiv).height((window.innerHeight - 60) * 0.98);
		this.canvas.width = $(this.masterDiv).width() * 0.98;
		this.canvas.height = $(this.masterDiv).height() * 0.95;
		//make a mock canvas context using canvas2svg. We use a C2S namespace for less typing:
		//this.ctx = new C2S(this.canvas.width,this.canvas.height); //width, height of your desired svg file

        
		var nameOffSet = 10 * cellWidth;
		var numOfBases = Math.min(Math.floor((this.canvas.width - cellWidth - nameOffSet)/cellWidth),this.seqData["maxLen"] );
	 	var numOfSeqs = Math.min(Math.floor((this.canvas.height - cellHeight)/cellHeight), this.seqData["seqs"].length);
	 	//console.log(numOfSeqs);
	 	//console.log(Math.floor((this.canvas.height - cellHeight)/cellHeight));
	 	//console.log(this.seqData["seqs"].length);
		this.painter = new SeqPainter(cellWidth, cellHeight, numOfSeqs, numOfBases, nameOffSet, baseColors);
		
		if(addQualChart){
			this.addedQualChart = true;
			this.chart = c3.generate({
				bindto: this.topDivName + " .qualChart",
			    data: {
			        json: {
			            qual: this.seqData["seqs"][this.currentSeq]["qual"]
			        }
			    }, 
				grid: {
			        y: {
			            lines: [{value: 20}]
			        }
			    },
			    axis: {
			    	y : {
			    		max: 50,
			            min: 0
			    	}
			    }
			});
		}
		
	};
	
	
	
	njhSeqView.prototype.setUpDefaultMenu = function(protein){
		var locSplit = window.location.toString().split(/[\/]+/);
		var rName = locSplit[2];
		var menuItems = {};
		var sortOptions = [];
		var self = this;
		sortOptions.push(new njhMenuItem("sortSeq", "Sequence",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost('/' + rName + '/sort/seq',postData, function(md){ mainData = md; });
		    //console.log(self.uid);
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortSeqCondensed", "Sequence Condensed",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost('/' + rName + '/sort/seqCondensed',postData, function(md){ mainData = md; });
	    	self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortTotalCount", "Total Read Count",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost('/' + rName + '/sort/totalCount',postData, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortSize", "Length",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost('/' + rName + '/sort/size',postData, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortName", "Name",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost( '/' + rName + '/sort/name',postData, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortName", "Reverse",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost( '/' + rName + '/sort/reverse',postData, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		menuItems["Sort"] = sortOptions;
		var alnOptions = [];
		alnOptions.push(new njhMenuItem("muscle", "muscle",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost( '/' + rName + '/muscle',postData, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		alnOptions.push(new njhMenuItem("removeGaps", "remove gaps",function(){
			var mainData;
			var postData = {"uid" : self.uid};
			if (self.selected.size > 0){
				postData["selected"] = setToArray(self.selected);
			}
			ajaxPost( '/' + rName + '/removeGaps',postData, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		menuItems["Aln"] = alnOptions;
		

		if(!protein){
			var editOptions = [];
			editOptions.push(new njhMenuItem("complement", "Reverse Complement",function(){
				var mainData;
				var postData = {"uid" : self.uid};
				if (self.selected.size > 0){
					postData["selected"] = setToArray(self.selected);
				}
				var ar = setToArray(self.selected);
				ajaxPost( '/' + rName + '/complement',postData, function(md){ mainData = md; });
			    self.updateData(mainData);
			}));
			menuItems["Edit"] = editOptions;
		}

		if(!protein){
			var translateOptions = [];
			translateOptions.push(new njhMenuItem("translate", "Translate",function(){
				var mainData;
				//console.log($("#startSiteInput", self.topDivName).val());
				var postData = {"uid" : self.uid, "start" : $("#startSiteInput", self.topDivName).val()};
				if (self.selected.size > 0){
					postData["selected"] = setToArray(self.selected);
				}
				var ar = setToArray(self.selected);
				ajaxPost( '/' + rName + '/translate',postData, function(md){ mainData = md; });
			    //console.log(mainData);
			    //console.log(self.topDivName.substring(1) + "_protein");
			    if($("#" + self.topDivName.substring(1) + "_protein").length){
			    	self.proteinViewer.updateData(mainData);
			    }else{
			    	var proteinColors = {};
					ajax('/' + rName + '/proteinColors', function(bc){ proteinColors = bc; });
			    	$( "<div id = \"" + self.topDivName.substring(1) + "_protein" +   "\"></div>" ).insertAfter( self.topDivName );
			    	self.proteinViewer = new njhSeqView("#" + self.topDivName.substring(1) + "_protein", mainData, self.painter.cw, self.painter.ch, proteinColors,false, "", true);
			    	self.proteinViewer.setUp();
			    	self.proteinViewer.paint();
			    }
			    $("#" + self.topDivName.substring(1) + "_protein").scrollView();
			    
			}));
			menuItems["Translate"] = translateOptions;
		}
		var windowOptions = [];
		if(!protein){
			windowOptions.push(new njhMenuItem("ShowQual", "Hide Qual Graph",function(){
				if( self.addedQualChart){
					self.addedQualChart = false;
					d3.select(self.topDivName +  " .njhSeqViewMenu #ShowQual").text("Show Qual Graph");
					d3.select(self.topDivName + " .qualChart").selectAll("*").remove();
					self.chart = null;
				}else{
					self.addedQualChart = true;
					d3.select(self.topDivName +  " .njhSeqViewMenu #ShowQual").text("Hide Qual Graph");
					self.chart = c3.generate({
						bindto: self.topDivName + " .qualChart",
					    data: {
					        json: {
					            qual: self.seqData["seqs"][self.currentSeq]["qual"]
					        }
					    }, 
						grid: {
					        y: {
					            lines: [{value: 20}]
					        }
					    },
					    axis: {
					    	y : {
					    		max: 50,
					            min: 0
					    	}
					    }
					});
					self.chart.xgrids([{value: self.currentBase, text:self.seqData["seqs"][self.currentSeq]["qual"][self.currentBase]}]);
				}
			}));
			windowOptions.push(new njhMenuItem("GenTree", "Gen Difference Graph",function(){
				if(!($(self.topDivName + " #minTreeChartTop #saveButton").length)){
					d3.select(self.topDivName + " #minTreeChartTop").append("button")
					.style("float", "top")
					.attr("class", "btn btn-success")
					.attr("id", "saveButton")
					.style("margin", "2px")
					.text("Save As Svg");
					addSvgSaveButton(self.topDivName + " #minTreeChartTop #saveButton", self.topDivName + " #minTreeChartTop #minTreeChart #chart", self.seqData["uid"])
				}
				if(!($(self.topDivName + " #minTreeChartTop #minTreeChart").length)){
					d3.select(self.topDivName + " #minTreeChartTop").append("svg").attr("id", "minTreeChart")
					.attr("width", "0px")
					.attr("height", "0px")
					.style("margin-left", "10px")
				}else{
					d3.select(self.topDivName + " #minTreeChart").selectAll("*").remove();
				}
				var jsonData;
				var postData = {"uid" : self.uid};
				if (self.selected.size > 0){
					postData["selected"] = setToArray(self.selected);
				}
				postData["numDiff"] = $("#numDiffInput", self.topDivName).val();
				var ar = setToArray(self.selected);
				ajaxPost( '/' + rName + '/minTreeDataDetailed',postData, function(md){ jsonData = md; });
				drawPsuedoMinTreeDetailed(jsonData, self.topDivName + " #minTreeChart", "minTreeChart",
						$("#treeWidthInput", self.topDivName).val(),$("#treeHeightInput", self.topDivName).val());
				$('#minTreeChart').scrollView();
			}));
			windowOptions.push(new njhMenuItem("HideTree", "Hide Difference Graph",function(){
				d3.select(self.topDivName + " #minTreeChartTop").selectAll("*").remove();
			}));

			menuItems["Graphs"] = windowOptions;
		}

		
		createSeqMenu(this.topDivName + " .njhSeqViewMenu", menuItems);
		
		if(!protein){
			var startSiteInput = d3.select(self.topDivName +  " .njhSeqViewMenu #TranslateDrops")
			.append("li")
				.append("div")
					.attr("style", "padding: 3px 20px;")
				.append("form")
					.attr("class", "form-inline")
					.attr("id", "startSiteForm");
			startSiteInput
				.append("label")
					.attr("id", "startSiteLabel")
					.attr("for","startSiteInput")
					.attr("class", "control-label")
					.text("Start")
					.style("margin-right", "5px");
			var divInputGroup = startSiteInput
				.append("div")
				.attr("class", "input-group");
			divInputGroup.append("input")
				.attr("type", "number")
				.attr("class", "form-control")
				.attr("id", "startSiteInput")
				.attr("step", "1")
				.attr("min", "0")
				.attr("max", "2")
				.attr("value", "0");
			$('#startSiteForm').submit(function(e){
		        e.preventDefault();
		        //console.log($("#startSiteInput").val());
		    });
			var treeWidthInput = d3.select(self.topDivName +  " .njhSeqViewMenu #GraphsDrops")
			.append("li")
				.append("div")
					.attr("style", "padding: 3px 20px;")
				.append("form")
					.attr("class", "form-inline")
					.attr("id", "treeWidthForm");
			treeWidthInput
				.append("label")
					.attr("id", "treeWidthLabel")
					.attr("for","treeWidthInput")
					.attr("class", "control-label")
					.text("Graph Window Width")
					.style("margin-right", "5px");;
			var divtreeWidthInputGroup = treeWidthInput
				.append("div")
				.attr("class", "input-group");
			divtreeWidthInputGroup.append("input")
				.attr("type", "number")
				.attr("class", "form-control")
				.attr("id", "treeWidthInput")
				.attr("step", "100")
				.attr("min", "0")
				.attr("value", "1000");
			$('#treeWidthForm').submit(function(e){
		        e.preventDefault();
		        //console.log($("#startSiteInput").val());
		    });
			var treeHeightInput = d3.select(self.topDivName +  " .njhSeqViewMenu #GraphsDrops")
			.append("li")
				.append("div")
					.attr("style", "padding: 3px 20px;")
				.append("form")
					.attr("class", "form-inline")
					.attr("id", "treeHeightForm");
			treeHeightInput
				.append("label")
					.attr("id", "treeHeightLabel")
					.attr("for","treeHeightInput")
					.attr("class", "control-label")
					.text("Graph Window Height")
					.style("margin-right", "5px");;
			var divTreeHeightInputGroup = treeHeightInput
				.append("div")
				.attr("class", "input-group");
			divTreeHeightInputGroup.append("input")
				.attr("type", "number")
				.attr("class", "form-control")
				.attr("id", "treeHeightInput")
				.attr("step", "100")
				.attr("min", "0")
				.attr("value", "1000");
			$('#treeHeightForm').submit(function(e){
		        e.preventDefault();
		        //console.log($("#startSiteInput").val());
		    });
			
			var numDiffInput = d3.select(self.topDivName +  " .njhSeqViewMenu #GraphsDrops")
			.append("li")
				.append("div")
					.attr("style", "padding: 3px 20px;")
				.append("form")
					.attr("class", "form-inline")
					.attr("id", "numDiffForm");
			numDiffInput
				.append("label")
					.attr("id", "numDiffLabel")
					.attr("for","numDiffInput")
					.attr("class", "control-label")
					.text("Num Diff\n0=min to connect all")
					.style("margin-right", "5px");;
			var divNumDiffInputGroup = numDiffInput
				.append("div")
				.attr("class", "input-group");
			divNumDiffInputGroup.append("input")
				.attr("type", "number")
				.attr("class", "form-control")
				.attr("id", "numDiffInput")
				.attr("step", "1")
				.attr("min", "0")
				.attr("value", "0");
			$('#numDiffForm').submit(function(e){
		        e.preventDefault();
		        //console.log($("#startSiteInput").val());
		    });
		}	
	};
	
	/*
	 * <div id= "cutOffDiv">
					<form id="fracCutOffForm" class="form-inline">
						<label id= "fracCutOffLabel" for="fracCutOffInput" class="control-label"></label>
						<div class = "input-group">
							<input type='number' class="form-control" id = "fracCutOffInput" step="0.01" min="0" max = "100"/>
							<div class="input-group-addon">%</div>
						</div>
						<button type="submit" class="btn btn-primary">Submit</button>
					</form>
				</div>
	 */
	njhSeqView.prototype.updateData = function(inputSeqData){
		this.seqData = inputSeqData;
		this.uid = inputSeqData["uid"];
		this.seqStart = 0;
		this.baseStart = 0;
		this.currentSeq = 0;
		this.currentBase = 0;
		this.painter.nSeqs = Math.min(Math.floor((this.canvas.height - this.painter.ch)/this.painter.ch), this.seqData["seqs"].length);
		this.painter.needToPaint = true;
		this.needToPaint = true;
		this.updateCanvas();
		this.setUpSliders();
		this.paint();
	};
	njhSeqView.prototype.setUp = function(){
		this.setUpCanvas();
		this.setUpSliders();
		this.setUpListeners();
		this.setSelector();
	};
	
	njhSeqView.prototype.updateSeqDims = function(){
		this.painter.nBases = Math.min(Math.floor((this.canvas.width - this.painter.cw - this.painter.nameOffSet)/this.painter.cw),this.seqData["maxLen"] );
		this.painter.nSeqs = Math.min(Math.floor((this.canvas.height - this.painter.ch)/this.painter.ch), this.seqData["seqs"].length);
	}
	
	njhSeqView.prototype.setUpCanvas = function(){
		$(this.masterDiv).width((window.innerWidth - 10) * 0.98);
		var maxPossHeight = this.painter.ch * (this.seqData["seqs"].length + 4);
		$(this.masterDiv).height(Math.min((window.innerHeight - 60) * 0.80, maxPossHeight));
		this.canvas.width = $(this.masterDiv).width() * 0.96;
		this.canvas.height = $(this.masterDiv).height() * 0.95;
		this.updateSeqDims();
		this.updateSelectors();
	};
	
	njhSeqView.prototype.updateCanvas = function(){
		var changingHeight = (window.innerHeight - 60) * 0.80;
		var changingWidth = (window.innerWidth - 10) * 0.98;
		$(this.masterDiv).width((window.innerWidth - 10) * 0.98);
		var maxPossHeight = this.painter.ch * (this.seqData["seqs"].length + 4);
		$(this.masterDiv).height(Math.min((window.innerHeight - 60) * 0.80, maxPossHeight));
		this.canvas.width = $(this.masterDiv).width() * 0.96;
		this.canvas.height = $(this.masterDiv).height() * 0.95;
		if(changingHeight > this.canvas.height){
			this.painter.needToPaint = true;
		}
		if(changingWidth > this.canvas.width){
			this.painter.needToPaint = true;
		}
		this.updateSeqDims();
		this.updateSelectors();
	};

	njhSeqView.prototype.paint = function(){

		this.painter.paintSeqs(this.context, this.seqData["seqs"], this.seqStart, this.baseStart);
		this.painter.placeBasePos(this.context, this.seqStart, this.baseStart);
		//this.painter.needToPaint = true;
		//this.painter.paintSeqs(this.ctx, this.seqData["seqs"], this.seqStart, this.baseStart);
		//this.painter.placeBasePos(this.ctx, this.seqStart, this.baseStart);
		

		this.paintSelectedSeq();
		this.setSelector();
		this.updateSelectors();
		//var myRectangle = this.ctx.getSerializedSvg(true); //true here will replace any named entities with numbered ones.
        //.
		/*d3.select(".SeqViewCanvasDiv")
        	.append("div")
        	.attr("id", "testCanvasToSvg")
        	.html(myRectangle);;*/
        //d3.select("#SeqViewCanvasDiv").select("#testCanvasToSvg").html(myRectangle);
       
	};
	
	njhSeqView.prototype.paintSelectedSeq = function(){
		this.painter.paintSelectedSeq(this.context, this.seqData["seqs"][this.currentSeq], this.currentBase );
		//this.painter.paintSelectedSeq(this.ctx, this.seqData["seqs"][this.currentSeq], this.currentBase );
	};
	
	njhSeqView.prototype.setSelector = function(){
		if(this.currentBase >= this.baseStart && this.currentSeq >=this.seqStart){
			$(this.sel).css('top', (this.currentSeq - this.seqStart) *this.painter.ch -1);
			$(this.sel).css('left', (this.currentBase - this.baseStart) *this.painter.cw + this.painter.nameOffSet -1);
			$(this.sel).show();
		}else{
			$(this.sel).hide();
		}
		/*console.log("setSelector");
		console.log(this.sel);
		console.log(this.currentBase);
		console.log(this.currentSeq );
		console.log(this.currentSeq *this.painter.ch);
		console.log(this.currentBase *this.painter.cw + this.painter.nameOffSet );*/
	};
	
	njhSeqView.prototype.updateSelectors = function(){
		var self = this;
		var selectors = d3.select(this.topDivName + " .SeqViewCanvasDiv").selectAll(".seqHighlight")
			.data(setToArray(this.selected));
		var lowerBound = this.seqStart;
		var upperBound = this.seqStart + this.painter.nSeqs;
		selectors.enter().append("div")
			.attr("class", "seqHighlight")
			.style("width", function(d){ 
				return self.painter.nameOffSet.toString() + "px";})
			.style("height",function(d){ return self.painter.ch.toString() + "px";});
		selectors.exit().remove();
		
		selectors.style("visibility",function(d){ 
				if(d >= lowerBound && d < upperBound){
					return "visible";
				}else{
					return "hidden";
				}})
			.style("top", function(d){ 
				if(d >= lowerBound && d < upperBound){
					return ((d - self.seqStart) *self.painter.ch).toString() + "px"; 
				}else{
					return 0;
				}})
			.style("left",function(d){ 
				//console.log(d);
				if(d >= lowerBound && d < upperBound){
					return 0;
				}else{
					return 0;
				}});
			
			
		
			

	};
	
    njhSeqView.prototype.mouseWheelUp = function(steps){
        if(this.seqStart > 0){
        	--this.seqStart;
        	$(this.rSlider).slider('value', this.seqData["numReads"] - this.seqStart - this.painter.nSeqs);
        	this.painter.needToPaint = true;
        	this.paint();
        }
    };

    njhSeqView.prototype.mouseWheelDown = function(steps){
        if(this.seqStart < Math.max(this.seqData["numReads"]- this.painter.nSeqs, 0)){
        	++this.seqStart;
        	$(this.rSlider).slider('value', this.seqData["numReads"] - this.seqStart - this.painter.nSeqs);
        	this.painter.needToPaint = true;
        	this.paint();
        }
    };
    
	njhSeqView.prototype.setUpSliders = function(){
    	$( this.bSlider ).css("left", this.painter.nameOffSet);
    	$( this.bSlider).css("width", this.painter.nBases * this.painter.cw);
    	$( this.rSlider ).css("height", this.painter.nSeqs * this.painter.ch);
	    $( this.bSlider).slider({
	      range: "min",
	      min: 0,
	      max: Math.max(this.seqData["maxLen"] - this.painter.nBases, 0),
	      value: 0,
	      slide :function(event, ui){
	      	this.baseStart = ui.value;
	      	this.painter.needToPaint = true;
	      	this.paint();
	      }.bind(this)
	      }).bind(this);
	    $( this.rSlider ).slider({
	      range: "max",
	      min: 0,
	      max: Math.max(this.seqData["numReads"]- this.painter.nSeqs, 0),
	      value: this.seqData["numReads"],
	      orientation: "vertical", slide :function(event, ui){
	      	this.painter.needToPaint = true;
	      	this.seqStart = this.seqData["numReads"] - this.painter.nSeqs - ui.value;
	      	this.paint();
	      }.bind(this)
	    }).bind(this);
   };
   njhSeqView.prototype.updateOnResize = function(){
	   this.updateCanvas();
	   this.setUpSliders();
	   this.paint();
	   if(this.proteinViewer){
		   this.proteinViewer.updateOnResize();
	   }
   }
   njhSeqView.prototype.clicked = function(e){
        var pt = getRelCursorPosition(e, this.canvas);
        if(pt[1] <= this.painter.nSeqs * this.painter.ch && 
        		pt[0] <= this.painter.nBases * this.painter.cw + this.painter.nameOffSet){
        	this.currentSeq = Math.ceil(pt[1]/this.painter.ch) + this.seqStart - 1;
        	if(pt[0] <= this.painter.nameOffSet){
                //console.log(pt);
                //console.log(this.currentSeq);
                //console.log(this.selected);
        		if(this.selected.has(this.currentSeq)){
        			this.selected.delete(this.currentSeq);
        		}else{
        			this.selected.add(this.currentSeq);
        		}
        		this.updateSelectors();
        		//console.log(this.selected);
        	}
            this.currentBase = Math.ceil(pt[0]/this.painter.cw) - this.painter.nameOffSet/this.painter.cw + this.baseStart -1;
            this.paintSelectedSeq();


            this.setSelector();
            var currentQual = this.seqData["seqs"][this.currentSeq]["qual"];
            if(this.chart){
        		this.chart.load({
        	        json: {
        	            qual: this.seqData["seqs"][this.currentSeq]["qual"]
        	        }
        	    });
        	    //this.chart.xgrids.remove();
        	    this.chart.xgrids([{value: this.currentBase, text:this.seqData["seqs"][this.currentSeq]["qual"][this.currentBase]}]);
            }

        }

    };
   njhSeqView.prototype.setUpListeners = function(){
	var self = this;
   	// add scrolling listener
   	addMouseScrollListener(this.canvas, this.mouseWheelUp.bind(this), this.mouseWheelDown.bind(this));
   	this.canvas.addEventListener("mousedown", this.clicked.bind(this), false);
	//add hover box listening 
	var moveLeft = 20;
    var moveDown = 10;
    //object.hover(function(e) {
    $(this.canvas).hover(function(e) {
      //fadeInBox.fadeIn(500);
      $(this.popUp).fadeIn(500);
      //.css('top', e.pageY + moveDown)
      //.css('left', e.pageX + moveLeft)
      //.appendTo('body');
    }, function() {
      //fadeInBox.hide();
      $(this.popUp).hide();
    }).bind(this);
    var popUpWindow = this.popUp;
    //var painter = this.painter;
    //var seqs = this.seqData["seqs"];
    //var seqStart = this.seqStart;
    //var baseStart = this.baseStart;
    $(this.canvas).mouseleave(function(e) {
    	$(popUpWindow).hide();
    });
    $(this.canvas).mousemove(function(e) {
    	var currentPoint = getRelCursorPosition(e, self.canvas);
        if(currentPoint[1] <= self.painter.nSeqs * self.painter.ch &&
        		currentPoint[0] <= self.painter.nBases * self.painter.cw + self.painter.nameOffSet){
        	$(popUpWindow).css('top', currentPoint[1] + moveDown).css('left', currentPoint[0] + moveLeft);
          	var currentBaseHover = Math.ceil(currentPoint[0]/self.painter.cw) - self.painter.nameOffSet/self.painter.cw + self.baseStart -1;
            var currentSeqHover = Math.ceil(currentPoint[1]/self.painter.ch) + self.seqStart - 1;
            var base = self.seqData["seqs"][currentSeqHover]["seq"][currentBaseHover];
            var qual = self.seqData["seqs"][currentSeqHover]["qual"][currentBaseHover];
			if(currentPoint[0] > self.painter.nameOffSet){
	        	//console.log($("#info", popUpWindow)[0]);
    	        $("#info", popUpWindow)[0].innerHTML = "name: " + 
    	        	self.seqData["seqs"][currentSeqHover]["name"]
    	        	+ "<br>base: "  + base
    	        	+ "<br>qual: " +  qual 
    	        	+ "<br>pos: " + currentBaseHover;
	        }else{
	        	$("#info", popUpWindow)[0].innerHTML = "name: " + 
	        	self.seqData["seqs"][currentSeqHover]["name"];
	        }
			$(popUpWindow).show();
        }else{
        	$(popUpWindow).hide();
        }
    });
   };

   
function initSeqViewer(SeqViewer){
	$(window).bind("resize", function(){
		SeqViewer.updateCanvas();
		SeqViewer.setUpSliders();
		SeqViewer.paint();
	});
	SeqViewer.setUp();
	SeqViewer.paint();
}
