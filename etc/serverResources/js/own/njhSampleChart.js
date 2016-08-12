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



function njhSampleChart(masterDivId, masterData, downLoadStubName, xCol, yCol,colorBy,hoverCols){
	//Create the margins and set the width and height for the sample chart
	this.xCol = xCol;
	this.yCol = yCol;
	this.colorBy = colorBy;
	this.hoverCols = hoverCols;
	this.masterDivId = masterDivId;
	this.masterData = masterData;
	this.downLoadStubName = downLoadStubName;
	this.margin = {top: 20, right: 30, bottom: 110, left: 40},
	this.width = $(window).width() -100 - this.margin.left - this.margin.right,
	this.height = 600 - this.margin.top - this.margin.bottom;
	var self = this;
	//create tootip for sample chart
	this.tooltip = d3.select("body")
				.append("div")
				.style("position", "absolute")
				.style("visibility", "hidden")
				.style("background-color", "#88aaaa")
				.attr("id", this.masterDivId.substr(1) + "_popHover");
	this.tooltipId = "#" + this.masterDivId.substr(1)+ "_popHover";
	//create the table for the tooltip
	this.hoverTab = createTable(this.tooltipId);
	this.hoverTab.style("border", "1px solid black");
	this.hoverTab.style("box-shadow", "3px 3px 1.5px rgba(0, 0, 0, 0.5)");

	this.saveButton = d3.select(this.masterDivId)
						.append("button")
							.attr("class", "njhSampleChartSave btn btn-success")
							.style("margin", "2px")
							.text("Save as Svg");
	this.svgDiv = d3.select(this.masterDivId)
						.append("div")
							.attr("class", "njhSampleChartSvgDiv");
							
	//Select the chart and update with width and height
	this.masterSvg = this.svgDiv.append("svg");;
	this.chart = this.masterSvg
			  			.append("g");
	//create the x and y axes 
	this.x = d3.scale.ordinal();
	
	this.y = d3.scale.linear();
			    
	this.xAxis = d3.svg.axis()
			    .scale(this.x)
			    .orient("bottom");
			
	this.yAxis = d3.svg.axis()
			    .scale(this.y)
			    .orient("left");
	//set up the colors and create the info needed for hover tooltip

	this.color = d3.scale.ordinal();
	//add the axes to the chart
	this.chart.append("g")
		.attr("class", "x axis");
	
	this.chart.append("g")
    	.attr("class", "y axis")

	//add legend
	this.legend = this.chart.append("g")
				  .attr("class","legend")
	//actually draw everything 
	this.draw();
	//added a save svg button for the sample chart
	var downLink = this.saveButton.append("a").attr("class", "njhSampleChartImgDownload");
	this.saveButton.on("click", function(){
	  	var html = self.masterSvg
	        .attr("version", 1.1)
	        .attr("xmlns", "http://www.w3.org/2000/svg")
	        .node().parentNode.innerHTML;
	  	//add the svg information to a and then click it to trigger the download
	  	var imgsrc = 'data:image/svg+xml;base64,'+ btoa(html);
	  	downLink.attr("download", self.downLoadStubName + ".svg");
	  	downLink.attr("href", imgsrc);
	  	downLink.node().click();
	});
	$(window).bind("resize", function(){
		self.draw();
	});
}

njhSampleChart.prototype.draw = function(){
	var self = this;
	this.width = $(window).width() -100 - this.margin.left - this.margin.right,
	this.height = 600 - this.margin.top - this.margin.bottom;
							
	//Select the chart and update with width and height
	this.masterSvg.attr("width", this.width + this.margin.left + this.margin.right)
			  			.attr("height", this.height + this.margin.top + this.margin.bottom);
	this.chart.attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");
	//create the x and y axes 
	this.x.rangeRoundBands([0, this.width - 150], .1);
	
	this.y.range([this.height, 0]);

	this.color.domain(this.masterData["tab"].map(function(d) { return d[self.colorBy]; }))
	    .range(this.masterData["popColors"]);
	//initiate the counts to zero and no samples have been added
	var sampleSum = {};
	var sampleMax = {};
	var samplePopNames = {};
	for(pos in this.masterData["tab"]){
		sampleSum[this.masterData["tab"][pos][self.xCol]] = 0;
		sampleMax[this.masterData["tab"][pos][self.xCol]] = 0;
		samplePopNames[this.masterData["tab"][pos][self.xCol]] = [];
		//this.masterData["tab"]["color"] = this.color(this.masterData["tab"][pos][self.colorBy])
	}

	for(pos in this.masterData["tab"]){
		var temp = {};
		for (sPos in self.hoverCols){
			temp[self.hoverCols[sPos]] = self.masterData["tab"][pos][self.hoverCols[sPos]];
		}
		sampleMax[self.masterData["tab"][pos][self.xCol]] = sampleMax[self.masterData["tab"][pos][self.xCol]] + self.masterData["tab"][pos][self.yCol];
		temp["color"] = self.color(self.masterData["tab"][pos][self.colorBy]);
		samplePopNames[self.masterData["tab"][pos][self.xCol]].push(temp);
		/*
		var actualPos = self.masterData["tab"].length - 1 - pos;
		var temp = {};
		for (sPos in self.hoverCols){
			temp[self.hoverCols[sPos]] = self.masterData["tab"][actualPos][self.hoverCols[sPos]];
		}
		sampleMax[self.masterData["tab"][pos][self.xCol]] = sampleMax[self.masterData["tab"][pos][self.xCol]] + self.masterData["tab"][pos][self.yCol];
		temp["color"] = self.color(self.masterData["tab"][actualPos][self.colorBy]);
		samplePopNames[self.masterData["tab"][actualPos][self.xCol]].push(temp);
		*/
	}
	var yMax = 0;
	for(sPos in samplePopNames){
		var sum = 0;
		for (pos in samplePopNames[sPos]){
			sum += samplePopNames[sPos][pos][self.yCol];
		}
		if(sum > yMax){
			yMax = sum;
		}
	}
	yMax = Math.round((yMax + 0.00001) * 1000) / 1000
		
	
	//console.log(yMax);
	//set up the x domain for the current sample names 
	this.x.domain(this.masterData["tab"].map(function(d) { return d[self.xCol]; }));
	this.y.domain([0,yMax]);
	//add the axes to the chart
	this.chart.select(".x.axis")
	      .attr("transform", "translate(0," + this.height + ")")
	      .call(this.xAxis)
	      .selectAll("text")
		    .attr("y", 5)
		    .attr("x", 9)
		    .attr("dy", ".35em")
		    .attr("transform", "rotate(45)")
		    .style("text-anchor", "start")
		    .style("font-family","\"Helvetica Neue\",Helvetica,Arial,sans-serif")
		    .style("font-size", "12px")
		    .style("font-weight", "bold");
	this.chart.selectAll(".x.axis path").style("display", "none");
	
	this.chart.select(".y.axis")
      .call(this.yAxis)
      .selectAll("text")
      	.style("font-family", "\"Helvetica Neue\",Helvetica,Arial,sans-serif")
	    .style("font-size", "12px")
		.style("font-weight", "bold");
	
	this.chart.selectAll(".axis path")
		.style("fill", "none")
		.style("stroke", "#000")
		.style("shape-rendering", "crispEdges");
	this.chart.selectAll(".axis line")
		.style("fill", "none")
		.style("stroke", "#000")
		.style("shape-rendering", "crispEdges");
			
	var bars = this.chart.selectAll(".bar")
	      .data(this.masterData["tab"]);

	//updateTableWithColors(self.hoverTab,samplePopNames[d[self.xCol]], [xCol, "h_popUID",self.colorBy, "c_barcodeFrac"]);
	bars.enter().append("rect")
	      .attr("class", "bar")
	bars.attr("x", function(d) { return self.x(d[self.xCol]); })
	      .attr("y", function(d) { 
	      	//adjust height for previously added samples 
	      		var ret = self.y(sampleMax[d[self.xCol]] - sampleSum[d[self.xCol]]);
	      		sampleSum[d[self.xCol]] = sampleSum[d[self.xCol]] + d[self.yCol];
	      		return ret; 
	      		})
	      .attr("height", function(d) { return self.height - self.y(d[self.yCol]); })
	      .attr("width", self.x.rangeBand())
	      .attr("fill", function(d){ return self.color(d[self.colorBy]);})
	      .style("stroke", "#000")
	      .attr("data-legend",function(d) { return d[self.colorBy];})
	      .attr("id", function(d){return d[self.colorBy];})
	   		.attr("data-legend",function(d) { return d[self.colorBy];})
	   		.on("mouseover", function(d){
	   			var currentId = d3.select(this).attr("id")
	   			self.chart.selectAll("rect").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			self.legend.selectAll("circle").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			self.legend.selectAll("text").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			updateTableWithColors(self.hoverTab,samplePopNames[d[self.xCol]], self.hoverCols);
	   			return self.tooltip.style("visibility", "visible");})
	   		.on("mousemove", function(){return self.tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
	   		.on("mouseout", function(d){
	   			var currentId = d3.select(this).attr("id")
	   			self.chart.selectAll("rect").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			self.legend.selectAll("circle").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			self.legend.selectAll("text").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			return self.tooltip.style("visibility", "hidden");});
	bars.exit()
      .remove();
	//add legend
	this.legend.attr("transform","translate(" + (this.width - 140) + ",10)");
	this.legend.call(d3.legend);
	this.legend.selectAll("text")
		.style("font-family","\"Helvetica Neue\",Helvetica,Arial,sans-serif")
	    .style("font-size", "12px")
	    .style("font-weight", "bold");
	this.legend.selectAll("circle").on("mouseover", function(d){
			var currentId = d3.select(this).attr("id")
   			self.chart.selectAll("rect").transition().style('opacity',function () {
   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
   		    });
   			self.legend.selectAll("circle").transition().style('opacity',function () {
   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
   		    });
   			self.legend.selectAll("text").transition().style('opacity',function () {
   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
   		    });
   			})
   		.on("mouseout", function(d){
   			self.chart.selectAll("rect").transition().style('opacity',function () {
   		        return 1.0;
   		    });
   			self.legend.selectAll("circle").transition().style('opacity',function () {
   		        return 1.0;
   		    });
   			self.legend.selectAll("text").transition().style('opacity',function () {
   		        return 1.0;
   		    });
   			});
	self.legend.selectAll("text").on("mouseover", function(d){
		var currentId = d3.select(this).attr("id");
			self.chart.selectAll("rect").transition().style('opacity',function () {
		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
		    });
			self.legend.selectAll("circle").transition().style('opacity',function () {
		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
		    });
			self.legend.selectAll("text").transition().style('opacity',function () {
		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
		    });
			})
		.on("mouseout", function(d){
			self.chart.selectAll("rect").transition().style('opacity',function () {
		        return 1.0;
		    });
			self.legend.selectAll("circle").transition().style('opacity',function () {
		        return 1.0;
		    });
			self.legend.selectAll("text").transition().style('opacity',function () {
		        return 1.0;
		    });
			});
	//add legend style directly so download will download it as well
	this.chart.select(".legend-box").style("fill", "#fff");
	this.chart.select(".legend-items").style("fill", "#000");


}

njhSampleChart.prototype.updateWithResize = function(){
	/**@todo add ability to have a set window size and resize proportionally to window change rather than just being the size of the window */
}



njhSampleChart.prototype.updateWithData = function(updatedDataTab){
	var self = this;

	
	var exiting = this.masterData["tab"].length > updatedDataTab["tab"].length;
	
	//update the data
	this.masterData = updatedDataTab;
	//update the data for the bars
	var bars = this.chart.selectAll(".bar").data(this.masterData["tab"]);
	self.x.domain(self.masterData["tab"].map(function(d) { return d[self.xCol]; }));

	var modTime = 0
	if(exiting){
	    //remove the bars no longer needed
		bars.exit()
	    .transition()
	    .duration(750)
	    
	    .attr("width", 1e-6)
	    .attr("height", 1e-6)
	    .remove();
		modTime = 750;
	}
	setTimeout(function(){
		self.color = d3.scale.ordinal()
		    .domain(self.masterData["tab"].map(function(d) { return d[self.colorBy]; }))
		    .range(self.masterData["popColors"]);
		var sampleSum = {};
		var sampleMax = {};
		var samplePopNames = {};
		for(pos in self.masterData["tab"]){
			sampleSum[self.masterData["tab"][pos][self.xCol]] = 0;
			sampleMax[self.masterData["tab"][pos][self.xCol]] = 0;
			samplePopNames[self.masterData["tab"][pos][self.xCol]] = [];
			//self.masterData["tab"]["color"] = self.color(self.masterData["tab"][pos][self.colorBy])
		}

		for(pos in self.masterData["tab"]){
			var temp = {};
			for (sPos in self.hoverCols){
				temp[self.hoverCols[sPos]] = self.masterData["tab"][pos][self.hoverCols[sPos]];
			}
			sampleMax[self.masterData["tab"][pos][self.xCol]] = sampleMax[self.masterData["tab"][pos][self.xCol]] + self.masterData["tab"][pos][self.yCol];
			temp["color"] = self.color(self.masterData["tab"][pos][self.colorBy]);
			samplePopNames[self.masterData["tab"][pos][self.xCol]].push(temp);
			/*
			var actualPos = self.masterData["tab"].length - 1 - pos;
			var temp = {};
			for (sPos in self.hoverCols){
				temp[self.hoverCols[sPos]] = self.masterData["tab"][actualPos][self.hoverCols[sPos]];
			}
			sampleMax[self.masterData["tab"][pos][self.xCol]] = sampleMax[self.masterData["tab"][pos][self.xCol]] + self.masterData["tab"][pos][self.yCol];
			temp["color"] = self.color(self.masterData["tab"][actualPos][self.colorBy]);
			samplePopNames[self.masterData["tab"][actualPos][self.xCol]].push(temp);
			*/
		}
	
		//create new bars as needed
		bars.enter().append("rect")
			.attr("height", 1)
			.attr("width", 1)
		    .attr("class", "bar");
		//update the x axis
		self.chart.select("g.x.axis")
			  .transition().duration(750).ease("sin-in-out")
		      .call(self.xAxis);
		self.chart.select("g.y.axis")
		  .transition().duration(750).ease("sin-in-out")
	      .call(self.yAxis);
		self.chart.select("g.x.axis").selectAll("text")
			.transition().duration(750).ease("sin-in-out")
		    .attr("y", 5)
		    .attr("x", 9)
		    .attr("dy", ".35em")
		    .attr("transform", "rotate(45)")
		    .style("text-anchor", "start")
		    .style("font-size", "12px")
		    .style("font-family", "\"Helvetica Neue\",Helvetica,Arial,sans-serif")
		    .style("font-weight", "bold")
		self.chart.selectAll(".axis line")
			.transition().duration(750).ease("sin-in-out")
			.style("fill", "none")
			.style("stroke", "#000")
			.style("shape-rendering", "crispEdges");
		bars.transition()
	    	.duration(750)
		  .attr("x", function(d) { return self.x(d[self.xCol]); })
	      .attr("y", function(d) { 
	      		var ret = self.y(sampleMax[d[self.xCol]] - sampleSum[d[self.xCol]]);
	      		sampleSum[d[self.xCol]] = sampleSum[d[self.xCol]] + d[self.yCol];
	      		return ret; 
	      		})
	      .attr("height", function(d) { return self.height - self.y(d[self.yCol]); })
	      .attr("width", self.x.rangeBand())
	      .attr("fill", function(d){ return self.color(d[self.colorBy]);})
	      .style("stroke", "#000");
	      
	   bars.attr("id", function(d){return d[self.colorBy];})
	   		.attr("data-legend",function(d) { return d[self.colorBy];})
	   		.on("mouseover", function(d){
	   			var currentId = d3.select(this).attr("id")
	   			self.chart.selectAll("rect").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			self.legend.selectAll("circle").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			self.legend.selectAll("text").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			updateTableWithColors(self.hoverTab,samplePopNames[d[self.xCol]], self.hoverCols);
	   			return self.tooltip.style("visibility", "visible");})
	   		.on("mousemove", function(){return self.tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
	   		.on("mouseout", function(d){
	   			var currentId = d3.select(this).attr("id")
	   			self.chart.selectAll("rect").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			self.legend.selectAll("circle").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			self.legend.selectAll("text").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			return self.tooltip.style("visibility", "hidden");});
	   setTimeout(function(){
			//update legend
			self.legend.call(d3.legend);
			self.legend.selectAll("text")
				.style("font-family","\"Helvetica Neue\",Helvetica,Arial,sans-serif")
			    .style("font-size", "12px")
			    .style("font-weight", "bold");
			self.legend.selectAll("circle").on("mouseover", function(d){
				var currentId = d3.select(this).attr("id")
	   			self.chart.selectAll("rect").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			self.legend.selectAll("circle").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			self.legend.selectAll("text").transition().style('opacity',function () {
	   		        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
	   		    });
	   			})
	   		.on("mouseout", function(d){
	   			self.chart.selectAll("rect").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			self.legend.selectAll("circle").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			self.legend.selectAll("text").transition().style('opacity',function () {
	   		        return 1.0;
	   		    });
	   			});
		self.legend.selectAll("text").on("mouseover", function(d){
			var currentId = d3.select(this).attr("id");
				self.chart.selectAll("rect").transition().style('opacity',function () {
			        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
			    });
				self.legend.selectAll("circle").transition().style('opacity',function () {
			        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
			    });
				self.legend.selectAll("text").transition().style('opacity',function () {
			        return (d3.select(this).attr("id") === currentId) ? 1.0 : 0.25;
			    });
				})
			.on("mouseout", function(d){
				self.chart.selectAll("rect").transition().style('opacity',function () {
			        return 1.0;
			    });
				self.legend.selectAll("circle").transition().style('opacity',function () {
			        return 1.0;
			    });
				self.legend.selectAll("text").transition().style('opacity',function () {
			        return 1.0;
			    });
				});
		}, 750);
	}, modTime);
};

