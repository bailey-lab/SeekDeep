var createMinTree = function(data, appendTo, name, width, height){
	// create a interconnected graph for a minimum spanning tree
	// data needs to have a "nodes" array and a "links" array 
	// nodes need to have at least the following variables, color (color of node), 
	// name (name of the node), and size (the size of the node)
	// links need to have at least color (the color of the link), target (the node position to connect to), 
	// source (the node position from where the connection is forming), and value (the value that controls the link distance)
	var svg = d3.select(appendTo).append("svg")
	    .attr("width", width)
	    .attr("height", height)
	    .attr("id", name);
	
	d3.json(data, function(error, graph)
	 {
	 	svg.append('rect')
	    .attr('width',width)
	    .attr('height', height)
	    .attr('fill',"#000");
	 var force = d3.layout.force()
	    .charge(-120)
	    .linkDistance(function(d, i){ return d.value * 10;})
	    .size([width, height]);
	    
	  force
	      .nodes(graph.nodes)
	      .links(graph.links)
	      .start();
	
	  var link = svg.selectAll(".link")
	  	  .data(graph.links)
	      .enter().append("line")
	      .attr("class", "link")
	      .style("stroke-width", function(d) { return Math.sqrt(d.value);})
	      .style("stroke", function(d) { return d.color;});
	 var node = svg.selectAll(".node")
	 	.data(graph.nodes)
	    .enter()
	    .append("circle")
	      .attr("class", "node")
	      .attr("r", function(d){ return d.size * 8;})
	      .style("fill", function(d) { return d.color; })
	      .call(force.drag);
	  node.append("title")
	      .text(function(d) { return d.name; });
	
	  force.on("tick", function() {
	    link.attr("x1", function(d) { return d.source.x; })
	        .attr("y1", function(d) { return d.source.y; })
	        .attr("x2", function(d) { return d.target.x; })
	        .attr("y2", function(d) { return d.target.y; });
	
	    node.attr("cx", function(d) { return d.x; })
	        .attr("cy", function(d) { return d.y; });

	  });
	});
	return svg;
};

var createMinTreeRawData = function(graph, appendTo, name, width, height){
	// create a interconnected graph for a minimum spanning tree
	// data needs to have a "nodes" array and a "links" array 
	// nodes need to have at least the following variables, color (color of node), 
	// name (name of the node), and size (the size of the node)
	// links need to have at least color (the color of the link), target (the node position to connect to), 
	// source (the node position from where the connection is forming), and value (the value that controls the link distance)
	var svg = d3.select(appendTo).append("svg:svg")
	    .attr("width", width)
	    .attr("height", height)
	    .attr("id", name);
	var vis = svg
  		.append('svg:g');/*
    	.call(d3.behavior.zoom().on("zoom", rescale))
    	.on("dblclick.zoom", null)
	  .append('svg:g')
	    .on("mousemove", mousemove)
	    .on("mousedown", mousedown)
	    .on("mouseup", mouseup);*/



	vis.append('svg:rect')
	    .attr('width', width)
	    .attr('height', height)
	    .attr('fill', '#FFF');
	 var force = d3.layout.force()
	    .charge(-120)
	    .linkDistance(function(d, i){ return d.value * 10;})
	    .size([width, height]);
	    
	  force
	      .nodes(graph.nodes)
	      .links(graph.links)
	      .start();
	
	  
	 var node = vis.selectAll(".node")
	 	.data(graph.nodes)
	    .enter().append("circle")
	      .attr("class",function(d) { return  "node " + d.name; })
	      .attr("r", function(d){ 
	      	if(d.size == 0.15){
	      		return 2.5;
	      	}else{
	      		return d.size * 8;
	      	}
	      	})
		  .style("fill", function(d) { 
		  	if(d.name == "00"){
		  		return "#489C28";
		  	}else if (d.name == "01"){
		  		return "#890309";
		  	}else if (d.name == "02"){
		  		return "#AED22D";
		  	}else if (d.name == "03"){
		  		return "#FEF035";
		  	}else if (d.name == "05"){
		  		return "#F68423";
		  	}else if (d.name == "06"){
		  		return "#BF4315";
		  	}else if (d.name == "08"){
		  		return "#1B5E44";
		  	}else if (d.name == "11"){
		  		return "#FDBB2C";
		  	}else {
		  		return d.color; 
		  	}
		  	})
	      .call(force.drag);
	  node.append("title")
	      .text(function(d) { return d.name; });
	  var link = vis.selectAll(".link").data(graph.links)
	      .enter().append("line")
	      .attr("class", "link")
	      .style("stroke-width", function(d) { return 1;})
	      .style("stroke", function(d) { return "#000";});
	  force.on("tick", function() {
	    link.attr("x1", function(d) { return d.source.x; })
	        .attr("y1", function(d) { return d.source.y; })
	        .attr("x2", function(d) { return d.target.x; })
	        .attr("y2", function(d) { return d.target.y; });
	    node.attr("cx", function(d) { return d.x; })
	        .attr("cy", function(d) { return d.y; });});
	return svg;
};

function drawPsuedoMinTreeDetailed(jsonData, addTo, hovIdStub,width,height){
	tooltipId = "#" + hovIdStub + "_popHover";
	if($(tooltipId).length){
		d3.select(tooltipId).remove();
	}
	var tooltip = d3.select("body")
				.append("div")
				.style("position", "absolute")
				.style("visibility", "hidden")
				.style("background-color", "#88aaaa")
				.attr("id", hovIdStub + "_popHover");
	
	d3.select(tooltipId).append("h3")
		.attr("id", hovIdStub + "_popHover_title")
		.style("margin", "5px");
	hoverTab = createTable(tooltipId);
	hoverTab.style("border", "1px solid black");
	hoverTab.style("box-shadow", "3px 3px 1.5px rgba(0, 0, 0, 0.5)");
	if(!width){
		width = 1000;
	}
	if(!height){
		height = 1000;
	}
	//var width = 1000, height = 1000;
	
	var force = d3.layout.force()
    	.charge(-120)
    	.linkDistance(30)
    	.size([width, height]);
	d3.select(addTo)
    	.attr("width", width)
    	.attr("height", height);
	var svg = d3.select(addTo).append("svg").attr("width", width).attr(
			"height", height).attr("id", "chart");

	svg.append('rect').attr('width', width).attr('height', height).attr('fill',
			jsonData.backgroundColor);
	force.nodes(jsonData.nodes).links(jsonData.links).on("tick", tick).start();

	var link = svg.selectAll(".link").data(jsonData.links).enter().append(
			"line").attr("class", "link").style("stroke", function(d) {
		return d.color;
	}).style("stroke-width", function(d) {
		return 2;
	});
	var node = svg.selectAll(".node").data(force.nodes()).enter().append("g")
			.attr("class", function(d) {
				return "node " + d.type;
			}).call(force.drag);
	// add the nodes
	node.append("circle")
      .attr("r", function(d) { return Math.sqrt(d.size/Math.PI); })
      .style("fill", function(d) { return d.color; })
      .style("stroke", "#fff")
      .style("stroke-width", "1.5px")
      .call(force.drag);

  var indels = svg.selectAll(".indel");
  var snps = svg.selectAll(".snp");
  indels.selectAll("circle").style("stroke", "#F00").style("stroke-width", 1);
  
  indels.selectAll("circle").on("mouseover", function(d){
  			d3.select("#" + hovIdStub + "_popHover_title").html("Indel")
  			updateTable(hoverTab, [{"SeqName":d.ref, "Pos":d.refPos, "GapSeq":d.refDisplay},{"SeqName":d.seq, "Pos":d.seqPos, "GapSeq":d.seqDisplay}], ["SeqName", "Pos", "GapSeq"]);
	      	 return tooltip.style("visibility", "visible");})
		  .on("mousemove", function(){return tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
		  .on("mouseout", function(){return tooltip.style("visibility", "hidden");});
  snps.selectAll("circle").style("stroke", "#00F").style("stroke-width", 1);
  snps.selectAll("circle").on("mouseover", function(d){
  			d3.select("#" + hovIdStub + "_popHover_title").html("Snp")
			 updateTable(hoverTab, [{"SeqName":d.ref, "Pos":d.refPos, "Base":d.refBase},{"SeqName":d.seq, "Pos":d.seqPos, "Base":d.seqBase}], ["SeqName", "Pos", "Base"]);
	      	 return tooltip.style("visibility", "visible");})
		  .on("mousemove", function(){return tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
		  .on("mouseout", function(){return tooltip.style("visibility", "hidden");});
		  
  var variants = svg.selectAll(".variant");
  variants.selectAll("circle").style("stroke", "#999");
  variants.selectAll("circle").on("mouseover", function(d){
  			 d3.select("#" + hovIdStub + "_popHover_title").html("Variant")
             updateTable(hoverTab, [{"SeqName":d.name, "ReadCnt":d.cnt, "Relative Abundance":d.frac}], ["SeqName", "ReadCnt", "Relative Abundance"]);
  			//tooltip.node().innerHTML = "<strong>Variant</strong>: " + d.name 
	      	 return tooltip.style("visibility", "visible");})
		  .on("mousemove", function(){return tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
		  .on("mouseout", function(){return tooltip.style("visibility", "hidden");});
  
  variants.append("text")
	  .attr("x", 12)
	  .attr("dy", ".35em")
	  .style("fill","#FFF")
	  .style("font-family", "\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica Neue\", Helvetica, Arial, \"Lucida Grande\", sans-serif")
	  .style("font-size", "12px")
	  .style("font-weight","900")
	  .style("pointer-events", "none")
	  .style("stroke", "#FFF")
	  .style("stroke-width", "1px")
	  .text(function(d) {return d.name;});
// add the curvy lines
function tick() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
};

}

function drawPsuedoMinTreeDetailedLink(jsonData, addTo){
	var tooltip = d3.select("body")
				.append("div")
				.style("position", "absolute")
				.style("visibility", "hidden")
				.style("background-color", "#88aaaa")
				.attr("id", addTo.substr(1) + "_popHover");
	tooltipId = "#" + addTo.substr(1) + "_popHover";
	d3.select(tooltipId).append("h3")
		.attr("id", addTo.substr(1) + "_popHover_title")
		.style("margin", "5px");
	hoverTab = createTable(tooltipId);
	hoverTab.style("border", "1px solid black");
	hoverTab.style("box-shadow", "3px 3px 1.5px rgba(0, 0, 0, 0.5)");
	var width = 2000,
    height = 2000;

var force = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);
	d3.select(addTo)
    	.attr("width", width)
    	.attr("height", height);
var svg = d3.select(addTo).append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("id", "chart");

d3.json(jsonData, function(error, graph) {
  svg.append('rect')
	    .attr('width',width)
	    .attr('height', height)
	    .attr('fill',graph.backgroundColor);
  
  force
      .nodes(graph.nodes)
      .links(graph.links)
      .on("tick",tick)
      .start();

  var link = svg.selectAll(".link")
      .data(graph.links)
    .enter().append("line")
      .attr("class", "link")
      .style("stroke", function(d) { return d.color; })
      .style("stroke-width", function(d) { return 2; });
var node = svg.selectAll(".node")
    .data(force.nodes())
  .enter().append("g")
    .attr("class", function(d) { return "node " + d.type; })
    .call(force.drag);

// add the nodes
	node.append("circle")
      .attr("r", function(d) { return Math.sqrt(d.size/Math.PI); })
      .style("fill", function(d) { return d.color; })
      .style("stroke", "#fff")
      .style("stroke-width", "1.5px")
      .call(force.drag);

  var indels = svg.selectAll(".indel");
  var snps = svg.selectAll(".snp");
  indels.selectAll("circle").style("stroke", "#F00").style("stroke-width", 1);
  
  indels.selectAll("circle").on("mouseover", function(d){
  			d3.select("#" + addTo.substr(1) + "_popHover_title").html("Indel")
  			updateTable(hoverTab, [{"SeqName":d.ref, "Pos":d.refPos, "GapSeq":d.refDisplay},{"SeqName":d.seq, "Pos":d.seqPos, "GapSeq":d.seqDisplay}], ["SeqName", "Pos", "GapSeq"]);
	      	 return tooltip.style("visibility", "visible");})
		  .on("mousemove", function(){return tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
		  .on("mouseout", function(){return tooltip.style("visibility", "hidden");});
  snps.selectAll("circle").style("stroke", "#00F").style("stroke-width", 1);
  snps.selectAll("circle").on("mouseover", function(d){
  			d3.select("#" + addTo.substr(1) + "_popHover_title").html("Snp")
			 updateTable(hoverTab, [{"SeqName":d.ref, "Pos":d.refPos, "Base":d.refBase},{"SeqName":d.seq, "Pos":d.seqPos, "Base":d.seqBase}], ["SeqName", "Pos", "Base"]);
	      	 return tooltip.style("visibility", "visible");})
		  .on("mousemove", function(){return tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
		  .on("mouseout", function(){return tooltip.style("visibility", "hidden");});
		  
  var variants = svg.selectAll(".variant");
  variants.selectAll("circle").style("stroke", "#999");
  variants.selectAll("circle").on("mouseover", function(d){
  			 d3.select("#" + addTo.substr(1) + "_popHover_title").html("Variant")
             updateTable(hoverTab, [{"SeqName":d.name, "ReadCnt":d.cnt, "Relative Abundance":d.frac}], ["SeqName", "ReadCnt", "Relative Abundance"]);
  			//tooltip.node().innerHTML = "<strong>Variant</strong>: " + d.name 
	      	 return tooltip.style("visibility", "visible");})
		  .on("mousemove", function(){return tooltip.style("top", (d3.event.layerY-10)+"px").style("left",(d3.event.layerX+10)+"px");})
		  .on("mouseout", function(){return tooltip.style("visibility", "hidden");});
  
  variants.append("text")
	  .attr("x", 12)
	  .attr("dy", ".35em")
	  .style("fill","#FFF")
	  .style("font-family", "\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica Neue\", Helvetica, Arial, \"Lucida Grande\", sans-serif")
	  .style("font-size", "12px")
	  .style("font-weight","900")
	  .style("pointer-events", "none")
	  .style("stroke", "#FFF")
	  .style("stroke-width", "1px")
	  .text(function(d) {return d.name;});
// add the curvy lines
function tick() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
};

});
}

function drawPsuedoMinTree(jsonData, addTo){
	var width = 2000,
    height = 2000;

var force = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);
	d3.select(addTo)
    	.attr("width", width)
    	.attr("height", height);
var svg = d3.select(addTo).append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("id", "chart");

d3.json(jsonData, function(error, graph) {
	//console.log(jsonData);
	//console.log(graph);
  svg.append('rect')
	    .attr('width',width)
	    .attr('height', height)
	    .attr('fill',graph.backgroundColor);
  
  force
      .nodes(graph.nodes)
      .links(graph.links)
      .on("tick",tick)
      .start();

  var link = svg.selectAll(".link")
      .data(graph.links)
    .enter().append("line")
      .attr("class", "link")
      .style("stroke", function(d) { return d.color; })
      .style("stroke-width", function(d) { return 2; });
var node = svg.selectAll(".node")
    .data(force.nodes())
  .enter().append("g")
    .attr("class", "node")
    .call(force.drag);

// add the nodes
	node.append("circle")
      .attr("r", function(d) { return Math.sqrt(d.size/Math.PI); })
      .style("fill", function(d) { return d.color; })
      .style("stroke", "#fff")
      .style("stroke-width", "1.5px")
      .call(force.drag);

  node.append("title")
      .text(function(d) { return d.name; });
      
  node.append("text")
	  .attr("x", 12)
	  .attr("dy", ".35em")
	  .style("fill","#FFF")
	  .style("font-family", "\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica Neue\", Helvetica, Arial, \"Lucida Grande\", sans-serif")
	  .style("font-size", "12px")
	  .style("font-weight","900")
	  .style("pointer-events", "none")
	  .text(function(d) { 
	  	if(d.size == 10){
	  		return "";
	  	}else{
	  		return d.name;
	  	}
	   });
// add the curvy lines
function tick() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
};

});
}