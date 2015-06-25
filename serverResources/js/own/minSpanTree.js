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