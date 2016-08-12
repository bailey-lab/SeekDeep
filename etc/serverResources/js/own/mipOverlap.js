mipOverlapWorker = function() {
  var mipOverlapWorker = {},
      nodeWidth = 24,
      nodePadding = 10,
      size = [1, 1],
      nodes = [],
      links = [];

  mipOverlapWorker.nodeWidth = function(_) {
    if (!arguments.length) return nodeWidth;
    nodeWidth = +_;
    return mipOverlapWorker;
  };

  mipOverlapWorker.nodePadding = function(_) {
    if (!arguments.length) return nodePadding;
    nodePadding = +_;
    return mipOverlapWorker;
  };

  mipOverlapWorker.nodes = function(_) {
    if (!arguments.length) return nodes;
    nodes = _;
    return mipOverlapWorker;
  };

  mipOverlapWorker.links = function(_) {
    if (!arguments.length) return links;
    links = _;
    return mipOverlapWorker;
  };

  mipOverlapWorker.size = function(_) {
    if (!arguments.length) return size;
    size = _;
    return mipOverlapWorker;
  };

  mipOverlapWorker.layout = function() {
    computeNodeLinks();
    computeNodeValues();
    computeNodeBreadths();
    //computeNodeDepths();
    computeNodeDepthsCalc(32)
    computeNodeColorFrac();
    computeLinkDepths();
    return mipOverlapWorker;
  };
  
  mipOverlapWorker.layoutByFrac = function() {
    computeNodeLinks();
    computeNodeValues();
    computeNodeBreadths();
    computeNodeDepths();
    computeNodeDepthsCalc(32)
    computeNodeColorFrac();
    computeLinkDepths();
    return mipOverlapWorker;
  };

  mipOverlapWorker.relayout = function() {
    computeLinkDepths();
    return mipOverlapWorker;
  };

  mipOverlapWorker.link = function() {
    var curvature = .5;

    function link(d) {
      var x0 = d.source.x + d.source.dx,
          x1 = d.target.x,
          xi = d3.interpolateNumber(x0, x1),
          x2 = xi(curvature),
          x3 = xi(1 - curvature),
          y0 = d.source.y + d.sy + d.dy / 2,
          y1 = d.target.y + d.ty + d.dy / 2;
      return "M" + x0 + "," + y0
           + "C" + x2 + "," + y0
           + " " + x3 + "," + y1
           + " " + x1 + "," + y1;
    }

    link.curvature = function(_) {
      if (!arguments.length) return curvature;
      curvature = +_;
      return link;
    };

    return link;
  };

  // Populate the sourceLinks and targetLinks for each node.
  // Also, if the source and target are not objects, assume they are indices.
  function computeNodeLinks() {
    nodes.forEach(function(node) {
      node.sourceLinks = [];
      node.targetLinks = [];
    });
    links.forEach(function(link) {
      var source = link.source,
          target = link.target;
      if (typeof source === "number") source = link.source = nodes[link.source];
      if (typeof target === "number") target = link.target = nodes[link.target];
      source.sourceLinks.push(link);
      target.targetLinks.push(link);
    });
  }

  // Compute the value (size) of each node by summing the associated links.
  function computeNodeValues() {
    //console.log(size[0]) width
    //console.log(size[1]) height
    nodes.forEach(function(node) {
      node.value = size[1] *.90 * node.frac
    });
  }
  
  function computeNodeColorByMipNum() {
	var maxMipNum = d3.max(nodes, function(node) {
        return node.mipNum;
      });
    var minMipNum = d3.min(nodes, function(node) {
        return node.mipNum;
      });
    var color = d3.scale.category20().domain([minMipNum, maxMipNum]);
    nodes.forEach(function(node) {
      if(!node.hasOwnProperty("color")){
      	node.color = color(node.mipNum)
      }
    });
  }
  function computeNodeColorFrac() {
	var maxNum = d3.max(nodes, function(node) {
        return node.yOrder;
      });
    var minNum = d3.min(nodes, function(node) {
        return node.yOrder;
      });
    var color = d3.scale.category20().domain([minNum, maxNum]);
    nodes.forEach(function(node) {
      if(!node.hasOwnProperty("color")){
      	node.color = color(node.yOrder)
      }
    });
  }


  function computeNodeBreadths() {
    x = 0;
  	nodes.forEach(function(node) {
        node.x = node.mipNum;
        if(node.mipNum > x){
        	x = node.mipNum
        }
        node.dx = nodeWidth;
      });
    ++x
    scaleNodeBreadths((size[0] - nodeWidth) / (x - 1));
  }

  function scaleNodeBreadths(kx) {
    nodes.forEach(function(node) {
      node.x *= kx;
    });
  }
  
  
  function computeNodeDepthsCalc(iterations) {
    var nodesByBreadth = d3.nest()
        .key(function(d) { return d.x; })
        .sortKeys(d3.ascending)
        .entries(nodes)
        .map(function(d) { return d.values; });

    //
    initializeNodeDepth();
    resolveCollisions();
    for (var alpha = 1; iterations > 0; --iterations) {
      relaxRightToLeft(alpha *= .99);
      resolveCollisions();
      relaxLeftToRight(alpha);
      resolveCollisions();
    }
	function nodeCompare(a,b) {
	  if (a.frac < b.frac)
	    return -1;
	  if (a.frac > b.frac)
	    return 1;
	  return 0;
	}
	
	
    function initializeNodeDepth() {
      var ky = d3.min(nodesByBreadth, function(nodes) {
        return (size[1] - (nodes.length - 1) * nodePadding) / d3.sum(nodes, value);
      });
      nodesByBreadth.forEach(function(nodes) {
        nodes.sort(nodeCompare);
        nodes.forEach(function(node, i) {
          node.y = i;
          node.yOrder = i
          node.dy = node.value * ky;
        });
      });
      nodesByBreadth.forEach(function(nodes) {
	      nodes.forEach(function(node, i) {
	        node.sourceLinks.forEach(function(link){
	        	link.dy = Math.min(link.value * node.dy, link.target.dy);
	        });
	      });
	    });
    }
/*
    function initializeNodeDepth() {
      var ky = d3.min(nodesByBreadth, function(nodes) {
        return (size[1] - (nodes.length - 1) * nodePadding) / d3.sum(nodes, value);
      });

      nodesByBreadth.forEach(function(nodes) {
        nodes.forEach(function(node, i) {
          node.y = i;
          node.dy = node.value * ky;
        });
      });

      links.forEach(function(link) {
        link.dy = link.value * ky;
      });
    }*/

    function relaxLeftToRight(alpha) {
      nodesByBreadth.forEach(function(nodes, breadth) {
        nodes.forEach(function(node) {
          if (node.targetLinks.length) {
            var y = d3.sum(node.targetLinks, weightedSource) / d3.sum(node.targetLinks, value);
            node.y += (y - center(node)) * alpha;
          }
        });
      });

      function weightedSource(link) {
        return center(link.source) * link.value;
      }
    }

    function relaxRightToLeft(alpha) {
      nodesByBreadth.slice().reverse().forEach(function(nodes) {
        nodes.forEach(function(node) {
          if (node.sourceLinks.length) {
            var y = d3.sum(node.sourceLinks, weightedTarget) / d3.sum(node.sourceLinks, value);
            node.y += (y - center(node)) * alpha;
          }
        });
      });

      function weightedTarget(link) {
        return center(link.target) * link.value;
      }
    }

    function resolveCollisions() {
      nodesByBreadth.forEach(function(nodes) {
        var node,
            dy,
            y0 = 0,
            n = nodes.length,
            i;

        // Push any overlapping nodes down.
        nodes.sort(ascendingDepth);
        for (i = 0; i < n; ++i) {
          node = nodes[i];
          dy = y0 - node.y;
          if (dy > 0) node.y += dy;
          y0 = node.y + node.dy + nodePadding;
        }

        // If the bottommost node goes outside the bounds, push it back up.
        dy = y0 - nodePadding - size[1];
        if (dy > 0) {
          y0 = node.y -= dy;

          // Push any overlapping nodes back up.
          for (i = n - 2; i >= 0; --i) {
            node = nodes[i];
            dy = node.y + node.dy + nodePadding - y0;
            if (dy > 0) node.y -= dy;
            y0 = node.y;
          }
        }
      });
    }

    function ascendingDepth(a, b) {
      return a.y - b.y;
    }
  }

  function computeNodeDepths() {
    var nodesByBreadth = d3.nest()
        .key(function(d) { return d.x; })
        .sortKeys(d3.ascending)
        .entries(nodes)
        .map(function(d) { return d.values; });

    initializeNodeDepth();
    resolveCollisions();

	function nodeCompare(a,b) {
	  if (a.frac < b.frac)
	    return -1;
	  if (a.frac > b.frac)
	    return 1;
	  return 0;
	}
	
	
    function initializeNodeDepth() {
      var ky = d3.min(nodesByBreadth, function(nodes) {
        return (size[1] - (nodes.length - 1) * nodePadding) / d3.sum(nodes, value);
      });
      nodesByBreadth.forEach(function(nodes) {
        nodes.sort(nodeCompare);
        nodes.forEach(function(node, i) {
          node.y = i;
          node.yOrder = i
          node.dy = node.value * ky;
        });
      });
      nodesByBreadth.forEach(function(nodes) {
	      nodes.forEach(function(node, i) {
	        node.sourceLinks.forEach(function(link){
	        	link.dy = Math.min(link.value * node.dy, link.target.dy);
	        });
	      });
	    });
    }


    function resolveCollisions() {
      nodesByBreadth.forEach(function(nodes) {
        var node,
            dy,
            y0 = 0,
            n = nodes.length,
            i;

        // Push any overlapping nodes down.
        nodes.sort(ascendingDepth);
        for (i = 0; i < n; ++i) {
          node = nodes[i];
          dy = y0 - node.y;
          if (dy > 0) node.y += dy;
          y0 = node.y + node.dy + nodePadding;
        }

        // If the bottommost node goes outside the bounds, push it back up.
        dy = y0 - nodePadding - size[1];
        if (dy > 0) {
          y0 = node.y -= dy;

          // Push any overlapping nodes back up.
          for (i = n - 2; i >= 0; --i) {
            node = nodes[i];
            dy = node.y + node.dy + nodePadding - y0;
            if (dy > 0) node.y -= dy;
            y0 = node.y;
          }
        }
      });
    }

    function ascendingDepth(a, b) {
      return a.y - b.y;
    }
  }

  function computeLinkDepths() {
    nodes.forEach(function(node) {
      node.sourceLinks.sort(ascendingTargetDepth);
      node.targetLinks.sort(ascendingSourceDepth);
    });
    nodes.forEach(function(node) {
      var sy = 0, ty = 0;
      node.sourceLinks.forEach(function(link) {
        link.sy = sy;
        sy += link.dy;
      });
      node.targetLinks.forEach(function(link) {
        link.ty = ty;
        ty += link.dy;
      });
    });

    function ascendingSourceDepth(a, b) {
      return a.source.y - b.source.y;
    }

    function ascendingTargetDepth(a, b) {
      return a.target.y - b.target.y;
    }
  }

  function center(node) {
    return node.y + node.dy / 2;
  }

  function value(link) {
    return link.value;
  }

  return mipOverlapWorker;
};



function MipOverlapper(addTo, graph, width, height, margin){
	this.initialWindowHeight = window.innerHeight
	this.initialWindowWidth = window.innerWidth
	this.margin = margin || {top: 10, right: 40, bottom: 6, left: 60};
	this.originalWidth =  width || window.innerWidth
	this.originalHeight =  height ||  window.innerHeight/2
	this.width = this.originalWidth - this.margin.left - this. margin.right;
	this.height = this.originalHeight - this.margin.top - this.margin.bottom;
	this.topId = addTo
	this.graph = graph
	this.mOverlap = mipOverlapWorker()
	    .nodeWidth(60)
	    .nodePadding(30)
	    .size([this.width, this.height])
	    .nodes(this.graph.nodes)
	    .links(this.graph.links)
	    .layout();
	this.path = this.mOverlap.link();
	this.svg = d3.select(this.topId).append("svg")
		.attr("class", "mipOverlapWorker")
		.attr("width", this.width + this.margin.left + this.margin.right)
	    .attr("height", this.height + this.margin.top + this.margin.bottom)
	  .append("g")
	    .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");

	
	var self = this;
	
	var link = this.svg.append("g")
	.attr("class", "mogLinks")
	  .selectAll(".link")
	  .data(this.graph.links)
	.enter().append("path")
	  .attr("class", "link")
	  .attr("d", this.path)
	  .style("stroke-width", function(d) { return Math.max(1, d.dy); })
	  .sort(function(a, b) { return b.dy - a.dy; });
	 link.append("title")
	  	.text(function(d) { return d.source.name + " → " + d.target.name + "\n" + d.value; });
	
	 function dragmove(d) {
	    d3.select(this).attr("transform", "translate(" + d.x + "," + (d.y = Math.max(0, Math.min(self.height - d.dy, d3.event.y))) + ")");
	    self.mOverlap.relayout();
	    link.attr("d", self.path);
	 }
	var node = this.svg.append("g").attr("class", "mogNodes").selectAll(".node")
	  .data(this.graph.nodes)
	.enter().append("g")
	  .attr("class", "node")
	  .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
	.call(d3.behavior.drag()
	  .origin(function(d) { return d; })
	  .on("dragstart", function() { this.parentNode.appendChild(this); })
	  .on("drag", dragmove));
	
	  node.append("rect")
	  .attr("height", function(d) { return d.dy; })
	  .attr("width", this.mOverlap.nodeWidth())
	  .style("fill", function(d) { return d.color })
	  .style("stroke", function(d) { return d3.rgb(d.color).darker(2); })
	.append("title")
	  .text(function(d) { return d.name + "\n" + d.value; });
	

}

MipOverlapper.prototype.updateSize = function(){
	var self = this;
	this.width = this.originalWidth * ( window.innerWidth  / this.initialWindowWidth) - this.margin.left - this. margin.right;
	this.height = this.originalHeight * ( window.innerHeight / this.initialWindowHeight) - this.margin.top - this.margin.bottom;
	this.mOverlap = mipOverlapWorker()
	    .nodeWidth(60)
	    .nodePadding(30)
	    .size([this.width, this.height])
	    .nodes(this.graph.nodes)
	    .links(this.graph.links)
	    .layout();
	this.path = this.mOverlap.link();
	d3.select(this.topId).select(".mipOverlapWorker")
		.attr("width", this.width + this.margin.left + this.margin.right)
	    .attr("height", this.height + this.margin.top + this.margin.bottom)
	this.redraw()
}

MipOverlapper.prototype.redraw = function(){
	var self = this;
	//select all link class objects and attach current data
	var links = this.svg.select(".mogLinks").selectAll(".link")
	  .data(this.graph.links)
	//enter any new links
	links.enter().append("path")
	  .attr("class", "link")
	  
	//exit any links that no longer have data
	links
	.exit().remove()
	
	//set up for the links
	links
	  .attr("d", this.path)
	  .style("stroke-width", function(d) { return Math.max(1, d.dy); })
	  .sort(function(a, b) { return b.dy - a.dy; });
	 links.append("title")
	  	.text(function(d) { return d.source.name + " → " + d.target.name + "\n" + d.value; });
	  
	 function dragmove(d) {
	    d3.select(this).attr("transform", "translate(" + d.x + "," + (d.y = Math.max(0, Math.min(self.height - d.dy, d3.event.y))) + ")");
	    self.mOverlap.relayout();
	    links.attr("d", self.path);
	 }
	//select all node class objects and attach current data
	var nodes = this.svg.select(".mogNodes").selectAll(".node")
	  .data(this.graph.nodes)
	//enter any new nodes
	nodes
	  .enter().append("g")
	  .attr("class", "node")
	  .append("rect")
	  .append("title")
	  
	//exit any nodes that no longer have data
	nodes
	.exit().remove()
	
	nodes.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
	.call(d3.behavior.drag()
	  .origin(function(d) { return d; })
	  .on("dragstart", function() { this.parentNode.appendChild(this); })
	  .on("drag", dragmove));
	
	 nodes.select("rect")
	  .attr("height", function(d) { return d.dy; })
	  .attr("width", this.mOverlap.nodeWidth())
	  .style("fill", function(d) { return d.color })
	  .style("stroke", function(d) { return d3.rgb(d.color).darker(2); })
	.select("title")
	  .text(function(d) { return d.name + "\n" + d.value; });
	
}

