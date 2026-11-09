//colors from RColorBrewer Spectral 11, the default ramp of HaplotypeRainbows
var rainbowSpectral = ["#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"];
//colorPalette_12 from HaplotypeRainbows, used to color by rank
var rainbowRankPalette = ["#9F0162", "#009F81", "#FF5AAF", "#00FCCF", "#8400CD", "#008DF9", "#00C2F9", "#FFB2FD", "#A40122", "#E20134", "#FF6E3A", "#FFC33B"];
var rainbowMissingColor = "#EDEDED";
var rainbowMaxClusterSamples = 3000;

function rainbowSaveBlob(blob, fileName){
	var url = URL.createObjectURL(blob);
	var link = document.createElement("a");
	link.href = url;
	link.download = fileName;
	document.body.appendChild(link);
	link.click();
	document.body.removeChild(link);
	setTimeout(function(){ URL.revokeObjectURL(url); }, 1000);
}

/**
 * Records the canvas calls the rainbow drawing uses as SVG elements, so the same drawing code makes the SVG download.
 * Clipping is ignored since the whole plot is drawn.
 * @param width width of the SVG
 * @param height height of the SVG
 */
function RainbowSvgContext(width, height){
	this.width = width;
	this.height = height;
	this.parts = [];
	this.lastFill = null;
	this.fillStyle = "#000000";
	this.strokeStyle = "#000000";
	this.lineWidth = 1;
	this.globalAlpha = 1;
	this.font = "11px sans-serif";
	this.textAlign = "left";
	this.textBaseline = "alphabetic";
	this.transforms = [];
	this.stack = [];
}

RainbowSvgContext.num = function(value){
	return String(Math.round(value * 100) / 100);
};

RainbowSvgContext.escape = function(str){
	return String(str).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
};

RainbowSvgContext.prototype.transformAttr = function(){
	return this.transforms.length > 0 ? ' transform="' + this.transforms.join(" ") + '"' : "";
};

RainbowSvgContext.prototype.alphaAttr = function(name){
	return this.globalAlpha < 1 ? " " + name + '="' + RainbowSvgContext.num(this.globalAlpha) + '"' : "";
};

RainbowSvgContext.prototype.fillRect = function(x, y, width, height){
	if(width <= 0 || height <= 0){
		this.lastFill = null;
		return;
	}
	var num = RainbowSvgContext.num;
	this.lastFill = {key: [x, y, width, height].join(","), pos: this.parts.length};
	this.parts.push('<rect x="' + num(x) + '" y="' + num(y) + '" width="' + num(width) + '" height="' + num(height) + '" fill="' + this.fillStyle + '"'
		+ this.alphaAttr("fill-opacity") + this.transformAttr() + "/>");
};

RainbowSvgContext.prototype.strokeRect = function(x, y, width, height){
	var num = RainbowSvgContext.num;
	var stroke = ' stroke="' + this.strokeStyle + '" stroke-width="' + num(this.lineWidth) + '"' + this.alphaAttr("stroke-opacity");
	//an outline right after the fill of the same rect becomes one element
	if(this.lastFill && this.lastFill.key === [x, y, width, height].join(",") && this.lastFill.pos === this.parts.length - 1){
		this.parts[this.parts.length - 1] = this.parts[this.parts.length - 1].replace(/\/>$/, stroke + "/>");
	} else {
		this.parts.push('<rect x="' + num(x) + '" y="' + num(y) + '" width="' + num(width) + '" height="' + num(height) + '" fill="none"' + stroke + this.transformAttr() + "/>");
	}
	this.lastFill = null;
};

RainbowSvgContext.prototype.fillText = function(text, x, y){
	var anchor = {left: "start", start: "start", right: "end", end: "end", center: "middle"}[this.textAlign] || "start";
	//plain attributes and a dy shift instead of the font shorthand and dominant-baseline, which some editors don't read
	var dy = {middle: "0.35em", top: "0.71em", hanging: "0.71em"}[this.textBaseline];
	var fontMatch = this.font.match(/([\d.]+)px\s+(.+)/);
	var fontAttrs = fontMatch ? ' font-size="' + fontMatch[1] + '" font-family="' + fontMatch[2] + '"' : "";
	this.parts.push('<text x="' + RainbowSvgContext.num(x) + '" y="' + RainbowSvgContext.num(y) + '"' + (dy ? ' dy="' + dy + '"' : "") + ' fill="' + this.fillStyle + '"' + fontAttrs
		+ ' text-anchor="' + anchor + '"' + this.alphaAttr("fill-opacity") + this.transformAttr() + ">" + RainbowSvgContext.escape(text) + "</text>");
	this.lastFill = null;
};

RainbowSvgContext.prototype.save = function(){
	this.stack.push({fillStyle: this.fillStyle, strokeStyle: this.strokeStyle, lineWidth: this.lineWidth, globalAlpha: this.globalAlpha,
		font: this.font, textAlign: this.textAlign, textBaseline: this.textBaseline, transforms: this.transforms.slice()});
};

RainbowSvgContext.prototype.restore = function(){
	var state = this.stack.pop();
	if(state){
		for(var key in state){
			this[key] = state[key];
		}
	}
};

RainbowSvgContext.prototype.translate = function(x, y){
	this.transforms.push("translate(" + RainbowSvgContext.num(x) + "," + RainbowSvgContext.num(y) + ")");
};

RainbowSvgContext.prototype.rotate = function(radians){
	this.transforms.push("rotate(" + RainbowSvgContext.num(radians * 180 / Math.PI) + ")");
};

RainbowSvgContext.prototype.clearRect = function(){};
RainbowSvgContext.prototype.beginPath = function(){};
RainbowSvgContext.prototype.rect = function(){};
RainbowSvgContext.prototype.clip = function(){};
RainbowSvgContext.prototype.setTransform = function(){};

RainbowSvgContext.prototype.toString = function(){
	return '<svg xmlns="http://www.w3.org/2000/svg" width="' + this.width + '" height="' + this.height + '" viewBox="0 0 ' + this.width + " " + this.height + '">\n'
		+ this.parts.join("\n") + "\n</svg>\n";
};

/**
 * The value HaplotypeRainbows maps onto the color ramp, the hue steps across the targets with a period so neighboring targets differ
 * @param rank the population rank of the haplotype in the target, 0 is the most common
 * @param hapCount the number of haplotypes in the target
 * @param targetPos the position of the target in the plot
 * @param period how many targets before the hues repeat
 * @param logScale spread the less common haplotypes out more
 */
function rainbowColorValue(rank, hapCount, targetPos, period, logScale){
	var frac = rank / hapCount;
	var perc = logScale ? Math.log(frac * 99 + 1) / Math.log(100) * 100 : frac * 100;
	var hueMod = (targetPos % period) + 1;
	return Math.round(Math.abs((perc + hueMod / period * 100) % 200 - 0.0001) % 100);
}

/**
 * Order of n items from average linkage clustering, built with the nearest neighbor chain
 * @param n the number of items
 * @param distFunc distance between items i and j
 * @return an array of item positions, items in the same cluster next to each other
 */
function averageLinkageOrder(n, distFunc){
	if(n < 3){
		return _.range(n);
	}
	var dist = new Float32Array(n * n);
	for(var i = 0; i < n; ++i){
		for(var j = i + 1; j < n; ++j){
			var d = distFunc(i, j);
			dist[i * n + j] = d;
			dist[j * n + i] = d;
		}
	}
	var size = new Int32Array(n).fill(1);
	var active = new Uint8Array(n).fill(1);
	var nodeOf = new Int32Array(n);
	for(var i = 0; i < n; ++i){
		nodeOf[i] = i;
	}
	var lefts = [];
	var rights = [];
	var chain = [];
	var remaining = n;
	var firstActive = 0;
	while(remaining > 1){
		if(chain.length === 0){
			while(!active[firstActive]){
				++firstActive;
			}
			chain.push(firstActive);
		}
		var a = chain[chain.length - 1];
		var prev = chain.length > 1 ? chain[chain.length - 2] : -1;
		var best = -1;
		var bestDist = Infinity;
		for(var k = 0; k < n; ++k){
			if(!active[k] || k === a){
				continue;
			}
			var d = dist[a * n + k];
			//ties go to the previous link so the chain can't cycle
			if(d < bestDist || (d === bestDist && k === prev)){
				bestDist = d;
				best = k;
			}
		}
		if(best === prev){
			chain.pop();
			chain.pop();
			for(var k = 0; k < n; ++k){
				if(!active[k] || k === a || k === prev){
					continue;
				}
				var merged = (size[a] * dist[a * n + k] + size[prev] * dist[prev * n + k]) / (size[a] + size[prev]);
				dist[a * n + k] = merged;
				dist[k * n + a] = merged;
			}
			lefts.push(nodeOf[a]);
			rights.push(nodeOf[prev]);
			nodeOf[a] = n + lefts.length - 1;
			size[a] += size[prev];
			active[prev] = 0;
			--remaining;
		} else {
			chain.push(best);
		}
	}
	var order = [];
	var stack = [n + lefts.length - 1];
	while(stack.length > 0){
		var node = stack.pop();
		if(node < n){
			order.push(node);
		} else {
			stack.push(rights[node - n]);
			stack.push(lefts[node - n]);
		}
	}
	return order;
}

/**
 * Interactive haplotype rainbow, columns are targets, rows are samples and each cell is split by the fraction of each haplotype in the sample
 * @param rName the root name of the viewer
 * @param data from /haplotypeRainbow
 * @param metaGroups sub group names by meta field, from /metaGroups
 * @param groupName the meta field the samples were taken from, empty for all samples
 * @param subGroupName the value of that field
 */
function HaplotypeRainbow(rName, data, metaGroups, groupName, subGroupName){
	var self = this;
	self.rName = rName;
	self.samples = data["samples"];
	self.targets = data["targets"];
	self.metaFields = data["metaFields"];
	self.sampleMeta = data["sampleMeta"];
	self.groupName = groupName;
	self.subGroupName = subGroupName;
	var sampleCount = self.samples.length;
	var targetCount = self.targets.length;
	self.sampleCount = sampleCount;
	self.targetCount = targetCount;

	//index the calls by cell, sample * targetCount + target, keeping the server order of most reads first
	var calls = data["calls"];
	var callCount = calls.length / 4;
	self.cellStart = new Int32Array(sampleCount * targetCount + 1);
	for(var c = 0; c < callCount; ++c){
		++self.cellStart[calls[4 * c] * targetCount + calls[4 * c + 1] + 1];
	}
	for(var cell = 1; cell <= sampleCount * targetCount; ++cell){
		self.cellStart[cell] += self.cellStart[cell - 1];
	}
	self.callRank = new Int32Array(callCount);
	self.callReads = new Float64Array(callCount);
	self.cellReads = new Float64Array(sampleCount * targetCount);
	var fillPos = self.cellStart.slice(0, sampleCount * targetCount);
	for(var c = 0; c < callCount; ++c){
		var cellPos = calls[4 * c] * targetCount + calls[4 * c + 1];
		var pos = fillPos[cellPos]++;
		self.callRank[pos] = calls[4 * c + 2];
		self.callReads[pos] = calls[4 * c + 3];
		self.cellReads[cellPos] += calls[4 * c + 3];
	}
	//the same calls ordered by rank within each cell
	self.rankOrder = new Int32Array(callCount);
	for(var cell = 0; cell < sampleCount * targetCount; ++cell){
		var start = self.cellStart[cell];
		var end = self.cellStart[cell + 1];
		var positions = [];
		for(var pos = start; pos < end; ++pos){
			positions.push(pos);
		}
		positions.sort(function(first, second){ return self.callRank[first] - self.callRank[second]; });
		for(var pos = start; pos < end; ++pos){
			self.rankOrder[pos] = positions[pos - start];
		}
	}
	self.targetsTyped = new Int32Array(sampleCount);
	for(var s = 0; s < sampleCount; ++s){
		for(var t = 0; t < targetCount; ++t){
			if(self.cellStart[s * targetCount + t + 1] > self.cellStart[s * targetCount + t]){
				++self.targetsTyped[s];
			}
		}
	}

	self.opts = {
		sampleSort: "alphabetical",
		targetSort: "alphabetical",
		cellOrder: "rank",
		layout: "stacked",
		colorMode: "rainbowLog",
		period: 11,
		rowHeight: Math.max(1, Math.min(14, Math.floor(1400 / Math.max(1, sampleCount)))),
		colWidth: 8,
		showSampleLabels: true,
		showTargetLabels: true
	};
	self.sampleOrder = _.range(sampleCount);
	self.targetOrder = _.range(targetCount);
	self.hover = null;
	self.drawPending = false;

	self.measureLabels();
	self.opts.colWidth = Math.max(3, Math.min(24, Math.floor((d34.select("#rainbowPlot").node().clientWidth - self.leftMargin() - 20) / Math.max(1, targetCount))));
	self.addSummary(data);
	self.addControls(metaGroups);
	self.setUpPlot();
	self.updateColors();
	self.updateLayout();
}

HaplotypeRainbow.prototype.addSummary = function(data){
	var self = this;
	var summary = d34.select("#rainbowSummary");
	var callCount = data["calls"].length / 4;
	summary.append("div").text(self.sampleCount + " samples with haplotypes, " + self.targetCount + " targets, " + callCount + " haplotype calls");
	var withoutData = data["samplesWithoutData"];
	if(withoutData.length > 0){
		summary.append("div").text(withoutData.length + " samples have no haplotypes in any target and aren't shown: " + withoutData.slice(0, 50).join(", ") + (withoutData.length > 50 ? ", ..." : ""));
	}
	summary.append("div")
		.style("color", "#555")
		.text("Each column is a target and each row a sample. Cells are split by the fraction of each haplotype in the sample, colors only mean something within a column: "
			+ "the same color in a column is the same haplotype, ranked by how many of these samples have it. Grey cells have no haplotypes. "
			+ "Click a cell for the sample's page for that target, a target name for the target, or a sample name for the sample across all targets.");
};

HaplotypeRainbow.prototype.addControls = function(metaGroups){
	var self = this;
	var controls = d34.select("#rainbowControls");
	var row = function(){
		return controls.append("div").style("margin-bottom", "8px");
	};
	var addSelect = function(parent, label, id, options, value, onChange){
		parent.append("label").attr("for", id).style("margin", "0 6px 0 12px").text(label);
		var sel = parent.append("select").attr("id", id);
		sel.selectAll("option").data(options).enter().append("option")
			.attr("value", function(d){ return d.value; })
			.property("disabled", function(d){ return !!d.disabled; })
			.text(function(d){ return d.text; });
		sel.property("value", value);
		sel.on("change", function(){ onChange(this.value); });
		return sel;
	};
	var addRange = function(parent, label, id, min, max, value, onChange){
		parent.append("label").attr("for", id).style("margin", "0 6px 0 12px").text(label);
		var input = parent.append("input").attr("type", "range").attr("id", id)
			.attr("min", min).attr("max", max).attr("value", value)
			.style("display", "inline-block").style("width", "140px").style("vertical-align", "middle");
		var shown = parent.append("span").style("margin-left", "4px");
		shown.text(value + "px");
		input.on("input", function(){
			shown.text(this.value + "px");
			onChange(parseInt(this.value));
		});
		return input;
	};
	var addCheck = function(parent, label, id, value, onChange){
		parent.append("input").attr("type", "checkbox").attr("id", id).style("margin-left", "12px").property("checked", value)
			.on("change", function(){ onChange(this.checked); });
		parent.append("label").attr("for", id).style("margin-left", "4px").text(label);
	};

	//which samples
	var groupOptions = [{value: "", text: "All samples"}];
	Object.keys(metaGroups).sort().forEach(function(group){
		metaGroups[group].forEach(function(subGroup){
			groupOptions.push({value: encodeURIComponent(group) + "/" + encodeURIComponent(subGroup), text: group + ": " + subGroup});
		});
	});
	var current = self.groupName ? encodeURIComponent(self.groupName) + "/" + encodeURIComponent(self.subGroupName) : "";
	addSelect(row(), "Samples", "rainbowGroup", groupOptions, current, function(value){
		window.location = "/" + self.rName + "/haplotypeRainbowPage" + (value ? "/" + value : "");
	});

	var orderRow = row();
	var tooManyToCluster = self.sampleCount > rainbowMaxClusterSamples;
	var clusterNote = tooManyToCluster ? " (" + rainbowMaxClusterSamples + " samples or fewer)" : "";
	var sampleSortOptions = [
		{value: "alphabetical", text: "name"},
		{value: "targetsTyped", text: "targets with haplotypes"},
		{value: "cluster:frequency", text: "clustering on all haplotypes, weighted by within-sample frequency" + clusterNote, disabled: tooManyToCluster},
		{value: "cluster:presence", text: "clustering on all haplotypes, presence only" + clusterNote, disabled: tooManyToCluster},
		{value: "cluster:major", text: "clustering on major haplotypes" + clusterNote, disabled: tooManyToCluster}
	];
	self.metaFields.forEach(function(field, fieldPos){
		sampleSortOptions.push({value: "meta:" + fieldPos, text: field});
	});
	addSelect(orderRow, "Order samples by", "rainbowSampleSort", sampleSortOptions, self.opts.sampleSort, function(value){
		self.opts.sampleSort = value;
		self.updateSampleOrder();
	});
	addSelect(orderRow, "Order targets by", "rainbowTargetSort", [
		{value: "alphabetical", text: "name"},
		{value: "haplotypes", text: "most haplotypes"},
		{value: "samples", text: "most samples"}
	], self.opts.targetSort, function(value){
		self.opts.targetSort = value;
		self.updateTargetOrder();
	});
	addSelect(orderRow, "Haplotypes in a cell by", "rainbowCellOrder", [
		{value: "rank", text: "population rank"},
		{value: "frac", text: "fraction in sample"}
	], self.opts.cellOrder, function(value){
		self.opts.cellOrder = value;
		self.requestDraw();
	});

	var lookRow = row();
	addSelect(lookRow, "Colors", "rainbowColorMode", [
		{value: "rainbowLog", text: "rainbow, log scaled rank"},
		{value: "rainbow", text: "rainbow"},
		{value: "rank", text: "by rank, same in every target"}
	], self.opts.colorMode, function(value){
		self.opts.colorMode = value;
		self.updateColors();
		self.requestDraw();
	});
	addSelect(lookRow, "Split cells", "rainbowLayout", [
		{value: "stacked", text: "stacked"},
		{value: "sideBySide", text: "side by side"}
	], self.opts.layout, function(value){
		self.opts.layout = value;
		self.requestDraw();
	});

	var sizeRow = row();
	addRange(sizeRow, "Row height", "rainbowRowHeight", 1, 30, self.opts.rowHeight, function(value){
		self.opts.rowHeight = value;
		self.updateLayout();
	});
	addRange(sizeRow, "Column width", "rainbowColWidth", 1, 40, self.opts.colWidth, function(value){
		self.opts.colWidth = value;
		self.updateLayout();
	});
	addCheck(sizeRow, "sample names", "rainbowSampleLabels", self.opts.showSampleLabels, function(value){
		self.opts.showSampleLabels = value;
		self.updateLayout();
	});
	addCheck(sizeRow, "target names", "rainbowTargetLabels", self.opts.showTargetLabels, function(value){
		self.opts.showTargetLabels = value;
		self.updateLayout();
	});

	var downloadRow = row();
	downloadRow.append("button").attr("class", "btn btn-default btn-sm").style("margin-left", "12px").text("Download PNG")
		.on("click", function(){ self.downloadPng(); });
	downloadRow.append("button").attr("class", "btn btn-default btn-sm").style("margin-left", "8px").text("Download SVG")
		.on("click", function(){ self.downloadSvg(); });
	downloadRow.append("button").attr("class", "btn btn-default btn-sm").style("margin-left", "8px").text("Download table")
		.on("click", function(){ self.downloadTable(); });
};

HaplotypeRainbow.prototype.fileStub = function(){
	var self = this;
	var stub = "haplotypeRainbow";
	if(self.groupName){
		stub += "_" + self.groupName + "_" + self.subGroupName;
	}
	return stub.replace(/[^A-Za-z0-9._-]/g, "_");
};

HaplotypeRainbow.prototype.setUpPlot = function(){
	var self = this;
	var plot = d34.select("#rainbowPlot");
	self.scroll = plot.append("div")
		.style("position", "relative")
		.style("overflow", "auto")
		.style("border", "1px solid #ccc")
		.style("max-height", "80vh");
	self.spacer = self.scroll.append("div").style("position", "relative");
	self.canvas = self.spacer.append("canvas")
		.style("position", "absolute")
		.style("left", "0px")
		.style("top", "0px");
	self.ctx = self.canvas.node().getContext("2d");
	self.tooltip = d34.select("body").append("div")
		.style("position", "fixed")
		.style("pointer-events", "none")
		.style("background", "#fff")
		.style("border", "1px solid #888")
		.style("border-radius", "4px")
		.style("padding", "6px 8px")
		.style("font-size", "12px")
		.style("max-width", "420px")
		.style("box-shadow", "2px 2px 6px rgba(0,0,0,0.25)")
		.style("z-index", "2000")
		.style("display", "none");
	self.metaLegend = d34.select("#rainbowLegend").style("margin", "6px 0");

	//what's under the mouse changes on scroll, so drop the hover until the mouse moves again
	self.scroll.on("scroll", function(){
		self.hover = null;
		self.tooltip.style("display", "none");
		self.requestDraw();
	});
	window.addEventListener("resize", function(){ self.updateLayout(); });
	self.canvas.on("mousemove", function(){
		var point = d34.mouse(this);
		self.onHover(point[0], point[1], d34.event);
	});
	self.canvas.on("mouseleave", function(){
		self.hover = null;
		self.tooltip.style("display", "none");
		self.requestDraw();
	});
	self.canvas.on("click", function(){
		var point = d34.mouse(this);
		self.onClick(point[0], point[1], d34.event);
	});
};

HaplotypeRainbow.prototype.measureLabels = function(){
	var self = this;
	var measure = document.createElement("canvas").getContext("2d");
	measure.font = "11px sans-serif";
	self.maxSampleLabel = 0;
	self.samples.forEach(function(sample){
		self.maxSampleLabel = Math.max(self.maxSampleLabel, measure.measureText(sample).width);
	});
	self.maxTargetLabel = 0;
	self.targets.forEach(function(target){
		self.maxTargetLabel = Math.max(self.maxTargetLabel, measure.measureText(target["name"]).width);
	});
};

HaplotypeRainbow.prototype.metaBandFieldPos = function(){
	var sort = this.opts.sampleSort;
	return 0 === sort.indexOf("meta:") ? parseInt(sort.substring(5)) : -1;
};

HaplotypeRainbow.prototype.sampleLabelsShown = function(){
	return this.opts.showSampleLabels && this.opts.rowHeight >= 8;
};

HaplotypeRainbow.prototype.targetLabelsShown = function(){
	return this.opts.showTargetLabels && this.opts.colWidth >= 8;
};

HaplotypeRainbow.prototype.metaBandWidth = function(){
	return this.metaBandFieldPos() >= 0 ? 14 : 0;
};

//room on the left for the sample names and meta band, kept even without names so rows can be hovered
HaplotypeRainbow.prototype.leftMargin = function(){
	return (this.sampleLabelsShown() ? Math.min(240, Math.ceil(this.maxSampleLabel)) + 8 : 12) + this.metaBandWidth();
};

//room on top for the target names, kept even without names so columns can be hovered
HaplotypeRainbow.prototype.topMargin = function(){
	return this.targetLabelsShown() ? Math.min(180, Math.ceil(this.maxTargetLabel)) + 8 : 14;
};

HaplotypeRainbow.prototype.updateLayout = function(){
	var self = this;
	var width = self.leftMargin() + self.targetCount * self.opts.colWidth;
	var height = self.topMargin() + self.sampleCount * self.opts.rowHeight;
	self.fullWidth = width;
	self.fullHeight = height;
	self.spacer.style("width", width + "px").style("height", height + "px");
	self.updateMetaLegend();
	self.requestDraw();
};

HaplotypeRainbow.prototype.updateColors = function(){
	var self = this;
	var colorMode = self.opts.colorMode;
	self.colors = new Array(self.targetCount);
	if("rank" === colorMode){
		self.targets.forEach(function(target, targetIdx){
			self.colors[targetIdx] = target["haps"].map(function(hap, rank){
				return rank < rainbowRankPalette.length ? rainbowRankPalette[rank] : "#AAAAAA";
			});
		});
		return;
	}
	//the ramp spans the values in use, like ggplot2's scale_fill_gradientn
	var values = new Array(self.targetCount);
	var minValue = Infinity;
	var maxValue = -Infinity;
	self.targetOrder.forEach(function(targetIdx, targetPos){
		var hapCount = self.targets[targetIdx]["haps"].length;
		values[targetIdx] = self.targets[targetIdx]["haps"].map(function(hap, rank){
			var value = rainbowColorValue(rank, hapCount, targetPos, self.opts.period, "rainbowLog" === colorMode);
			minValue = Math.min(minValue, value);
			maxValue = Math.max(maxValue, value);
			return value;
		});
	});
	if(minValue === maxValue){
		maxValue = minValue + 1;
	}
	var stops = rainbowSpectral.map(function(color, pos){
		return minValue + (maxValue - minValue) * pos / (rainbowSpectral.length - 1);
	});
	var scale = d34.scaleLinear().domain(stops).range(rainbowSpectral).interpolate(d34.interpolateLab);
	self.targets.forEach(function(target, targetIdx){
		self.colors[targetIdx] = values[targetIdx].map(function(value){ return scale(value); });
	});
};

HaplotypeRainbow.prototype.updateMetaLegend = function(){
	var self = this;
	self.metaLegend.selectAll("*").remove();
	var fieldPos = self.metaBandFieldPos();
	if(fieldPos < 0){
		self.metaColors = null;
		return;
	}
	var values = _.uniq(self.sampleMeta.map(function(meta){ return meta[fieldPos]; })).sort();
	var palette = values.length <= 10 ? d34.schemeCategory10 : d34.schemeCategory20;
	self.metaColors = {};
	values.forEach(function(value, pos){
		self.metaColors[value] = "NA" === value ? "#CCCCCC" : palette[pos % palette.length];
	});
	self.metaLegend.append("span").style("margin-right", "8px").style("font-weight", "bold").text(self.metaFields[fieldPos] + ":");
	var shownValues = values.slice(0, 60);
	var items = self.metaLegend.selectAll(".rainbowMetaItem").data(shownValues).enter().append("span")
		.attr("class", "rainbowMetaItem")
		.style("margin-right", "12px")
		.style("white-space", "nowrap");
	items.append("span")
		.style("display", "inline-block").style("width", "12px").style("height", "12px").style("vertical-align", "middle").style("margin-right", "4px")
		.style("background", function(d){ return self.metaColors[d]; });
	items.append("span").text(function(d){ return d; });
	if(values.length > shownValues.length){
		self.metaLegend.append("span").text("and " + (values.length - shownValues.length) + " more");
	}
};

HaplotypeRainbow.prototype.updateSampleOrder = function(){
	var self = this;
	var sort = self.opts.sampleSort;
	var byName = _.range(self.sampleCount);
	if("alphabetical" === sort){
		self.sampleOrder = byName;
	} else if("targetsTyped" === sort){
		self.sampleOrder = byName.slice().sort(function(first, second){
			return self.targetsTyped[second] - self.targetsTyped[first] || first - second;
		});
	} else if(0 === sort.indexOf("cluster:")){
		var gif = prsentDivGifLoading();
		//let the loading screen show before the clustering blocks
		setTimeout(function(){
			self.sampleOrder = self.clusterSampleOrder(sort.substring(8));
			removeAllDivGifLoading();
			self.updateLayout();
		}, 30);
		return;
	} else {
		var fieldPos = self.metaBandFieldPos();
		self.sampleOrder = byName.slice().sort(function(first, second){
			var firstValue = self.sampleMeta[first][fieldPos];
			var secondValue = self.sampleMeta[second][fieldPos];
			return firstValue < secondValue ? -1 : firstValue > secondValue ? 1 : first - second;
		});
	}
	self.updateLayout();
};

/**
 * Sample order from average linkage clustering, the distance between two samples is the mean of a per target distance over the targets both have haplotypes in,
 * 1 when they share none
 * @param mode how targets are compared:
 *   major: 0 if the haplotypes with the most reads match, otherwise 1
 *   frequency: 1 - the sum over haplotypes of the smaller within-sample fraction of the two samples (Bray-Curtis)
 *   presence: 1 - shared haplotypes / haplotypes in either sample (Jaccard)
 */
HaplotypeRainbow.prototype.clusterSampleOrder = function(mode){
	var self = this;
	var targetCount = self.targetCount;
	var cellStart = self.cellStart;
	var callRank = self.callRank;
	//within-sample fractions, split evenly when a cell has no reads as it's drawn
	var callFrac = new Float64Array(callRank.length);
	for(var cell = 0; cell < self.sampleCount * targetCount; ++cell){
		var count = cellStart[cell + 1] - cellStart[cell];
		for(var pos = cellStart[cell]; pos < cellStart[cell + 1]; ++pos){
			callFrac[pos] = self.cellReads[cell] > 0 ? self.callReads[pos] / self.cellReads[cell] : 1 / count;
		}
	}
	var major = "major" === mode;
	var presence = "presence" === mode;
	return averageLinkageOrder(self.sampleCount, function(first, second){
		var shared = 0;
		var total = 0;
		var firstOff = first * targetCount;
		var secondOff = second * targetCount;
		for(var t = 0; t < targetCount; ++t){
			var firstStart = cellStart[firstOff + t];
			var firstEnd = cellStart[firstOff + t + 1];
			var secondStart = cellStart[secondOff + t];
			var secondEnd = cellStart[secondOff + t + 1];
			if(firstEnd === firstStart || secondEnd === secondStart){
				continue;
			}
			++shared;
			//calls in a cell are kept with the most reads first, and most cells have one haplotype
			if(major || (1 === firstEnd - firstStart && 1 === secondEnd - secondStart)){
				if(callRank[firstStart] !== callRank[secondStart]){
					total += 1;
				}
				continue;
			}
			var matched = 0;
			for(var i = firstStart; i < firstEnd; ++i){
				for(var j = secondStart; j < secondEnd; ++j){
					if(callRank[i] === callRank[j]){
						matched += presence ? 1 : Math.min(callFrac[i], callFrac[j]);
						break;
					}
				}
			}
			total += presence ? 1 - matched / (firstEnd - firstStart + secondEnd - secondStart - matched) : Math.max(0, 1 - matched);
		}
		return 0 === shared ? 1 : total / shared;
	});
};

HaplotypeRainbow.prototype.updateTargetOrder = function(){
	var self = this;
	var sort = self.opts.targetSort;
	var byName = _.range(self.targetCount);
	if("haplotypes" === sort){
		self.targetOrder = byName.slice().sort(function(first, second){
			return self.targets[second]["haps"].length - self.targets[first]["haps"].length || first - second;
		});
	} else if("samples" === sort){
		self.targetOrder = byName.slice().sort(function(first, second){
			return self.targets[second]["samples"] - self.targets[first]["samples"] || first - second;
		});
	} else {
		self.targetOrder = byName;
	}
	//the hues step along the plotted order
	self.updateColors();
	self.updateLayout();
};

HaplotypeRainbow.prototype.requestDraw = function(){
	var self = this;
	if(self.drawPending){
		return;
	}
	self.drawPending = true;
	window.requestAnimationFrame(function(){
		self.drawPending = false;
		self.drawViewport();
	});
};

HaplotypeRainbow.prototype.drawViewport = function(){
	var self = this;
	var scrollNode = self.scroll.node();
	var viewWidth = Math.min(scrollNode.clientWidth, self.fullWidth);
	var viewHeight = Math.min(scrollNode.clientHeight, self.fullHeight);
	var scrollLeft = scrollNode.scrollLeft;
	var scrollTop = scrollNode.scrollTop;
	var ratio = window.devicePixelRatio || 1;
	var canvasNode = self.canvas.node();
	if(canvasNode.width !== Math.round(viewWidth * ratio) || canvasNode.height !== Math.round(viewHeight * ratio)){
		canvasNode.width = Math.round(viewWidth * ratio);
		canvasNode.height = Math.round(viewHeight * ratio);
		self.canvas.style("width", viewWidth + "px").style("height", viewHeight + "px");
	}
	//the canvas follows the scroll so only what's in view is drawn
	self.canvas.style("left", scrollLeft + "px").style("top", scrollTop + "px");
	self.ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
	self.draw(self.ctx, scrollLeft, scrollTop, viewWidth, viewHeight, true);
};

/**
 * Draw the plot as seen from a scroll position, the target and sample names stay at the top and left
 */
HaplotypeRainbow.prototype.draw = function(ctx, scrollLeft, scrollTop, viewWidth, viewHeight, withHover){
	var self = this;
	var opts = self.opts;
	var rowHeight = opts.rowHeight;
	var colWidth = opts.colWidth;
	var left = self.leftMargin();
	var top = self.topMargin();
	var targetCount = self.targetCount;
	ctx.clearRect(0, 0, viewWidth, viewHeight);
	ctx.fillStyle = "#FFFFFF";
	ctx.fillRect(0, 0, viewWidth, viewHeight);

	var firstCol = Math.max(0, Math.floor(scrollLeft / colWidth));
	var lastCol = Math.min(targetCount - 1, Math.floor((scrollLeft + viewWidth - left) / colWidth));
	var firstRow = Math.max(0, Math.floor(scrollTop / rowHeight));
	var lastRow = Math.min(self.sampleCount - 1, Math.floor((scrollTop + viewHeight - top) / rowHeight));
	var colX = function(col){ return left + col * colWidth - scrollLeft; };
	var rowY = function(row){ return top + row * rowHeight - scrollTop; };

	//leave a gap between rows and outline the haplotypes once they're big enough to see
	var barHeight = rowHeight >= 6 ? rowHeight * 0.8 : rowHeight;
	var barOffset = (rowHeight - barHeight) / 2;
	//black outlines once bars are 3 px, thinner on small bars so the colors still show
	var smallestSide = Math.min(colWidth, barHeight);
	var outline = smallestSide >= 3;
	var order = "rank" === opts.cellOrder ? self.rankOrder : null;
	var stacked = "stacked" === opts.layout;
	//while a haplotype is hovered the rest of its column fades so the samples sharing it stand out
	var hoverCol = withHover && self.hover && self.hover.rank >= 0 ? self.hover.col : -1;
	var hoverRank = hoverCol >= 0 ? self.hover.rank : -1;
	ctx.lineWidth = smallestSide >= 6 ? 1 : 0.5;
	ctx.strokeStyle = "#000000";

	ctx.save();
	ctx.beginPath();
	ctx.rect(left, top, viewWidth - left, viewHeight - top);
	ctx.clip();
	for(var row = firstRow; row <= lastRow; ++row){
		var sampleIdx = self.sampleOrder[row];
		var y = rowY(row) + barOffset;
		for(var col = firstCol; col <= lastCol; ++col){
			var targetIdx = self.targetOrder[col];
			var cell = sampleIdx * targetCount + targetIdx;
			var start = self.cellStart[cell];
			var end = self.cellStart[cell + 1];
			var x = colX(col);
			if(start === end){
				ctx.fillStyle = rainbowMissingColor;
				ctx.fillRect(x, y, colWidth, barHeight);
				continue;
			}
			var total = self.cellReads[cell];
			var targetColors = self.colors[targetIdx];
			var offset = 0;
			var fadeOthers = col === hoverCol;
			for(var pos = start; pos < end; ++pos){
				var callPos = order ? order[pos] : pos;
				var frac = 0 === total ? 1 / (end - start) : self.callReads[callPos] / total;
				ctx.globalAlpha = fadeOthers && self.callRank[callPos] !== hoverRank ? 0.25 : 1.0;
				ctx.fillStyle = targetColors[self.callRank[callPos]];
				if(stacked){
					ctx.fillRect(x, y + offset * barHeight, colWidth, frac * barHeight);
					if(outline){
						ctx.strokeRect(x, y + offset * barHeight, colWidth, frac * barHeight);
					}
				} else {
					ctx.fillRect(x + offset * colWidth, y, frac * colWidth, barHeight);
					if(outline){
						ctx.strokeRect(x + offset * colWidth, y, frac * colWidth, barHeight);
					}
				}
				offset += frac;
			}
			ctx.globalAlpha = 1.0;
		}
	}
	if(withHover && self.hover){
		ctx.lineWidth = 2;
		ctx.strokeStyle = "rgba(0, 0, 0, 0.85)";
		if(self.hover.col >= 0){
			ctx.strokeRect(colX(self.hover.col), top - 1, colWidth, self.sampleCount * rowHeight + 2);
		}
		if(self.hover.row >= 0){
			ctx.strokeRect(left - 1, rowY(self.hover.row), targetCount * colWidth + 2, rowHeight);
		}
	}
	ctx.restore();

	//target names along the top
	ctx.fillStyle = "#FFFFFF";
	ctx.fillRect(0, 0, viewWidth, top);
	ctx.font = "11px sans-serif";
	if(self.targetLabelsShown()){
		ctx.save();
		ctx.beginPath();
		ctx.rect(left, 0, viewWidth - left, top);
		ctx.clip();
		ctx.textAlign = "left";
		ctx.textBaseline = "middle";
		for(var col = firstCol; col <= lastCol; ++col){
			var hovered = withHover && self.hover && self.hover.col === col;
			ctx.fillStyle = hovered ? "#C00000" : "#000000";
			ctx.save();
			ctx.translate(colX(col) + colWidth / 2, top - 4);
			ctx.rotate(-Math.PI / 2);
			ctx.fillText(self.targets[self.targetOrder[col]]["name"], 0, 0);
			ctx.restore();
		}
		ctx.restore();
	} else if(withHover && self.hover && self.hover.col >= 0){
		ctx.fillStyle = "#C00000";
		ctx.fillRect(colX(self.hover.col), 2, colWidth, top - 4);
	}

	//sample names and the meta band down the left
	ctx.fillStyle = "#FFFFFF";
	ctx.fillRect(0, top, left, viewHeight - top);
	ctx.save();
	ctx.beginPath();
	ctx.rect(0, top, left, viewHeight - top);
	ctx.clip();
	var bandWidth = self.metaBandWidth();
	var fieldPos = self.metaBandFieldPos();
	ctx.textAlign = "right";
	ctx.textBaseline = "middle";
	for(var row = firstRow; row <= lastRow; ++row){
		var sampleIdx = self.sampleOrder[row];
		var hovered = withHover && self.hover && self.hover.row === row;
		if(bandWidth > 0 && self.metaColors){
			ctx.fillStyle = self.metaColors[self.sampleMeta[sampleIdx][fieldPos]];
			ctx.fillRect(left - bandWidth - 2, rowY(row), bandWidth, rowHeight);
		}
		if(self.sampleLabelsShown()){
			ctx.fillStyle = hovered ? "#C00000" : "#000000";
			ctx.fillText(self.samples[sampleIdx], left - bandWidth - 6, rowY(row) + rowHeight / 2);
		} else if(hovered){
			ctx.fillStyle = "#C00000";
			ctx.fillRect(2, rowY(row), left - bandWidth - 6, rowHeight);
		}
	}
	ctx.restore();
	ctx.fillStyle = "#FFFFFF";
	ctx.fillRect(0, 0, left, top);
};

//what's under a point of the canvas, row or col is -1 when the point is in the names along the top or left
HaplotypeRainbow.prototype.locate = function(x, y){
	var self = this;
	var scrollNode = self.scroll.node();
	var left = self.leftMargin();
	var top = self.topMargin();
	var col = x >= left ? Math.floor((x - left + scrollNode.scrollLeft) / self.opts.colWidth) : -1;
	var row = y >= top ? Math.floor((y - top + scrollNode.scrollTop) / self.opts.rowHeight) : -1;
	if(col >= self.targetCount || row >= self.sampleCount || (col < 0 && row < 0)){
		return null;
	}
	var ret = {col: col, row: row, rank: -1};
	if(col >= 0 && row >= 0){
		//how far into the cell the point is along the direction the haplotypes are split
		var opts = self.opts;
		var barHeight = opts.rowHeight >= 6 ? opts.rowHeight * 0.8 : opts.rowHeight;
		var into = "stacked" === opts.layout
			? (y - top + scrollNode.scrollTop - row * opts.rowHeight - (opts.rowHeight - barHeight) / 2) / barHeight
			: (x - left + scrollNode.scrollLeft - col * opts.colWidth) / opts.colWidth;
		into = Math.min(Math.max(into, 0), 0.999999);
		var cellCalls = self.cellCalls(self.sampleOrder[row], self.targetOrder[col]);
		var offset = 0;
		for(var pos = 0; pos < cellCalls.length; ++pos){
			offset += cellCalls[pos].frac;
			if(into < offset || pos === cellCalls.length - 1){
				ret.rank = cellCalls[pos].rank;
				break;
			}
		}
	}
	return ret;
};

HaplotypeRainbow.prototype.cellCalls = function(sampleIdx, targetIdx){
	var self = this;
	var cell = sampleIdx * self.targetCount + targetIdx;
	var ret = [];
	var order = "rank" === self.opts.cellOrder ? self.rankOrder : null;
	var count = self.cellStart[cell + 1] - self.cellStart[cell];
	for(var pos = self.cellStart[cell]; pos < self.cellStart[cell + 1]; ++pos){
		var callPos = order ? order[pos] : pos;
		//split evenly when there are no reads, as the cells are drawn
		ret.push({rank: self.callRank[callPos], reads: self.callReads[callPos], frac: self.cellReads[cell] > 0 ? self.callReads[callPos] / self.cellReads[cell] : 1 / count});
	}
	return ret;
};

HaplotypeRainbow.prototype.onHover = function(x, y, event){
	var self = this;
	var spot = self.locate(x, y);
	var spotKey = function(hover){ return hover ? hover.col + ":" + hover.row + ":" + hover.rank : ""; };
	var changed = spotKey(spot) !== spotKey(self.hover);
	self.hover = spot;
	if(changed){
		self.requestDraw();
	}
	if(!spot){
		self.tooltip.style("display", "none");
		self.canvas.style("cursor", "default");
		return;
	}
	self.canvas.style("cursor", "pointer");
	var esc = function(str){ return _.escape(String(str)); };
	var html = "";
	if(spot.row < 0){
		html = self.targetHoverHtml(self.targetOrder[spot.col], esc);
	} else if(spot.col < 0){
		html = self.sampleHoverHtml(self.sampleOrder[spot.row], esc);
	} else {
		var sampleIdx = self.sampleOrder[spot.row];
		var targetIdx = self.targetOrder[spot.col];
		var target = self.targets[targetIdx];
		html = "<b>" + esc(self.samples[sampleIdx]) + "</b> in <b>" + esc(target["name"]) + "</b><br>";
		var cellCalls = self.cellCalls(sampleIdx, targetIdx);
		if(0 === cellCalls.length){
			html += "no haplotypes";
		} else {
			html += "<table style='margin-top:4px'><tr><th></th><th style='padding-right:8px'>haplotype</th><th style='padding-right:8px'>fraction</th><th style='padding-right:8px'>reads</th><th>samples with it</th></tr>";
			cellCalls.forEach(function(call){
				html += "<tr" + (call.rank === spot.rank ? " style='font-weight:bold'" : "") + "><td><span style='display:inline-block;width:12px;height:12px;margin-right:4px;border:1px solid #000;background:" + self.colors[targetIdx][call.rank] + "'></span></td>"
					+ "<td style='padding-right:8px'>" + esc(target["haps"][call.rank]) + "</td>"
					+ "<td style='padding-right:8px'>" + call.frac.toFixed(3) + "</td>"
					+ "<td style='padding-right:8px'>" + call.reads + "</td>"
					+ "<td>" + target["hapSampleCounts"][call.rank] + " of " + target["samples"] + "</td></tr>";
			});
			html += "</table>";
		}
	}
	self.tooltip.html(html).style("display", "block");
	//keep the tooltip on screen
	var node = self.tooltip.node();
	var tipLeft = event.clientX + 16;
	var tipTop = event.clientY + 16;
	if(tipLeft + node.offsetWidth > window.innerWidth - 8){
		tipLeft = event.clientX - node.offsetWidth - 16;
	}
	if(tipTop + node.offsetHeight > window.innerHeight - 8){
		tipTop = event.clientY - node.offsetHeight - 16;
	}
	self.tooltip.style("left", Math.max(4, tipLeft) + "px").style("top", Math.max(4, tipTop) + "px");
};

HaplotypeRainbow.prototype.targetHoverHtml = function(targetIdx, esc){
	var self = this;
	var target = self.targets[targetIdx];
	var html = "<b>" + esc(target["name"]) + "</b><br>"
		+ target["samples"] + " of " + self.sampleCount + " samples with haplotypes, " + target["haps"].length + " haplotypes";
	var shown = Math.min(10, target["haps"].length);
	if(shown > 0){
		html += "<table style='margin-top:4px'><tr><th></th><th style='padding-right:8px'>haplotype</th><th>samples with it</th></tr>";
		for(var rank = 0; rank < shown; ++rank){
			html += "<tr><td><span style='display:inline-block;width:12px;height:12px;margin-right:4px;border:1px solid #000;background:" + self.colors[targetIdx][rank] + "'></span></td>"
				+ "<td style='padding-right:8px'>" + esc(target["haps"][rank]) + "</td><td>" + target["hapSampleCounts"][rank] + "</td></tr>";
		}
		html += "</table>";
		if(target["haps"].length > shown){
			html += "and " + (target["haps"].length - shown) + " more";
		}
	}
	return html;
};

HaplotypeRainbow.prototype.sampleHoverHtml = function(sampleIdx, esc){
	var self = this;
	var html = "<b>" + esc(self.samples[sampleIdx]) + "</b><br>"
		+ "haplotypes in " + self.targetsTyped[sampleIdx] + " of " + self.targetCount + " targets";
	self.metaFields.forEach(function(field, fieldPos){
		html += "<br>" + esc(field) + ": " + esc(self.sampleMeta[sampleIdx][fieldPos]);
	});
	return html;
};

HaplotypeRainbow.prototype.onClick = function(x, y, event){
	var self = this;
	var spot = self.locate(x, y);
	if(!spot){
		return;
	}
	var url = null;
	if(spot.row < 0){
		url = "/" + self.rName + "/mainProjectPage/" + encodeURIComponent(self.targets[self.targetOrder[spot.col]]["name"]);
	} else if(spot.col < 0){
		url = "/" + self.rName + "/sampleAcrossTargetsPage/" + encodeURIComponent(self.samples[self.sampleOrder[spot.row]]);
	} else {
		var sampleIdx = self.sampleOrder[spot.row];
		var targetIdx = self.targetOrder[spot.col];
		if(self.cellStart[sampleIdx * self.targetCount + targetIdx + 1] === self.cellStart[sampleIdx * self.targetCount + targetIdx]){
			return;
		}
		url = "/" + self.rName + "/sampleMainPage/" + encodeURIComponent(self.targets[targetIdx]["name"]) + "/" + encodeURIComponent(self.samples[sampleIdx]);
	}
	if(event.metaKey || event.ctrlKey){
		window.open(url, "_blank");
	} else {
		window.location = url;
	}
};

HaplotypeRainbow.prototype.downloadPng = function(){
	var self = this;
	//browsers refuse canvases much bigger than this
	if(self.fullWidth > 16000 || self.fullHeight > 16000 || self.fullWidth * self.fullHeight > 120000000){
		alert("The plot is " + self.fullWidth + " by " + self.fullHeight + " pixels, too big for one image, make the rows or columns smaller first");
		return;
	}
	var canvas = document.createElement("canvas");
	canvas.width = self.fullWidth;
	canvas.height = self.fullHeight;
	self.draw(canvas.getContext("2d"), 0, 0, self.fullWidth, self.fullHeight, false);
	canvas.toBlob(function(blob){
		rainbowSaveBlob(blob, self.fileStub() + ".png");
	});
};

HaplotypeRainbow.prototype.svgString = function(){
	var self = this;
	var svg = new RainbowSvgContext(self.fullWidth, self.fullHeight);
	self.draw(svg, 0, 0, self.fullWidth, self.fullHeight, false);
	return svg.toString();
};

HaplotypeRainbow.prototype.downloadSvg = function(){
	var self = this;
	//one element per haplotype and empty cell, past this the file is too big for most programs to open
	var elements = self.callRank.length + self.sampleCount * self.targetCount;
	if(elements > 2000000){
		alert("The SVG would have about " + elements + " elements, too many to open well, use a group of samples or the PNG instead");
		return;
	}
	rainbowSaveBlob(new Blob([self.svgString()], {type: "image/svg+xml"}), self.fileStub() + ".svg");
};

HaplotypeRainbow.prototype.downloadTable = function(){
	var self = this;
	var lines = [["library_sample_name", "target_name", "microhaplotype_name", "population_rank", "reads", "within_sample_freq"].concat(self.metaFields).join("\t")];
	self.sampleOrder.forEach(function(sampleIdx){
		self.targetOrder.forEach(function(targetIdx){
			var target = self.targets[targetIdx];
			self.cellCalls(sampleIdx, targetIdx).forEach(function(call){
				lines.push([self.samples[sampleIdx], target["name"], target["haps"][call.rank], call.rank + 1, call.reads, call.frac].concat(self.sampleMeta[sampleIdx]).join("\t"));
			});
		});
	});
	rainbowSaveBlob(new Blob([lines.join("\n") + "\n"], {type: "text/tab-separated-values"}), self.fileStub() + ".tsv");
};

$(document).ready(function(){
	var rName = getRootName();
	//the path is the root name, the page, then optionally the meta field and its value
	var pathSegs = window.location.pathname.split("/").filter(function(seg){ return seg.length > 0; });
	var groupUrl = pathSegs.length >= 4 ? pathSegs[2] + "/" + pathSegs[3] : "";
	var groupName = pathSegs.length >= 4 ? decodeURIComponent(pathSegs[2]) : "";
	var subGroupName = pathSegs.length >= 4 ? decodeURIComponent(pathSegs[3]) : "";
	var title = "Haplotype Rainbow" + (groupUrl ? " " + groupName + ": " + subGroupName : "");
	addDiv("body", "topNav");
	//full width so wide panels fit more targets
	addMainDiv("body", true).attr("class", "container-fluid theme-showcase");
	setHeadTitle(title);
	$("#jumboTitle").text(title);
	addPanelWithDiv("#mainContent", "rainbowSummary", "Summary");
	addPanelWithDiv("#mainContent", "rainbowControls", "Options");
	addPanelOnlyHead("#mainContent", "Haplotypes Across Samples and Targets");
	addDiv("#mainContent", "rainbowLegend");
	addDiv("#mainContent", "rainbowPlot");
	getJSON('/' + rName + '/projectsNames').then(function(projectsNames){
		createMainPageNavBar("#topNav", projectsNames["projects"]);
	}).catch(logRequestError);
	prsentDivGifLoading();
	Promise.all([
		getJSON('/' + rName + '/haplotypeRainbow' + (groupUrl ? "/" + groupUrl : "")),
		getJSON('/' + rName + '/metaGroups')
	]).then(function(results){
		new HaplotypeRainbow(rName, results[0], results[1]["groups"], groupName, subGroupName);
		//only PMOs with protein variants have amino acid calls
		if(!results[0]["hasAminoAcidCalls"]){
			return;
		}
		addPanelOnlyHead("#mainContent", "Amino Acids in These Samples");
		addDiv("#mainContent", "aminoAcidNote");
		d3.select("#aminoAcidNote").append("p")
			.text("Every amino acid at each typed position, the reference first. AlleleFrequency is the mean within sample read frequency over the samples typed at the position (SampleTotal), "
				+ "SamplePrevalence is the fraction of those samples with the amino acid, AlleleCount and AlleleTotal count sample haplotypes.");
		addDiv("#mainContent", "aminoAcidTable");
		return getJSON('/' + rName + '/aminoAcidChanges' + (groupUrl ? "/" + groupUrl : "")).then(function(aminoAcids){
			var stub = "aminoAcids" + (groupUrl ? "_" + groupName + "_" + subGroupName : "_allSamples");
			new njhTable("#aminoAcidTable", aminoAcids, stub.replace(/[^A-Za-z0-9_.-]/g, "_"), false);
		});
	}).catch(logRequestError).then(function(){
		removeAllDivGifLoading();
	});
});
