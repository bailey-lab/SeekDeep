
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
		var currentMenuItemOptions = currentMenuItem.append("ul").attr("class", "dropdown-menu");
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
   
   SeqPainter.prototype.paintSelectedSeq = function(seqContext,seq, currentBase){
   		seqContext.font = "bold 15px Arial, sans-serif";
   		seqContext.textAlign = "left";
    	seqContext.textBaseline = "middle";
   		var logInfo = "name: " + 
        	seq["name"]
        	+ " base: "  + seq["seq"][currentBase] 
        	+ " qual: " +  seq["qual"][currentBase]
        	+ " pos: " + currentBase;
        var tWidth = seqContext.measureText(logInfo).width;
        seqContext.fillStyle = "#FFFFFF";
	    seqContext.fillRect(this.nameOffSet + this.cw + 2, (this.nSeqs)*this.ch + 2 , this.nBases * this.cw - (3 *this.cw), this.ch);
        seqContext.fillStyle = "#EEEEEE";
	    seqContext.fillRect(this.nameOffSet + this.cw + 2, (this.nSeqs)*this.ch + 2 , tWidth, this.ch);
	    seqContext.fillStyle = "#000000";
	    //console.log(tWidth);
	    seqContext.fillText(logInfo,this.nameOffSet + this.cw + 2, (this.nSeqs)*this.ch + 2  + this.ch/2 );
   };

	function njhSeqView(viewName, seqData, cellWidth, cellHeight, baseColors,addQualChart, minTreeLink){
		//need to add style and html, currently just there
		//retrieve html elements 
		this.topDivName = viewName;
		this.topDiv = document.getElementById(viewName);
		$(this.topDiv).addClass("SeqView");
		this.uid = seqData["uid"];
		this.menuDiv = d3.select(viewName).append("div")
				.attr("class", "njhSeqViewMenu");
		this.setUpDefaultMenu();
		if(minTreeLink != undefined){
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
			.style("border", "2px solid black")
			.style("padding", "2px")
			.style("margin", "5px")
			.style("width", "100px")
			.style("float", "left");
		this.masterDivD3 = d3.select(viewName).append("div").attr("class", "SeqViewCanvasDiv");
		this.masterDivD3.append("div").attr("class", "rightSlider");
		this.masterDivD3.append("canvas").attr("class", "canvas");
		this.masterDivD3.append("div").attr("class", "bottomSlider");
		this.masterDivD3.append("div").attr("class", "pop-up").append("p").attr("id", "info");
		this.masterDivD3.append("div").attr("class", "select");
		d3.select(viewName).append("div").attr("class", "qualChart");
		var self = this;
		var linkButton = d3.select(this.topDivName + " .downFastaDiv")
			.append("button")
			.text("Download Fasta")
			.attr("class", "fastaSaveButton");
		linkButton.append("a")
				.attr("class", "fastaDownLink");
		
		linkButton.on("click", function(){
				    var mainTable = [];
				    //
				    for (i = 0; i <self.seqData["seqs"].length ; i++) { 
						mainTable.push([">" + self.seqData["seqs"][i]["name"]]);
						mainTable.push([self.seqData["seqs"][i]["seq"]]);
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
		var nameOffSet = 10 * cellWidth;
		var numOfBases = Math.floor((this.canvas.width - cellWidth - nameOffSet)/cellWidth);
	 	var numOfSeqs = Math.min(Math.floor((this.canvas.height - cellHeight)/cellHeight), this.seqData["seqs"].length);
	 	//console.log(numOfSeqs);
	 	//console.log(Math.floor((this.canvas.height - cellHeight)/cellHeight));
	 	//console.log(this.seqData["seqs"].length);
		this.painter = new SeqPainter(cellWidth, cellHeight, numOfSeqs, numOfBases, nameOffSet, baseColors);
		//this.seqs = seqs;
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
			    }
			});
		}

		//
	};
	
	
	
	njhSeqView.prototype.setUpDefaultMenu = function(){
		var locSplit = window.location.toString().split(/[\/]+/);
		var rName = locSplit[2];
		var menuItems = {};
		var sortOptions = [];
		var self = this;
		sortOptions.push(new njhMenuItem("sortSeq", "Sequence",function(){
			var mainData;
		    ajax('/' + rName + '/sort/' + self.uid+'/seq', function(md){ mainData = md; });
		    console.log(self.uid);
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortSeqCondensed", "Sequence Condensed",function(){
			var mainData;
	   		ajax('/' + rName + '/sort/' + self.uid+'/seqCondensed', function(md){ mainData = md; });
	    	self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortTotalCount", "Total Read Count",function(){
			var mainData;
		    ajax('/' + rName + '/sort/' + self.uid+'/totalCount', function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortSize", "Length",function(){
			var mainData;
		    ajax('/' + rName + '/sort/' + self.uid+'/size', function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		sortOptions.push(new njhMenuItem("sortName", "Name",function(){
			var mainData;
		    ajax('/' + rName + '/sort/' + self.uid+'/name', function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		menuItems["Sort"] = sortOptions;
		var alnOptions = [];
		alnOptions.push(new njhMenuItem("muscle", "muscle",function(){
			var mainData;
		    ajax('/' + rName + '/muscle/' + self.uid, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		alnOptions.push(new njhMenuItem("removeGaps", "remove gaps",function(){
			var mainData;
		    ajax('/' + rName + '/removeGaps/' + self.uid, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		
		menuItems["Aln"] = alnOptions;
		var editOptions = [];
		editOptions.push(new njhMenuItem("complement", "Reverse Complement",function(){
			var mainData;
		    ajax('/' + rName + '/complement/' + self.uid, function(md){ mainData = md; });
		    self.updateData(mainData);
		}));
		
		menuItems["Edit"] = editOptions;
	
		createSeqMenu(this.topDivName + " .njhSeqViewMenu", menuItems);
	};
	
	njhSeqView.prototype.updateData = function(inputSeqData){
		this.seqData = inputSeqData;
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
	
	njhSeqView.prototype.setUpCanvas = function(){
		$(this.masterDiv).width((window.innerWidth - 10) * 0.98);
		var maxPossHeight = this.painter.ch * (this.seqData["seqs"].length + 4);
		$(this.masterDiv).height(Math.min((window.innerHeight - 60) * 0.80, maxPossHeight));
		this.canvas.width = $(this.masterDiv).width() * 0.96;
		this.canvas.height = $(this.masterDiv).height() * 0.95;
		this.painter.nBases = Math.floor((this.canvas.width - this.painter.cw - this.painter.nameOffSet)/this.painter.cw);
	 	this.painter.nSeqs = Math.min(Math.floor((this.canvas.height - this.painter.ch)/this.painter.ch),this.seqData["seqs"].length);
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
		this.painter.nBases = Math.floor((this.canvas.width - this.painter.cw - this.painter.nameOffSet)/this.painter.cw);
	 	this.painter.nSeqs = Math.floor((this.canvas.height - this.painter.ch)/this.painter.ch);
	};

	njhSeqView.prototype.paint = function(){
		this.painter.paintSeqs(this.context, this.seqData["seqs"], this.seqStart, this.baseStart);
		this.painter.placeBasePos(this.context, this.seqStart, this.baseStart);
		this.paintSelectedSeq();
		this.setSelector();
	};
	
	njhSeqView.prototype.paintSelectedSeq = function(){
		this.painter.paintSelectedSeq(this.context, this.seqData["seqs"][this.currentSeq], this.currentBase );
	};
	
	njhSeqView.prototype.setSelector = function(){
		if(this.currentBase >= this.baseStart && this.currentSeq >=this.seqStart){
			$(this.sel).css('top', (this.currentSeq -this.seqStart) *this.painter.ch -1);
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
   njhSeqView.prototype.clicked = function(e){
        var pt = getRelCursorPosition(e, this.canvas);
        this.currentBase = Math.ceil(pt[0]/this.painter.cw) - this.painter.nameOffSet/this.painter.cw + this.baseStart -1;
        this.currentSeq = Math.ceil(pt[1]/this.painter.ch) + this.seqStart - 1;
        this.paintSelectedSeq();
        //console.log(pt);
        //console.log(this.currentSeq);
        //console.log(this.currentBase);
        //console.log(this.seqData["seqs"][this.currentSeq]["name"]);
        //console.log(this.seqData["seqs"][this.currentSeq]["seq"][this.currentBase]);
        //console.log(this.seqData["seqs"][this.currentSeq]["qual"][this.currentBase]);
        this.setSelector();
        var currentQual = this.seqData["seqs"][this.currentSeq]["qual"];
		this.chart.load({
	        json: {
	            qual: this.seqData["seqs"][this.currentSeq]["qual"]
	        }
	    });
	    //this.chart.xgrids.remove();
	    this.chart.xgrids([{value: this.currentBase, text:this.seqData["seqs"][this.currentSeq]["qual"][this.currentBase]}]);

    };
   njhSeqView.prototype.setUpListeners = function(){
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
	    $(popUpWindow).hide();
	    //console.log(popUpWindow);
	    var rect = this.canvas.getBoundingClientRect();
		//console.log(rect.left, rect.top, rect.right, rect.bottom );
		var currentPoint = getRelCursorPosition(e, this.canvas);
		//console.log("WindowX:" + (currentPoint[0]));
		//console.log("WindowY:" + (currentPoint[1]));
		//console.log("AdjustX:" + (currentPoint[0] - rect.left));
		//console.log("AdjustY:" + (currentPoint[1] - rect.top));
    	$(popUpWindow).css('top', currentPoint[1] + moveDown).css('left', currentPoint[0] + moveLeft);
    	
      	var currentBaseHover = Math.ceil(currentPoint[0]/this.painter.cw) - this.painter.nameOffSet/this.painter.cw + this.baseStart -1;
        var currentSeqHover = Math.ceil(currentPoint[1]/this.painter.ch) + this.seqStart - 1;
        var base = "";
        var qual = "";
        if(this.seqData["seqs"][currentSeqHover]["seq"][currentBaseHover]){
        	base = this.seqData["seqs"][currentSeqHover]["seq"][currentBaseHover];
        	qual = this.seqData["seqs"][currentSeqHover]["qual"][currentBaseHover];
        }

		var maxHeight = Math.min(this.painter.nSeqs, this.seqData["numReads"]) * this.painter.ch;
		var maxWidth = Math.min(this.painter.nBases, this.seqData["maxLen"]) * this.painter.ch + this.painter.nameOffSet;
		if(currentPoint[1] <= maxHeight && 
			currentPoint[0] <= maxWidth){
			        if(currentPoint[0] > this.painter.nameOffSet){
	        	//console.log($("#info", popUpWindow)[0]);
	        $("#info", popUpWindow)[0].innerHTML = "name: " + 
	        	this.seqData["seqs"][currentSeqHover]["name"]
	        	+ "<br>base: "  + base
	        	+ "<br>qual: " +  qual 
	        	+ "<br>pos: " + currentBaseHover;
	        }else{
	        	$("#info", popUpWindow)[0].innerHTML = "name: " + 
	        	this.seqData["seqs"][currentSeqHover]["name"];
	        }
			//$(popUpWindow).fadeIn(500);
			$(popUpWindow).show();
		}else{
			$(popUpWindow).hide();
		}
    }.bind(this)).bind(this);
   };

