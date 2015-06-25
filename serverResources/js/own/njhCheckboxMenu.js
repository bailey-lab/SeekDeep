
function njhCheckboxMenu(divSelector, names, updateFunction){
	this.divSelector = divSelector;
	this.names = names;
	this.updateFunc = updateFunction;
	//create the inputs checkboxes for the sample name 
	var menu = d3.select(this.divSelector);
	var inputs = menu.selectAll("input")
		.data(this.names)
		.enter()
		.append("label")
			.html(function(d){return d;})
			.attr("style", "border:2px solid black; padding: 2px;margin: 1px")
		.append("input")
			.attr("type","checkbox")
			.attr("name", "check")
			.attr("checked", true)
			.attr("value",function(d){return d;})
			.attr("id", function(d){return d;});
	//Add check all and uncheck all buttons on their own lines
	menu.append("br");
	menu.append("a")
		.attr("href", "javascript:void(0);")
		.attr("class", "njhCheckAll")
		.text("Check All");
	menu.append("br");
	menu.append("a")
		.attr("href", "javascript:void(0);")
		.attr("class", "njhUncheckAll")
		.text("Uncheck All");

	//add the un-check all and check all functions
	var self = this;
	menu.select(".njhCheckAll").on("click", function(){
		menu.selectAll("input").property("checked", true);
		self.updateFunc();
	});
	menu.select(".njhUncheckAll").on("click", function(){
		menu.selectAll("input").property("checked", false);
		self.updateFunc();
	});
	//add update function to the checkboxes
	menu.selectAll("input").on("click", this.updateFunc);
};

