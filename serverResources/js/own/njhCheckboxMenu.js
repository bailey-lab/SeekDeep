// Allow Bootstrap dropdown menus to have forms/checkboxes inside, 
// and when clicking on a dropdown item, the menu doesn't disappear.
$(document).on('click', '.dropdown-menu.dropdown-menu-form', function(e) {
  e.stopPropagation();
});

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


function njhCheckboxMenuOrganized(divSelector, names, updateFunction){
	this.divSelector = divSelector;
	this.names = names;
	this.updateFunc = updateFunction;
	var self = this;
	var menuCats = {"misc": []}
	for (obj in names){
		
		if(names[obj].indexOf("_") > -1){
			var pre = names[obj].substr(0, names[obj].indexOf("_") + 1);
			if(pre in menuCats){
				menuCats[pre].push(names[obj]);
			}else{
				menuCats[pre] = [names[obj]];
			}
		}else if (names[obj].indexOf(".") > -1){
			var pre = names[obj].substr(0, names[obj].indexOf(".") + 1);
			if(pre in menuCats){
				menuCats[pre].push(names[obj]);
			}else{
				menuCats[pre] = [names[obj]];
			}
		}else{
			menuCats["misc"].push(names[obj]);
		}
	}
	var menuCatsData = [];
	keys = Object.keys(menuCats),
	keys.sort();
	for (i = 0; i < keys.length; i++){
		cat = keys[i];
		var currentCat = {name : cat, subNames : menuCats[cat]};
		menuCatsData.push(currentCat);
	}
	//console.log(menuCatsData);
	//console.log(menuCats);
	var menu = d3.select(this.divSelector);
	var nav = menu.append("ul").attr("class", "nav nav-pills");
	var cats = menu.select("ul").selectAll("li")
		.data(menuCatsData)
		.enter().append("li")
			.attr("class", "dropdown active")
		
			
	cats.append("a")
		.attr("class", "dropdown-toggle")
		.attr("data-toggle", "dropdown")
		.attr("href", "#")
		.html(function(d){ return d.name + "<span class=\"caret\"></span>"});
	var labs = cats.append("ul")
		.attr("class", "dropdown-menu dropdown-menu-form")
		.attr("id", function(d){ return "menu_" + d.name.replace(".", "");})
		.selectAll("li")
		.data(function(d){ return d.subNames;})
		.enter().append("li")
			.append("div")
				.attr("class", "checkbox")
				.attr("style", "margin-left: 2px")
			.append("label");
	labs.html(function(d){return "<input type=\"checkbox\" name=\"check\" checked=\"true\" value=" + d + " id=" + d + ">" + d;})
		;//.attr("style", "border:2px solid black; padding: 2px;margin: 1px");
	/*labs.append("input")
		.attr("type","checkbox")
		.attr("name", "check")
		.attr("checked", true)
		.attr("value",function(d){return d;})
		.attr("id", function(d){return d;});*/

		
	var drops = cats.select("ul");
	drops.append("li").append("a")
		.attr("href", "javascript:void(0);")
		.text("Check All").on("click", function(d){
			cats.select("#menu_" + d.name.replace(".", "")).selectAll("input").property("checked", true);
			self.updateFunc();
		});
	drops.append("li").append("a")
		.attr("href", "javascript:void(0);")
		.text("Uncheck All").on("click", function(d){
			cats.select("#menu_" + d.name.replace(".", "")).selectAll("input").property("checked", false);
			self.updateFunc();
		});
	//Add check all and uncheck all buttons on their own lines
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

