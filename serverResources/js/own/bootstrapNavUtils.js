function addNavDrop(selector, title, id){
	var top = d3.select(selector)
		.append("li")
			.attr("class", "dropdown");
	top.append("a")
		.attr("href", "#")
		.attr("class", "dropdown-toggle")
		.attr("data-toggle", "dropdown")
		.attr("role", "button")
		.attr("aria-haspopup", "true")
		.attr("aria-expanded", "false")
		.text(function(){ return title + " "})
		.append("span")
			.attr("class", "caret");
	top.append("ul")
		.attr("id", function(){return id;})
		.attr("class", "dropdown-menu scrollable-menu");
} 

function addNavLink(selector, title, link, id, active){
	var top = d3.select(selector)
	top.append("li")
		.attr("class", function(){
			if(active){
				return "active";
			}else{
				return "";
			}
		})
		.append("a")
			.attr("id", id)
			.attr("href", link)
			.text(title);
}

function addFixedTopNavSkeleton(selector, mainTitle, leftNavId, rightNavId){
	var top = d3.select(selector).append("nav")
		.attr("class", "navbar navbar-default navbar-fixed-top")
		.append("div")
			.attr("class", "container");
	var header = top.append("div")
		.attr("class", "navbar-header");
	var mobileNavBut = header.append("button")
			.attr("class", "navbar-toggle collapsed")
			.attr("data-toggle", "collapse")
			.attr("data-target", "#navbar")
			.attr("aria-expanded", "false")
			.attr("aria-controls", "navbar");
	mobileNavBut.append("span")
		.attr("class", "sr-only")
		.text("Toggle navigation");
	mobileNavBut.append("span")
		.attr("class", "icon-bar");
	mobileNavBut.append("span")
		.attr("class", "icon-bar");
	mobileNavBut.append("span")
		.attr("class", "icon-bar");
	header.append("a")
		.attr("class", "navbar-brand")
		.style("pointer-events", "none")
		.attr("href", "")
		.text(mainTitle);
	var navBar = top.append("div")
		.attr("id", "navbar")
		.attr("class", "navbar-collapse collapse");
	navBar.append("ul")
		.attr("class", "nav navbar-nav")
		.attr("id", leftNavId);
	navBar.append("ul")
		.attr("class", "nav navbar-nav navbar-right")
		.attr("id", rightNavId);
	
}




