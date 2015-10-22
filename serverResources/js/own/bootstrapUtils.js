function addPanelWithDiv(addToSelector, divId, panelTitle, panelType){
	var top = d3.select(addToSelector)
				.append("div")
					.attr("class", "panel " + (panelType || "panel-primary"));
	top.append("div")
		.attr("class", "panel-heading")
			.append("h1")
				.text(panelTitle);
	top.append("div")
		.attr("class", "panel-body")
			.append("div")
				.attr("id", divId);
}

function addPanelOnlyHead(addToSelector, panelTitle, panelType){
	var top = d3.select(addToSelector)
				.append("div")
					.attr("class", "panel " + (panelType || "panel-primary"));
	top.append("div")
		.attr("class", "panel-heading")
			.append("h1")
				.text(panelTitle);
}

function addMainDiv(addToSelector, addJumbo){
	var main = d3.select(addToSelector).append("div")
		.attr("class", "container theme-showcase")
		.attr("id", "mainContent")
		.attr("role", "main");
	if(addJumbo){
		main.append("div")
			.attr("class", "jumbotron")
			.append("h1")
				.attr("id", "jumboTitle");
	}
	return main;
}

