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




