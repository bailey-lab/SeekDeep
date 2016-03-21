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

$.fn.scrollView = function () {
    return this.each(function () {
        $('html, body').animate({
            scrollTop: $(this).offset().top
        }, 1000);
    });
}


function drawCircle(x, y, radius, color, borderColor) {
	context.beginPath();
	context.arc(x, y, radius, 0, 2 * Math.PI, false);
	context.fillStyle = color;
	context.fill();
	context.lineWidth = 5;
	context.strokeStyle = borderColor;
	context.stroke();
}

function drawLine(sx, sy, ex, ey, width) {
	context.beginPath();
	context.moveTo(sx, sy);
	context.lineTo(ex, ey, width);
	context.stroke();
}


function setToArray(set) {
	var array_var = [];
	set.forEach(function(element){array_var.push(element)});
	return array_var;
}

function ajax(url, func) {
	$.ajax({
		url : url,
		dataType : 'json',
		async : false,
		success : function(ct) {
			func(ct);
		}
	});
}

function ajaxRet(url){
	var data;
	ajax(url, function(d){ data = d});
	return data;
}

function ajaxPost(url, data, func) {
	$.ajax({
		type: 'POST',
        url: url,
        datatype : 'json',
        async : false,
        data: data,
		success: function(ct) {
			func(ct);
		}
    });
}


function ajaxAsync(url, func) {
	$.ajax({
		url : url,
		dataType : 'json',
		async : true,
		success : function(ct) {
			func(ct);
		}
	});
}

var range = function(start, end, step) {
	// from http://stackoverflow.com/questions/3895478/does-javascript-have-a-range-equivalent
	var range = [];
	var typeofStart = typeof start;
	var typeofEnd = typeof end;

	if (step === 0) {
		throw TypeError("Step cannot be zero.");
	}

	if (typeofStart == "undefined" || typeofEnd == "undefined") {
		throw TypeError("Must pass start and end arguments.");
	} else if (typeofStart != typeofEnd) {
		throw TypeError("Start and end arguments must be of same type.");
	}

	typeof step == "undefined" && ( step = 1);

	if (end < start) {
		step = -step;
	}

	if (typeofStart == "number") {

		while (step > 0 ? end >= start : end <= start) {
			range.push(start);
			start += step;
		}

	} else if (typeofStart == "string") {

		if (start.length != 1 || end.length != 1) {
			throw TypeError("Only strings with one character are supported.");
		}

		start = start.charCodeAt(0);
		end = end.charCodeAt(0);

		while (step > 0 ? end >= start : end <= start) {
			range.push(String.fromCharCode(start));
			start += step;
		}

	} else {
		throw TypeError("Only string and number types are supported");
	}

	return range;

}; 



function addMouseScrollListener(obj, up, down) {
	// from http://www.sitepoint.com/html5-javascript-mouse-wheel/
	this.handler = function(e) {

		var e = window.event || e;
		// old IE support

		var delta = Math.max(-1, Math.min(1, (e.wheelDelta || -e.detail)));
		if (delta > 0) {
			up(delta);
		} else {
			down(delta);
		}
		e.preventDefault();
	};
	if (obj.addEventListener) {
		// IE9, Chrome, Safari, Opera
		obj.addEventListener("mousewheel", handler, false);
		// Firefox
		obj.addEventListener("DOMMouseScroll", handler, false);
	} else {
		// IE 6/7/8
		obj.attachEvent("onmousewheel", handler);
	}
}


var getRelCursorPosition = function(event, obj) {
	// from http://stackoverflow.com/a/5417934
	var objOffset = $(obj).offset();
	var x = event.clientX + document.body.scrollLeft + document.documentElement.scrollLeft - Math.floor(objOffset.left);
	var y = event.clientY + document.body.scrollTop + document.documentElement.scrollTop - Math.floor(objOffset.top) + 1;
	return [x, y];
}; 

    

function addDiv(parentId, childName) {
	$("<div id =\"" + childName + "\"></div>").appendTo(parentId);
};

function addH1(selector, text){
	d3.select(selector).append("h1").text(text);
}

var sort_by = function(field, reverse, primer){
	//from http://stackoverflow.com/questions/979256/sorting-an-array-of-javascript-objects
   var key = function (x) {return primer ? primer(x[field]) : x[field];};

   return function (a,b) {
	  var A = key(a), B = key(b);
	  return ( (A < B) ? -1 : ((A > B) ? 1 : 0) ) * [-1,1][+!!reverse];                  
   };
};

function arrayContains(arr, val){
	return (arr.indexOf(val) > -1);
}


function addSvgSaveButton(buttonId, topSvg) {
	d3.select(buttonId).append("a").attr("id", "imgDownload");
	d3.select(buttonId).on(
			"click",
			function() {
				var html = $(
						d3.select(topSvg).attr("version", 1.1).attr("xmlns",
								"http://www.w3.org/2000/svg").node()).clone()
						.wrap('<p/>').parent().html();
				;
				// add the svg information to a and then click it to trigger the
				// download
				var imgsrc = 'data:image/svg+xml;base64,' + btoa(html);
				d3.select("#imgDownload").attr("download", "graph.svg");
				d3.select("#imgDownload").attr("href", imgsrc);
				var a = $("#imgDownload")[0];
				a.click();
			});
}

function setHeadTitle(title){
	if(!$("title", "head").length){
		$("head").append("<title></title>");
	}
	$("title", "head").html(title);
}
