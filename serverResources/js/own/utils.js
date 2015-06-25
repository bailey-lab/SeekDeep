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


var sort_by = function(field, reverse, primer){
	//from http://stackoverflow.com/questions/979256/sorting-an-array-of-javascript-objects
   var key = function (x) {return primer ? primer(x[field]) : x[field];};

   return function (a,b) {
	  var A = key(a), B = key(b);
	  return ( (A < B) ? -1 : ((A > B) ? 1 : 0) ) * [-1,1][+!!reverse];                  
   };
};

