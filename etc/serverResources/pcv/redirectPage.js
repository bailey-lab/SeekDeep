$(document).ready(function(){
	var rName = getRootName();
	addDiv("body", "topNav");
	getJSON('/' + rName + '/projectsNames').then(function (projectsNames) {
		createMainPageNavBar("#topNav", projectsNames);
		addMainDiv("body", true);
		setHeadTitle(rName);
		$("#jumboTitle").html("Ruh Roh");
		$(".jumbotron", "#mainContent").append("<h3>You attempted to go to a page that doesn't exist, or it could be there is no data for where you attempted to go.</h3>");
	}).catch(logRequestError);

});