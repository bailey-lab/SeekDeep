    function tabulate(data, columns, divId) {
		/** Add json data organized row wise into a table in the div with the divId with only the columns given
		 *
		 */
		//adapted from http://stackoverflow.com/questions/9268645/creating-a-table-linked-to-a-csv-file
	    var table = d3.select(divId).append("table"),
	        thead = table.append("thead"),
	        tbody = table.append("tbody");
	
	    // append the header row
	    thead.append("tr")
	        .selectAll("th")
	        .data(columns)
	        .enter()
	        .append("th")
	            .attr("style", "font-weight: bold; padding: 2px 4px;")
	            .text(function(column) { return column; });
	
	    // create a row for each object in the data
	    var rows = tbody.selectAll("tr")
	        .data(data)
	        .enter()
	        .append("tr");
	
	    // create a cell in each row for each column
	    var cells = rows.selectAll("td")
	        .data(function(row) {
	            return columns.map(function(column) {
	                return {column: column, value: row[column]};
	            });
	        })
	        .enter()
	        .append("td")
	            .attr("style", "padding: 2px 4px;")
	            .text(function(d) { return d.value; });
	    
	    return table;
	}
	
	function createTable(divId) {
		/**create the table and return it to be manipulate*/
		var table = d3.select(divId).append("table");
		table.append("thead").append("tr");
		table.append("tbody");
		return table;
	}
	
	function updateTable(tab, data, columns){
		//console.log(data);
		//ensure header row
		var headerRow = tab.select("thead")
			.selectAll("tr")
			.data([1])
			.enter();
		//attach column name data to header
		var header = tab.select("thead").select("tr")
	        .selectAll("th")
	        .data(columns).text(function(column) { return column; });
	    header
	        .enter()
			.append("th")
				.attr("style", "font-weight: bold; padding: 2px 4px;")
	            .text(function(column) { return column; });
	   //create headers as needed and add bolding 
	  /*header.enter()
	        .append("th")*/
	            
	  //remove any headers that don't have data attached to them
	  //console.log(columns);
	  header.exit()
        	.remove();
		
	    // create a row for each object in the data
	    var newRows = tab.select("tbody").selectAll("tr")
	        .data(data)
	        .enter()
	        .append("tr");
	    tab.select("tbody").selectAll("tr")
	        .data(data).exit().remove();
	    var currentColor = "#e9e9e9";
		var rows = tab.select("tbody").selectAll("tr").style("background-color",function(d,i){
	        		if(i == 0){
	        			return currentColor;
	        		}else {
	        			if (d[columns[0]] != data[i-1][columns[0]]){
	        				if(currentColor == "#e9e9e9"){
	        					currentColor = "#c9c9c9";
	        				}else{
	        					currentColor = "#e9e9e9";
	        				}
	        			}
	        		}
	        		return currentColor;
	        		});;
	    //create a cell in each row for each column
	    //console.log(rows);
	    var cells = rows.selectAll("td")
	        .data(function(row) {
	        	var ret = columns.map(function(column) {
	                return {column: column, value: row[column]};
	            });         
	            return ret;
	        });
	   	cells.enter()
	        .append("td")
	            .attr("style", "padding: 2px 4px;")
	            .text(function(d) { return d.value; });
	    cells.text(function(d) { return d.value; });
	    //remove cells as needed
	    cells.exit()
        	 .remove();
       
	}
	function updateTableWithColors(tab, data, columns){
		//console.log(data);
		//ensure header row
		var headerRow = tab.select("thead")
			.selectAll("tr")
			.data([1])
			.enter();
		//attach column name data to header
		var header = tab.select("thead").select("tr")
	        .selectAll("th")
	        .data(columns).text(function(column) { return column; });
	    header
	        .enter()
			.append("th")
				.attr("style", "font-weight: bold; padding: 2px 4px;")
	            .text(function(column) { return column; });
	   //create headers as needed and add bolding 
	  /*header.enter()
	        .append("th")*/
	            
	  //remove any headers that don't have data attached to them
	  //console.log(columns);
	  header.exit()
        	.remove();
		
	    // create a row for each object in the data
	    var newRows = tab.select("tbody").selectAll("tr")
	        .data(data)
	        .enter()
	        .append("tr");
	    tab.select("tbody").selectAll("tr")
	        .data(data).exit().remove();
	    var currentColor = "#e9e9e9";
		var rows = tab.select("tbody")
		.selectAll("tr")
		.style("background-color",function(d){
	        		return d.color;});
	    //create a cell in each row for each column
	    //console.log(rows);
	    var cells = rows.selectAll("td")
	        .data(function(row) {
	        	var ret = columns.map(function(column) {
	                return {column: column, value: row[column]};
	            });         
	            return ret;
	        });
	   	cells.enter()
	        .append("td")
	            .attr("style", "padding: 2px 4px;")
	            .text(function(d) { return d.value; });
	    cells.text(function(d) { return d.value; });
	    //remove cells as needed
	    cells.exit()
        	 .remove();
       
	}
	
			function createLinksTable(addToSelector,linkPrefix, links, colNum, mouseOverColor, mouseLeaveColor){
				var dataset = [],
				tmpDataset = [],
				i, j;
				var rowNum = Math.ceil(links.length/colNum);
				for (i = 0; i < rowNum; ++i) {
				    for (j = 0, tmpDataset = []; j < colNum && i * colNum + j < links.length; ++j) {
				        tmpDataset.push({text: links[i * colNum + j],link: linkPrefix + links[i * colNum + j] });
				    }
				    dataset.push(tmpDataset);
				}
				var tab = d3.select(addToSelector)
				    .append("table")
				    .style("border-collapse", "collapse")
				    .style("border", "2px black solid")
				    
				    .selectAll("addToSelector tr")
				    .data(dataset)
				    .enter().append("tr")
				    
				    .selectAll("td")
				    .data(function(d){return d;})
				    .enter().append("td")
				    .style("border", "1px black solid")
				    .style("padding", "10px")
				    .on("mouseover", function(){d3.select(this).style("background-color", mouseOverColor);}) 
				    .on("mouseout", function(){d3.select(this).style("background-color", mouseLeaveColor);}) 
				    .append("a")
		  				.attr("href", function(d){return d.link;})
		  				.html(function(d){return d.text;});
		  		return tab;
			}
			
			function createLinksTableDifferentLinks(addToSelector,linkPrefix, links, encodedLinks, colNum, mouseOverColor, mouseLeaveColor){
				var dataset = [],
				tmpDataset = [],
				i, j;
				var rowNum = Math.ceil(links.length/colNum);
				for (i = 0; i < rowNum; ++i) {
				    for (j = 0, tmpDataset = []; j < colNum && i * colNum + j < links.length; ++j) {
				        tmpDataset.push({text: links[i * colNum + j],link: linkPrefix + encodedLinks[i * colNum + j] });
				    }
				    dataset.push(tmpDataset);
				}
				var tab = d3.select(addToSelector)
				    .append("table")
				    .style("border-collapse", "collapse")
				    .style("border", "2px black solid")
				    
				    .selectAll("addToSelector tr")
				    .data(dataset)
				    .enter().append("tr")
				    
				    .selectAll("td")
				    .data(function(d){return d;})
				    .enter().append("td")
				    .style("border", "1px black solid")
				    .style("padding", "10px")
				    .on("mouseover", function(){d3.select(this).style("background-color", mouseOverColor);}) 
				    .on("mouseout", function(){d3.select(this).style("background-color", mouseLeaveColor);}) 
				    .append("a")
		  				.attr("href", function(d){return d.link;})
		  				.html(function(d){return d.text;});
		  		return tab;
			}