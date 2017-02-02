function foodweblayout() {
	
	//----------------------------
	// Default parameters
	//----------------------------
	
	var paddingLeft = 50,
		paddingRight = 10, 
		paddingTop = 10, 
		paddingBottom = 10,
		tllim = [NaN, NaN],
		rlim = [1, 50],
		blim = [NaN, NaN],
		totwidth = 960,
		totheight = 600,
		colfun = NaN,
		strengthflxlink = 0.0001,
		strengthgrplink = 0.1,
		strengthtrophic = 0.9,
		strengthxcenter = 0.05,
		nodepad = 10,
		animating = true,
		animationStep = 0;
		
	//----------------------------
	// Main plotting routine
	//----------------------------

	function chart(selection) {
		selection.each(function(data, i) {
			
			// Add SVG canvas
			
			var svg = d3.select(this).append("svg")
					.attr("width", totwidth)
					.attr("height", totheight);
			 
			// Draw rectangle in axis plotting area, for reference

			svg.append("rect")
				.attr("x", paddingLeft)
				.attr("y", paddingTop)
				.attr("width", totwidth - paddingLeft - paddingRight)
				.attr("height", totheight - paddingTop - paddingBottom)
				.attr("fill", "white");
				
			// Set up force simulation
	
			var simulation = d3.forceSimulation();
			
			// Extract data

			nodedata = data.nodes[0]
			flxdata = data.links1[0]
			grpdata = data.links2[0]

			// Scales

			tllimdefault = d3.extent(nodedata, function(d) {return d.TL});
			
			if (isNaN(tllim[0])) {
				tllim[0] = tllimdefault[0]
			}
			if (isNaN(tllim[1])) {
				tllim[1] = tllimdefault[1]
			}

			blimdefault = d3.extent(nodedata, function(d) {return d.B});
			if (isNaN(blim[0])) {
				blim[0] = blimdefault[0]
			}
			if (isNaN(blim[1])) {
				blim[1] = blimdefault[1]
			}

			var tl2y = d3.scaleLinear()
				.domain([tllim[0], tllim[1]])
				.range([totheight-paddingBottom, paddingTop]);

			var yAxis = d3.axisLeft(tl2y);

			var b2r = d3.scaleSqrt()
				.domain(blim)
				.range(rlim)
				.clamp(true);
				
			// Color scale: If nodes have a cval property, use that with a 
			// continuous color scale (by default, blue-pink-yellow one on the 
			// domain of 0-1).	Otherwise, use node type for color.	 Flux lines 
			// without a cval property match the value of their source node.

			var coltype = d3.scaleOrdinal()
				.domain([0,1,2,3])
				.range(["#b3cde3", "#ccebc5", "#fbb4ae", "#decbe4"])
				.unknown("#ffffff");

			var colgrp = d3.scaleOrdinal()
				.domain([1,2,3,4,5,6,7,8,9,10])
				.range(d3.schemeCategory10);

			var colval = d3.scaleLinear()
				.domain([0,1])
				.interpolate(d3.interpolateHcl)
				.range([d3.rgb("#007AFF"), d3.rgb('#FFF500')]);
			
			if ('cval' in nodedata[0]) {
				if (typeof colfun === "number" & isNaN(colfun)) {
					colfun = colval;
				}
			} else {
				for (var ii = 0; ii < nodedata.length; ii++) {
					nodedata[ii].cval = nodedata[ii].type;
				}
				if (typeof colfun === "number" & isNaN(colfun)) {
					colfun = coltype;
				}
			}
				
			// Add axis
	
			svg.append("g")
				 .attr("class", "y-axis")
				 .attr("transform", "translate(" + (paddingLeft) + ", 0)")
				 .call(yAxis); 
				 
			// Set up node circles and link lines

			var flink = svg.append("g")
				.attr("class", "flinks")
				.selectAll("line")
				.data(flxdata)
				.enter().append("line");
	
			var glink = svg.append("g")
				.attr("class", "glinks")
				.selectAll("line")
				.data(grpdata)
				.enter().append("line")
					.attr("stroke", function(d) {return colgrp(d.TG-1); })
					.attr("stroke-width", 2);

			var node = svg.append("g")
				.attr("class", "nodes")
				.selectAll("circle")
				.data(nodedata)
				.enter().append("circle")
					.attr("class", "nodecircle")
					.attr("r", function(d) {return b2r(d.B); })
					.style("fill", function(d) {return colfun(d.cval); })
					.style("stroke", "white")
					.style("stroke-width", 0.5)
					.call(d3.drag()
						.on("start", dragstarted)
						.on("drag", dragged)
						.on("end", dragended));

			node.append("title")
				.text(function(d) { return d.id; });
				 
			// Start nodes with y = trophic level, x = center

			tgmax = d3.max(nodedata, function(d) {return d.TG});
			dx = (totwidth-paddingLeft-paddingRight)/(tgmax+1);
			for (var i = 0; i < nodedata.length;	i++) {
				if (!('x' in nodedata[i])) {
					nodedata[i].x = nodedata[i].TG*dx + paddingLeft;
					nodedata[i].y = tl2y(nodedata[i].TL);
				}
			}
		
			// Add forces to force simulation (in order of priority):
			// trophic = nudge all nodes to y position based on trophic axis
			// linkgrp = links between trophic groups, strong
			// linkflx = food web biomass flux links, weak
			// collide = prevent overlap of node circles
			// charge	 = default repelling force between nodes 
			// centerx = weak centering force to keep nodes on svg canvas
		
			simulation		
				.force("linkflx", d3.forceLink().id(function(d) { return d.id; }).strength(strengthflxlink))
				.force("linkgrp", d3.forceLink().id(function(d) { return d.id; }).strength(strengthgrplink))
				.force("charge",	d3.forceManyBody())
				.force("centerx", d3.forceX().strength(strengthxcenter).x(totwidth/2))
				.force("trophic", d3.forceY().strength(strengthtrophic).y(function (d) {return d.type >= 4 ? totheight/2 : tl2y(d.TL); } ))
				.force("collide", d3.forceCollide().radius(function(d) {return b2r(d.B) + nodepad; }));
		
			simulation
				.nodes(nodedata)
				.on("tick", ticked);

			simulation.force("linkflx")
				.links(flxdata);
	
			simulation.force("linkgrp")
				.links(grpdata);
				
			// On tick, update lines and circles

			function ticked() {
				simulation.stop();
				flink
						.attr("x1", function(d) { return d.source.x; })
						.attr("y1", function(d) { return d.source.y; })
						.attr("x2", function(d) { return d.target.x; })
						.attr("y2", function(d) { return d.target.y; });
		
				glink
						.attr("x1", function(d) { return d.source.x; })
						.attr("y1", function(d) { return d.source.y; })
						.attr("x2", function(d) { return d.target.x; })
						.attr("y2", function(d) { return d.target.y; });

				node
						.attr("cx", function(d) { return d.x; })
						.attr("cy", function(d) { return d.y; });
				
				if (animating) {
							setTimeout(function() { simulation.restart(); },animationStep);
				}
			}
			
			function dragstarted(d) {
				if (!d3.event.active) simulation.alphaTarget(0.3).restart();
				d.fx = d.x;
				d.fy = d.y;
			}

			function dragged(d) {
				d.fx = d3.event.x;
				d.fy = d3.event.y;
			}

			function dragended(d) {
				if (!d3.event.active) simulation.alphaTarget(0);
				d.fx = null;
				d.fy = null;
			}
			
		
		}); // end selection.each
	} // end chart
	
	
	
	//----------------------------
	// Fill in user-set variables
	//----------------------------
		
	chart.paddingLeft = function(value) {
			if (!arguments.length) return paddingLeft;
			paddingLeft = value;
			return chart;
	};
	chart.paddingRight = function(value) {
			if (!arguments.length) return paddingRight;
			paddingRight = value;
			return chart;
	};
	chart.paddingTop = function(value) {
			if (!arguments.length) return paddingTop;
			paddingTop = value;
			return chart;
	};
	chart.paddingBottom = function(value) {
			if (!arguments.length) return paddingBottom;
			paddingBottom = value;
			return chart;
	};
	chart.tllim = function(value) {
			if (!arguments.length) return tllim;
			tllim = value;
			return chart;
	};
	chart.rlim = function(value) {
			if (!arguments.length) return rlim;
			rlim = value;
			return chart;
	};
	chart.blim = function(value) {
			if (!arguments.length) return blim;
			blim = value;
			return chart;
	};
	chart.totwidth = function(value) {
			if (!arguments.length) return totwidth;
			totwidth = value;
			return chart;
	};
	chart.totheight = function(value) {
			if (!arguments.length) return totheight;
			totheight = value;
			return chart;
	};
	chart.colfun = function(value) {
			if (!arguments.length) return colfun;
			colfun = value;
			return chart;
	};
	
	chart.strengthflxlink = function(value) {
			if (!arguments.length) return strengthflxlink;
			strengthflxlink = value;
			return chart;
	};
	chart.strengthgrplink = function(value) {
			if (!arguments.length) return strengthgrplink;
			strengthgrplink = value;
			return chart;
	};
	chart.strengthtrophic = function(value) {
			if (!arguments.length) return strengthtrophic;
			strengthtrophic = value;
			return chart;
	};
	chart.strengthxcenter = function(value) {
			if (!arguments.length) return strengthxcenter;
			strengthxcenter = value;
			return chart;
	};
	chart.nodepad = function(value) {
			if (!arguments.length) return nodepad;
			nodepad = value;
			return chart;
	};
	chart.animating = function(value) {
			if (!arguments.length) return animating;
			animating = value;
			return chart;
	};
	chart.animationStep = function(value) {
			if (!arguments.length) return animationStep;
			animationStep = value;
			return chart;
	};

	return chart;
	
	
	
} // end foodweblayout