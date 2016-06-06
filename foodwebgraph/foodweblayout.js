function foodweblayout() {

	var marl      = 50,      // left margin for main grouping
		marr      = 50,      // right margin
		mart      = 50,      // top margin
		marb      = 50,      // bottom margin
        totwidth  = 800,     // Width of the svg canvas
        totheight = 600,     // Height of the svg canvas
        padding   = 10,      // separation between nodes
        tlnudge   = 2.0,     // strength of trophic level nudging
        linkdist  = 5,       // link distance
        charge    = -100,    // charge
	    hlineop   = 0.5,     // opacity for temporary trophic group lines
	    flineop   = 0.3,     // opacity for flux lines
        fontsz    = 8,       // font size for labels
	    seed      = 'hello', // string used to initialize random number generator
	    cpfac     = 0.8,     // scaling factor for box used to pack circles, compared to svg canvas 
		drawmode  = 'cpack', // drawing mode: 'cpack' or 'force'
	    init      = 'cpack', // node initial poisition mode
		tllim     = [NaN,NaN]; // trophic axis limits

    function chart(selection) {
        selection.each(function(data, i) {

			//--------------------
			// Setup
			//--------------------

			// Seed the random number generator.  Collision-avoidance 
			// includes a bit of randomness, so this is necessary for 
			// reproducability

			Math.seedrandom(seed);

			// Width of grouping object to which all plotted things will 
			// be attached
			
	        var width = totwidth - marr - marl,
	            height = totheight - mart - marb;

            // First, reformat node data from flat input to properly-nested
            // hierarchical format

            var dataMap = data.nodes.reduce(function(map, node) {
                map[node.Name] = node;
                return map;
            	}, {});

            var treeData = [];
            data.nodes.forEach(function(node) {
                var parent = dataMap[node.parent];
                if (parent) {
                    (parent.children || (parent.children = []))
                    .push(node);
                } else {
                    treeData.push(node);
                }
            })

		    // Append SVG to selected DOM object, and append a <g> to
			// that with a bit of a margin

		    var main = d3.select(this).append("svg")
		        .attr("width", totwidth)
			    .attr("height", totheight);
			main.append("rect")
				.attr("width", "100%")
				.attr("height", "100%")
				.attr("fill", "#F8F9F9");
			
			var svg = main.append("g") // misleading variable name, sorry
				.attr("id", "g_svg")
	            .attr("transform", "translate(" + marl + "," + mart + ")");

	        // Use circle packing to calculate hierarchy details and
	        // initialize positions

	        var tree = d3.layout.pack()
	            .size([width*cpfac,height*cpfac])
	            .value(function(d) { return Math.max((d.type >= 4 && d.type < 7 ? 0 : d.B), 0.001); })
				.sort(function(a, b){ return a.TL - b.TL; });
				
            // Assign data to the pack layout, and save tree-ified nodes
		    // for force layout use

            var tnodes = tree.nodes(treeData[0]),
                tlinks = tree.links(tnodes);
				
			// Custom initial positions, if provided
				
			if (init == 'coords') {
				tnodes.forEach(function(d) {
					d.x = d.xinit;
					d.y = d.yinit;
					d.r = d.rinit;
				})
			}

			// Color

            var c10 = d3.scale.category20();
            var col = d3.scale.ordinal()
				      .domain([0,1,2,3])
				      .range(["#AEC7E8", "#98DF8A", "#C49C94", "#FF9896"]);
					  // .range(["rgb(236,237,135)", "rgb(186,222,131)", "rgb(240,180,139)", "rgb(152,222,224)"]);

			switch (drawmode) {
				
			//--------------------
			// Plot initial layout
			//--------------------
					
			case 'cpack': // Only draw the initial positions
				
		        svg.selectAll('circles')
		            .data(tnodes)
		          .enter().append('svg:circle')
		  	  	    .attr("class", function(d) { return d.children ? "node" : "leaf node"; })
		            .attr('cx', function(d) { return d.x; })
		            .attr('cy', function(d) { return d.y; })
		  	        .attr('r', function(d) { return d.r; })
		  		    .style("fill", function(d) {return d.type >= 4 ? "white" : col(d.type); })
		  		    .style("stroke", function(d) {return c10(d.TG); })
					.style("stroke-width", 2)
				
				// Add corner reference points (used by 
				// waitforfoodweb to check that everything is 
				// finished)
					
       			svg.append("circle").attr("cx",0).attr("cy",0).attr("r",1).attr("class","ref")
        		svg.append("circle").attr("cx",width).attr("cy",height).attr("r",1).attr("class","ref")
				
				
				break;
				
			//--------------------
			// Run force layout
			//--------------------	
				
			case 'force': // Run the force layout
				
	            // Eliminate central-web-links for layout purposes
                
	            var tlinks2 = [];
	            tlinks.forEach(function(lnk) {
	                if (!(lnk.source.Name == "web" || lnk.target.Name == "web")) {
	                    tlinks2.push(lnk)
	                }
	            })
				
	            // Initialize layout for nodes and labels

	            var force = d3.layout.force()  // Nodes positioning
	                .nodes(tnodes)
	                .links(tlinks2)
	                .size([width, height])
	                .linkDistance(function(d) { return  d.source.Name.substring(0,3)=="lev1" ? linkdist*5 : linkdist; })
					.charge(charge)
					.gravity(0.1)
	                .start();
					
			    // Create hierarchical edge lines (for testing) 
        
			    var hlink = svg.selectAll(".link")
			                    .data(tlinks2)
			                  .enter().append("line")
			                    .attr("class", "link")
			                    .style("stroke", "rgb(150,150,150)")
								.style("opacity", hlineop); 
								

			    // Container for flux lines, so under nodes
						
				fluxdata = data.edges[0];				
				for(var i = 0, j = fluxdata.length; i < j; i++) {
				    for(var k = 0, l = tnodes.length; k < l; k++) {
				        if( fluxdata[i].source == tnodes[k].Name)  {
				            fluxdata[i].sidx = k;
				        }
				        if( fluxdata[i].target == tnodes[k].Name) {
				            fluxdata[i].tidx = k;
				        }
				    }
				}
				
				// Add flux lines
			
				var fluxlines = svg.selectAll(".fluxline")
		            .data(fluxdata)
	            .enter().append("line")
		            .attr("class", "fluxline")
		            .attr("x1", function(d) { return (tnodes[d.sidx].x); })
		            .attr("y1", function(d) { return (tnodes[d.sidx].y); })
		            .attr("x2", function(d) { return (tnodes[d.tidx].x); })
		            .attr("y2", function(d) { return (tnodes[d.tidx].y); })
		            .attr("stroke-width", 1)
		            .attr("stroke", "#D1D0CE")
					.attr("opacity", flineop);
				
				var fluxg = svg.append("g").selectAll(".fluxline");
								
	            // Create container for label-connector edges 
			    // now so they're layered under the nodes and labels
                            
				var lline = svg.append("g").selectAll(".labelline")	
								
	            // Add the nodes, with circles
                
	            var node = svg.selectAll(".node")
	                .data(tnodes)
	                .enter().append("g");
					
				var node_drag = d3.behavior.drag()
				        .on("dragstart", dragstart)
				        .on("drag", dragmove)
				        .on("dragend", dragend);
						
				function dragstart(d, i) {
				        force.stop() // stops the force auto positioning before you start dragging
				    }

			    function dragmove(d, i) {
			        d.px += d3.event.dx;
			        // d.py += d3.event.dy;
			        d.x += d3.event.dx;
			        // d.y += d3.event.dy;
			        tick({type: "tick", alpha: 0.00}); // this is the key to make it work together with updating both px,py,x,y on d !
			    }

			    function dragend(d, i) {
			        d.fixed = true;
			        tick({type: "tick", alpha: 0.05});
			        force.resume();
			    }

			    node.call(node_drag);

	            var circle = node.append("circle")
	                .attr("class", "node")
	                .attr("r", function(d) {return d.radius = (d.type >= 4 ? 0 : d.r); })
	                .style("fill", function(d) {return d.type >= 4 ? "white" : col(d.type); }) 
	                .style("stroke", function(d) {return d.TG == 0 ? "white" : c10(d.TG-1); })
	                .style("opacity", 0.9)
				    .style("stroke-width", 2);

	            circle.append("title")
	                .text(function(d) { return d.Name; });  
    
	            var maxRadius = d3.max(tnodes, function(d) {return d.depth == 0 ? 0 : d.r});
				
	            // Scales for node position: 
	            // y direction corresponds to trophic level
	            // x is arbitrarty, not currently used
                
				if (isNaN(tllim[0])) {
	            	tllim = d3.extent(tnodes, function(x) {return isNaN(x.TL) ? NaN : x.TL});
				}
                
	            var tlscale = d3.scale.linear()
	                .range([height, 0])
	                .domain(tllim); 
                
	            var xscale = d3.scale.linear()
	                .range([0, width])
	                .domain([0, 1]);
					
	            // Set up node/edge repositioners
            
	            var updateLink = function() {
	                this.attr("x1", function(d) {
	                    return d.source.x;
	                }).attr("y1", function(d) {
	                    return d.source.y;
	                }).attr("x2", function(d) {
	                    return d.target.x;
	                }).attr("y2", function(d) {
	                    return d.target.y;
	                });

	            }
				
	            var updateNode = function() {
	                this.attr("transform", function(d) {
	                    d.x = Math.max(0, Math.min(width, d.x))
	                    return "translate(" + d.x + "," + d.y + ")";
	                });

	            }
				
				var updateFLink = function () {
	                this.attr("x1", function(d) {
	                    return tnodes[d.sidx].x;
	                }).attr("y1", function(d) {
	                    return tnodes[d.sidx].y;
	                }).attr("x2", function(d) {
	                    return tnodes[d.tidx].x;
	                }).attr("y2", function(d) {
	                    return tnodes[d.tidx].y;
	                });
				}
				
				
				node.call(updateNode);
				
				
	            // While force layout runs...
                
				force.on("tick", function(e) {tick(e)});
				
				function tick(e) {
					
					// Release previously-fixed nodes (custom drags) when nearly cooled
					
					if (e.alpha < 0.005) {
						node.each(function(d) {d.fixed = false}); 
					}	
					
	                // Nudge nodes towards prescribed x/y position
                
	                var k = tlnudge * e.alpha;
                
	                node.each(function(d,i) {
	                    if (!isNaN(d.TL)) {
	                        d.y += (tlscale(d.TL) - d.y) * k;
	                    }
	                    if (!isNaN(d.xpos)) {
	                        d.x += (xscale(d.xpos) - d.x) * k;
	                    }
	                })
                
	                // Update positions of circles and hierarchical lines
                
	                node.call(updateNode);
	                hlink.call(updateLink);
					fluxlines.call(updateFLink);
					
                
	                // Collision detection
                
	                circle.each(collide(0.5))
			      
			    };				
				
				break;
			default: console.log('Unrecognized mode');
			};
			
		    // *****************
		    // Subfunctions
		    // *****************
      
		    function collide(alpha) {
		      var quadtree = d3.geom.quadtree(tnodes);
		      return function(d) {
		        var r = d.radius + maxRadius + padding,
		            nx1 = d.x - r,
		            nx2 = d.x + r,
		            ny1 = d.y - r,
		            ny2 = d.y + r;
		        quadtree.visit(function(quad, x1, y1, x2, y2) {
		          if (quad.point && (quad.point !== d)) {
		            var x = d.x - quad.point.x,
		                y = d.y - quad.point.y,
		                l = Math.sqrt(x * x + y * y),
		                r = d.radius + quad.point.radius + padding;
		            if (l < r) {
		              l = (l - r) / l * alpha;
		              d.x -= x *= l;
		              d.y -= y *= l;
		              quad.point.x += x;
		              quad.point.y += y;
		            }
		          }
		          return x1 > nx2 || x2 < nx1 || y1 > ny2 || y2 < ny1;
		        });
		      };
		    }

		}); // end selection.each

	}; // end chart
	
	
    //****************************
    // Fill in user-set variables
    //****************************
    
    chart.marl = function(value) {
        if (!arguments.length) return marl;
        marl = value;
        return chart;
    };
    chart.marr = function(value) {
        if (!arguments.length) return marr;
        marr = value;
        return chart;
    };
    chart.mart = function(value) {
        if (!arguments.length) return mart;
        mart = value;
        return chart;
    };
    chart.marb = function(value) {
        if (!arguments.length) return marb;
        marb = value;
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
    
    chart.padding = function(value) {
        if (!arguments.length) return padding;
        padding = value;
        return chart;
    };
    
    chart.tlnudge = function(value) {
        if (!arguments.length) return tlnudge;
        tlnudge = value;
        return chart;
    };
    
    chart.linkdist = function(value) {
        if (!arguments.length) return linkdist;
        linkdist = value;
        return chart;
    };
    
    chart.charge = function(value) {
        if (!arguments.length) return charge;
        charge = value;
        return chart;
    };
    
    chart.hlineop = function(value) {
        if (!arguments.length) return hlineop;
        hlineop = value;
        return chart;
    };
	
    chart.flineop = function(value) {
        if (!arguments.length) return flineop;
        flineop = value;
        return chart;
    };
    
    chart.fontsz = function(value) {
        if (!arguments.length) return fontsz;
        fontsz = value;
        return chart;
    };
	
    chart.seed = function(value) {
        if (!arguments.length) return seed;
        seed = value;
        return chart;
    };
	
    chart.cpfac = function(value) {
        if (!arguments.length) return cpfac;
        cpfac = value;
        return chart;
    };
	
    chart.drawmode = function(value) {
        if (!arguments.length) return drawmode;
        drawmode = value;
        return chart;
    };
	
    chart.tllim = function(value) {
        if (!arguments.length) return tllim;
        tllim = value;
        return chart;
    };
	
    chart.init = function(value) {
        if (!arguments.length) return init;
        init = value;
        return chart;
    };
	
    return chart;
}; // end foodweblayout