<!DOCTYPE html>
<html>

  <head>
    <meta name="description" content="math.js | basic usage">
    <title>math.js | basic usage</title>
    <script src="https://unpkg.com/mathjs/lib/browser/math.js"></script>
    
    <script src="https://d3js.org/d3.v4.min.js"></script>
    <script src="https://d3js.org/d3-hsv.v0.1.min.js"></script>
    <script src="https://d3js.org/d3-contour.v1.min.js"></script>
    <script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>
    <style>

      /* Make the chart container fill the page using CSS. */
      #chart {
        position: fixed;
        left: 0px;
        right: 0px;
        top: 0px;
        bottom: 0px;
      }
    </style>
  </head>

<!-- Create a div where the graph will take place -->
<div id="my_dataviz">
  <svg id="click" xmlns="http://www.w3.org/2000/svg">
      <defs>
          <g id="pointer" transform="scale(0.5)">
              <circle cx="0" cy="0" r="20" id="dragcircle" />
          </g>
      </defs>
  </svg>
</div>



  <body style='overflow:hidden'>
  
    
    <script>
    
      const vw = Math.max(document.documentElement.clientWidth || 0, window.innerWidth || 0)
      const vh = Math.max(document.documentElement.clientHeight || 0, window.innerHeight || 0)
      
      var height 	= math.min([vw,vh]); // 450
      var width 	= math.min([vw,vh]);
    
      var p1 				= [0.25, 0.25];
      var p2 				= [0.75, 0.75];
      
      var w1 				= [0.6, 0.4];
      var rw1 			= 0.01*height/450; // Well screen radius
      var sw1 			= -1; // Well strength
      
      var w2 				= [0.4, 0.6];
      var rw2 			= 0.01*height/450; // Well screen radius
      var sw2 			= 1; // Well strength
      
      var vres 			= 11;
      
      var vertices 	= linspace_complex(p1, p2, vres);
      var zc 			 	= linspace_complex_segments(p1, p2, vres);
      var nvec 			= normal_vector(p1,p2);
      var A 				= noflowboundary_gradients(vertices,zc,nvec);
      var solvec 		= solution_vector(zc,nvec);
      var s 				= solve_linear_system(A,solvec);
      var x_res			= 51;
      var y_res 		= 51;
      var grid 			= construct_grid([0,1],[0,1],x_res,y_res);
      var output 		= evaluate(grid, vertices, s);


      


      // helper function to output formatted results.
      function print(value) {
        var precision = 14;
        document.write(math.format(value, precision) + '<br>');
      }

      // This function creates 11 points equally spaced between a tuple p1 and a tuple p2
      function linspace(start,end,resolution) {
        var spacing = [];
        // Go through a for-loop
        var i;
        for (i = 0; i < resolution; i++) {
        	spacing.push(start + (end-start)*i/(resolution-1))
        }
        return spacing; // The function returns the linspace between p1 and p2
      }

      // This function creates 11 points equally spaced between a tuple p1 and a tuple p2
      function linspace_complex(p1, p2, resolution) {
        var spacing = [];
        // Go through a for-loop
        var i;
        for (i = 0; i < resolution; i++) {
        	spacing.push(math.complex(
          	p1[0] + (p2[0] -p1[0] )*i/(resolution-1),
            p1[1]  + (p2[1]-p1[1])*i/(resolution-1)))
        }
        return spacing; // The function returns the linspace between p1 and p2
      }

      // This function creates 11 points equally spaced between a tuple p1 and a tuple p2
      function linspace_complex_segments(p1, p2, resolution) {
        var spacing = [];
        // Go through a for-loop
        xdif = (p2[0]-p1[0])/(resolution-1)/2
        ydif = (p2[1]-p1[1])/(resolution-1)/2
        var i;
        for (i = 0; i < resolution-1; i++) {
        	spacing.push(math.complex(
          	p1[0] + (p2[0] -p1[0] )*i/(resolution-1) + xdif,
            p1[1]  + (p2[1]-p1[1])*i/(resolution-1) + ydif))
        }
        return spacing; // The function returns the linspace between p1 and p2
      }
      
			// =============================================================================
      // EVALUATE THE MODEL
      // =============================================================================
      
      function evaluate(x, vertices, s) {
      
      	// Segments are between the vertices, so we have one less than vertices
        var segments = vertices.length-1
        
        // Go through a for-loop
        var i;
        for (i = 0; i < segments; i++) {
        
        	// Project the control points on the real line between -1 and 1
          Z = cdivide(
          	math.subtract(
            	math.dotMultiply(x, 2),
              cadd(vertices[i],vertices[i+1])),
            csub(vertices[i+1],vertices[i]));
              
          // Calculate term 1
          term1 	= math.dotMultiply(
          	math.add(
            	Z,
              1),
           	math.log(
            	math.dotDivide(
                math.subtract(
                  Z,
                  1),
                math.add(
                  Z,
                  1)) ) )
          
          // Calculate term 2
          term2 	= math.dotMultiply(
          	math.subtract(
            	Z,
              1),
           	math.log(
            	math.dotDivide(
                math.subtract(
                  Z,
                  1),
                math.add(
                  Z,
                  1))))
          
          numer = math.dotMultiply(
              s[i][0],
              math.subtract(
                  term1,
                  term2))
          
          denom = math.dotMultiply(
          	4*math.PI,
            math.complex(0,1) )
          
          temp 	= math.dotDivide(
          	numer,
            denom)

          // Add the results to the output vector
          if (i == 0) {
          	// In the first iteration, initialize the output vector
          	var results = temp
          } else {
          	results 		 = math.add(
            	results,
              temp)
          }
        }
        
        // Add the effect of the well
        
        //# Find the indices of wells      
        //dist            = np.abs(z-self.zc)
        //idx_inside      = np.where(dist < self.rw)[0]
        
        //# Correct the coordinates ---------------------------------------------
        //# Center the evaluation points on the well
        //zs  = z.copy()-self.zc
        
        //# Snap points inside the well to the well edge
        //zs[idx_inside] = self.rw + 0j    
        
        //# Calculate the complex potential
        //res = -self.strength/(2*np.pi)*np.log(zs) - temp
        
				
        // Get the local coordinates
        zs 	= math.subtract(
          x,
          math.complex(w1[0],w1[1]))
          
        //print(zs)

        // Truncate any coordinates inside the well screen
        for (j = 0; j < zs.length; j++) {
          radius = math.abs(zs[j])
          if (radius < rw1) {
            zs[j] = math.dotMultiply(zs[j],rw1/radius);
          }}

        // Calculate the gradient denominator
        temp = math.dotDivide(
        	math.dotMultiply(
          	sw1,
            math.log(zs)),
         	math.dotMultiply(
          	2,
            math.PI))
        

        // Add this to the stack (actually subtract a negative value)
        results = math.subtract(
          results,
          temp)
          
          
          
        // Get the local coordinates
        zs 	= math.subtract(
          x,
          math.complex(w2[0],w2[1]))
          
        //print(zs)

        // Truncate any coordinates inside the well screen
        for (j = 0; j < zs.length; j++) {
          radius = math.abs(zs[j])
          if (radius < rw2) {
            zs[j] = math.dotMultiply(zs[j],rw2/radius);
          }}

        // Calculate the gradient denominator
        temp = math.dotDivide(
        	math.dotMultiply(
          	sw2,
            math.log(zs)),
         	math.dotMultiply(
          	2,
            math.PI))
        

        // Add this to the stack (actually subtract a negative value)
        results = math.subtract(
          results,
          temp)
          
          
          
          
          
          
          

				
        // Also add the uniform flow
        results 		 = math.add(
          results,
          x)
          
        return creal(results);
      }
     
      
			// =============================================================================
      // CONSTRUCT GRID
      // =============================================================================
      
      function construct_grid(x_limits,y_limits,x_res,y_res) {
      
      	// This function creates a rectangular grid
        var grid 	= [];
        var ydif  = y_limits[1]-y_limits[0];
        var xdif 	= x_limits[1]-x_limits[0];
        
        // Go through a for-loop
        var i;
        var j;
        for (j = 0; j < y_res; j++) {
          for (i = 0; i < x_res; i++) {
						grid.push(math.complex(
            	x_limits[0] + xdif*i/(x_res-1)+0.0000000001,
              y_limits[0] + ydif*j/(y_res-1)+0.0000000001));
          }
        }
        
				//print(grid)
        return grid; 
      }
      
      
			// =============================================================================
      // SOLVE SYSTEM OF EQUATIONS
      // =============================================================================
      
      function solve_linear_system(A,c) {
      
      	// In this case, the solution vector is really simple: Just standard uniform
        // flow from east to west with a gradient of one. As such:

        return math.lusolve(A, c); 
      }
      
      
			// =============================================================================
      // FIND NORMAL VECTOR
      // =============================================================================
      
      function normal_vector(p1,p2) {
      
      	// Define the vector between the endpoints, get its orthogonal part
      	nvec 	= [
        	p1[1] - p2[1],
          p2[0] - p1[0]];

        // Normalize this vector
        nvec = math.dotDivide(
        	nvec,
          math.sqrt(math.pow(nvec[0],2) + math.pow(nvec[1],2)));

        return nvec;
      }
      
			// =============================================================================
      // SOLUTION VECTOR
      // =============================================================================
      
      function solution_vector(zc) {
      
      	// In this case, the solution vector is really simple: Just standard uniform
        // flow from east to west with a gradient of one. As such:
        gradient 	= 1;
        solvec 		= [];
        
        // Go through a for-loop
        var i;
        for (i = 0; i < zc.length; i++) {
          
          // Initialize the background gradient
          temp 	= math.complex(-gradient,0);
          
          // Add the well gradient
          // grad[idx_valid] = -self.strength/(2*np.pi)/zs[idx_valid]
          
          // Get the local coordinates
          zs 	= math.subtract(
            zc[i],
            math.complex(w1[0],w1[1]))
            
          // Truncate any coordinates inside the well screen
          radius = math.abs(zs)
          if (radius < rw1) {
            zs 	= math.dotMultiply(
            	zs,rw1/radius);
          }
          
          // Calculate the gradient denominator
          denom = math.dotMultiply(
          	2,
            math.dotMultiply(
            	math.PI,
              zs))
          
          // Calculate the gradient numerator
          numer = sw1
          
          // Calculate the well's induced complex potential gradient
          well_effect = math.dotDivide(
            	numer,
              denom)
              
          // Convert it to the partial derivatives of the hydraulic potential
          well_effect = math.conj(well_effect)
              
          // Add this to the stack (actually subtract a negative value)
          temp 	= math.add(
          	temp,
            well_effect)
            
            
            
            
          // Get the local coordinates
          zs 	= math.subtract(
            zc[i],
            math.complex(w2[0],w2[1]))
            
          // Truncate any coordinates inside the well screen
          radius = math.abs(zs)
          if (radius < rw2) {
            zs 	= math.dotMultiply(
            	zs,rw2/radius);
          }
          
          // Calculate the gradient denominator
          denom = math.dotMultiply(
          	2,
            math.dotMultiply(
            	math.PI,
              zs))
          
          // Calculate the gradient numerator
          numer = sw2
          
          // Calculate the well's induced complex potential gradient
          well_effect = math.dotDivide(
            	numer,
              denom)
              
          // Convert it to the partial derivatives of the hydraulic potential
          well_effect = math.conj(well_effect)
              
          // Add this to the stack (actually subtract a negative value)
          temp 	= math.add(
          	temp,
            well_effect)
            
            
            
            
            
            

          
          // Now calculate the inner product between the partial derivatives 
          // and the normal vector
          temp = math.add(
            math.dotMultiply(
              nvec[0],
              temp.re),
            math.dotMultiply(
              nvec[1],
              temp.im));
              
          
          
          // Append the results to the matrix
					solvec.push(temp);
          
        } 
        
        return solvec; // The function returns the linspace between p1 and p2
      }

			// =============================================================================
      // INFLUENCE MATRIX
      // =============================================================================
      
      function noflowboundary_gradients(vertices,zc, nvec) {
      
      	// Segments are between the vertices, so we have one less than vertices
        var segments = vertices.length-1

				// Initialize the output array
        var results = [];
        
        // Go through a for-loop
        var i;
        for (i = 0; i < segments; i++) {
        
        	// Project the control points on the real line between -1 and 1
          var Z = cdivide(
          	math.subtract(
            	math.dotMultiply(zc, 2),
              cadd(vertices[i],vertices[i+1])),
            csub(vertices[i+1],vertices[i]));

          // Convert to dOmega(Z)/dZ
          temp = math.dotDivide(
            math.complex(0,1),
            math.subtract(
            	math.PI,
              math.dotMultiply(
              	math.PI,
                cpow(
                	Z,
                  2))));
                  
          // Multiply with dZ/dz to obtain dOmega(Z)/dz
          temp = math.dotDivide(
            	math.dotMultiply(
          			temp,
                2),
              csub(vertices[i+1],vertices[i]) )
          
					// We require the partial derivatives of the real component, get the conjugate
          temp = math.conj(temp)
          
          // Now calculate the inner product between the partial derivatives 
          // and the normal vector
          temp = math.add(
          	math.dotMultiply(
            	nvec[0],
              creal(temp)),
            math.dotMultiply(
            	nvec[1],
              cimag(temp)));
              
          // Append the results to the matrix
          results.push(temp);
        } 

        return results; // Take its transpose: math.transpose(results)
      }
      
			// =============================================================================
      // FUNCTION CMULTIPLY
      // =============================================================================
      
      function cmultiply(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        
        var c1_len = c1.length;
        if (c1_len == undefined) {c1_len = 1;}
        
        var c2_len = c2.length;
        if (c2_len == undefined) {c2_len = 1;}
        
      	if (c1_len == 1 && c2_len == 1) { // Both variables are scalars
          var result 		= math.complex(
          	num1.re*num2.re-num1.im*num2.im, 
            num1.re*num2.im+num1.im*num2.re);
        } else if (c1_len > 1 && c2_len > 1) { // Both variables are vectors
        	var result 		= []
          for (i = 0; i < num1.length; i++) {
            result.push(math.complex(
              num1[i].re*num2[i].re-num1[i].im*num2[i].im, 
              num1[i].re*num2[i].im+num1[i].im*num2[i].re));
          }
        } else if (c1_len > 1 && c2_len == 1) { // The first variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num1.length; i++) {
            result.push(math.complex(
              num1[i].re*num2.re-num1[i].im*num2.im, 
              num1[i].re*num2.im+num1[i].im*num2.re));
          } 
        } else if (c1_len == 1 && c2_len > 1) { // The second variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num2.length; i++) {
            result.push(math.complex(
              num1.re*num2[i].re-num1.im*num2[i].im, 
              num1.re*num2[i].im+num1.im*num2[i].re));
          } 
        } else {
          var result 		= False;
          }

      return result;   
      }
      
			// =============================================================================
      // FUNCTION CDIVIDE
      // =============================================================================
      
      function cdivide(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        
        var c1_len = c1.length;
        if (c1_len == undefined) {c1_len = 1;}
        
        var c2_len = c2.length;
        if (c2_len == undefined) {c2_len = 1;}
        
      	if (c1_len == 1 && c2_len == 1) { // Both variables are scalars
          var denom 		= num2.im * num2.im + num2.re * num2.re;
          var real 			= (num1.re * num2.re + num1.im * num2.im) /denom;
          var imaginary = (num2.re * num1.im - num1.re * num2.im) /denom; 
          var result 		= math.complex(real, imaginary);
        } else if (c1_len > 1 && c2_len > 1) { // Both variables are vectors
        	var result 		= []
          for (i = 0; i < num1.length; i++) {
            var denom 		= num2[i].im * num2[i].im + num2[i].re * num2[i].re;
            var real 			= (num1[i].re * num2[i].re + num1[i].im * num2[i].im) /denom;
            var imaginary = (num2[i].re * num1[i].im - num1[i].re * num2[i].im) /denom; 
            result.push(math.complex(real, imaginary));
          }
        } else if (c1_len > 1 && c2_len == 1) { // The first variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num1.length; i++) {
            var denom 		= num2.im * num2.im + num2.re * num2.re;
            var real 			= (num1[i].re * num2.re + num1[i].im * num2.im) /denom;
            var imaginary = (num2.re * num1[i].im - num1[i].re * num2.im) /denom; 
            result.push(math.complex(real, imaginary));
          } 
        } else if (c1_len == 1 && c2_len > 1) { // The second variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num2.length; i++) {
            var denom 		= num2[i].im * num2[i].im + num2[i].re * num2[i].re;
            var real 			= (num1.re * num2[i].re + num1.im * num2[i].im) /denom;
            var imaginary = (num2[i].re * num1.im - num1.re * num2[i].im) /denom; 
            result.push(math.complex(real, imaginary));
          } 
        } else {
          var result 		= False;
          }

      return result;   
      }
      
      function cadd(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        var real = (num1.re + num2.re);
        var imaginary = (num1.im + num2.im); 
      return math.complex(real, imaginary);   
      }
      
      function csub(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        var real = (num1.re - num2.re);
        var imaginary = (num1.im - num2.im); 
      return math.complex(real, imaginary);   
      }
      
      function cpow(c, exp) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push(math.pow(c[i],exp))
        }
      return res;   
      }
      
      function creal(c) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push(c[i].re)
        }
      return res;   
      }
      
      function cimag(c) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push(c[i].im)
        }
      return res;   
      }
      
      function rvec_to_cvec(c) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push([c[i]])
        }
      return res;   
      }
      
      
			// =============================================================================
      // UPDATE MODEL AND FIGURE
      // =============================================================================
      
      function update() {
      
      	p1 				= [d3.select("#point1").attr("x")/width,d3.select("#point1").attr("y")/height];
        p2 				= [d3.select("#point2").attr("x")/width,d3.select("#point2").attr("y")/height];
        w1 				= [d3.select("#well1").attr("x")/width,d3.select("#well1").attr("y")/height];
        w2 				= [d3.select("#well2").attr("x")/width,d3.select("#well2").attr("y")/height];
        vertices 	= linspace_complex(p1, p2, vres);
        zc 			 	= linspace_complex_segments(p1, p2, vres);
        nvec 			= normal_vector(p1,p2);
        A 				= noflowboundary_gradients(vertices,zc,nvec);
        solvec 		= solution_vector(zc,nvec);
        s 				= solve_linear_system(A,solvec);
        output 		= evaluate(grid, vertices, s);
        
        
        // Remove the previous elements
        d3.select("#line_element").remove();
        svg.selectAll("path").remove();


        // array of threshold values 
        
        //thresholds = d3.range(
        //d3.min(output) - (d3.max(output)- d3.min(output))/21,
        //d3.max(output) + (d3.max(output)- d3.min(output))/21*2,
        //(d3.max(output)- d3.min(output))/21);
        
        output_dif = math.max(output)-math.min(output)
        thresholds 	= linspace(
        	math.min(output),
          math.max(output),
          21)

				color = d3.scaleLinear()
          .domain(d3.extent(thresholds))
          .interpolate(function() { return d3.interpolateRgbBasis(["#152d3b","#295c79","#4794c1","#c3e7f9"])});  

        // initialise contours
        contours = d3.contours()
            .size([y_res, x_res])
            .thresholds(thresholds)
            (output);

        // make and project the contours
        svg.selectAll("path")
            .data(contours)
            .enter().append("path")
                .attr("d", d3.geoPath(projection))
                .attr("fill", function(d) { return color(d.value); })

        svg
            .append("line")
            .attr("x1",p1[0]*width)
            .attr("x2",p2[0]*width)
            .attr("y1",p1[1]*height)
            .attr("y2",p2[1]*height)
            .attr("stroke", "#666666")
            .attr("stroke-width", 10)
            .attr("stroke-linecap","round")
            .attr("id","line_element");
                    
        point1.raise()
        point2.raise()
        well1.raise()
        well2.raise()
        
      return;   
      }
      

			// =============================================================================
      // PLOT THE CONTOURS
      // =============================================================================

      // set x and y scale to maintain 1:1 aspect ratio  
      // Extract the width and height that was computed by CSS.

      
      var scaling 	= math.min([width/x_res,height/y_res]);
      //print(scaling)
      
      var projection = d3.geoTransform({
          point: function(px, py) {
              this.stream.point(px*scaling, py*scaling);
          }
      });
      
			
        
      var svg = d3.select("#click") // This selects the div
          .attr("width", width) // This defines the canvas' width
          .attr("height", height) // This defines the canvas' height
      
      // array of threshold values 
      var thresholds = d3.range(
      	d3.min(output),
        d3.max(output),
        (d3.max(output)- d3.min(output))/21);

      // color scale  
      var color = d3.scaleLinear()
          .domain(d3.extent(thresholds))
          .interpolate(function() { return d3.interpolateRgbBasis(["#152d3b","#295c79","#4794c1","#c3e7f9"])});  
          //.interpolate(function() { return d3.interpolateRdBu; });  
          
      // initialise contours
      var contours = d3.contours()
          .size([y_res, x_res])
          .thresholds(thresholds)
          (output);

      // Resize the drag circles
      d3.select("#dragcircle").attr("r",20*height/450)
      
      // make and project the contours
      svg.selectAll("path")
          .data(contours)
          .enter().append("path")
              .attr("d", d3.geoPath(projection))
              .attr("fill", function(d) { return color(d.value); })
              
			svg
      		.append("line")
          .attr("x1",p1[0]*width)
          .attr("x2",p2[0]*width)
          .attr("y1",p1[1]*height)
          .attr("y2",p2[1]*height)
          .attr("stroke", "#666666")
          .attr("stroke-width", 10*height/450)
          .attr("stroke-linecap","round")
          .attr("id","line_element");
          

          
			// =============================================================================
      // CREATE THE POINTER
      // =============================================================================
          
      var point1 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", p1[0]*width)
          .attr("y", p1[1]*height)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#666666")
          .attr("stroke-width", 5*height/450)
          .attr("id","point1");
          
      var point2 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", p2[0]*width)
          .attr("y", p2[1]*height)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#666666")
          .attr("stroke-width", 5*height/450)
          .attr("id","point2");
          
      var well1 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", w1[0]*width)
          .attr("y", w1[1]*height)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#666666")
          .attr("stroke-width", 5*height/450)
          .attr("id","well1");
          
      var well2 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", w2[0]*width)
          .attr("y", w2[1]*height)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#666666")
          .attr("stroke-width", 5*height/450)
          .attr("id","well2");
          
          //.attr("transform","scale("+string(0.5*height/450)+")");

      var deltaX, deltaY;

      var dragHandler = d3.drag()
          .on("start", function () {
              var current = d3.select(this);
              deltaX = current.attr("x") - d3.event.x;
              deltaY = current.attr("y") - d3.event.y;
          })
          .on("drag", function () {
          		var movex = d3.event.x + deltaX;
              var movey = d3.event.y + deltaY;
              if (movex < 0.05*width) {
              	movex = 0.05*width
              } else if (movex > width*0.95) {
              	movex = width*0.95
              }
              if (movey < 0.05*height) {
              	movey = 0.05*height
              } else if (movey > height*0.95) {
              	movey = height*0.95
              }
              d3.select(this)
                  .attr("x", movex)
                  .attr("y", movey);
          })
          .on("end", function () {
              update();
          });


      dragHandler(svg.selectAll("use"));




			
      // ===========================================================================
      // This function shifts the contours
      // ===========================================================================
      
      var increment = 0;
      
      function shift_contours() {


				threshdif = thresholds[1]-thresholds[0];
				increment += threshdif/10
        
        
        if (increment >= threshdif){
        	increment = 0;
        }

        // initialise contours
        contours = d3.contours()
            .size([y_res, x_res])
            .thresholds(math.subtract(thresholds,increment))
            (output);

        
				svg.selectAll("path").remove();
        
        // make and project the contours
        svg.selectAll("path")
            .data(contours)
            .enter().append("path")
                .attr("d", d3.geoPath(projection))
                .attr("fill", function(d) { return color(d.value); })

				d3.select("#line_element").raise()
        point1.raise()
        point2.raise()
				well1.raise()
        well2.raise()
				return }


  		// Make the contour flow
      var repater = setInterval(function() {
        shift_contours();
      }, 25);
      
  
  
  
    </script>
  </body>

</html>
