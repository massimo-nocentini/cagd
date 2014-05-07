
require("deCasteljau.jl")

function exercise_one ()
    controlPoints = [
                     [1.0 1 0];
                     [3.0 4 2];
                     [5.0 6 1];
                     [7.0 8 9];
                     [10.0 2 8];
                     [1.0 1 0];
                     ]
    
    bezierCurvePoints = drawCurve(controlPoints,
				  linspace(0,1,200))

    writeArrayForGnuplot(bezierCurvePoints, "bezier-exercise-one.coordinates")
end

function exercise_three()
    x = t -> 1 + t + t^2
    y = t -> t^3
    points = zeros(200, 2)
    params = linspace(0,1,200)
    for i=1:200
        t = params[i]
        points[i,:] =
            [x(t) y(t)]
    end
    writeArrayForGnuplot(points, "test-exercise-three.coordinates")
end


function exercise_two ()

    maximaControlPoints = [
                           [1.0 0];
                           [4/3 0];
                           [2.0 0];
                           [3.0 1];
                           ]
    
    maximaBezierCurvePoints = drawCurve(maximaControlPoints,
				  linspace(0,1,200))

    writeArrayForGnuplot(maximaControlPoints,
                         "control-poly-exercise-three-maxima.coordinates")

    writeArrayForGnuplot(maximaBezierCurvePoints,
                         "bezier-exercise-three-maxima.coordinates")
    
    controlPoints = parametricCoordinates([
			                   x -> 1 + x + x^2;
			                   y -> y^3
		                           ], 200)

    ## bezierCurvePoints = drawCurve(controlPoints,
    ##     			  linspace(0,1,200))

    ## writeArrayForGnuplot(bezierCurvePoints, "bezier-exercise-two.coordinates")
    writeArrayForGnuplot(controlPoints, "parametric-curve-exercise-three.coordinates")

    controlPoints = parametricCoordinates([
        		                   x -> 1 + x + x^2;
        		                   y -> y^3;
        		                   z -> z^2 - pi/2
        	                           ], 6)

    bezierCurvePoints = drawCurve(controlPoints,
        			  linspace(0,1,200))

    writeArrayForGnuplot(bezierCurvePoints, "bezier-three-axes-exercise-two.coordinates")
    writeArrayForGnuplot(controlPoints, "control-three-axes-poly-exercise-two.coordinates")
end

function exercise_four ()
    outerControlPoints = [
                          [0.0 0];
                          [1.0 3];
                          [2.5 4];
                          [5 4.5];
                          [4.5 2.5];
                          [6.0 1];
                          ]
    writeArrayForGnuplot(outerControlPoints, "control-poly-exercise-four.coordinates")
    t_star = .25
    q_value, Diagonal, Last = decasteljau(outerControlPoints, t_star)
    writeArrayForGnuplot(Diagonal, "cyan-control-poly-exercise-four.coordinates")
    writeArrayForGnuplot(Last, "magenta-control-poly-exercise-four.coordinates")

    #diagonalBezierCurvePoints = drawCurve(Diagonal, [t/t_star for t in linspace(0,1,100)])
    diagonalBezierCurvePoints = drawCurve(Diagonal, linspace(0,1,100))
    lastBezierCurvePoints = drawCurve(Last, linspace(0,1,100))
    writeArrayForGnuplot(diagonalBezierCurvePoints, "bezier-diagonal-exercise-four.coordinates")
    writeArrayForGnuplot(lastBezierCurvePoints, "bezier-last-exercise-four.coordinates")
end

function exercise_six()
    originalControlPoints = [
                             [1.0 1];
                             [3.0 4];
                             [5.0 6];
                             [7.0 1.5];
                             [10 2]
                             ]
    writeArrayForGnuplot(originalControlPoints, "control-poly-exercise-six.coordinates")
    originalBezierCurvePoints = drawCurve(originalControlPoints, linspace(0,1,200))
    writeArrayForGnuplot(originalBezierCurvePoints, "bezier-original-exercise-six.coordinates")

    oneMoreDegreeControlPoints = increaseDegree(originalControlPoints)
    oneMoreDegreeBezierCurvePoints = drawCurve(oneMoreDegreeControlPoints, linspace(0,1,200))
    writeArrayForGnuplot(oneMoreDegreeControlPoints, "control-poly-exercise-six-one-more-degree.coordinates")
    writeArrayForGnuplot(oneMoreDegreeBezierCurvePoints, "bezier-one-more-degree-exercise-six.coordinates")

    twoMoreDegreeControlPoints = increaseDegree(oneMoreDegreeControlPoints)
    twoMoreDegreeBezierCurvePoints = drawCurve(twoMoreDegreeControlPoints, linspace(0,1,200))
    writeArrayForGnuplot(twoMoreDegreeControlPoints, "control-poly-exercise-six-two-more-degree.coordinates")
    writeArrayForGnuplot(twoMoreDegreeBezierCurvePoints, "bezier-two-more-degree-exercise-six.coordinates")

    threeMoreDegreeControlPoints = increaseDegree(twoMoreDegreeControlPoints)
    threeMoreDegreeBezierCurvePoints = drawCurve(threeMoreDegreeControlPoints, linspace(0,1,200))
    writeArrayForGnuplot(threeMoreDegreeControlPoints, "control-poly-exercise-six-three-more-degree.coordinates")
    writeArrayForGnuplot(threeMoreDegreeBezierCurvePoints, "bezier-three-more-degree-exercise-six.coordinates")
end

function exercise_five()
    PointToRepeat = Float64[10 1]

    originalControlPoints = [ 	
		             [2.0 4]; 
		             [6 12]; 
		             PointToRepeat; 
		             [12 12] 
	                     ]

    repeatedControlPoints = [ 
		             [2.0 4]; 
		             [6 12]; 
		             PointToRepeat; 
		             PointToRepeat; 
		             PointToRepeat; 
		             PointToRepeat; 
		             [12 12] 
	                     ]

    originalBezierCurvePoints = drawCurve(originalControlPoints, linspace(0,1,200))
    repeatedBezierCurvePoints = drawCurve(repeatedControlPoints, linspace(0,1,200))
    writeArrayForGnuplot(originalBezierCurvePoints, "bezier-original-exercise-five.coordinates")
    writeArrayForGnuplot(repeatedBezierCurvePoints, "bezier-repeated-exercise-five.coordinates")
end

function exercise_seven()

    controlPoints = [
		     [0.0 0 0];
		     [1.0 2 0];
		     [3.0 2 0];
		     [6.0 -1 0];
	             ]

    continuityPoint, tangentPoint, obsculatorPoint, a_succ_i =
        junction_right(controlPoints, 1, 1)

    restControlPoints = [
                         [15.0 3 0];
		         [15.0 3 0];
		         [15.0 3 0];
		         [15.0 3 0];
		         [23.0 -4 0];
		         [24.0 0 0];
	                 ]

    continuityControlPoints = [
                               continuityPoint;
                               restControlPoints;
                               ]

    tangentControlPoints = [
                            continuityPoint;
                            tangentPoint;
                            restControlPoints;
                            ]

    obsculatingControlPoints = [
                                continuityPoint;
                                tangentPoint;
                                obsculatorPoint;
                                restControlPoints;
                                ]
    
    baseBezierCurvePoints = drawCurve(controlPoints,
                                      linspace(0,1,200))
    continuityBezierCurvePoints = drawCurve(continuityControlPoints,
                                            linspace(0,1,200))
    tangentBezierCurvePoints = drawCurve(tangentControlPoints,
                                         linspace(0,1,200))
    obsculatingBezierCurvePoints = drawCurve(obsculatingControlPoints,
                                             linspace(0,1,200))

    writeArrayForGnuplot(controlPoints,
                         "control-poly-base-exercise-seven-exercise-seven.coordinates")
    writeArrayForGnuplot(continuityControlPoints, 
		         "control-poly-continuity-exercise-seven.coordinates")
    writeArrayForGnuplot(tangentControlPoints, 
		         "control-poly-tangent-exercise-seven.coordinates")
    writeArrayForGnuplot(obsculatingControlPoints, 
		         "control-poly-obsculating-exercise-seven.coordinates")

    writeArrayForGnuplot(baseBezierCurvePoints, "bezier-base-exercise-seven.coordinates")
    writeArrayForGnuplot(continuityBezierCurvePoints, "bezier-continuity-exercise-seven.coordinates")
    writeArrayForGnuplot(tangentBezierCurvePoints, "bezier-tangent-exercise-seven.coordinates")
    writeArrayForGnuplot(obsculatingBezierCurvePoints, "bezier-obsculating-exercise-seven.coordinates")

    writeArrayForGnuplot([
                          controlPoints;
                          obsculatingControlPoints;
                          ], 
		         "control-poly-base-obsculating-exercise-seven.coordinates")

    writeArrayForGnuplot([
                          controlPoints[end-1,:];
                          a_succ_i;                      
                          obsculatingControlPoints[2,:];
                          ], 
		         "control-poly-base-obsculating-with-a-succ-i-exercise-seven.coordinates")


end

function exercise_seven_left()

    controlPoints = [
		     [10.0 12 0];
		     [13.0 8.5 0];
		     [16.0 -2 0];
		     [26.0 0 0];
	             ]

    continuityPoint, tangentPoint, obsculatorPoint, a_i =
        junction_left(controlPoints, 1, 1)

    restControlPoints = [
		         [0.0 3 0];
		         [1.0 -4 0];
		         [5.0 0 0];
	                 ]

    continuityControlPoints = [
                               restControlPoints;
                               continuityPoint;
                               ]

    tangentControlPoints = [
                            restControlPoints;
                            tangentPoint;
                            continuityPoint;
                            ]

    obsculatingControlPoints = [
                                restControlPoints;
                                obsculatorPoint;
                                tangentPoint;
                                continuityPoint;
                                ]
    
    baseBezierCurvePoints = drawCurve(controlPoints,
                                      linspace(0,1,200))
    continuityBezierCurvePoints = drawCurve(continuityControlPoints,
                                            linspace(0,1,200))
    tangentBezierCurvePoints = drawCurve(tangentControlPoints,
                                         linspace(0,1,200))
    obsculatingBezierCurvePoints = drawCurve(obsculatingControlPoints,
                                             linspace(0,1,200))

    writeArrayForGnuplot(controlPoints,
                         "control-poly-base-left-exercise-seven-exercise-seven.coordinates")
    writeArrayForGnuplot(continuityControlPoints, 
		         "control-poly-continuity-left-exercise-seven.coordinates")
    writeArrayForGnuplot(tangentControlPoints, 
		         "control-poly-tangent-left-exercise-seven.coordinates")
    writeArrayForGnuplot(obsculatingControlPoints, 
		         "control-poly-obsculating-left-exercise-seven.coordinates")

    writeArrayForGnuplot(baseBezierCurvePoints, "bezier-base-left-exercise-seven.coordinates")
    writeArrayForGnuplot(continuityBezierCurvePoints, "bezier-continuity-left-exercise-seven.coordinates")
    writeArrayForGnuplot(tangentBezierCurvePoints, "bezier-tangent-left-exercise-seven.coordinates")
    writeArrayForGnuplot(obsculatingBezierCurvePoints, "bezier-obsculating-left-exercise-seven.coordinates")

    writeArrayForGnuplot([
                          obsculatingControlPoints;
                          controlPoints;
                          ], 
		         "control-poly-base-obsculating-left-exercise-seven.coordinates")

    writeArrayForGnuplot([
                          obsculatingControlPoints[end-1,:];
                          a_i;
                          controlPoints[2,:];                          
                          ], 
		         "control-poly-base-obsculating-with-a-i-left-exercise-seven.coordinates")


end

function exercise_derivative()
    controlPoints = [
                     [1.0 1 0];
                     [3.0 4 2];
                     [5.0 6 1];
                     [7.0 8 9];
                     [10.0 2 8];
                     [1.0 1 0]
                     ]

    params = 200
    differencePoints = difference(controlPoints)
    derivativePoints = zeros(params, 3)

    j = 1
    for i=linspace(0,1,200)
        derivativePoints[j,:] = derivative(differencePoints,i)
        j += 1
    end

    writeArrayForGnuplot(drawCurve([
                                    [0.0 0 0];
                                    differencePoints;
                                    [0.0 0 0];
                                    ],
                                   linspace(0,1,200)),
		         "differences-of-exercise-one.coordinates")
    
    writeArrayForGnuplot(derivativePoints,
		         "derivative-of-exercise-one.coordinates")

end

function exercise_polar()
    controlPoints = [
                     [0.0 0.0];
                     [0.0 4.0];
                     [2.0 4.0];
                     [1.0 2.0];
                     ]

    n = length(controlPoints[:,1])
    params = 200
    differencePoints = difference(controlPoints)
    derivativePoints = zeros(params, 2)
    polarPoints = zeros(params, 2)
    bezierPoints = zeros(params, 2)

    j = 1
    for i=linspace(0,1,params)
        derivativePoints[j,:] = derivative(differencePoints,i)
	aPoint, Diagonal, Last = decasteljau(controlPoints, i)
	bezierPoints[j,:] = aPoint
        ## polarPoints[j,:] = aPoint + ((0.4 - i)/n)*derivativePoints[j,:]
        polarPoints[j,:] = polar(controlPoints, i, 0.4)
        j += 1
    end

    writeArrayForGnuplot(controlPoints,
                         "control-poly-polar.coordinates")

    ## bezierCurvePoints = drawCurve(controlPoints,
    ##     			  linspace(0,1,200))

    writeArrayForGnuplot(bezierPoints, "bezier-polar.coordinates")

    
    ## writeArrayForGnuplot(drawCurve([
    ##                                 [0.0 0 0];
    ##                                 differencePoints;
    ##                                 [0.0 0 0];
    ##                                 ],
    ##                                linspace(0,1,200)),
    ##     	         "differences-of-exercise-one.coordinates")
    
    writeArrayForGnuplot(polarPoints,
		         "derivative-polar.coordinates")

end
