function drawCurve(aMatrix, aParameterVector)
    
    n = length(aParameterVector)
    x = zeros(n, length(aMatrix[1,:]))
    for i = 1:n
	param = aParameterVector[i]
	aPoint, Diagonal, Last = decasteljau(aMatrix, param)
	x[i,:] = aPoint
    end
    
    x
end

function writeArrayForGnuplot(anArray, aFilename)
    writedlm(aFilename, anArray, '\t')
    anArray
end

function parametricCoordinates(functionsVector, pointsNumber)

    parameterVector = linspace(0, 1, pointsNumber)	

    numberOfCoordinates = length(functionsVector)

    coordinates = zeros(pointsNumber, numberOfCoordinates)

    for i = 1:pointsNumber
	t = parameterVector[i]

	coord = zeros(numberOfCoordinates)'

	for j = 1:numberOfCoordinates
	    func = functionsVector[j]
	    coord[j] = func(t)
	end

	coordinates[i,:] = coord
    end

    coordinates

end

function decasteljau(controlPointsMatrix, t)
    n = length(controlPointsMatrix[:,1])
    Q = copy(controlPointsMatrix)

    Diagonal = copy(controlPointsMatrix)
    Last = copy(controlPointsMatrix)

    Diagonal[1,:] = Q[1,:]
    Last[1,:] = Q[n,:]
    
    for k = 1:n
	for i = 1:n-k
	    Q_i = Q[i,:]
	    Q_isucc = Q[i+1,:]	
	    Q[i,:] = (1 - t)*Q_i + t*Q_isucc
	end
	
	if k < n	
	    Diagonal[k+1,:] = Q[1,:]
	    Last[k+1,:] = Q[n-k,:]
	end
    end
    
    Q[1,:], Diagonal, Last
end

function increaseDegree(originalControlPoints)
    columns = length(originalControlPoints[1,:])
    sandbox = [
	       zeros(columns)'; 
	       copy(originalControlPoints); 
	       zeros(columns)'; 
	       ]
    
    increasing = [
                  copy(originalControlPoints); 
		  zeros(columns)' 
	          ]
    
    n = length(originalControlPoints[:,1])
    
    for i=1:n+1
	ratio = (i-1)/(n)
	increasing[i,:] = ratio*sandbox[i,:] + (1 - ratio)*sandbox[i+1,:] 
    end

    increasing
end

function junction_right(controlPoints, h_i, h_succ_i)

    lastControlPoint = controlPoints[end,:]
    lastButOneControlPoint = controlPoints[end - 1,:]
    lastButTwoControlPoint = controlPoints[end-2,:]
    
    intervalLengthSum = h_i + h_succ_i
    
    # in order to get continuity of degree zero
    continuityPoint = lastControlPoint
    
    # in order to get continuity of degree one
    tangentPoint = (intervalLengthSum*lastControlPoint - h_succ_i * lastButOneControlPoint) / h_i

    # in order to get continuity of degree two
    obsculatorPoint = tangentPoint - lastButOneControlPoint -
    (h_succ_i / h_i) * (lastButOneControlPoint - lastButTwoControlPoint) +
    (h_i / h_succ_i) * tangentPoint

    obsculatorPoint = (h_succ_i / h_i) * obsculatorPoint

    a_succ_i = tangentPoint + (h_i / h_succ_i) * (tangentPoint - obsculatorPoint)
    
    continuityPoint, tangentPoint, obsculatorPoint, a_succ_i
end


function junction_left(controlPoints, h_i, h_succ_i)

    firstControlPoint = controlPoints[1,:]
    secondControlPoint = controlPoints[2,:]
    thirdControlPoint = controlPoints[3,:]
    
    intervalLengthSum = h_i + h_succ_i
    
    # in order to get continuity of degree zero
    continuityPoint = firstControlPoint
    
    # in order to get continuity of degree one
    tangentPoint = (firstControlPoint - (h_i/ intervalLengthSum)*secondControlPoint) *
    (intervalLengthSum / h_succ_i)

    # in order to get continuity of degree two
    obsculatorPoint = tangentPoint *(1 + (h_succ_i / h_i)) - secondControlPoint -
    (h_i / h_succ_i) * (secondControlPoint - thirdControlPoint)

    obsculatorPoint = (h_i / h_succ_i) * obsculatorPoint

    a_i = tangentPoint + (h_succ_i / h_i) * (tangentPoint - obsculatorPoint)
    
    continuityPoint, tangentPoint, obsculatorPoint, a_i
end

function difference(controlPoints)
    n = length(controlPoints[:,1])
    derivativePoints = zeros(n-1, length(controlPoints[1,:]))

    for i=1:n-1
        derivativePoints[i,:] = controlPoints[i+1,:] - controlPoints[i,:]        
    end

    derivativePoints
end

function derivative(differencePoints, param)
    n = length(differencePoints[:,1])
    value = zeros(differencePoints[1,:])
    B = i -> binomial(n, i) * param^i * (1 - param)^(n - i); 
    for i=1:n
        value = value + differencePoints[i,:]*B(i)
    end
    value * (n+1)
end

function readSimplePolygon()
    readdlm("simple-control-polygon.dat")
end

function exercise_one ()
    controlPoints = readSimplePolygon()
    
    bezierCurvePoints = drawCurve(controlPoints,
				  linspace(0,1,200))

    writeArrayForGnuplot(bezierCurvePoints, "bezier-exercise-one.coordinates")
end

function exercise_two ()
    controlPoints = parametricCoordinates([
			                   x -> 1 + x + x^2;
			                   y -> y^3
		                           ], 6)

    bezierCurvePoints = drawCurve(controlPoints,
				  linspace(0,1,200))

    writeArrayForGnuplot(bezierCurvePoints, "bezier-exercise-two.coordinates")
    writeArrayForGnuplot(controlPoints, "control-poly-exercise-two.coordinates")

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
    outerControlPoints = readdlm("control-poly-exercise-four.coordinates")
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
    originalControlPoints = readdlm("control-poly-exercise-six.coordinates")
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
