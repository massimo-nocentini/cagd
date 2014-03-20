function drawCurve(aMatrix, aParameterVector)
	
	n = length(aParameterVector)
	x = zeros(n, 2)
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

function parametricCoordinates(xFn, yFn, pointsNumber)

	parameterVector = linspace(0, 1, pointsNumber)	

	coordinates = zeros(pointsNumber, 2)

	for i = 1:pointsNumber
		t = parameterVector[i]
		coordinates[i,:] = [xFn(t) yFn(t)]		
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
	sandbox = [ [0 0]; copy(originalControlPoints); [0 0] ]
	increasing = [copy(originalControlPoints); [0 0] ]
	n = length(originalControlPoints[:,1])
	
	for i=1:n+1
#		if i == 1
#			increasing[1,:] = sandbox[2,:] 
#		elseif i == n+1
#			increasing[n+1,:] = sandbox[n+1,:] 
#		else
			ratio = (i-1)/(n)
			increasing[i,:] = ratio*sandbox[i,:] + (1 - ratio)*sandbox[i+1,:] 
#		end
	end

	increasing
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
	controlPoints = parametricCoordinates(x -> 1 + x + x^2, y -> y^3, 6)

	bezierCurvePoints = drawCurve(controlPoints,
					linspace(0,1,200))

	writeArrayForGnuplot(bezierCurvePoints, "bezier-exercise-two.coordinates")
	writeArrayForGnuplot(controlPoints, "control-poly-exercise-two.coordinates")
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


end
