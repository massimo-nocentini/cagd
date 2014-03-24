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
	controlPoints = parametricCoordinates(
		[
			x -> 1 + x + x^2;
			y -> y^3
		], 6)

	bezierCurvePoints = drawCurve(controlPoints,
					linspace(0,1,200))

	writeArrayForGnuplot(bezierCurvePoints, "bezier-exercise-two.coordinates")
	writeArrayForGnuplot(controlPoints, "control-poly-exercise-two.coordinates")

	controlPoints = parametricCoordinates(
		[
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
