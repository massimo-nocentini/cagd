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
    value = zeros( length(differencePoints[1,:]))'
    B = i -> binomial(n, i) * (param^i) * ((1 - param)^(n - i))
    for i=1:n
        value = value + differencePoints[i,:]*B(i)
    end
    value * (n+1)
end


function polar(controlPoints, param, respect_to)

    n = length(controlPoints[:,1]) - 1
    
    polarpoint = zeros(length(controlPoints[1,:]))'

    B = i -> binomial(n, i) * (param^i) * ((1 - param)^(n - i))

    for i=1:n
        polarpoint += ((1-respect_to)*controlPoints[i,:] + respect_to*controlPoints[i+1,:])*B(i)
    end

    polarpoint
end
