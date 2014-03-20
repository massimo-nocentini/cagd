function drawCurve(aMatrix, aParameterVector)
	
	n = length(aParameterVector)
	x = zeros(n, 2)
	for i = 1:n
		param = aParameterVector[i]
		aPoint = decasteljau(aMatrix, param)
		x[i,:] = aPoint
	end
	writedlm("simple.dat", x, '\t')
end


function decasteljau(controlPointsMatrix, t)
	n = length(controlPointsMatrix[:,1])
	Q = copy(controlPointsMatrix)

	for k = 1:n
		for i = 1:n-k
			Q_i = Q[i,:]
			Q_isucc = Q[i+1,:]	
			Q[i,:] = (1 - t)*Q_i + t*Q_isucc
		end
	end
	return Q[1,:]
end

function readSimplePolygon()
	 readdlm("simple-control-polygon.dat")
end
