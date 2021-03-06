
function writeArrayForGnuplot(anArray, aFilename)
    writedlm(aFilename, anArray, '\t')
    anArray
end

type Clumped end
type Uniformed end

function drawBSpline(k, controlPoints, fittingPoints)

    if controlPoints[1,:] == controlPoints[end,:]

        m = length(controlPoints[:,1])
        newPartition = [i for i in (-k/(m-1):1/(m-1):(k+m-1)/(m-1))]

        M = sum(ones(m))
        matrix = buildFunctionsMatrix(k, M, newPartition)

        for i=1:k
            controlPoints = [controlPoints; controlPoints[i+1,:]]
        end
        
        return drawCurve(newPartition, k, M, controlPoints, matrix, fittingPoints)
    end

end

function drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)

    
    ## paramRange = linspace(minimum(extendedPartition),
    ##                       maximum(extendedPartition),
    ##                       fittingPoints)

    ## fittingPoints += 1
    paramRange = linspace(0,
                          1-eps(Float64),
                          fittingPoints)

    bspline = zeros(fittingPoints, length(controlPoints[1,:]))
    for i=1:fittingPoints
        bspline[i,:] = bsplineAtParam(paramRange[i],
                                      k,
                                      M,
                                      controlPoints,
                                      matrix)
    end
    bspline
    ## bspline[1:end-1, :]
end

function bsplineAtParam(t, k, M, controlPoints, matrix)
    value = zeros(length(controlPoints[1,:]))'
    ## if length(controlPoints[:,1]) != (k+M)
    ##     error("Control points aren't enough")
    ## end
    for i=1:(k+M)
        v = controlPoints[i,:]
        N_ik = matrix[i,k]
        value += v*N_ik(t)
    end
    value
end

function extendedPartitionAttempt(k, partition, spacing::Clumped)
    previous_zeros = zeros(k)
    following_ones = ones(k)
    M = length(partition)
    middle = zeros(M)
    range = linspace(1/(M+1), M/(M+1), M)
    for i=1:length(range) 
        middle[i] = range[i]
    end
    [previous_zeros; middle; following_ones;]
end

function extendedPartitionAttempt(k, partition, spacing::Uniformed)
    M = length(partition)
    range = (1-k)/(M+1):1/(M+1):(1+ (k-1)/(M+1))
    nel = length(range)
    extended = zeros(nel)
    for i=1:nel
        extended[i] = range[i]
    end
    extended
end


function buildFunctionsMatrix(k, M, extendedPartition)
    if (2*k + M) != length(extendedPartition)
        error("partition length mismatch dimension")
    end
    firstColumn = []
    for i=1:(k+M+1)
        fn = (t ->
              if (extendedPartition[i] <= t && t < extendedPartition[i+1])
                  1
              else
                  0
              end)
        firstColumn = [firstColumn; fn;]
    end
    ## firstColumn = [firstColumn; (t -> 0);]

    matrix = [firstColumn]
    for h=2:k
        column = []
        for i=1:(k+M)

            augend_coeff = 0
            augend_denominator = extendedPartition[i+h-1] - extendedPartition[i]
            if augend_denominator != 0
                augend_coeff =  1 / augend_denominator
            end

            addend_coeff = 0
            addend_denominator = extendedPartition[i+h] - extendedPartition[i+1]
            if addend_denominator != 0
                addend_coeff =  1 / addend_denominator
            end
            
            fn = t -> ((t - extendedPartition[i]) * augend_coeff * matrix[i,h-1](t) +
                       (extendedPartition[i+h] - t) * addend_coeff * matrix[i+1,h-1](t))
            
            column = [column; fn;]
        end
        lastFn = (t ->
              if (extendedPartition[k+M+1] <= t && t < extendedPartition[k+M+2])
                  0
              else
                  0
              end)

        column = [column; lastFn;]
        matrix = [matrix column]
    end
    
    matrix
    
end

type Maximum end

function normalizeRespect(vector, _::Maximum, k = 1)
    m = maximum(vector[1:end - (k-1)])
    return map(t -> t/m, vector)
end

function expandPartition(partition, multiplicityVector)

    extended = []
    
    for i=1:length(multiplicityVector)
        augmenting = [partition[i] for _=1:multiplicityVector[i] ]
        extended = [extended; augmenting;]
    end

    extended
end

function knotsInsertion_head(controlPoints, partition, n, L)

    fittingPoints = 100
    working_partition = partition[n:n+L]
    paramRange = linspace(working_partition[1],
                          working_partition[end],
                          fittingPoints)
    
    ## u = working_partition[1] + rand()*(working_partition[end]-working_partition[1])
    
    results = zeros(fittingPoints, length(controlPoints[1,:]))
    for j=1:fittingPoints-1

        u = paramRange[j]

        workingControlPoints = controlPoints
        workingPartition = partition
        workingL = L

        for i=1:n
            workingControlPoints, workingPartition, workingL =
                knotsInsertion(workingControlPoints, workingPartition, n, workingL, u)
        end

        ## for i=1:length(workingControlPoints[:,1])-1
        ##     if(workingControlPoints[i,:] == workingControlPoints[i+1,:])
                k = findfirst(x -> x == u, workingPartition) 
                results[j,:] = workingControlPoints[k,:]
        ##         break
        ##     end
        ## end
    end

    results[1:end-1,:]
end
    
function knotsInsertion(controlPoints, partition, n, L, u)
    ## working_partition = partition[n:n+L]
    ## ## u = working_partition[1] + rand()*(working_partition[end]-working_partition[1])
    ## u = 2.5

    ## if (count(x -> x == u, partition) == n)
    ##     return controlPoints, partition, n, L
    ## end

    I = -1
    for i = length(partition):-1:(n+1)
        if(partition[i] > u && partition[i-1] <= u)
            I = i-1
            break
        end
    end

    last_index = L+n+1
    refined_partition = zeros(last_index+1)
    
    for i=1:I-n+2
        refined_partition[i] = partition[i]
    end

    for i=I-n+3:I+1
        refined_partition[i] = (1/n)*(sum(partition[i:i+n-2]) + u)
    end

    for i=I+2:last_index+1
        refined_partition[i] = partition[i-1]
    end

    refined_controlPoints = zeros(last_index+1, length(controlPoints[1,:]))

    for i=1:I-n+2
        refined_controlPoints[i,:] = controlPoints[i,:]
    end

    for i=I-n+3:I+1
        diff = partition[i+n-1]-partition[i-1]
        coeff = if (diff > 0) 1/diff else 0 end
        first_coeff = partition[i+n-1] - u
        second_coeff = u - partition[i-1]
        refined_controlPoints[i,:] = coeff*(first_coeff*controlPoints[i-1,:] +
                                            second_coeff*controlPoints[i,:])
    end

    for i=I+2:last_index+1
        refined_controlPoints[i,:] = controlPoints[i-1,:]
    end

    partition = [partition; u]
    refined_controlPoints, sort(partition), L+1
end

function deBoor (extendedBreakpointArray,
                 controlPointsMatrix,
                 tabulationArray)

    dimension, numberOfControlPoints = size (controlPointsMatrix)
    k = length (extendedBreakpointArray) - numberOfControlPoints

    C = zeros (dimension, length (tabulationArray))
    ind = 0

    for r = k:numberOfControlPoints

        comparison = r < n ? .< : .<= ;
        isd = (tabulationArray .>= extendedBreakpointArray [r]) & (comparison(tabulationArray, extendedBreakpointArray [r+1]))

        nloc = sum (isd)

        if nloc < 1
            continue
        end

        tloc = tabulationArray[isd]
        Qloc = zeros (k, nloc, dimension)

        ## the following is a one line in substitution of the next line of code
        ## broadcast(+, [[1 2 5]; [3 4 6];], zeros(2, 3, 2))
        [Qloc[:,:,i] = Q[i, r-k+1:r]' * ones(1,nloc) for i=1:dimension]
        
        for j=1:k-1
            alfa = zeros(k-1,nloc)
            for i=1:k-j
                alfa[i,:] = (tloc - t[i+r-k+j])/ (t[i+r] - t[i+r-k+j])
                [Qloc[i,:,s] = (1-alfa[i,:]) .* Qloc[i,:,s] + alfa[i,:] .* Qloc[i+1,:,s] for s=1:dimension]
            end
        end            

        [C[i, ind+1:ind+nloc] = Qloc[1, 1:nloc, i] for i=1:dimension]
        ind += nloc
        
    end

end
