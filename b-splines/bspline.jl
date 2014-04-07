
function writeArrayForGnuplot(anArray, aFilename)
    writedlm(aFilename, anArray, '\t')
    anArray
end

type Clumped end
type Uniformed end

function exercise_zero()
    controlPoints = [
                     [-0.2 2];
                     [-0.3 6.2];
                     [-1.2 4.8];
                     [-2.8 8.8];
                     [-0.7 14];
                     [1.4 14.7];
                     [3.6 10.2];
                     [3.2 5.1];
                     [1.5 6.2];
                     [1.4 2];
                     ]


    ## continuityVector = ones(length(controlPoints[:,1]))
    k = 4
    extendedPartition = extendedPartitionAttempt(k, ones(6), Clumped())
    M = sum(ones(6))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)
   
    writeArrayForGnuplot(controlPoints,
                         "exercise-zero-clamped-control-poly.coordinates")

    writeArrayForGnuplot(bspline,
                         "exercise-zero-clamped-bspline.coordinates")

    extendedPartition = extendedPartitionAttempt(k, ones(6), Uniformed())
    M = sum(ones(6))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)
   
    writeArrayForGnuplot(controlPoints,
                         "exercise-zero-uniformed-control-poly.coordinates")

    writeArrayForGnuplot(bspline,
                         "exercise-zero-uniformed-bspline.coordinates")

end

function drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)
    paramRange = linspace(minimum(extendedPartition),
                          maximum(extendedPartition),
                          fittingPoints)
    bspline = zeros(fittingPoints, length(controlPoints[1,:]))
    for i=1:fittingPoints
        bspline[i,:] = bsplineAtParam(paramRange[i], k, M, controlPoints, matrix)
    end
    bspline
end

function bsplineAtParam(t, k, M, controlPoints, matrix)
    value = zeros(length(controlPoints[1,:]))'
    if length(controlPoints[:,1]) != (k+M)
        error("Control points aren't enough")
    end
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

function extendPartition(k, partition, continuityVector, spacing::Clumped)

    extended = []
    extended = [extended; [partition[1,:] for j=1:k]]
        print(typeof(extended))
        print(typeof(      extended[:,1]))
    
    for i=1:length(continuityVector)        
        augmenting = [partition[i + 1,:] for q=1:(continuityVector[i])]
        extended = [extended; augmenting;]
    end

    tail = [partition[end,:] for j=1:k]
    extended = [extended; tail;]
    extended
end

        

