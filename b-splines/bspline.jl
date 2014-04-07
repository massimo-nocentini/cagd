
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


    continuityVector = ones(length(controlPoints[:,1]))
    k = 4
    hukairs = -3
    extendedPartition = extendPartition(k,
                                        [
                                         [-3 hukairs];
                                         controlPoints;
                                         [4 hukairs];
                                         ],
                                        continuityVector,
                                        Clumped())

    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, extendedPartition[:,1])

    writeArrayForGnuplot(controlPoints,
                         "exercise-zero-clamped-control-poly.coordinates")

    
    
    extendedPartition, matrix
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

function buildFunctionsMatrix(k, M, extendedPartition)
    firstColumn = []
    for i=1:(k+M)
        fn = (t ->
              if (extendedPartition[i] <= t && t < extendedPartition[i+1])
                  1
              else
                  0
              end)
        firstColumn = [firstColumn; fn;]
    end
    firstColumn = [firstColumn; (t -> 0);]

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
            
            fn = t -> ((t - extendedPartition[i]) * augend_coeff * matrix[i,end](t) +
                       (extendedPartition[i+h] - t) * addend_coeff * matrix[i+1,end](t))
            
            column = [column; fn;]
        end
        column = [column; (t -> 0);]
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

        

