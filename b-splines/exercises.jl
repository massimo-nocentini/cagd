require("bspline.jl")

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
                         "exercise-zero-clumped-control-poly.coordinates")

    writeArrayForGnuplot(bspline,
                         "exercise-zero-clumped-bspline.coordinates")

    extendedPartition = extendedPartitionAttempt(k, ones(6), Uniformed())
    M = sum(ones(6))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)
   
    writeArrayForGnuplot(controlPoints,
                         "exercise-zero-uniformed-control-poly.coordinates")

    writeArrayForGnuplot(bspline,
                         "exercise-zero-uniformed-bspline.coordinates")


    closedControlPoints = [controlPoints; controlPoints[1,:];]
    closedBSpline = drawBSpline(k, closedControlPoints, 200)
    writeArrayForGnuplot(closedControlPoints,
                         "exercise-zero-closed-control-poly.coordinates")

    writeArrayForGnuplot(closedBSpline,
                         "exercise-zero-closed-bspline.coordinates")

    
end

function exercise_one()
    controlPoints = [
                     [-2.5 -1];
                     [-3 2];
                     [0.3 2.0];
                     [-0.3 -2.0];
                     [4 -2.0];
                     [3 1.0];
                     ]

    writeArrayForGnuplot(controlPoints,
                         "exercise-one-control-poly.coordinates")
    
    fittingPoints = 200
    
    k = 2
    continuityVector = ones(4)
    k_two_extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, k_two_extendedPartition)
    k_two_bspline = drawCurve(k_two_extendedPartition,
                              k,
                              M,
                              controlPoints,
                              matrix,
                              fittingPoints)
    writeArrayForGnuplot(k_two_bspline,
                         "exercise-one-clumped-bspline-k-two.coordinates")

    k = 3
    continuityVector = ones(3)
    k_three_extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, k_three_extendedPartition)
    k_three_bspline = drawCurve(k_three_extendedPartition,
                              k,
                              M,
                              controlPoints,
                              matrix,
                              fittingPoints)
    writeArrayForGnuplot(k_three_bspline,
                         "exercise-one-clumped-bspline-k-three.coordinates")

    k = 4
    continuityVector = ones(2)
    k_four_extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, k_four_extendedPartition)
    k_four_bspline = drawCurve(k_four_extendedPartition,
                              k,
                              M,
                              controlPoints,
                              matrix,
                              fittingPoints)
    writeArrayForGnuplot(k_four_bspline,
                         "exercise-one-clumped-bspline-k-four.coordinates")

    k = 5
    continuityVector = ones(1)
    k_five_extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, k_five_extendedPartition)
    k_five_bspline = drawCurve(k_five_extendedPartition,
                              k,
                              M,
                              controlPoints,
                              matrix,
                              fittingPoints)
    writeArrayForGnuplot(k_five_bspline,
                         "exercise-one-clumped-bspline-k-five.coordinates")

    k = 6
    continuityVector = ones(0)
    k_six_extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, k_six_extendedPartition)
    k_six_bspline = drawCurve(k_six_extendedPartition,
                              k,
                              M,
                              controlPoints,
                              matrix,
                              fittingPoints)
    writeArrayForGnuplot(k_six_bspline,
                         "exercise-one-clumped-bspline-k-six.coordinates")

end

function exercise_two()
    controlPoints = [
                     [-0.7 -0.4];
                     [-1 -1];
                     [-1.3 -0.3];
                     [-1 1.2];
                     [-0.3 1.4];
                     [-0.15 0];
                     [0.15 0];
                     [0.3 1.6];
                     [2 0.6];
                     [2.6 0.2];
                     ]

    writeArrayForGnuplot(controlPoints,
                         "exercise-two-control-poly.coordinates")

    k = 4
    extendedPartition = expandPartition(0:10, [4 [1 for _ in 1:10]'])
    extendedPartition = normalizeRespect(extendedPartition, Maximum(), k)
    M = sum(ones(6))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)   

    writeArrayForGnuplot(bspline,
                         "exercise-two-bspline-first.coordinates")

    k = 4
    extendedPartition = expandPartition(0:4, [4 2 2 2 4])
    extendedPartition = normalizeRespect(extendedPartition, Maximum())
    M = sum(ones(6))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)   

    writeArrayForGnuplot(bspline,
                         "exercise-two-bspline-second.coordinates")

    k = 6
    extendedPartition = expandPartition(0:5, [6 1 1 1 1 6])
    extendedPartition = normalizeRespect(extendedPartition, Maximum())
    M = sum(ones(4))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)   

    writeArrayForGnuplot(bspline,
                         "exercise-two-bspline-three.coordinates")

    k = 6
    extendedPartition = expandPartition(0:4, [6 1 1 2 6])
    extendedPartition = normalizeRespect(extendedPartition, Maximum())
    M = sum(ones(4))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)   

    writeArrayForGnuplot(bspline,
                         "exercise-two-bspline-four.coordinates")

    k = 8
    extendedPartition = expandPartition(0:3, [8 1 1 8])
    extendedPartition = normalizeRespect(extendedPartition, Maximum())
    M = sum(ones(2))
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    fittingPoints = 200
    bspline = drawCurve(extendedPartition, k, M, controlPoints, matrix, fittingPoints)   

    writeArrayForGnuplot(bspline,
                         "exercise-two-bspline-five.coordinates")

    
end

function exercise_three()
    controlPoints = [
                     [-3.0 2];
                     [-1.0 -2];
                     [1 1.0];
                     [1 1.0];
                     [3.0 -2];
                     [5 3.0];
                     ]

    writeArrayForGnuplot(controlPoints,
                         "exercise-three-control-poly.coordinates")
    
    fittingPoints = 200
    
    k = 3
    continuityVector = ones(3)
    extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)
    writeArrayForGnuplot(bspline,
                         "exercise-three-bspline.coordinates")
end

function exercise_four()
    controlPoints = [
                     [-2.0 2];
                     [-3.0 4];
                     [0 6];
                     [0 6];
                     [3 4];
                     [2.0 2];
                     ]

    writeArrayForGnuplot(controlPoints,
                         "exercise-four-control-poly.coordinates")
    
    fittingPoints = 200
    
    k = 4
    ## continuityVector = ones(3)
    ## extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    extendedPartition = [0 0 0 0 1/4 3/4 1 1 1 1]
    M = 2
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)
    writeArrayForGnuplot(bspline,
                         "exercise-four-bspline-first-partition.coordinates")

    extendedPartition = [0 0 0 0 1/2 1/2 1 1 1 1]
    M = 2
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)
    writeArrayForGnuplot(bspline,
                         "exercise-four-bspline-second-partition.coordinates")

end


function exercise_five()
    controlPoints = [
                     [-2.0 1];
                     [-1.0 -1];
                     [0 1];
                     [1 -1];
                     [2.0 1];
                     ]

    writeArrayForGnuplot(controlPoints,
                         "exercise-five-control-poly.coordinates")
    
    fittingPoints = 200
    
    k = 4
    continuityVector = ones(1)
    extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)
    writeArrayForGnuplot(bspline,
                         "exercise-five-bspline-single.coordinates")

    controlPoints = [
                     [-2.0 1];
                     [-1.0 -1];
                     [0 1];
                     [0 1];
                     [1 -1];
                     [2.0 1];
                     ]

    fittingPoints = 200
    
    k = 4
    continuityVector = ones(2)
    extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)
    writeArrayForGnuplot(bspline,
                         "exercise-five-bspline-doubled.coordinates")

    
    controlPoints = [
                     [-2.0 1];
                     [-1.0 -1];
                     [0 1];
                     [0 1];
                     [0 1];
                     [1 -1];
                     [2.0 1];
                     ]

    fittingPoints = 200
    
    k = 4
    continuityVector = ones(3)
    extendedPartition = extendedPartitionAttempt(k, continuityVector, Clumped())
    M = sum(continuityVector)
    matrix = buildFunctionsMatrix(k, M, extendedPartition)
    bspline = drawCurve(extendedPartition,
                        k,
                        M,
                        controlPoints,
                        matrix,
                        fittingPoints)
    writeArrayForGnuplot(bspline,
                         "exercise-five-bspline-tripled.coordinates")

end

function exercise_six()

    k = 4
    
    firstClosedControlPoints = [
                                [1.0 4];
                                [-0.8 4];
                                [-1.9 8.2];
                                [0.3 13.2];
                                [2.8 4];
                                [1.0 4];
                                ]

    firstClosedBSpline = drawBSpline(k, firstClosedControlPoints, 200)
    writeArrayForGnuplot(firstClosedControlPoints,
                         "exercise-six-first-closed-control-poly.coordinates")

    writeArrayForGnuplot(firstClosedBSpline,
                         "exercise-six-first-closed-bspline.coordinates")

    secondClosedControlPoints = [
                                 [0.4 6];
                                 [-1.8 3];
                                 [-2 14.8];
                                 [0.4 11];
                                 [3.0 14];
                                 [3.4 2.2];
                                 [0.4 6];
                                 ]
   
    secondClosedBSpline = drawBSpline(k, secondClosedControlPoints, 200)
    writeArrayForGnuplot(secondClosedControlPoints,
                         "exercise-six-second-closed-control-poly.coordinates")

    writeArrayForGnuplot(secondClosedBSpline,
                         "exercise-six-second-closed-bspline.coordinates")

    
end

function exercise_knotsInsertion()
    controlPoints = [
                     [0.0 0.0];
                     [1.0 2.0];
                     [2.0 2.5];
                     [3.0 2.2];
                     [4.0 0.5];
                     ]

    partition = [-1 0.5 1.8 2.8 3.2 4]

    k = 3
    refined_controlPoints, refined_partition, n, L =
        knotsInsertion(controlPoints, partition, k, 1)
    
    writeArrayForGnuplot(controlPoints,
                         "exercise-knots-insertion-original-control-points.coordinates")

    writeArrayForGnuplot(refined_controlPoints,
                         "exercise-knots-insertion.coordinates")
    
    ## deBoorBSpline = drawBSpline(k, controlPoints, 200)
    ## writeArrayForGnuplot(deBoorBSpline,
    ##                      "exercise-knots-insertion-deBoor.coordinates")

    
end
