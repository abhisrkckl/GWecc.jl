export mismatch

using LinearAlgebra

function mismatch(a, b)
    aa = dot(a, a)
    bb = dot(b, b)
    ab = dot(a, b)
    match = ab / sqrt(aa * bb)
    return 1.0 - match
end

function mismatch(K, a, b)
    aa = dot(a, K, a)
    bb = dot(b, K, b)
    ab = dot(a, K, b)
    match = ab / sqrt(aa * bb)
    return 1.0 - match
end