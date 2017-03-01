% frep: Adds a singleton first dimension, then concatenates n copies
%        of M along the first dimension.

function repM = frep(M, n)
    sz = size(M);
    M = reshape(M, [1 sz]);
    repM = M(ones(n,1), :);
    repM = reshape(repM, [n sz]);
end