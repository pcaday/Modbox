function rA = repmattosize(A, szB)
    szA = size(A);
    dB = length(szB);

    assert(dB >= 2, 'The size argument must have at least two dimensions.');

    if isequal(szA, szB)
        rA = A;
    elseif isempty(A),
        assert(prod(szB) == 0, 'If A is empty then B must be also.');
        rA = reshape(A, szB);
    else
        dA = ndims(A);
        assert(dA <= dB, 'A has more dimensions than B.');
        szA(dA+1:dB) = 1;

        repA = szB ./ szA;
        assert(all(fix(repA) == repA), 'A does not fit evenly into B');
        rA = repmat(A, repA);
    end
end