function B = fsqueeze(A)
    assert(size(A,1) == 1, 'fsqueeze:firstDimNot1', ...
        'A has more than 1 row');
    B = shiftdim(A, 1);
end