function reshapeVals = preshape(vals, inputs)
    sz = size(inputs);
    reshapeVals = reshape(vals, [sz(2:end) 1]);
end