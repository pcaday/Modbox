function fMasked = ApplyWedgeCreatorMask(f, wcr)
    global GlobalDFTRoutine
    global GlobalIDFTRoutine

    
    fhat = GlobalDFTRoutine(f);
    assert(all(fhat.grid.ns == wcr.grid.ns), 'The function and WedgeCreator have incompatible grids.');
    fhat.f = fhat.f .* wcr.Sum();
    fMasked = GlobalIDFTRoutine(fhat, f.grid);
end