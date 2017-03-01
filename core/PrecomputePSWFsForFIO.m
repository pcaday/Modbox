function PrecomputePSWFsForFIO(pswfCache, fio)
    n = fio.dims;
    % Construct a list of all c values used.
    allCs = [];
    for patch = fio.patches(:)
        allCs = [allCs patch.expCs]; %#ok<AGROW>
    end
    allCs = unique(allCs);
    % Now call GetCachedPSWF to fetch them all into the cache.
    for c = allCs(:)
        [~] = pswfCache.GetCachedPSWF(n, c);
    end
end