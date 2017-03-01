function subs = DFTCenterSubs(dftGrid)
    subs = ceil((dftGrid.ns+1) / 2);
    subs = num2cell(subs);
end