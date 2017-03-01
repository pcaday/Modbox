classdef SinglePatchFIO < FIO
    methods
        function fio = SinglePatchFIO(patch, wrapFlag, outLogicalGrid)
            if nargin < 2, wrapFlag = 0; end
            if nargin < 3,
                % Default logical grid is nonperiodic.
                e = ones(patch.dims, 1);
                outLogicalGrid = Grid(-inf*e, inf*e, 2*e);
            end
            
            fio.dims = patch.dims;
            fio.patches = patch;
            fio.wrapFlag = wrapFlag;
            fio.outLogicalGrid = outLogicalGrid;
            
            for p = patch(:).'
                p.parent = fio;
            end
        end
    end
end