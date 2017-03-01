classdef SARReflectFIO < FIO
    methods
        function fio = SARReflectFIO(curveInfo, maskInside)
            if nargin < 2,
                maskInside = false;
            end
            
            fio.dims = 2;
            fio.patches = [SARReflectFIOPatch(curveInfo,+1,maskInside,fio) ...
                           SARReflectFIOPatch(curveInfo,-1,maskInside,fio)];
            fio.wrapFlag = 0;
            
            % Output is a nonperiodic grid.
            fio.outLogicalGrid = Grid([0 0], [1 1], [2 2], false);
        end
    end
end