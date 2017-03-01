classdef Radon2DFIO < FIO
    methods
        function fio = Radon2DFIO()
            fio.dims = 2;
            fio.patches = [Radon2DPatch(+1,fio)...
                           Radon2DPatch(-1,fio)];
            fio.wrapFlag = 2;       % Really should only wrap in the theta
                                    % variable, not s, but this isn't yet
                                    % implemented.
            fio.outLogicalGrid = Grid([-inf,0], [inf,2*pi], [2 2], [false true]);
        end
    end
end