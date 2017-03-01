classdef IdentityFIO < FIO
    methods
        function fio = IdentityFIO(n)
            fio.dims = n;
            fio.patches = IdentityFIOPatch(fio);
            fio.wrapFlag = 0;
        end
    end
end