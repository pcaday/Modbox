classdef ProcessedFIO < dynamicprops
    properties
        fio                 % Parent FIO object
        patches             % Array of ProcessedLocalFIOPatch objects
        inGrid              % Grid object for domain
    end
    
    
    methods
        % Create a ProcessedFIO object, which contains the precomputations
        %  necessary for applying the modified box algorithm to the FIO.
        %
        %  pFIO = ProcessedFIO(fio, ig, og, diffeos, wcc, zp, class)
        %
        % Inputs:
        %     fio: FIO object.
        %  ig, og: Grid objects representing the grids for the domain
        %           and codomain of the FIO. The codomain grid's size is
        %           not important, just whether it's periodic or not.
        %     wcc: ("WedgeCreator creator") function handle which
        %           takes a Grid object for a frequency grid and returns a
        %           WedgeCreator object on that grid.
        %      zp: zero-padding factor (default 1)
        %   class: 'single' or 'double', specifying which class
        %           to use for storing the data.
        function pFIO = ProcessedFIO(fio, inGrid, ~, diffeos, wcc, zp, valueClass)
            if nargin < 7, valueClass = 'single'; end
            if nargin < 6, zp = 1; end
            
            pFIO.fio = fio;
            pFIO.inGrid = inGrid;
            
            nPatches = length(fio.patches);
            patchList = cell(nPatches,1);
            
            for i = 1:nPatches
                lPatches = LocalFIOPatch(fio.patches(i), diffeos); 
                patchList{i} = ProcessedLocalFIOPatch(lPatches, inGrid, wcc, zp, valueClass);
            end
            pFIO.patches = [patchList{:}];
        end
    end
end