classdef ProcessedLocalFIOPatch < matlab.mixin.Copyable
    properties
        localPatch          % LocalFIOPatch that we are derived from.
        warps               % Cell array of warp functions
        expFs               % Cell array of \tilde{f} functions
        expCs               % Array of PSWF c values
        cExtras             % Square root of the ratio of the rounded c
                            %  to the computed c.
        amplitude           % Amplitude function (stored as a matrix,
                            %  with x's as a single linear row index, and
                            %  xi as the column index.)
        inGrid              % Grid object for domain
        wedgeCreator        % WedgeCreator object
    end
    
    
    methods
        function pPatches = ProcessedLocalFIOPatch(localPatches, inGrid, wcc, zeroPadding, valueClass)
            % 
            % pPatches = ProcessedLocalFIOPatch(locals, inGrid, wcc, zp, class)
            %
            % Perform precomputations for the modified box algorithm for
            %  the given local patch(es). The precomputations are stored
            %  in an array of ProcessedLocalFIOPatch objects.
            %
            %   locals: array of LocalFIOPatch objects
            %   inGrid: x grid
            %      wcc: function handle which takes a frequency grid and
            %             zero-padding factor, and returns a WedgeCreator
            %             object for the grid.
            %       zp: zero-padding factor.
            %    class: 'single' or 'double' 
            %
            % pPatches: array of ProcessedLocalFIOPatch objects, with the
            %            same size as the local patches arrays.
            
            global GlobalAmplitudeDoXCutoff
            global GlobalAmplitudeSuppThreshold
            global GlobalPSWFCache
            global GlobalUsePseudopolarFFT
            global GlobalDFTGridRoutine
            global GlobalCWarningThreshold
            global GlobalNewProcessVerbosity
            
            % No inputs: do nothing and let MATLAB make empty object.
            if nargin == 0, return; end
            
            % Otherwise: allocate pPatches to be an array of the same
            %  size as localPatches by creating the last item in the array.
            sz = num2cell(size(localPatches));
            pPatches(sz{:}) = ProcessedLocalFIOPatch;
            
            % Process the patch.
            for k = 1:numel(localPatches)
                pPatches(k).localPatch = localPatches(k);
                diffeo = localPatches(k).diffeo;
                n = localPatches(k).globalPatch.dims;
            
                % Find an x~ grid (inGridT) covering the x grid.
                % For pseudopolar DFT, make sure the grid has the same
                %  number of grid points in each coordinate, and that this
                %  is an even number.
                inGridT = diffeo.TransformGrid(inGrid);
                if GlobalUsePseudopolarFFT,
                    inGridT.ns(:) = ceil(max(inGridT.ns) / 2) * 2;
                end
                pPatches(k).inGrid = inGridT;

                % Create the WedgeCreator.
                inDFTGridT = GlobalDFTGridRoutine(inGridT, zeroPadding);               
                wcrT = wcc(inDFTGridT);
                pPatches(k).wedgeCreator = wcrT;
                
                % Get some information from it.
                nBoxes = wcrT.Count;
                nus = wcrT.Centers;
                
                % Get grid points we'll need.
                edgeGrid = inGridT.PixelEdgeGrid;
                xeT = edgeGrid.AllPoints;
                xT = inGridT(:);
                xiT = inDFTGridT(:);
                f = zeros(inGridT.NPoints, n*(n-1)/2);
                
                % Allocate variables
                pPatches(k).expFs = zeros([size(f) nBoxes]);
                pPatches(k).expCs = zeros(nBoxes, 1);
                pPatches(k).cExtras = zeros(nBoxes, 1);
                pPatches(k).warps = cell(nBoxes, 1);
                pPatches(k).amplitude = cell(nBoxes, 1);
                
                
                for l = 1:nBoxes
                    % Compute the warps for each wedge: i.e., the images of
                    %  each (x~,nu~), where x~ is a corner point of a
                    %  pixel in the grid (i.e., a grid point shifted by 1/2)
                    pPatches(k).warps{l} = localPatches(k).CanonicalTransformationY(xeT,nus(:,l));
                    
                    % Compute amplitudes at each (x~,nu), where x~ is a
                    %  point on the grid.
                    [amp, dxtdxit] = localPatches(k).Amplitude(xT,nus(:,l));

                    if isempty(dxtdxit),
                        % If the amplitude is identically zero, dx~/dxi~
                        %  will be empty, and we can skip all the other
                        %  calculations.
                        roundedC = 0;
                        cExtra = 1;
                        pPatches(k).expFs(:,:,l) = 0;
                    else
                        % Cut off the amplitude outside the image of the
                        %  (x,xi) domain, if configured.
                        if GlobalAmplitudeDoXCutoff,
                            x = diffeo.InverseTransform(xT, nus(:,l));
                            amp(~inGrid.IsInside(x)) = 0.;
                        end
                        
                        % Convert any NaN's or infinities to zero. Warn about
                        %  them.
                        nonfinite = ~isfinite(amp);
                        if any(nonfinite),
                            warning('ProcessedLocalFIOPatch:ctor:nonfiniteamp', 'Nonfinite amplitudes set to zero.');
                            amp(nonfinite) = 0;
                        end
                        
                        % Compute unscaled f = xi''^T (dx~/dxi~) xi'' / xi'.
                        nuTBasis = cast(wcrT.NuBasis(l), valueClass);
                        idx = 1;
                        for j = 2:n
                            dxtdxit_nuj = ParallelMatrixMultiply(dxtdxit, nuTBasis(:,j));
                            for i = 2:j
                                f(:,idx) = ParallelMatrixMultiply(nuTBasis(:,i).', dxtdxit_nuj);
                                idx = idx + 1;
                            end
                        end
                        
                        % Get the amplitude's support.
                        %  (not the precise mathematical support, just
                        %   where the amplitude is "numerically significant").
                        ampSupp = abs(amp) > GlobalAmplitudeSuppThreshold;
                        
                        % Find the maximum of f on the amplitude's support.
                        fmax = sqrt(max(sum(f.*f, 2) .* ampSupp));
                        
                        % Compute the sup of g's 2-norm (unscaled) over the
                        %  wedge's support.
                        [~,gmax] = ComputeUnscaledG(xiT, wcrT, l, valueClass);
                        
                        % Compute c and round to nearest cached c value.
                        % Record the extra.
                        c = fmax * gmax;
                        
                        if c > 0,
                            roundedC = GlobalPSWFCache.RoundC(c);
                            cExtra = sqrt(roundedC / c);
                        else
                            roundedC = 0;
                            cExtra = 1;
                        end
                        
                        % Divide f's by c_extra * sup |f|_2 and record.
                        pPatches(k).expFs(:,:,l) = f * (1 / (fmax * cExtra));

                    end

                    % Store amplitude, c, and extra c factors
                    pPatches(k).expCs(l) = roundedC;
                    pPatches(k).cExtras(l) = cExtra;
                    pPatches(k).amplitude{l} = amp;
                end
                
                % Check for high c values, and warn the user if they occur.
                cmax = max(pPatches(k).expCs);
                if GlobalNewProcessVerbosity > 0,
                    fprintf('Coordinate choice #%d: max c = %g\n', k, cmax);
                end
                if cmax > GlobalCWarningThreshold,
                    warning('ProcessedLocalFIOPatch:highC', 'Large c value(s) have been found for coordinate choice #%d. Caustics may be present and not handled correctly.', k);
                end
            end            
        end
    end
end