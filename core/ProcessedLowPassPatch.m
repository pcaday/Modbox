classdef ProcessedLowPassPatch < matlab.mixin.Copyable
    properties
        nAngles         % Number of angles in frequency domain to consider
        nRadii          % Number of radii to consider
        nFreqs          % Total number of frequencies = nAngles * nRadii;
        
        wedgeCreator    % Associated WedgeCreator object.
        
        freqs           % n x nRadii x nAngles array of frequencies
        warps           % length-nAngles cell array of n x i1 x ... x in
                        %  arrays representing the map given by the
                        %  canonical transformation restricted to the xi
                        %  value at the center of the array.
                        % As with ProcessedFIOPatch, these are the images
                        %  of the pixel corner grid.
        amps            % length-nAngles cell array of values of the
                        %  amplitude.
        pPatch          % ProcessedFIOPatch that we are piggybacking on.
    end
    
    methods
        function lp = ProcessedLowPassPatch(pPatch, angleOS, radOS)
            % Do precomputation for the low-pass portion of the algorithm,
            %  constructing a ProcessedLowPassPatch object.
            %
            % p = ProcessedLowPassPatch(pPatch)
            % p = ProcessedLowPassPatch(pPatch, angleOS, radOS)
            %
            % Inputs
            %   pPatch: ProcessedFIOPatch array, or ProcessedFIO
            %  angleOS: oversampling factor in the angular direction
            %             (default: 1, no oversampling)
            %    radOS: oversampling factor in the radial direction
            %             (default: 1, no oversampling)
            if nargin == 0,
                % No-argument constructor: do nothing.
                return;
            elseif isa(pPatch, 'ProcessedFIO'),
                lp = ProcessedLowPassPatch(pPatch.patches, angleOS, radOS);
                return;
            elseif length(pPatch) > 1,
                % Initialize ProcessedLowPassPatch array
                sz = num2cell(size(pPatch));
                lp(sz{:}) = ProcessedLowPassPatch();
                
                for i = 1:numel(lp)
                    lp(i) = ProcessedLowPassPatch(pPatch(i), angleOS, radOS);
                end
                return;
            end
            
            % Substitute default arguments if not given
            if nargin < 2 || isempty(angleOS),  angleOS = 1;    end
            if nargin < 3 || isempty(radOS),    radOS = 1;      end

            % Get the WedgeCreator and frequency grid.
            wedgeCreator_ = pPatch.wedgeCreator;
            freqGrid = wedgeCreator_.grid;
            
            % Get number of dimensions
            dims = pPatch.localPatch.globalPatch.dims;
            
            % Get the highest frequency not covered by the WedgeCreator
            lpCeil = wedgeCreator_.LowPassCeiling();
            
            % How many discrete frequencies is that in the cartesian DFT?
            if isa(freqGrid, 'PseudopolarGrid'),
                ds = repmat(freqGrid.dr, dims, 1);
            else
                ds = freqGrid.ds;
            end
            lpCeil = ceil(lpCeil ./ ds);
            R = max(lpCeil);

            % Get the number of radii and angles to use (angles should
            % always be a multiple of 4)
            nRadii_ = ceil(radOS * R);
            nAngles_ = ceil(angleOS * 2 * R) * 4;
            
            % Store properties
            lp.wedgeCreator = wedgeCreator_;
            lp.nRadii = nRadii_;
            lp.nAngles = nAngles_;
            lp.nFreqs = nRadii_ * nAngles_;
            lp.pPatch = pPatch;
            
            % Create a pseudopolar grid.
            M = nAngles_ / 4;
            N = nRadii_;
            maxFreq = lpCeil .* ds;
            lp.freqs = PseudopolarGrid.GeneratePseudopolarGrid(maxFreq,M,N);
            
            % Remove the zero frequency.
            lp.freqs(:,1,:) = [];
            
            % Calculate the warps and amplitudes
            xis = squeeze(lp.freqs(:,end,:));
            xis = bsxfun(@rdivide, xis, sum(xis.^2, 1));
            localPatch = pPatch.localPatch;
            inGrid = pPatch.inGrid;
            edgeGrid = inGrid.PixelEdgeGrid;
            xeT = edgeGrid.AllPoints;
            xT = inGrid.AllPoints;
            
            lp.warps = cell(nAngles_,1);
            lp.amps = cell(nAngles_,1);
            
            haveNonfinite = false;
            for i = 1:nAngles_
                lp.warps{i} = localPatch.CanonicalTransformationY(xeT, xis(:,i));
                a = localPatch.Amplitude(xT, xis(:,i));
                
                nonfinite = ~isfinite(a);
                if any(nonfinite(:)),
                    haveNonfinite = true;
                    a(nonfinite) = 0;
                end
                
                lp.amps{i} = a;
%                if GlobalAmplitudeDoXCutoff,
%                    x = localPatch.diffeo.InverseTransform(xT, xis(:,i));
%                    lp.amps{i}(~inGrid
%                end
            end
            
            if haveNonfinite,
                warning('ProcessedLowPassPatch:ctor:nonfiniteamp', 'Nonfinite amplitudes set to zero.');
            end
        end
    end
end
    
    