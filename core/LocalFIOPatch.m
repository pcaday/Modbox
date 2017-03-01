classdef LocalFIOPatch < matlab.mixin.Copyable
    properties
        globalPatch             % The FIOPatch we are derived from
        diffeo                  % The Diffeo object representing our local
                                %  coordinates.
        allDiffeos              % Array of all Diffeos used on the global
                                %  patch.
        diffeoIndex             % Index of this diffeo in the allDiffeos
                                %  array.
    end
    
    methods
        function lPatch = LocalFIOPatch(gPatch, allDiffeos, diffeoIndex)
        % Create a local representation(s) of an FIO patch.
        %
        %     lPatches = LocalFIOPatch(gPatch, allDiffeos)
        %       lPatch = LocalFIOPatch(gPatch, allDiffeos, i)
        %
        % Inputs
        %     gPatch: FIOPatch object
        % allDiffeos: cell array of Diffeo objects for the coordinate
        %              changes to use on the global patch.
        %          i: which coordinate change to use for this local 
        %              representation of the FIO.
        %             This is an index into the allDiffeos array.
        %             If omitted, an array of LocalFIOPatches will be
        %              created, one for each coordinate change.
            if nargin >= 3,
                lPatch.globalPatch = gPatch;
                lPatch.diffeo = allDiffeos{diffeoIndex};
                lPatch.allDiffeos = allDiffeos;
                lPatch.diffeoIndex = diffeoIndex;
            elseif nargin > 0,
                % Preallocate...
                nDiffeos = numel(allDiffeos);
                lPatch(nDiffeos) = LocalFIOPatch();
                % and create each LocalFIOPatch.
                for i = 1:nDiffeos
                    lPatch(i) = LocalFIOPatch(gPatch, allDiffeos, i);
                end
            end
        end
        
        
        function [a,dXtdXit] = Amplitude(patch, xt, xit)
        % Calculate the amplitude of the local representation at given
        %  (x~,xi~) values.
        %
        % a = localPatch.Amplitude(xt, xit)
        % 
        % xt and xit are n x anything arrays of x~ and xi~ coordinates
        %  at which to compute the amplitude. They must either have the
        %  same size, or one of them can be a single vector, in which
        %  case the same x or xi value is repeated for every value of the
        %  other variable.
        %
        % [a,dXtdXit] = localPatch.Amplitude(...)
        %
        % Also returns in dXtdXit the derivative dx~/dxi~ of the canonical
        %  transformation, used by the second-order PSWF approximation
        %  algorithm.
        
            global GlobalProcessNoDyDxFactor
            global GlobalProcessMaslovSign
            global GlobalProcessNewMaslov
            global GlobalDebugMaslov
            global GlobalDiffeoThreshLow
            global GlobalDiffeoThreshHigh
            global GlobalSignatureTol
        
            % Set dXtdXit in case we exit early
            dXtdXit = [];
            
            % Get global patch
            gPatch = patch.globalPatch;
            
            % Change coordinates to (x,xi)
            [x,xi] = patch.diffeo.InverseTransform(xt,xit);
        
            % Start with principal symbol
            a = gPatch.PrincipalSymbol(x,xi);

            % Shortcut: if the principal symbol's identically zero here,
            %  we're already done.
            if ~any(a(:)), return; end
            
            % Calculate the first derivatives dy/dx~, dy/dxi~.
            % We also need to calculate the smallest singular values
            %  of dy/dx~ for all coordinate choices on the global patch.
            % We'll do both tasks at the "same time."
            
            % First dy/dx, dy/dxi.
            [dYdX, dYdXi] = gPatch.DCanonicalTransformationY(x,xi);
            assert(~isempty(dYdX), 'LocalFIOPatch:Amplitude:needDerivs', ...
                'The first derivatives dy/dx, dy/dxi of the canonical transformation must be implemented.');
            
            % [dYdX, dYdXi] = gPatch.DCTYQuickFiniteDifferences(x,xi);

            nDiffeos = numel(patch.allDiffeos);
            allMinSVs = cell(nDiffeos, 1);
            maxMinSVs = zeros(size(a), class(a));
            for k = 1:nDiffeos
                diffeoK = patch.allDiffeos{k};
                if isa(diffeoK, 'IdentityDiffeo'),
                    % Identity diffeomorphism:
                    %   dy/dx~ = dy/dx, dy/dxi~ = dy/dxi.
                    
                    % Remember dy/dx temporarily.
                    thisDyDxt = dYdX;
                    
                    % Only remember dy/dx, dy/dxi~ in our coordinates.
                    if k == patch.diffeoIndex
                        dYdXit = dYdXi;
                        dYdXt = dYdX;
                    end
                else
                    % Get dy/dx~ and dy/dxi~ with the chain rule.
                    % First, get dx/dx~, dxi/dx~, dxi/dxi~.
                    [dXdXt,dXidXt,dXidXit] = diffeoK.DInverseCotangent(x,xi);
                    
                    % Remember dy/dx~ temporarily
                    thisDyDxt = ParallelMatrixMultiply(dYdX,dXdXt)...
                              + ParallelMatrixMultiply(dYdXi,dXidXt);
                            
                    % Only remember dy/dx, dy/dxi~ in our coordinates.
                    if k == patch.diffeoIndex
                        dYdXit = ParallelMatrixMultiply(dYdXi,dXidXit);
                        dYdXt = thisDyDxt;
                    end
                end
                
                % Compute the minimum singular values of dy/dx~
                %  for these coordinates, and store.
                minSVs = min(abs(ParallelSingularValues(thisDyDxt)), [], 1);
                minSVs = reshape(minSVs, size(a));
                
                allMinSVs{k} = minSVs;
                maxMinSVs = max(maxMinSVs, minSVs);
            end

            % Now, compute the cutoff to apply to the principal symbol
            %  on our local patch.
            
            % First, make a smooth cut off from the minimum SVs which is:
            % * 1, if this minimum SV is the maximum of all minimum SVs for 
            %      all patches.
            % * 0, if this minimum SV is less than GlobalBadDiffeoThreshold
            %      times the maximum of all patches.
            % * smoothly varying in between.
            %
            % (Everything is done pointwise)
            %
            % Then total up the resulting cut-off minimum SVs (totalCOMinSVs)
            %  and save our cut-off minimum SVs.
            totalCOMinSVs = zeros(size(maxMinSVs), class(a));
            
            for k = 1:nDiffeos
                coMinSVs = allMinSVs{k} .* ...
                                HalfWindow(allMinSVs{k} ./ maxMinSVs, ...
                                      GlobalDiffeoThreshLow, ...
                                      GlobalDiffeoThreshHigh, 'smooth');
                if k == patch.diffeoIndex,
                    ourCOMinSVs = coMinSVs;
                end
                totalCOMinSVs = totalCOMinSVs + coMinSVs;
            end
            
            % Shortcut: if we are completely outside the cutoff, a is zero,
            %  and we can return.
            if ~any(ourCOMinSVs(:)),
                a(:) = 0;      % Some quick tests indicated this is
                return;        %  just marginally faster than 
            end                %  a = zeros(size(a));
            
            % All right, divide our cutoffs by the total to get
            %  the scaled, final cutoff for this local representation.
            a = a .* (ourCOMinSVs ./ totalCOMinSVs);
            
            % Multiply the principal symbol by the half-density factor.
            if ~GlobalProcessNoDyDxFactor
                % Compute |det (dy/dx~)|
                adetDxty = abs(ParallelDet(dYdXt));                
                a = a ./ sqrt(adetDxty);
            end

            % Multiply by appropriate Maslov factor.
            
            % First, compute dx~/dxi~ = -(dy/dx~)^(-1)(dy/dxi~) in our
            %  coordinates.
            dXtdXit = -ParallelMatrixMultiply(ParallelMatrixInverse(dYdXt), dYdXit);

            if GlobalProcessNewMaslov
                % Count positive eigenvalues of dy/dxi~ (kappa+
                %  in Safarov's notation)
                %
                % More precisely, count eigenvalues greater than
                %  GlobalSignatureTol so we don't get false
                %  positives from theoretically zero eigenvalues
                %  (of which there is always at least one).
                kappaplus = ParallelSignature(dXtdXit, GlobalSignatureTol, 1);
                % Maslov factor is i^(-kappa+(dx~/dxi~))
                maslovFactor = 1i.^(-GlobalProcessMaslovSign*kappaplus);
            else
                % Get signature of dx~/dxi~.
                sgn = ParallelSignature(dXtdXit, GlobalSignatureTol);
                % Maslov factor is exp(pi*i/4*sgn(dx~/dxi~))
                maslovFactor = exp(-GlobalProcessMaslovSign*pi*0.25i*sgn);
            end
            if GlobalDebugMaslov,
                compleximagesc(maslovFactor.');
                title('Maslov factor')
                drawnow
            end
            a = a .* maslovFactor;
        end
        
        function y = CanonicalTransformationY(patch, xt, xit)
        % Calculate the canonical transformation for the local
        %  representation at given (x~,xi~) values.
        
        % All we need to do is change coordinates to (x,xi) then apply the
        %  global patch's canonical transformation.
            [x,xi] = patch.diffeo.InverseTransform(xt,xit);
            y = patch.globalPatch.CanonicalTransformationY(x,xi);
        end
    end
end