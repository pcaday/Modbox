% "Reference" SAR reflection operator computation using the 
%   wave solver.
%
% f and Uf are Function objects, while curve has to be a Curve
%  object (not a curveInfo structure...)
function Uf = ReferenceSARReflect(f, curveInfo, curve, maskInside, showprogress)
    if nargin < 5, showprogress = true; end
    
    % Get the grid and add a PML to it.
    xyg = curve.grid;
    xyg.pml = PML(xyg, min(f.grid.maxes-f.grid.mins) / 8, 10);
    
    % Translate f to an XYFunc.
    xyf = XYFunc(xyg);
    xyf.f = real(f.f).';      % who cares about the complex part? :)
    
    % Run the reflection algorithm
    TWParams.curve = curve;
    TWParams.ic = xyf;
    TWParams.progressFlag = showprogress;
    TWParams.tCFunc = @tLinearCutoff;
    
    if showprogress, figure; end
    Uxyf = waveReconstruct(TWParams);
    if showprogress, close; end
    
    % Convert the result from an XYFunc to a Function.
    Uf = Function.WithValues(f.grid, Uxyf.f.');
    
    % If the curve has an inside, set the result to zero inside
    %  the curve
    if maskInside && curveInfo.hasInside,
        x = f.grid.AllPoints;
        if maskInside == 1,
            Uf.f = Uf.f .* (~curveInfo.lrTestF(x));
        elseif maskInside == -1,
            Uf.f = Uf.f .* (curveInfo.lrTestF(x));
        end            
    end


% tLinearCutoff: Linear time cutoff function
%
% Inputs:
%      t: input time (can be a vector)
%   Tmax: max time on the grid
%
% Outputs:
%      c: cutoff value (between 0 and 1)
%
function cutoff = tLinearCutoff(t, Tmax)
	cutoff = HalfWindow(t, Tmax, Tmax/2, 'smooth');
end

end