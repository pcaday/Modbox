
function [Aus,errs] = TestParameterVariants(baseArgs,vdArgs,variant,relError)

if nargin < 2, vdArgs.dummy = []; end
if nargin < 3, variant = 'grabbag'; end
if nargin < 4, relError = false; end

switch variant
    case 'grabbag'
        args = repmat({baseArgs}, 5, 4);

        % Row 1: zero-padding
        args{1,1}.zeroPad = 1;
        args{1,2}.zeroPad = 2;
        args{1,3}.zeroPad = 3;
        args{1,4}.zeroPad = 4;

        % Row 2: cone count
        args{2,1}.coneRatio = 1;
        args{2,2}.coneRatio = 2;
        args{2,3}.coneRatio = 4;
        args{2,4}.coneRatio = 16;

        % Row 3: cone smoothness
        args{3,1}.wcrSmoothness = -1;
        args{3,2}.wcrSmoothness = 0.05;
        args{3,3}.wcrSmoothness = 0.25;
        args{3,4}.wcrSmoothness = 0.5;

        % Row 4: coordinate bending
        args{4,1}.bendFactor = 0.25;
        args{4,2}.bendFactor = 0.5;
        args{4,3}.bendFactor = 2;
        args{4,4}.bendFactor = 3;

        % Row 5: weighting smoothness
        args{5,1}.GlobalBadDiffeoThreshold = 0.1;
        args{5,2}.GlobalBadDiffeoThreshold = 0.25;
        args{5,3}.GlobalBadDiffeoThreshold = 0.5;
        args{5,4}.GlobalBadDiffeoThreshold = 0.9;

    case 'grabbag2'
        if isfield(baseArgs, 'pseudopolar')
            pp = baseArgs.pseudopolar;
        else
            pp = false;
        end
        
        args = repmat({baseArgs}, ~pp + 2, 4);

        if ~pp,
            % Row 1: zero-padding
            args{1,1}.zeroPad = 1;
            args{1,2}.zeroPad = 2;
            args{1,3}.zeroPad = 3;
            args{1,4}.zeroPad = 4;
        end
        
        % Row 2: cone count
        args{~pp+1,1}.coneRatio = 1;
        args{~pp+1,2}.coneRatio = 2;
        args{~pp+1,3}.coneRatio = 4;
        args{~pp+1,4}.coneRatio = 16;

        % Row 3: coordinate bending
        args{~pp+2,1}.bendFactor = 0.25;
        args{~pp+2,2}.bendFactor = 0.5;
        args{~pp+2,3}.bendFactor = 2;
        args{~pp+2,4}.bendFactor = 3;
        
    case 'zeropad'
        args = repmat({baseArgs}, 1, 4);

        args{1,1}.zeroPad = 1;
        args{1,2}.zeroPad = 2;
        args{1,3}.zeroPad = 3;
        args{1,4}.zeroPad = 4;
        
    case 'lowpass'
        args = repmat({baseArgs}, 1, 4);

        baseArgsLP = baseArgs;
        baseArgsLP.lowPass = true;
        baseArgsLP.maskMode = 'separate';
        
        % PPFFT / lowpass
        %
        % In order: Cartesian DFT, no low pass
        %           Cartesian DFT, low pass
        %           Pseudopolar DFT, no low pass
        %           Pseudopolar DFT, low pass
        args{1,1}.pseudopolar = false;
        args{1,2} = baseArgsLP;
        args{1,2}.pseudopolar = false;
        args{1,3}.pseudopolar = true;
        args{1,3}.zeroPad = 1;
        args{1,3}.wcrSmoothness = 0.25;
        args{1,4} = baseArgsLP;
        args{1,4}.pseudopolar = true;
        args{1,4}.zeroPad = 1;
        args{1,4}.wcrSmoothness = 0.25;
        
    case 'diffeo-thresh'
        args = repmat({baseArgs}, 3, 3);
        
        bads = [0.75 0.5 0.25];
        fracs = [0.5 0.25 0.1];
        
        for i = 1:3
            for j = 1:3
                args{i,j}.GlobalDiffeoThreshLow = 1 - bads(i);
                args{i,j}.GlobalDiffeoThreshHigh = 1 - (bads(i) * fracs(j));
            end
        end
end

if relError,
    fprintf('Computing base variant\n');
    vdArgs.ref = ApplyFIO(baseArgs);
end

[Aus,errs] = ApplyFIOVariantDisplay(args,vdArgs);

end