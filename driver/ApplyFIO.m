% Apply an FIO to a test function.
%
%   Au = ApplyFIO(args)
%   [Au,vars] = ApplyFIO(args)
%
% The arguments are specified as fields in a structure. Possible fields are
%             'n': grid points per dimension
%         'image': string specifying the test image (see TestImage.m)
%        'inGrid': grid for the input
%           'fio': string specifying the FIO to use
%                        'identity': the identity
%                           'riesz': a Riesz operator
%                        'waveprop': constant-speed wave propagator
%                      'sarreflect': the SAR reflection operator
%                      'sarinverse': the (one-sided) inverse circular Radon
%                                     transform
%                           'radon': the 2D Radon transform
%                            'xray': the 2D geodesic X-ray transform
%                        'halfwave': variable-speed half-wave equation
%                                     propagator.
%       'diffeos': string specifying which diffeos to use
%                        'identity': just the identity
%                             'one': a single quadratic change of variables
%                             'two': two quadratic changes of variables
%               'two-plus-identity': same as 'two', plus the identity
%                            'four': four quadratic changes of variables
%                                      (two in each variable, bending in
%                                      opposite directions)
%              'four-plus-identity': same as 'four', plus the identity.
%       'lowPass': true/false: true to apply low-pass algorithm
%                   (set maskMode to 'separate' if low-pass algorithm
%                    desired.)
%                   Default: false.
%      'highPass': true/false: true to apply high-pass (standard) algorithm
%                   Default: true.
%      'maskMode': masking mode; controls the handling of low-frequency
%                  components
%                        'separate': use this with the low-pass algorithm
%                         'no-mask': use this for handling the high-pass
%   'pseudopolar': true/false: true to use pseudopolar DFT; false for
%                   Cartesian DFT.
%         'class': data storage class, 'single' or 'double'
%     'coneRatio': controls the size of cones (see SqrtWedgeCreator)
%       'zeroPad': zero-padding factor
% 'wcrSmoothness': smoothness for cones (see SqrtWedgeCreator)
%     'angleOSLP': angular oversampling factor for low-pass algorithm
%    'radialOSLP': radial oversampling factor for low-pass algorithm
%  'adjustDegree': If set, the given amount will be added to the degree of
%                   the FIO during the algorithm, and 
%                   (will be countered by a Riesz operator afterwards)
%    'bendFactor': amount of coordinate bending
%         'alpha': exponent of Riesz potential (fio = 'riesz')
%             'T': propagation time (fio = 'halfwave', 'waveprop')

function [Au, vars] = ApplyFIO(args)
%#ok<*NUSED>

% Anything declared as global here can be modified by the arguments
%  structure.
global GlobalUsePseudopolarFFT
global GlobalDFTGridRoutine
global GlobalBadDiffeoThreshold 
global GlobalDiffeoThreshLow
global GlobalDiffeoThreshHigh

% Initialize global variables if necessary
if isempty(GlobalUsePseudopolarFFT),
    init
end


% Convert parameter-value pairs to a structure.
%assert(mod(nargin,2) == 0, 'Arguments to ApplyFIO are expected to be parameter-value pairs.');
%args = cell2struct(varargin(1:2:end), varargin(2:2:end), 2);

% Default values for arguments
n = 64;
inGrid = [];
image = 'square-and-circle';
fio = 'radon';
diffeos = 'two';
lowPass = false;
highPass = true;
maskMode = 'no-mask';
class = 'single';
coneRatio = [];                % [] for default
curve = 'default';
metric = 'soundspeed';
soundSpeed = 'deep-well';
bendFactor = 1;                 % Amount of coordinate bending.
pseudopolar = false;
zeroPad = 2;
wcrSmoothness = 0.25;
alpha = -1;                     % Exponent for 'riesz'
T = 0.2;                        
periodic = false;
maskInside = 0;                 % Masking mode for 'sarreflect'

% Oversampling factors for low-pass algorithm
angleOSLP = 1;                  % angular oversampling
radOSLP = 1;                    % radial oversampling

% For testing, override degree of FIO by setting
%  this variable to a scalar value. The value will be added
%  to the degrees of all FIOs.
adjustDegree = [];


assignall(args);

% Pseudopolar DFT support.
% The zero-padding factor must be 1 for pseudopolar (can't remember why)
UsePPDFT(pseudopolar);
if pseudopolar,
    zeroPad = 1;
end

% Create default grid if none specified
if isempty(inGrid),
    inGrid = Grid([-1 -1], [1 1], [n n]);
end

% Make sure grid is periodic if specified and not otherwise.
inGrid.periodic(:) = periodic;

%%%
%%% Create the FIO and specify test image and wedge creators
%%%

inGrid = convertClass(inGrid, class);

refFunc = [];
bigTitle = [];
extraDrawFunc = @(subplot_no) void() ;
testImageFunc = @TestImage;               % By default: can be overriden as needed.


switch fio
    case 'identity'
        A = IdentityFIO(2);
        refFunc = @(f) f;
        outGrid = inGrid;
    case 'riesz'
        A = SinglePatchFIO(RieszPotentialPatch(alpha));
        refFunc = @(f) ReferenceRieszPotential(f, alpha);
        outGrid = inGrid;
    case 'waveprop' 
        A = SinglePatchFIO(WavePropFIOPatch(T), 1);
        refFunc = @(f) ReferenceWaveProp(f, T);
        outGrid = inGrid;
        bigTitle = ['Propagating ' image ' for time ' num2str(T)];
        % It's periodic, so...
        %  zeroPadding = 1;
    case 'sarreflect'
        if strcmp(curve, 'default'), curve = 'parabola'; end
        
        curveObj = [];
        switch curve
            case 'vline'
                x_intercept = 0;
                curveInfo = VerticalLineCurveInfo(x_intercept);
                refFunc = @(f) Function.WithValues(f.grid, flipdim(f.f, 1));
            case 'circle'
                % CurveInfo structure for the new algorithm
                rad = 0.55;
                curveInfo = CircleCurveInfo(rad);
                % Corresponding Curve structure for the old reference
                %  reflection algorithm
                curveObj = Curve.Circle(GridToXYTGrid(inGrid), rad);
            case 'lower-half-circle'
                rad = 0.55;
                thetas = [pi 2*pi];
                curveInfo = CircleCurveInfo(rad, thetas);
                curveObj = Curve.Circle(GridToXYTGrid(inGrid), rad, thetas);
            case 'parabola'
                inGrid.mins = [-3 -3];
                inGrid.maxes = [3 3];
                a = 1;
                curveInfo = ParabolaCurveInfo(a);
                curveObj = Curve.Parabola(GridToXYTGrid(inGrid), -3, 3);
            case 'shallow-parabola'
                inGrid.mins = [-0.5 -0.5];
                inGrid.maxes = [0.5 0.5];
                a = 1;
                curveInfo = ParabolaCurveInfo(a);
                curveObj = Curve.Parabola(GridToXYTGrid(inGrid), -0.5, 0.5);                
        end
        bigTitle = [image ' reflected across ' curve];
        
        A = SARReflectFIO(curveInfo, maskInside);
        outGrid = inGrid;
        
        if ~isempty(curveObj),
            refFunc = @(f) ReferenceSARReflect(f, curveInfo, curveObj, maskInside, false);
            extraDrawFunc = @(sub) curveObj.drawDiscretized();
        end
    case 'sarinverse'
        if strcmp(curve,'default'), curve = 'tallparabola'; end
        
        side = +1;                      % +1 = left, -1 = right
        switch (curve)
            case 'vline'
                x_intercept = 0;
                curveInfo = VerticalLineCurveInfo(x_intercept);
            case 'hline'
                y_intercept = 0;
                curveInfo = LineCurveInfo([0 y_intercept], [1 0]);
            case 'circle'
                rad = 0.55;
                curveInfo = CircleCurveInfo(rad);
                image = 'center-dot';
            case 'parabola'
                a = 1;
                curveInfo = ParabolaCurveInfo(a);
            case 'tallparabola'
                a = 2;
                curveInfo = ParabolaCurveInfo(a);
        end
        A = SinglePatchFIO(SARInverseFIOPatch(curveInfo, side),0,...
                           Grid([-inf 0], [inf inf], [2 2]));
        
        origGrid = inGrid;
        origImage = TestImage(origGrid, image);
        
        tmax = norm(origGrid.maxes-origGrid.mins);
        smin = curveInfo.sRange(1);
        smax = curveInfo.sRange(2);
        if ~isfinite(smin), smin = -1; end  % Substitute defaults for infinite curves
        if ~isfinite(smax), smax =  1; end
        if smin >= smax, smax = smin + 2; end   % Avoid empty s range after substituting defaults.
        
        % Create (s,t) and output grids.
        stGrid = Grid([smin 0], [smax tmax], [n n]);
        outGrid = origGrid;
        inGrid = stGrid;
        
        % Apply forward SAR transform to original image.
        fprintf('Applying forward SAR transform... ');
        uOrig = TestImage(origGrid, image);
        Ru = ReferenceSARForward(uOrig, stGrid, curveInfo);
        fprintf('done.\n');
        
        % Pass the SAR-transformed image forward in the "test image"
        %  function handle, and the original image to the reference
        %  function.
        testImageFunc = @(varargin) Ru;
        refFunc       = @(varargin) uOrig; 
        extraDrawFunc = @(sub)      (sub ~= 1) && DrawCurve(curveInfo);   % Equivalent to: if sub ~= 1, DrawCurve(curveInfo); end
                                                                            % control statements aren't allowed in anonymous functions.        
        % Use more wedges.
        coneRatio = 16;
    case 'radon'
        A = Radon2DFIO();
        outGrid = Grid([-1.6 0], [+1.6 2*pi], [n n], [false true]);
        refFunc = @(f) ReferenceRadon2D(f, outGrid);        
    case 'oldradon'
        A = SinglePatchFIO(Radon2DPatch(+1), 2,...
                            Grid([-inf 0], [inf 2*pi], [2 2], [false true]));
        outGrid = Grid([-1.6 0], [+1.6 2*pi], [n n], [false true]);
        refFunc = @(f) ReferenceRadon2D(f, outGrid);
    case 'xray'
        if strcmp(curve, 'default'), curve = 'circle'; end
        kappa = 0;
        
        switch curve
            case 'circle'
                rad = 1;
                curveInfo = CircleCurveInfo(rad);
        end
        
        m = ConstantCurvature2DMetric(kappa);
        
        outGrid = Grid([0 0], [2*pi pi], [n n], [true false]);
        
        A = SinglePatchFIO(GeodesicXRay2DPatch(m,curveInfo), 2, outGrid);
        
        if kappa == 0, A.patches.srds(1) = []; end
        
        refFunc = [];
    case 'halfwave'
        dims = 2;
        c = [];
        dt = [];
        
        switch metric
            case 'flat'
                kappa = 0;
                m = ConstantCurvature2DMetric(kappa);
            case 'constant'
                speed = 2;
                s = ConstantSoundSpeed(speed);
                m = EuclideanConformalMetric(s);
            case 'soundspeed'
                switch soundSpeed
                    case 'well'
                        s = SoundSpeedWell([],[],[],[],[],inGrid);
                    case 'deep-well'
                        s = SoundSpeedWell(2,[],-1.5,[],[],inGrid);
                end
                
                m = EuclideanConformalMetric(s);
            case 'soundspeed-discrete'
                if strcmp(soundSpeed, 'well'),
                    soundSpeed = 'soundspeed-well';
                end
        
                % Expand the grid to create the sound speed grid
                cGrid = inGrid.copy();
                cGrid.mins = cGrid.mins * 2;
                cGrid.maxes = cGrid.maxes * 2;
                cGrid.ns = cGrid.ns * 2;
                c = TestImage(cGrid, soundSpeed);
                c.f = 1./c.f;
                m = EuclideanConformalMetric.DiscreteMultiplier(c);
            case 'constantcurvature'
                kappa = -4;
                m = ConstantCurvature2DMetric(kappa);
        end
        
        outGrid = inGrid;
        
        A = SinglePatchFIO(VarHalfWavePatch(dims,m,T), periodic*2, outGrid);
        
        if ~isempty(dt), A.patches.dt = dt; end
        
        if isempty(c), c = m; end
        
        refFunc = @(f) ReferenceVarHalfWave(f,T,c);
        
        bigTitle = ['Variable-speed wave propagation for time ' num2str(T)];
end

% Choose the list of coordinate choices.
diffList = {};

switch diffeos
    case 'identity'
        diffList = {IdentityDiffeo(2)};
    case 'one'
        diffList = {Bend2DDiffeo(1,bendFactor)};
    case 'one-x1'       % Same as 'one'
        diffList = {Bend2DDiffeo(1,bendFactor)};
    case 'one-x2'       % Bend x_2 only
        diffList = {Bend2DDiffeo(2,bendFactor)};
    case 'two'
        diffList = {Bend2DDiffeo(1,bendFactor)...
            Bend2DDiffeo(2,bendFactor)};
    case 'two-plus-identity'
        diffList = {IdentityDiffeo(2)...
            Bend2DDiffeo(1,bendFactor)...
            Bend2DDiffeo(2,bendFactor)};
    case 'four'
        diffList = {Bend2DDiffeo(1,+bendFactor)...
            Bend2DDiffeo(1,-bendFactor)...
            Bend2DDiffeo(2,+bendFactor)...
            Bend2DDiffeo(2,-bendFactor)};
    case 'four-plus-identity'
        diffList = {IdentityDiffeo(2)...
            Bend2DDiffeo(1,+bendFactor)...
            Bend2DDiffeo(1,-bendFactor)...
            Bend2DDiffeo(2,+bendFactor)...
            Bend2DDiffeo(2,-bendFactor)};
    otherwise
        error('Unknown choice of SRDs.');
end


switch maskMode
    case 'no-mask'
        SWCParams = {coneRatio [] []};
    case 'mask-input'
        SWCParams = {coneRatio [] []};
    case 'wedge-creator-mask'
        SWCParams = {coneRatio 0.08 inf};
    case 'wcr-mask-clear'
        SWCParams = {coneRatio 0.08 inf};
    case 'separate'
        SWCParams = {coneRatio 0.08 inf};        
    otherwise
        error('Unknown mask mode');
end

% Convert grid again in case it was changed by the FIO-specific code.
inGrid = convertClass(inGrid, class);

if GlobalUsePseudopolarFFT,
    wedgeCreatorCreator = @(grid) PPSqrtWedgeCreator(grid, SWCParams{:}, wcrSmoothness, class);
else
    wedgeCreatorCreator = @(grid) SqrtWedgeCreator2D(grid, SWCParams{:}, wcrSmoothness, class);
end

wcr = wedgeCreatorCreator(GlobalDFTGridRoutine(inGrid));

% Override degree if requested for debugging
if ~isempty(adjustDegree)
    for p = A.patches(:).'
        p.degree = p.degree + adjustDegree;
    end
end

%%% Convert to local patches and process the FIO


tic
Ap = ProcessedFIO(A, inGrid, outGrid, diffList, wedgeCreatorCreator, zeroPad, class);
T = toc;
fprintf('\n** FIO processed in %f seconds\n\n', T);

if lowPass,
    ApLP = ProcessedLowPassPatch(Ap, angleOSLP, radOSLP);
end

%%% Run the modified box algorithm

u = testImageFunc(inGrid, image);

switch maskMode
    case 'mask-input'
        wedgeCreatorCreator_masked = @(grid) SqrtWedgeCreator2D(grid,[],0.08,1);
        wcr_masking = wedgeCreatorCreator_masked(GlobalDFTGridRoutine(u.grid));
        
        uHP = ApplyWedgeCreatorMask(u, wcr_masking);
    case 'wedge-creator-mask'
        uHP = ApplyWedgeCreatorMask(u, wcr);
    case 'wcr-mask-clear'
        for i = 1:length(Ap.patches)
            Ap.patches(i).wedgeCreator.mask(:) = 1;
        end
        wcr.mask(:) = 1;
        uHP = u;
    case 'separate'
        uHP = u;
    otherwise
        uHP = u;
end


if highPass,
    tic
    [Au,Au_patches] = Modified2Box(Ap, uHP, outGrid, class, zeroPad);
    T = toc;
    fprintf('\n** FIO applied in %f seconds\n\n', T);
else
    Au = Function.Zeros(outGrid, class);
    Au_patches = {};
end

% Compute the lowpass parts

if lowPass,
    tic
    [AuLP, AuLP_patches] = EvaluateLowPassFIO(ApLP, u, outGrid, class, zeroPad);
    T = toc;
    Au_patches = [Au_patches(:); AuLP_patches(:)];
    Au.f = Au.f + AuLP.f;
    fprintf('\n** Low pass computations in %f seconds\n\n', T);
end




% Compute 'effective input'
switch maskMode
    case 'no-mask'
        uEffective = u;
    case 'mask-input'
        uEffective = uHP;
    case 'wedge-creator-mask'
        % If lowpass algorithm is used, effective input is same as input.
        % If not, it's the high-frequency portions.
        if lowPass
            uEffective = u;
        else
            uEffective = uHP;
        end
        %uEffective = ApplyWedgeCreatorMask(uHP, wcr);
    case 'wcr-mask-clear'
        uEffective = u;
    case 'separate'
        uEffective = u;
end

% Also compute reference output if we know how.
if ~isempty(refFunc)
    Au_ref = refFunc(uEffective);
end

% If degree overridden, undo by applying appropriate Riesz operator
if ~isempty(adjustDegree),
    Au = ReferenceRieszPotential(Au, adjustDegree);
end

% Save all our local variables for the caller if requested.
if nargout > 1,
    vars = saveall(who);
end

end





% For each field in the structure vars, create a variable in the
%  caller's workspace.
function assignall(vars)
    for fn = fieldnames(vars).'
        assignin('caller', fn{1}, vars.(fn{1}));
    end
end

% Save all variables in the caller's workspace listed by name in w
%  to a structure vars.
% e.g. saveall(who) to retrieve all variables.
function vars = saveall(w)
    for fn = w.'
        vars.(fn{1}) = evalin('caller', fn{1});
    end
end

% Empty function
function void
end