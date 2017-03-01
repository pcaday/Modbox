%#ok<*UNRCH>

%
%   Initialize global variables used by the algorithm.
%

global GlobalPSWFCache              % used in Modified2Box
global GlobalProcessVerbosity       % controls debug output for ProcessedFIOPatch
global GlobalNewProcessVerbosity    % controls debug output for ProcessedLocalFIOPatch
global GlobalProcessNoDyDxFactor    % used in ProcessedFIOPatch
global GlobalProcessMaslovSign      % used in ProcessedFIOPatch
global GlobalProcessNewMaslov       % used in ProcessedFIOPatch
global GlobalDebugMaslov            % used in ProcessedFIOPatch
global GlobalM2BVerbosity           % controls debug output for Modified2Box
global GlobalWeightDiffeosVerbosity % controls debug output for WeightDiffeos
global GlobalWeightDiffeosFakeSVs   % used in WeightDiffeos
global GlobalGeoXRayVerbosity       % controls debug output for GeodesicXRay2DPatch
global GlobalSARReflectVerbosity    % controls debug output for SARReflectFIOPatch
global GlobalHalfWaveVerbosity      % controls debug output for VarHalfWavePatch
global GlobalRefVarHalfWaveVerbosity% controls debug output for ReferenceVarHalfWave
global GlobalRefVarHalfWaveDt       % used in ReferenceVarHalfWave
global GlobalEvalLPVerbosity        % used in EvaluateLowPassFIO
global GlobalAlphaThreshold         % used in Modified2Box
global GlobalUkHatThreshold         % currently unused
global GlobalCWarningThreshold      % used in ProcessedFIOPatch
global GlobalCutoffFeather          % currently unused (was used in old WeightDiffeos)
global GlobalBadDiffeoThreshold     % used in WeightDiffeos
global GlobalDiffeoThreshLow        % used in LocalFIOPatch
global GlobalDiffeoThreshHigh       % used in LocalFIOPatch
global GlobalWeightingType          % used in WeightDiffeos
global GlobalGmaxMinXiRadius        % used in ComputeGTilde
global GlobalFreqNoiseFloor         % used in Modified2Box
global GlobalUkHatBug               % used in Modified2Box
global GlobalSignatureTol           % tolerance for signature calculations
global GlobalUseMEX                 % used for wave-based reflection algorithm
                                    %  (ReferenceSARReflection)
global GlobalFunctionCxPlot         % controls whether we display complex
                                    %  functions in color (used in
                                    %  Function.plot)
global GlobalWedgeThreshold         % used in SqrtWedgeCreator2D
global GlobalAmplitudeDoXCutoff     % used in ProcessedLocalFIOPatch
global GlobalAmplitudeSuppThreshold % used in ProcessedLocalFIOPatch
                                    
GlobalPSWFCache = PSWFCache();      % Cache for prolate spheroidal wave functions
GlobalPSWFCache.pswfTol = 1e-5;
GlobalPSWFCache.pswfK = 40;

GlobalAlphaThreshold = 1e-3;
GlobalUkHatThreshold = 1e-5;
GlobalCWarningThreshold = 512;      % c values larger than this trigger a warning.

% Weighting parameters:


% Weighting method to use, 'multipullback' or 'onesource'
GlobalWeightingType = 'multipullback';
GlobalCutoffFeather = 0.2;
GlobalBadDiffeoThreshold = 0.5;

% New thresholding:
tauW = 0.5;
GlobalDiffeoThreshLow = 1 - tauW;
GlobalDiffeoThreshHigh = 1 - tauW * 0.25;

GlobalGmaxMinXiRadius = 0.08;

GlobalUkHatBug = false;             % Emulate a bug that was in Modified2Box
GlobalFreqNoiseFloor = 0;           % Where |\hat u| is less than this value, it is considered zero.

GlobalRefVarHalfWaveDt = 0.002;

% The following two flags are for testing the components of
%  transforming principal symbol into amplitude.
% Set the first flag to true and the second to 0 to effectively
%  use the principal symbol function as the amplitude.
% For normal use, the first should be false, and the second should
%  be 1.
% 
GlobalProcessNoDyDxFactor = false;  % Don't include |dy/dx| factor when computing amplitude from principal symbol.
GlobalProcessMaslovSign = 1;        % Sign in exponent of Maslov factor when computing amplitude from principal symbol.
                                    %  Should be 1.
                                    
GlobalProcessNewMaslov = false;     % Set to true to use Safarov's singular
                                    %  symbol; false to use the old Maslov
                                    %  code.

GlobalSignatureTol = 1e-4;

% The verbosity flags control how much debug output
%  is emitted by various functions.
% 1 is typically some status messages, while higher
%  values are used for more intensive debugging.
GlobalProcessVerbosity = 1;
GlobalNewProcessVerbosity = 1;
GlobalM2BVerbosity = 1;
GlobalWeightDiffeosVerbosity = 1;
GlobalGeoXRayVerbosity = 1;
GlobalSARReflectVerbosity = 1;
GlobalHalfWaveVerbosity = 1;
GlobalRefVarHalfWaveVerbosity = 1;
GlobalDebugMaslov = 0;
GlobalEvalLPVerbosity = 1;

% Used by the old (wave-based) reflection code
GlobalUseMEX = false;

% Set to 'abrupt-stripes' or 'smooth' for debugging to fake the minimum
%  singular values for WeightDiffeos (to control the cutoffs we get).
% Set to empty to use the real minimum singular values.
GlobalWeightDiffeosFakeSVs = '';

% Pseudopolar support
UsePPDFT(false);

% Enable/disable color plotting for complex values
GlobalFunctionCxPlot = true;

% Threshold for wedges, to eliminate small values.
GlobalWedgeThreshold = 0;

% True to cut off the principal symbol outside of the domain in x.
GlobalAmplitudeDoXCutoff = true;

% The amplitude at a point must be at least this large for the point to
%  count as being in the amplitude's support.
GlobalAmplitudeSuppThreshold = 1e-4;
