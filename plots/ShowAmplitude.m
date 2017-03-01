% Show the amplitude function for a processed FIO patch,
%  using ShowCotangentFun
%
%  SHOWAMPLITUDE(patch, ...)
%  SHOWAMPLITUDE(fio, ...)
%
% The first argument can be a ProcessedFIOPatch object,
%  or a ProcessedFIO object (in which case the amplitude for
%  each patch will be shown in separate figures, numbered
%  41, 42, ...).
% Any extra arguments (optional) are passed to ShowCotangentFun.

function ShowAmplitude(obj, varargin)
    if isa(obj, 'ProcessedFIOPatch') || isa(obj, 'ProcessedLocalFIOPatch')
        ShowCotangentFun([obj.amplitude{:}], obj.wedgeCreator, obj.inGrid, varargin{:});
    elseif isa(obj, 'ProcessedFIO')
        for k = 1:length(obj.patches)
            figure(40+k);
            ShowCotangentFun([obj.patches(k).amplitude{:}], obj.patches(k).wedgeCreator, obj.patches(k).inGrid, varargin{:});
        end
    end
end