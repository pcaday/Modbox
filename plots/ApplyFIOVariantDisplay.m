% Run ApplyFIO with different arguments, and display the results
%  in a grid.
%
%   ApplyFIOVariantDisplay(args)
%   ApplyFIOVariantDisplay(args, captions)
%
% args is a 2D cell array of argument structures to pass to ApplyFIO.
% captions (optional) is a 2D cell array of caption strings for the plots.
%

function [Aus,errs] = ApplyFIOVariantDisplay(args, vdArgs)

% Substitute (almost) empty parameter arrays if omitted.
if nargin < 1, args.dummy = []; end
if nargin < 2, vdArgs.dummy = []; end

global GlobalFunctionCxPlot

cxPlot = ~isempty(GlobalFunctionCxPlot) && GlobalFunctionCxPlot;
doColorbars = ~cxPlot ...
    && (~isfield(vdArgs, 'doColorbar') || ~vdArgs.doColorbar);
doErrors = isfield(vdArgs, 'doErrors') && vdArgs.doErrors;
showArgs = isfield(vdArgs, 'showArgs') && vdArgs.showArgs;
if isfield(vdArgs, 'ref'),
    ref = vdArgs.ref;
else
    ref = [];
end

assert(iscell(args), 'Expected cell array.');
[m,n] = size(args);


fOutput = figure;
if doErrors, fError = figure; end

sa = nan(m,n);
Aus = cell(m,n);
errs = cell(m,n);
e = [];

for i = 1:m
    for j = 1:n
        fprintf('Plot (%d,%d)...\n', i, j);
        if showArgs, disp(args{i,j}); end

        [Au,vars] = ApplyFIO(args{i,j});

        figure(fOutput);
        sa(i,j) = subplot(m,n,j+(i-1)*n);
        Au.plot
        
       
        hold on
        vars.extraDrawFunc(2);
        hold off
        
        colormap copper;
        
        if doErrors,
            figure(fError);
        
            e = Au.copy;
            if ~isempty(ref)
                e.f = Au.f - ref.f;
            else
                e.f = Au.f - vars.Au_ref.f;
            end
            
            % Take absolute value of error if not making a complex plot.
            %
            % If GlobalFunctionCxPlot is empty, that's equivalent to it being
            % false.
            if isempty(GlobalFunctionCxPlot) || ~GlobalFunctionCxPlot,
                e.f = abs(e.f);
            end
            subplot(m,n,j+(i-1)*n);
            e.plot;
            
            normtitle('Error', e);
            
            figure(fOutput);
        end
        
        Aus{i,j} = Au;
        errs{i,j} = e;
    end
end


% Add title over all the plots if one is provided.
if ~isempty(vars.bigTitle),
    suptitle(vars.bigTitle);
end

% Then, add colorbars (except for complex plots) as needed
if doColorbars,
    for i = 1:m
        for j = 1:n
            if ~isnan(sa(i,j))
                axes(sa(i,j)); %#ok<LAXES>
            end
            colorbar
        end
    end
end


end
