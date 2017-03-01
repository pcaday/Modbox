function ReferenceComparisonPlot(vars)

global GlobalFunctionCxPlot

assignall(vars);

%%% Create a summary display
clf

sa = nan(4,1);
sa(1) = subplot(2,2,1);

uEffective.plot
normtitle('Effective input', uEffective)

hold on
extraDrawFunc(1);
hold off

sa(2) = subplot(2,2,2);
Au.plot
normtitle('Computed FIO', Au)

hold on
extraDrawFunc(2);
hold off

colormap copper

drawnow;

    
if ~isempty(Au_ref)
    % Compute error.
    e = Au.copy;
    e.f = Au.f - Au_ref.f;
    
    % Take absolute value of error if not making a complex plot.
    %
    % If GlobalFunctionCxPlot is empty, that's equivalent to it being
    % false.
    if isempty(GlobalFunctionCxPlot) || ~GlobalFunctionCxPlot,
        e.f = abs(e.f);
    end
    
    sa(3) = subplot(2,2,3);
    e.plot
    normtitle('Error', e)
    
    hold on
    extraDrawFunc(3);
    hold off

    sa(4) = subplot(2,2,4);
    Au_ref.plot
    normtitle('Reference FIO', Au_ref)
    
    hold on
    extraDrawFunc(4);
    hold off
end


% Add title over all the plots if one is provided.
if ~isempty(bigTitle),
    suptitle(bigTitle);
end

% Then, add colorbars.
for i = 1:4
    if ~isnan(sa(i))
        axes(sa(i)); %#ok<LAXES>
    end
    colorbar
end

% Link the axes.
%linkaxes(sa(~isnan(sa)), 'xy');

%
%  Display the contributions from each patch separately, if requested.
if isfield(args, 'doPieces') && args.doPieces,
    figure;
    colormap copper
    PlotGallery(Au_patches, true);
end

end





% For each field in the structure vars, create a variable in the
%  caller's workspace.
function assignall(vars)
    for fn = fieldnames(vars).'
        assignin('caller', fn{1}, vars.(fn{1}));
    end
end
