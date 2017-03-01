% Simple slider-creating function.
%
%    slider(fun, [min max])
%
%  fun is a function handle taking a single argument (the slider value)
%   which should update the plot appropriately for that value.
%  [min max] is the range of values for the slider.
%
%   h = slider(...)
%
%  Returns the handle to the created slider uicontrol.
%
function hh = slider(drawfun, range)

% create slider
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20],'max',range(2),'min',range(1),'Value',range(1));

% add a function that is called when the slider is moved
addlistener(h, 'ContinuousValueChange', ...
    @(hObject, event) sliderUpdate(hObject, event, drawfun));

% This function is called when the slider is changed I think...
%  the ContinuousValueChange doesn't get called when the user
%  releases the mouse button, so this callback catches when that
%  happens.
set(h, 'Callback', @(hObject, event) sliderUpdate(hObject, event, drawfun));

if nargout > 0, hh = h; end

sliderUpdate(h, [], drawfun);

end


% Function that gets called when slider moved
function sliderUpdate(hObject,~,drawfun)
% Get the current value of the slider
val = get(hObject,'Value');

% Call the user's function and draw.
drawfun(val);
drawnow;

end