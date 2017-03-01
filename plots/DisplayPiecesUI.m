% Create a UI where different pieces (from the modified box algorithm)
%  can be overlaid, enabled, or disabled.
function DisplayPiecesUI(Au_pieces)
assert(iscell(Au_pieces));
npieces = length(Au_pieces);

f = figure;

set(f, 'Toolbar', 'none', 'WindowStyle', 'normal');

iconmask = zeros(16);
toggles = zeros(npieces, 1);
for i = 1:npieces
    a = mod(i-1,4);
    b = (i-1-a)/4;
    iconmask(b*4+1:b*4+3, a*4+1:a*4+3) = 1;
    
    icon = bsxfun(@times, iconmask, reshape(Color(i), [1 1 3]));
    
    toggles(i) = uitoggletool('CData', icon, 'ClickedCallback', @RefreshCallback);
end

% Get the sum
piecesum = 0;
for i = 1:npieces
    piecesum = piecesum + Au_pieces{i}.f;
end

global GlobalFunctionCxPlot


% Get the mins and maxes of
%   (real plots) the sum's real part
%   (complex plots) the sum's absolute value.
if GlobalFunctionCxPlot,
    piecesum = abs(piecesum);
else
    piecesum = real(piecesum);
end
max_sum = max(piecesum(:));
min_sum = min(piecesum(:));

% Renormalize so images between 0 and 1.
%     for i = 1:npieces
%         Au_pieces{i} = (Au_pieces{i} - min_sc) / (max_sc - min_sc);
%     end




colormap copper

Refresh








    function RefreshCallback(~,~)
        Refresh;
    end

    function Refresh
        im = Function.Zeros(Au_pieces{1}.grid);
        for j = 1:npieces
            if isequal(get(toggles(j),'State'), 'on')
                %                color_image = bsxfun(@times, Au_pieces{j},...
                %                    reshape(Color(j), [1 1 3]));
                %                im = im + color_image;
                im.f = im.f + Au_pieces{j}.f;
            end
        end
        im.plot([min_sum max_sum]);
        colorbar
    end
end


% Get the i'th color.
function c = Color(i)
colors =    [1 0 0; ...
    0 0 1;  ...
    0 1 0;  ...
    1 1 0;  ...
    1 0 1;  ...
    1 0.5 0;...
    0.5 0 1;...
    0 1 1;  ...
    0.5 1 0];
c = colors(mod(i-1,9)+1,:).';
end

