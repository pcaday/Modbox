% Plot complex-valued images in color.
% 
% The phase is mapped to color/hue (0 degrees to red, 60 degrees to yellow,
%  120 degrees to green, and so on), and the magnitude is mapped to
%  brightness.
%
% The magnitude range is determined automatically from the range of
%  magnitudes in the image, or can be specified manually by including an
%  rlim parameter at the end of the argument list.
%
% COMPLEXIMAGESC supports all the syntax supported by imagesc:
%
% h = COMPLEXIMAGESC(C)
% h = COMPLEXIMAGESC(x,y,C)
% h = COMPLEXIMAGESC(..., rlim)
% h = COMPLEXIMAGESC(..., 'Property1', Value1, ... 'PropertyN', ValueN)
%
% If no image C is provided, COMPLEXIMAGESC is the same as image.

function ih = compleximagesc(varargin)
    rlim = [];
    switch nargin,
        case 0,
            C_index = [];
        case 1,
            C_index = 1;
        otherwise,
            % Determine if last input is a radius limit
            % Taken from imagesc ver. 5.11.4.5, R2012b.
            if isequal(size(varargin{end}),[1 2])
                str = false(length(varargin),1);
                for n=1:length(varargin)
                    str(n) = ischar(varargin{n});
                end
                str = find(str);
                if isempty(str) || (rem(length(varargin)-min(str),2)==0),
                    rlim = varargin{end};
                    varargin(end) = []; % Remove last cell
                else
                    rlim = [];
                end
            else
                rlim = [];
            end
            
            % Determine whether we have x,y arguments
            if length(varargin) > 2 && ...
                    ~ischar(varargin{1}) && ...
                    ~ischar(varargin{2}),
                % Yes, C is argument 3
                C_index = 3;
            elseif ~ischar(varargin{1}) 
                % No x,y arguments and C is argument 1
                C_index = 1;
            else
                % No x,y or C arguments.
                C_index = [];
            end
    end

    if ~isempty(C_index),
        % Convert complex C to RGB image.
        C = varargin{C_index};
        
        if ~ismatrix(C), error('The image must be a matrix'); end

        % Convert complex image C to argument and magnitude.
        arg = angle(C);
        r = abs(C);

        if isempty(rlim),
            % Find magnitude limits if not specified.
            rmin = 0;
            rmax = max(r(isfinite(r)));
            if isempty(rmax), rmax = 1; end
        else
            rmin = rlim(1); rmax = rlim(2);
        end
        
        % Argument becomes hue
        h = mod(arg/(2*pi), 1);
        % Saturation is always 1.
        s = ones(size(h));
        % Magnitude becomes value
        v = (r - rmin) / (rmax - rmin);
        % Clamp value between 0 and 1 outside the magnitude range.
        v(r >= rmax) = 1.;
        v(r <= rmin) = 0.;
        
        hsv = cat(3, h, s, v);

        varargin{C_index} = hsv2rgb(hsv);
    end
    
    hh = image(varargin{:});
    
    if nargout > 0,
        ih = hh;
    end    
end