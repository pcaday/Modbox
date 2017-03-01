function [f,colorimage] = TestImage(gres, kind)
if nargin < 2, kind = 'square-and-circle'; end

% Get the grid
if isa(gres, 'Grid'),
    xg = gres;
    assert(ndims(xg) == 2 && xg.ns(1) == xg.ns(2), 'Input grid must be square and 2D.'); %#ok<ISMAT>
else
    xg = Grid([-1 -1], [1 1], [gres gres]);
end
res = xg.ns(1);
cf = xg.coordFuns;
x = cf{1}; y = cf{2};
nx = (x - xg.mins(1)) / (xg.maxes(1) - xg.mins(1));
ny = (y - xg.mins(2)) / (xg.maxes(2) - xg.mins(2));

i = zeros(res);
colorimage = [];

switch kind
    case 'square-and-circle'
        i = LoadBWImage('TestImage-SquareAndCircle.png');
    case 'text'
        i = LoadBWImage('TestImage-Text.png');
    case 'world'
        i = LoadColorImage('TestImage-World.jpg');
        
        colorimage = i;
        i(:,:,[2 3]) = [];
    case 'iguassu'
        i = LoadColorImage('TestImage-Iguassu.jpg');
        
        colorimage = i;
        i(:,:,[2 3]) = [];
    case 'chameleon'
        i = LoadColorImage('TestImage-Chameleon.jpg');
        
        colorimage = i;
        i = i(:,:,1)+i(:,:,2);
    case 'circle'        
        i = HalfWindow(x.^2+y.^2, (0.27)^2, (0.23)^2, 'smooth');
    case 'offset-circle'        
        i = HalfWindow(x.^2+(y-0.5).^2, (0.27)^2, (0.23)^2, 'smooth');
    case 'square'
        i(max(abs(x),abs(y)) < 0.3001) = 1;
    case 'diamond'
        i((abs(x)+abs(y)) < 0.4001) = 1;
    case 'low-diamond'
        i((abs(x)+abs(y+0.5)) < 0.4001) = 1;
    case 'high-diamond'
        i((abs(x)+abs(y-0.5)) < 0.4001) = 1;
    case 'dash'
        i(abs(x) < 0.75 & abs(y+0.1) < 0.05) = 1;
    case 'slight-low-dash'
        i(abs(x) < 0.75 & abs(y+0.3) < 0.05) = 1;
    case 'big-dash'
        i(abs(x) < 1.75 & abs(y+0.5) < 0.07) = 1;
    case 'low-dash'
        i(abs(x) < 0.75 & abs(y+0.8) < 0.03) = 1;
    case 'tilted-dash'
        angle = 0.5;
        x0 = -0.05;
        y0 = -0.1;
        xt = (x-x0) .* cos(angle) - (y-y0) .* sin(angle);
        yt = (x-x0) .* sin(angle) + (y-y0) .* cos(angle);
        i(abs(xt) < 0.75 & abs(yt+0.1) < 0.05) = 1;
    case 'center-dot'
        i = HalfWindow(x.^2+y.^2, (0.14)^2, (0.11)^2, 'smooth');        
    case 'tiny-dot'
        i = HalfWindow(x.^2+y.^2, (0.05)^2, (0.03)^2, 'smooth');        
    case 'offcenter-dot'
        i = HalfWindow((x-0.3).^2+y.^2, (0.14)^2, (0.11)^2, 'smooth');        
    case 'hollow-ellipse'
        i = FullWindow(4*(y-0.2).^2+12*(x-0.1).^2, 0.8, 0.9, 1.1, 1.2, 'smooth');
    case 'wave-packet'
        iHat = Function.Zeros(xg.DFTGrid);
        wcr = SqrtWedgeCreator2D(iHat.grid,2,0.2,0.5);
        iHat.f = wcr.Wedge(1);
        ifun = iHat.IDFT(xg);
        i = ifun.f;
        i = fftshift(i);
    case 'sound-speed-well'
        c0 = 2;
        kappa = -0.4;
        sigma = 0.5;
        i = c0 + kappa*exp(-(x.^2+y.^2)/sigma^2);
    case 'all-ones'
        i(:) = 1;
    case 'plane-wave'
        i = cos(x * 2 * pi - y * 3 * pi);
    case 'very-smooth'
        iHat = Function.Zeros(xg.DFTGrid);
        c = ceil((size(xg.DFTGrid)+1)/2);
        iHat.f(c(1)+(-2:2), c(2)+(-2:2)) = ...
            [1  3 4 3 1;
             2  0 4 0 2;
             -1 2 3 5 0;
             0 -1 0 3 1;
             2  1 0 2 2];
        i = iHat.IDFT(xg);
        i = fftshift(i.f);
    case 'quarter-circle'
        i = double(nx.^2 + ny.^2 < 0.625^2);
    otherwise
        error('TestImage: I don''t know that image.');
end


% Create the Function object.
f = Function.WithValues(xg, i);





    function i = LoadBWImage(name)
        ihigh = imread(name);
        ihigh(:,:,2:end) = [];
        ihigh = ihigh.';
        ihigh = max(ihigh(:)) - ihigh;
        ihigh = double(ihigh) / double(max(ihigh(:)));
        
        i = Downsize(ihigh);        
    end

    function i = LoadColorImage(name)
        ihigh = imread(name);
        ihigh = permute(ihigh, [2 1 3]);            % Swap x and y in each color plane
        ihigh = flipdim(ihigh,2);
        ihigh = double(ihigh) / double(max(ihigh(:)));

        for k = 1:3
            i(:,:,k) = Downsize(ihigh(:,:,k));
        end        
    end

    % New (slow but nice) downsampling using RenderQuadMesh.
    function i = Downsize(ihigh)
        [m,n] = size(ihigh);
        x1gv = linspace(0.5,double(res)+0.5,m+1);
        x2gv = linspace(0.5,double(res)+0.5,n+1);
        [x1,x2] = ndgrid(x1gv,x2gv);
        i = zeros(res,res);
        i = RenderQuadMeshInternal(i, x1, x2, ihigh, 0);
    end

    % Old downsampling; not very good.
    function i = OldDownsize(ihigh)
        m = size(ihigh, 1);
        n = size(ihigh, 2);
        i = interp2(ihigh, linspace(1,n,res), linspace(1,m,res).', 'cubic');
        i = max(min(i,1),0);        
    end
end