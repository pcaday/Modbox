% Compute a pushforward.
%
%   Pushforward(O, f, images, wrap)
%
% Input/Output
%    O: A Function object. The pushforward of f (multiplied by a Jacobian
%        factor) will be added to O.
%
%       The Jacobian factor is such that calling Pushforward with an
%        invertible map should be equivalent to pulling back by the map's
%        inverse.
%  
% Inputs
%      f: A Function object representing the function to be pushed forward.
% images: An n x i1 x ... x in array, the images of the grid points of the
%          the grid of f in the grid of O.
%   wrap: A flag specifying whether the grid of O should be treated as
%          periodic, and if the coordinates images in 'images' should be
%          unwrapped (see RenderQuadMesh).
%         Defaults to 0 (no wrapping) if O.grid is not periodic in any
%          direction, and 2 (periodic + unwrapping)
function Pushforward(O, f, gridImages, wrapFlag)

if nargin < 4, wrapFlag = 2 * any(O.grid.periodic); end

switch O.grid.dims
    case 2
        RenderQuadMesh(O, squeeze(gridImages(1,:,:)), squeeze(gridImages(2,:,:)), f.f, wrapFlag);
    otherwise
        error('Pushforward only implemented for 2D currently.');
end

end