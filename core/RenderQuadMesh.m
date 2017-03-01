% Render an (anti-aliased) quadrilateral and add it to the image
%
%  RenderQuadMesh(O, x1, x2, v, p)
%  RenderQuadMesh(O, x1, x2, v, p, w)
%
%      O: image (m x n array)
% x1, x2: (k+1) by (l+1) matrices with the coordinates of the mesh points.
%      v: k by l matrix with the values to assign to each quadrilateral
%      p: Flag controlling wraparound:
%          p = 0: no wraparound (nonperiodic grid)
%          p = 1: wraparound (periodic grid)
%          p = 2: wraparound + unwrap coordinates of each quad
%      w: optional m x n weight function.
%          If w is given, the rendered mesh is multiplied by w before
%          being added to the image O.
%
% O, v, and w may be complex.
%
% The real work is done by the MEX function RenderQuadMeshInternal.

function RenderQuadMesh(O, x1, x2, v, p, w)

x1 = O.grid.toGridSpace(x1,1);
x2 = O.grid.toGridSpace(x2,2);
if (nargin < 6)
    O.f = RenderQuadMeshInternal(O.f, x1, x2, v, p);
else
    O.f = RenderQuadMeshInternal(O.f, x1, x2, v, p, w);
end

end