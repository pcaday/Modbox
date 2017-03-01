% SurfaceAreaOfNSphere
%
% Return the surface area of the unit n-sphere.
%
function SA = SurfaceAreaOfNSphere(n)
    SA = (n+1) * pi^((n+1)/2) / (gamma((n+1)/2+1));
end