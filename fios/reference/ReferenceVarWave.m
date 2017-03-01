% Reference solver for wave propagation using the variable-speed
%  wave equation
%
%  d^2u/dt^2 - c^2\Delta u = 0
%
% and Dirichlet boundary conditions.
%
% Inputs
%
% u0,v0: Cauchy data at t = 0:
%          u(x,t) = u0(x),   du/dt(x,t) = v0(x).
%     T: final time.
%    dt: time step.
% dx,dy: grid spacing in x and y.
%     c: sound speed array.
%  verb: Boolean, true to show progress.
%
% Outputs
%
%   u,v: Values of u and du/dt at final time T.
%
%
% Note: assumes that u0(i,j), v0(i,j) represent u0(x_i,y_j), etc.
%  (i.e. the first index on u0,v0,... is x and the second is y)
%
function [u,v] = ReferenceVarWave(u0,v0,T,dx,dy,dt,c,verbose)
    % Round dt down to the nearest fraction of the total time t
    nt = ceil(T / dt);
    dt = T / nt;
    
    % u and v hold current values of u and du/dt
    u = u0;
    v = v0;
    
    % Compute c^2
    c2 = c .* c;
    
    % March forward in time.
    for ts = 1:nt
        % Forward Euler:
        %  u = u + dt * ComputeDuDt(u);
        
        % RK4 time stepping:
        [uq1,vq1] = ComputeDuAndDv(u,         v,        dx,dy,dt);
        [uq2,vq2] = ComputeDuAndDv(u + uq1/2, v + vq1/2,dx,dy,dt);
        [uq3,vq3] = ComputeDuAndDv(u + uq2/2, v + vq2/2,dx,dy,dt);
        [uq4,vq4] = ComputeDuAndDv(u + uq3,   v + vq3,  dx,dy,dt);
        
        % Update u and v for next step.
        u = u + (uq1 + 2*uq2 + 2*uq3 + uq4) / 6;
        v = v + (vq1 + 2*vq2 + 2*vq3 + vq4) / 6;
        
        if verbose,
            subplot(1,2,1,'replace')
            imagesc(real(u).'); axis xy; colorbar; title('u');
            subplot(1,2,2,'replace')
            imagesc(real(v).'); axis xy; colorbar;
            title(sprintf('v, time = %.2f', ts*dt));
            drawnow;
        end
    end
    
    function [du,dv] = ComputeDuAndDv(u,v,dx,dy,dt)
        % du/dt is just v
        dudt = v;
        % dv/dt is c^2 times the Laplacian of u.
        dvdt = c2 .* LaplacianDirichlet(u,dx,dy);
        % Multiply by dt.
        du = dudt * dt;
        dv = dvdt * dt;
    end
end

% Compute Laplacian of a function.
%
% This is where the boundary conditions come in; we use Dirichlet here.
function Lu = LaplacianDirichlet(u,dx,dy)
    odx2 = 1 / (dx * dx);
    ody2 = 1 / (dy * dy);
    
    sz = size(u);
    zr = zeros(1,sz(2));
    zc = zeros(sz(1),1);
    
    Lu =  (-2*u + [u(2:end,:); zr] + [zr; u(1:end-1,:)]) * odx2 ...
        + (-2*u + [u(:,2:end)  zc] + [zc  u(:,1:end-1)]) * ody2;
end