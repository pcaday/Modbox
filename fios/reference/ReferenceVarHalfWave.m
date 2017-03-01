% Reference solver for wave propagation using the half-wave equation
%
%  du/dt + i Pu = 0
%
% where P is the pseudodifferential operator with symbol
%  |c(x)||xi|.
%
% Inputs:
%     f: input function (Function object).
%     t: total time to propagate the waves
%     c: function handle for the wave speed (c) function, or Function
%         object on f's grid, or a EuclideanConformalMetric, 
%         or a plain array of c values.
function Af = ReferenceVarHalfWave(f,t,c)
    global GlobalRefVarHalfWaveDt
    global GlobalRefVarHalfWaveVerbosity
    
    
    % Convert the wave speed to a raw array
    if isa(c, 'Function')
        c = c.f;
    elseif isa(c, 'function_handle')
        c = c(f.grid.AllPoints());
    elseif isa(c, 'EuclideanConformalMetric')
        c = exp(-c.lambda(f.grid.AllPoints()));
    end
    
    % Get the dt to use.
    dt = GlobalRefVarHalfWaveDt;
    
    % Round dt down to the nearest fraction of the total time t
    nt = ceil(t / dt);
    dt = t / nt;
    
    % Get an array of |xi| values.
    xi = f.grid.DFTGrid.AllPoints();
    xiR = sqrt(fsqueeze(sum(xi.^2, 1)));
    
    % Use ifftshift to move from the DFT grid's coordinate system
    %  (which has the zero frequency in the center)
    %  to the frequency coordinate system arrangement of ifftn/fftn.
    xiR = ifftshift(xiR);
    
    % March forward in time.
    u = f.f;
    for ts = 1:nt
        % Forward Euler:
        %  u = u + dt * ComputeDuDt(u);
        
        % RK4 time stepping:
        q1 = dt*ComputeDuDt(u);
        q2 = dt*ComputeDuDt(u + q1/2);
        q3 = dt*ComputeDuDt(u + q2/2);
        q4 = dt*ComputeDuDt(u + q3);
        
        u = u + (q1 + 2*q2 + 2*q3 + q4) / 6;
        
        if GlobalRefVarHalfWaveVerbosity > 1,
            subplot(1,2,1,'replace')
            compleximagesc(q1.'); axis xy; colorbar; title('Pu');
            subplot(1,2,2,'replace')
            compleximagesc(u.'); axis xy; colorbar;
            uf = Function.WithValues(f.grid,u);
            title(sprintf('u, time = %.2f, norm = %g', ts*dt, uf.norm));
            drawnow;
        end
    end
    
    Af = Function.WithValues(f.grid, u);
    
    
    
    
    function dudt = ComputeDuDt(u)
        % Apply the pseudodifferential operator P
        dudt = -1i * ifftn(xiR .* fftn(u .* c));
    end
end