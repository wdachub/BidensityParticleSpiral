function [z,sol,f1,nIter] = monodisperseSolver(absTol,Re,Ri,rho,alpha,R,r,Kv,Kc,phim,hr,phitotal)
% monodisperseSolver: Computes the equilibrium particle volume fraction 
% in a spiral separator using the shooting method.
%
% This function solves a nonlinear boundary value problem using a shooting 
% method combined with the secant method to determine the correct initial 
% condition. The ODE solver terminates when the particle volume fraction 
% reaches zero or exceeds the maximum packing fraction.
%
% Inputs:
%   absTol   - Absolute tolerance for the secant method convergence.
%   Re       - Reynolds number.
%   Ri       - Richardson number.
%   rho      - Relative density of the particles.
%   alpha    - Pitch parameter (2*pi*alpha is the pitch of the spiral separator).
%   R        - Channel width of the spiral separator.
%   r        - Radial variable.
%   Kv, Kc   - Parameters governing shear-induced migration.
%   phim     - Maximum packing fraction of particles.
%   hr       - Height of the domain.
%   phitotal - Total particle volume within the domain.
%
% Outputs:
%   z        - Grid points along the height.
%   sol      - Particle volume fraction at each grid point.
%   f1       - Error in satisfying the boundary condition.
%   nIter    - Number of iterations used in the secant method.

% Set up ODE solver termination conditions: stop if particle volume fraction 
% reaches zero or exceeds the maximum packing fraction.
options = odeset('Events', @(t,y) EventFunction(t,y,phim));

% Compute the critical volume fraction phic based on physical parameters.
temp = 2*r*R/(9*alpha*Kc) + 1/(rho - 1);
phic = min(phim, 0.5 * (sqrt(temp^2 + (8*r*R) / (9*alpha*Kc)) - temp));

% Initialize iteration counter
nIter = 0;

% Define the spatial domain [0, hr]
zspan = [0 hr];

% Initial guess for the shooting method
x0 = phitotal / hr;

% Compute the error for the first initial guess
f0 = shooting(Re, Ri, rho, alpha, R, r, Kv, Kc, phim, zspan, x0, phitotal, options);

% Generate a second initial guess based on phic
if phic > (phitotal / hr)
    x1 = phitotal / hr * 0.99; % Slightly reduce the guess if below critical
else
    x1 = min(phitotal / hr * 1.01, phim * 0.99); % Slightly increase the guess
end

% Compute the error for the second initial guess
f1 = shooting(Re, Ri, rho, alpha, R, r, Kv, Kc, phim, zspan, x1, phitotal, options);

% Secant method loop: adjust x1 until f1 is within the tolerance
while abs(f1) > absTol
    % Compute the next approximation using the secant method
    xnew = (x0 * f1 - x1 * f0) / (f1 - f0);
    
    % Update variables for the next iteration
    x0 = x1;
    f0 = f1;
    x1 = xnew;
    
    % Compute the new function value and solution profile
    [f1, z, sol] = shooting(Re, Ri, rho, alpha, R, r, Kv, Kc, phim, zspan, x1, phitotal, options);
    
    % Increment iteration counter
    nIter = nIter + 1;
end
end


function [f1, z, sol] = shooting(Re, Ri, rho, alpha, R, r, Kv, Kc, phim, zspan, x0, phitotal, options)
% shooting: Auxiliary function for the shooting method.
%
% This function integrates the governing equation using an initial guess 
% and computes the deviation from the desired boundary condition.
%
% Inputs:
%   Re, Ri, rho, alpha, R, r, Kv, Kc, phim - Physical and model parameters.
%   zspan     - Interval [0, hr] defining the domain.
%   x0        - Initial guess for the volume fraction at z = 0.
%   phitotal  - Total particle volume fraction in the system.
%   options   - Options for the ODE solver, including termination criteria.
%
% Outputs:
%   f1        - Residual error at the boundary.
%   z         - Grid points along the height.
%   sol       - Particle volume fraction at each grid point.

% Extract height limit
hr = zspan(2);

% Define the initial condition for the ODE solver
% Initial condition includes volume fraction and an estimated derivative
initialValue = [x0, -alpha * Re * Ri / (r * R) * (hr + (rho - 1) * phitotal)];

% Solve the ODE system using ode45
[z, sol] = ode45(@(z, phiSigma) particleEquation(z, phiSigma, Re, Ri, rho, alpha, R, r, Kv, Kc, phim), ...
                 zspan, initialValue, options);

% If the solver terminates before reaching the right boundary (hr), 
% extrapolate the final point using the analytical formula of the exact solution of the governing equation.
if z(end) < hr
    sol = [sol;
           sol(end,1), sol(end,2) + alpha * Re * Ri / (r * R) * (1 + (rho - 1) * sol(end,1)) * (hr - z(end))];
    z = [z; hr];
end

% Extract the residual error at the boundary
f1 = sol(end,2);
end
