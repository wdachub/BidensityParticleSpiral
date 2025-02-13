function [z,sol,f1,x1,Niter] = bidensitySolver(Re,Ri,rhop,alpha,R,r,Kv,Kc,phim, hr, phitotal, absTol)
% BIDENSITYSOLVER Solves for the bidisperse particle-laden flow using a shooting method.
%
% This function determines the steady-state solution of a bidisperse suspension
% flow problem by solving a boundary value problem using a shooting method.
% Since two unknown parameters need to be determined, the Broyden method is used
% instead of the secant method for faster convergence.
%
% INPUTS:
%   Re       - Reynolds number
%   Ri       - Richardson number
%   rhop     - Density ratios of the two particle species relative to the fluid [rho1, rho2]
%   alpha    - Spiral pitch parameter (2*pi*alpha is the pitch)
%   R        - Characteristic length scale (e.g., channel width)
%   r        - Radial position in the channel
%   Kv       - Shear-induced migration parameter
%   Kc       - Shear-induced migration parameter
%   phim     - Maximum packing fraction
%   hr       - Domain length in the z-direction
%   phitotal - Total initial volume fractions of the two species [phi1, phi2]
%   Tol      - Tolerance for convergence of the shooting method
%
% OUTPUTS:
%   z     - Computed z-coordinates
%   sol   - Solution matrix containing [phi, logchi, sigma, integral of phi]
%   f1    - Residual function values at the final iteration
%   x1    - Final values of the shooting parameters
%   Niter - Number of iterations required for convergence

%% defulat value
if ~exist('absTol','var')
    absTol=1e-3;
end

%% Compute critical volume fractions for each species
temp = 2*r*R/(9*alpha*Kc) + 1/(rhop(1)-1);
phic1 = min(phim, 0.5 * (sqrt(temp^2 + (8*r*R)/(9*alpha*Kc)) - temp));
temp = 2*r*R/(9*alpha*Kc) + 1/(rhop(2)-1);
phic2 = min(phim, 0.5 * (sqrt(temp^2 + (8*r*R)/(9*alpha*Kc)) - temp));

%% Set up ODE solver options
ODEoptions = odeset('Events', @(t,y) EventFunction(t,y,phim), 'RelTol', 1e-5, 'AbsTol', 1e-7);
zspan = [0 hr];

% Initial guess for sigma based on force balance
sigma0 = -alpha * Re * Ri / (r * R) * (hr + (rhop(1)-1)*phitotal(1) + (rhop(2)-1)*phitotal(2));

% Convert phitotal into a total volume fraction and ratio of species concentrations
ratio = log(phitotal(1) / (phitotal(1) + phitotal(2)));
phitotal = sum(phitotal);

% Set initial guesses for the shooting method
if phitotal < max(phic2, phic1)
    x0 = [phitotal/hr; ratio];
    x1 = [phitotal/hr * 0.95; ratio * 1.1];
else
    x0 = [phitotal/hr; ratio];
    x1 = [min(phitotal/hr * 1.05, phim); ratio * 1.1];
end

%% Perform the first shooting method evaluations
initialValue = [x0; sigma0; 0];
[z, sol] = shotting(zspan, initialValue, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, ODEoptions);
f0 = [sol(end,3); sol(end,4) - phitotal]; % Residuals from first guess

initialValue = [x1; sigma0; 0];
[z, sol] = shotting(zspan, initialValue, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, ODEoptions);
f1 = [sol(end,3); sol(end,4) - phitotal]; % Residuals from second guess

%% Iterate using the Broyden method to adjust shooting parameters
Niter = 0;
invJ0 = eye(2); % Initial inverse Jacobian estimate

while max(abs(f1)) > absTol  % Convergence check
    deltax = x1 - x0;
    deltaf = f1 - f0;

    % Update inverse Jacobian using Broyden's method
    invJ1 = invJ0 + ((deltax - invJ0 * deltaf) / (deltax' * invJ0 * deltaf)) * (deltax' * invJ0);
    x2 = x1 - invJ1 * f1; % New estimate of shooting parameters

    % Ensure the new solution remains within the valid domain
    while (x2(1) > 1 || x2(1) < 0 || x2(2) > 0 || norm(deltax,2) > 0.2)
        deltax = deltax / 2;
        x1 = x0 + deltax;
        initialValue = [x1; sigma0; 0];
        [z, sol] = shotting(zspan, initialValue, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, ODEoptions);
        f1 = [sol(end,3); sol(end,4) - phitotal];
        deltaf = f1 - f0;
        invJ1 = invJ0 + ((deltax - invJ0 * deltaf) / (deltax' * invJ0 * deltaf)) * (deltax' * invJ0);
        x2 = x1 - invJ1 * f1;
    end

    % Update values for next iteration
    x0 = x1;
    f0 = f1;
    x1 = x2;
    invJ0 = invJ1;

    % Compute new residuals
    initialValue = [x1; sigma0; 0];
    [z, sol] = shotting(zspan, initialValue, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, ODEoptions);
    f1 = [sol(end,3); sol(end,4) - phitotal];

    Niter = Niter + 1;
end

end

%% Helper function: Shooting method
function [z,sol] = shotting(zspan, initialValue, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, ODEoptions)
% SHOTTING Solves the ODE system using the shooting method for a given initial guess.
%
% This function integrates the ODE system for a given initial guess of
% the shooting parameters and returns the solution.
%
% INPUTS:
%   zspan        - Interval for integration
%   initialValue - Initial conditions including guessed parameters
%   Re, Ri       - Reynolds and Richardson numbers
%   rhop         - Density ratios of the two species
%   alpha, R, r  - Geometric parameters
%   Kv, Kc       - Shear-induced migration parameters
%   phim         - Maximum packing fraction
%   hr           - Height of the domain
%   ODEoptions   - Options for the ODE solver
%
% OUTPUTS:
%   z   - Computed z-coordinates
%   sol - Solution matrix containing [phi, logchi, sigma, integral of phi]

[z, sol] = ode15s(@(z, phiSigma) particleEquation(z, phiSigma, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim), zspan, initialValue, ODEoptions);

% If the solver terminates before reaching hr, extrapolate the final values
if z(end) < hr
    sol = [sol;
        sol(end,1), sol(end,2), ...
        sol(end,3) + alpha * Re * Ri / (r * R) * (1 + (rhop(1) - 1) * sol(end,1) + rhop(2) * sol(end,2)) * (hr - z(end)), ...
        sol(end,4)];
    z = [z; hr];
end

end
