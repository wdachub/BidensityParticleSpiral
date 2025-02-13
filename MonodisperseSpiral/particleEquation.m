function dphiSigma = particleEquation(z, phiSigma, Re, Ri, rho, alpha, R, r, Kv, Kc, phim)
% PARTICLEEQUATION Defines the right-hand side of the nonlinear ODE system
% 
% This function computes the derivatives of the volume fraction (phi) and 
% the shear stress (sigma) for use in an ODE solver.
%
% INPUTS:
%   z        - Independent variable (typically height or axial position)
%   phiSigma - Vector containing [phi; sigma]
%   Re       - Reynolds number
%   Ri       - Richardson number
%   rho      - Particle density relative to fluid density
%   alpha    - Spiral pitch parameter (2*pi*alpha is the pitch)
%   R        - Characteristic length scale (e.g., channel width)
%   r        - Radial position in the channel
%   Kv       - Shear-induced migration parameter
%   Kc       - Shear-induced migration parameter
%   phim     - Maximum packing fraction
%
% OUTPUT:
%   dphiSigma - Column vector [dphi/dz; dsigma/dz], representing the 
%               rate of change of phi and sigma with respect to z

% Extract variables from input vector
phi = phiSigma(1, :);  % Volume fraction of particles
sigma = phiSigma(2, :); % Logarithmic concentration difference

% Compute the buoyancy-adjusted density difference
rhod = rho - 1;

% Compute the system of differential equations
dphiSigma = Re * Ri * [...
    % Equation for dphi/dz: balance of migration and buoyancy effects
    (-alpha * phi .* (1 + rhod * phi) / (r * R) + 2 * (1 - phi) * rhod / (9 * Kc)) ./ ...
    (1 + (Kv - Kc) / Kc * 2 * phi ./ (phim - phi)) ./ sigma;  
    
    % Equation for dsigma/dz: strain rate term
    alpha * (1 + rhod * phi) / (r * R)
];

end
