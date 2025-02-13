function dphiSigma = particleEquation(z, phiSigma, Re, Ri, rhop, alpha, R, r, Kv, Kc, phim)
% PARTICLEEQUATION Defines the right-hand side of the nonlinear ODE system
% 
% This function computes the derivatives of the volume fraction (phi), 
% the logarithm of the concentration ratio between two species (logchi), 
% and the shear stress (sigma) for use in an ODE solver.
%
% Since the concentration ratio between two species (chi) can vary by several 
% orders of magnitude, to ensure numerical stability, it is preferable to solve 
% for logchi instead of chi directly.
%
% INPUTS:
%   z        - Independent variable (typically height or axial position)
%   phiSigma - Vector containing [phi; logchi; sigma]
%   Re       - Reynolds number
%   Ri       - Richardson number
%   rhop     - Density of the two particle species relative to the fluid [rho1, rho2]
%   alpha    - Spiral pitch parameter (2*pi*alpha is the pitch)
%   R        - Characteristic length scale (e.g., channel width)
%   r        - Radial position in the channel
%   Kv       - Shear-induced migration parameter
%   Kc       - Shear-induced migration parameter
%   phim     - Maximum packing fraction
%
% OUTPUT:
%   dphiSigma - Column vector [dphi/dz; dlogchi/dz; dsigma/dz; phi], representing 
%               the rate of change of phi, logchi, sigma, and the integral of phi 
%               with respect to z

% Extract variables from input vector
phi = phiSigma(1, :);       % Volume fraction of particles
logchi = phiSigma(2, :);    % Logarithm of concentration ratio between two species
sigma = phiSigma(3, :);     % Shear stress
chi = exp(logchi);          % Convert logchi back to chi

% Compute the effective density variation
rhod = rhop - 1;  % Relative density difference
rho = 1 + phi .* (rhod(1) * chi + rhod(2) * (1 - chi)); % Effective density

% Shear-induced diffusion coefficient
Kt = 0.35;   
phia = 0.4;  % Cutoff volume fraction beyond which Dtr is constant
Dtr = (phi <= phia) .* (Kt * phi.^2) + (phi > phia) .* (Kt * phia^2);

% Compute the system of differential equations
dphiSigma = Re * Ri * [
    % Equation for dphi/dz: Balance of migration and buoyancy effects
    (-alpha * rho .* phi / (r * R) + 2 * (1 - phi) .* (rhod(1) * chi + rhod(2) * (1 - chi)) / (9 * Kc)) ./ ...
    (1 + (Kv - Kc) / Kc * 2 * phi ./ (phim - phi)) ./ sigma; 
    
    % Equation for dlogchi/dz: Evolution of species concentration ratio
    2 * (1 - chi) * (rhop(1) - rhop(2)) / (9 * sigma * Dtr * (1 - phi / phim));
    
    % Equation for dsigma/dz: Shear stress evolution
    alpha * rho / (r * R);
    
    % Integral of phi for diagnostics
    phi
];

end


% backup, the following code advancing chi directly, which suffers
% numerical problem
% phi=phiSigma(1,:);
% chi=phiSigma(2,:);
% sigma=phiSigma(3,:);
% 
% rhod=rhop-1;
% 
% rho=1+phi*(rhod(1)*chi+rhod(2)*(1-chi));
% 
% dphiSigma=Re*Ri*...
% [(-alpha*rho*phi/r^2/R+2*(1-phi)*(rhod(1)*chi+rhod(2)*(1-chi))/(9*Kc))./...
% (1+(Kv-Kc)/Kc*2*phi/(phim-phi))/sigma;%phi 
% 2*chi*(1-chi)*(rhop(1)-rhop(2))/(9*sigma*0.5*phi^2*(1-phi/phim));%chi
% alpha*rho/(r^2*R);%sigma
% phi;%the integral of phi
% ];




