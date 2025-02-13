% Test the performance of the function monodisperseSolver

clear
close all


%% Define Physical and Model Parameters
Re = 1;         % Reynolds number
Ri = 10;        % Richardson number
rho = 2.5;      % Relative particle density
alpha = 1;      % Spiral pitch parameter (2*pi*alpha is the pitch)
R = 1;          % Channel width of the spiral separator
r = 5;          % Radial position
Kv = 0.62;      % Shear-induced migration parameter
Kc = 0.41;      % Shear-induced migration parameter
phim = 0.68;    % Maximum packing fraction
hr = 1;         % Height of the domain

%% Numerical Tolerance for the Solver
absTol = 1e-4;  % Absolute tolerance for secant method convergence

%% Compute Critical Volume Fraction
temp = 2*r*R/(9*alpha*Kc) + 1/(rho - 1);
phic = min(phim, 0.5 * (sqrt(temp^2 + (8*r*R) / (9*alpha*Kc)) - temp));

% Set total particle volume fraction
phitotal = 0.2;

%% Solve for Equilibrium Particle Distribution
[z, sol, f1, nIter] = monodisperseSolver(absTol, Re, Ri, rho, alpha, R, r, Kv, Kc, phim, hr, phitotal);

% Display iteration count and residual error
disp('Number of iterations:'); disp(nIter)
disp('Residual error:'); disp(f1)

%% Plot Results
figure;
movegui('west')
yyaxis left
plot(z, sol(:,1), '-b', 'LineWidth', 1.5) % Plot particle volume fraction
hold on 
plot([0, hr], phic * [1, 1], '--k')       % Plot critical volume fraction
hold off
ylabel('Particle Volume Fraction')

yyaxis right
plot(z, sol(:,2), '--r', 'LineWidth', 1.5) % Plot auxiliary variable (e.g., flux)
ylabel('Auxiliary Variable')

xlabel('Height (z)')
legend({'$\phi$', '$\phi_{max}$', '$\sigma$'}, 'Interpreter', 'latex', 'Location', 'northwest')
grid on

%% Compute and Display Additional Metrics
particleVolume = trapz(z, sol(:,1));  % Compute total particle volume using trapezoidal integration
finalFlux = sol(end,2);               % Extract final flux value

disp('Total integrated particle volume:'); disp(particleVolume)
disp('Final flux value:'); disp(finalFlux)
