% Test the performance of the function monodisperseSolver
clear  % Clear all variables from the workspace
close all  % Close all open figure windows
movegui('west')  % Move the figure window to the left side of the screen

%% Define problem parameters
Re = 1;          % Reynolds number
Ri = 1;          % Richardson number
rhop = [2.5, 3.8]; % Relative density of two particle species
alpha = 1;       % Spiral pitch parameter (2*pi*alpha is the pitch)
R = 1;           % Channel width of the spiral separator
r = 1;           % Radial position
Kv = 0.62;       % Shear-induced migration parameter
Kc = 0.41;       % Shear-induced migration parameter
phim = 0.61;     % Maximum packing fraction
hr = 1;          % Domain height

%% Compute critical volume fractions for both particle species
temp = 2 * r^2 * R / (9 * alpha * Kc) + 1 / (rhop(1) - 1);
phic1 = min(phim, 0.5 * (sqrt(temp^2 + (8 * r^2 * R) / (9 * alpha * Kc)) - temp));

temp = 2 * r^2 * R / (9 * alpha * Kc) + 1 / (rhop(2) - 1);
phic2 = min(phim, 0.5 * (sqrt(temp^2 + (8 * r^2 * R) / (9 * alpha * Kc)) - temp));

%% Define total particle volume fraction for the two species
phitotal = [0.05, 0.4];  % Initial volume fractions of two particle species

%% Solve the bidensity particle-laden flow problem using the shooting method
[z, sol, f1, x1, Niter] = bidensitySolver(Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, phitotal);

%% Display results
f1  % Display final function values (should be close to zero if solution converged)
x1  % Display final solution for the volume fraction and ratio parameter
Niter  % Display number of iterations taken by the Broyden method

%% Compute and display integral quantities
% Compute the integral of the first component weighted by exp(sol(:,2))
trapz(z, sol(:,1) .* exp(sol(:,2)))

% Compute the integral of the first component weighted by (1 - exp(sol(:,2)))
trapz(z, sol(:,1) .* (1 - exp(sol(:,2))))

%% Plot results
yyaxis left  % Set left y-axis
plot(z, sol(:,1) .* exp(sol(:,2)), '-b', 'linewidth', 1.5)  % Plot first species volume fraction
hold on  
plot(z, sol(:,1) .* (1 - exp(sol(:,2))), '-k', 'linewidth', 1.5)  % Plot second species volume fraction
plot([0, hr], phic1 * [1,1], '--b')  % Plot critical volume fraction for species 1
plot([0, hr], phic2 * [1,1], '--k')  % Plot critical volume fraction for species 2
hold off

yyaxis right  % Switch to right y-axis
plot(z, sol(:,3), '--r', 'linewidth', 1.5)  % Plot the stress function

%% Add legend
legend('$\phi$', '$\sigma$', 'Interpreter', 'latex', 'location', 'northwest')
