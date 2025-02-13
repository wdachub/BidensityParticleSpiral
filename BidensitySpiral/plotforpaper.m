% This script produces Figure 4(a) and 4(b) in the paper.
% It solves the bidensity particle transport problem and plots the
% equilibrium volume fractions of two particle species.

clear
close all

% Initialize figure properties: 3.5'' width, aspect ratio 0.618 (Golden Ratio)
figs = figure('Units', 'inch', 'Position', [0 0 3.5 3.5 * 0.618]);
movegui(figs, 'west'); % Move figure to the left side of the screen

% Define axes position for plotting
leftMargin = 0.12;
BottomMargin = 0.15;
width = 0.85;
height = 0.82;

ax = axes('Position', [leftMargin, BottomMargin, width, height]);

%% Define Physical and Model Parameters
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

%% Compute Critical Volume Fractions for Two Densities
temp = 2 * r^2 * R / (9 * alpha * Kc) + 1 / (rhop(1) - 1);
phic1 = min(phim, 0.5 * (sqrt(temp^2 + (8 * r^2 * R) / (9 * alpha * Kc)) - temp));

temp = 2 * r^2 * R / (9 * alpha * Kc) + 1 / (rhop(2) - 1);
phic2 = min(phim, 0.5 * (sqrt(temp^2 + (8 * r^2 * R) / (9 * alpha * Kc)) - temp));

%% Solve for Equilibrium Particle Distribution - Figure 4(b)
phitotal = [0.05, 0.4]; % Total volume fractions of species 1 and 2

[z, sol, f1, x1, Niter] = bidensitySolver(Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, phitotal);

% Plot the equilibrium volume fractions for both species
plot(z, sol(:,1) .* exp(sol(:,2)), '-r', 'LineWidth', 1.5) % Species 1
hold on 
plot(z, sol(:,1) .* (1 - exp(sol(:,2))), '-b', 'LineWidth', 1.5) % Species 2
plot([0, hr], phic1 * [1, 1], '--r') % Critical volume fraction for species 1
plot([0, hr], phic2 * [1, 1], '--b') % Critical volume fraction for species 2
hold off

% Format plot labels and legend
xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'middle')
ylabel('$\phi$', 'Interpreter', 'latex', 'Rotation', 0, 'FontSize', 11, 'HorizontalAlignment', 'right')
axis tight
legend('$\phi_1$', '$\phi_2$', 'Interpreter', 'latex', 'Location', 'southwest')

% Save the figure as a vector graphic
exportgraphics(figs, 'criticalVolumeFractionBi1.pdf', 'ContentType', 'Vector')

%% Solve for Equilibrium Particle Distribution - Figure 4(1)
clf % Clear the figure for the next plot
ax = axes('Position', [leftMargin, BottomMargin, width, height]);

% New total volume fractions for species 1 and 2
phitotal = [0.30, 0.05];

[z, sol, f1, x1, Niter] = bidensitySolver(Re, Ri, rhop, alpha, R, r, Kv, Kc, phim, hr, phitotal);

% Plot the equilibrium volume fractions for both species
plot(z, sol(:,1) .* exp(sol(:,2)), '-r', 'LineWidth', 1.5) % Species 1
hold on 
plot(z, sol(:,1) .* (1 - exp(sol(:,2))), '-b', 'LineWidth', 1.5) % Species 2
plot([0, hr], phic1 * [1, 1], '--r') % Critical volume fraction for species 1
plot([0, hr], phic2 * [1, 1], '--b') % Critical volume fraction for species 2
hold off

% Format plot labels and legend
xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'middle')
ylabel('$\phi$', 'Interpreter', 'latex', 'Rotation', 0, 'FontSize', 11, 'HorizontalAlignment', 'right')
axis tight
ylim([0, 0.5]) % Adjust y-axis limits for better visualization
legend('$\phi_1$', '$\phi_2$', 'Interpreter', 'latex', 'Location', 'east')

% Save the figure as a vector graphic
exportgraphics(figs, 'criticalVolumeFractionBi2.pdf', 'ContentType', 'Vector')
