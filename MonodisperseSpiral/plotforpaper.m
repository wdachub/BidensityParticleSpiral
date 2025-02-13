% This script produces Figure 3(a) and 3(b) in the paper.

clear  % Clear workspace
close all  % Close all figure windows

% Create figure with specific dimensions
figs = figure('Units','inch','Position',[0 0 3.5 3.5*0.618]); % 3.5'' width
movegui(figs,'west'); % Move figure to the left side of the screen

% Define margins and axis positioning
leftMargin = 0.12;
BottomMargin = 0.15;
width = 0.85;
hight = 0.82;
wDistant = 0.49;
hDistant = 0.49;

% Create axes for the plot
ax = axes('Position',[leftMargin, BottomMargin, width, hight]);

%% Define physical parameters
Re = 1;    % Reynolds number
Ri = 10;   % Richardson number
rho = 2.5; % Density ratio
alpha = 1; % Diffusivity parameter
R = 1;     % Outer radius
r = 1;     % Inner radius
Kv = 0.62; % Vertical permeability
Kc = 0.41; % Horizontal permeability
phim = 0.61; % Maximum volume fraction
hr = 1;      % Height ratio

%% Set absolute tolerance for numerical solver
absTol = 1e-4;

%% Generate Figure 3(a)
% Compute the critical volume fraction phic
temp = 2 * r * R / (9 * alpha * Kc) + 1 / (rho - 1);
phic = min(phim, 0.5 * (sqrt(temp^2 + (8 * r * R) / (9 * alpha * Kc)) - temp));

% Define offset values for different cases in the plot
offset = 0.03:0.03:0.06;
linesytle = {'-r', '-b'}; % Line styles for each offset value

% Loop over different offset values and compute solutions
for indi = 1:length(offset)
    phitotal = offset(indi);
    
    % Solve for concentration profile using numerical solver
    [z, sol, f1, nIter] = monodisperseSolver(absTol, Re, Ri, rho, alpha, R, r, Kv, Kc, phim, hr, phitotal);
    
    % Plot solution
    plot(z, sol(:,1), linesytle{indi}, 'linewidth', 1.5);
    hold on;

    % Store legend text
    legendtext{indi} = ['$\bar{\phi}=$', num2str(offset(indi))];
end

% Adjust axis limits
axis tight;

% Create legend
hleg = legend(legendtext, 'Interpreter', 'latex', 'location', 'southwest', 'FontSize', 10);
hleg.NumColumns = 2;
hleg.ItemTokenSize = [14,9];

% Plot critical volume fraction line
plot([0,1], phic * [1,1], '--k', 'linewidth', 1.5);

% Add legend entry for critical volume fraction
legendtext{end+1} = ['$\bar{\phi}=\phi_c$'];

% Update legend location
hleg = legend(legendtext, 'Interpreter', 'latex', 'location', 'northeast');
hleg.NumColumns = 3;

% Label axes
xlabel('$z$', 'Interpreter', 'latex', 'fontsize', 12, 'verticalalignment', 'middle');
ylabel('$\phi$', 'Interpreter', 'latex', 'rotation', 0, 'fontsize', 11, 'horizontalalignment', 'right');

% Set y-axis limits
ylim([0,0.68]);

% Export figure as a vector graphic
exportgraphics(figs, 'criticalVolumeFractionMonodisperse1.pdf', 'ContentType', 'Vector');

% Clear current figure for the next plot
clf;

%% Generate Figure 3(b)
% Recompute critical volume fraction phic
temp = 2 * r * R / (9 * alpha * Kc) + 1 / (rho - 1);
phic = min(phim, 0.5 * (sqrt(temp^2 + (8 * r * R) / (9 * alpha * Kc)) - temp));

% Define offset values relative to phic
offset = 1 + (-0.05:0.05:0.05);
linesytle = {'-r', '--k', '-b'}; % Line styles for different cases

% Loop over offset values and compute solutions
for indi = 1:length(offset)
    phitotal = phic * offset(indi);
    
    % Solve for concentration profile
    [z, sol, f1, nIter] = monodisperseSolver(absTol, Re, Ri, rho, alpha, R, r, Kv, Kc, phim, hr, phitotal);
    
    % Plot solution
    plot(z, sol(:,1), linesytle{indi}, 'linewidth', 1.5);
    hold on;

    % Store legend text
    if abs(offset(indi) - 1) < eps
        legendtext{indi} = ['$\phi_c$'];
    else
        legendtext{indi} = [num2str(offset(indi)), '$\phi_c$'];
    end
end

% Adjust axis limits
axis tight;

% Create legend
hleg = legend(legendtext, 'Interpreter', 'latex', 'location', 'southwest', 'FontSize', 10);
hleg.ItemTokenSize = [14,9];

% Label axes
xlabel('$z$', 'Interpreter', 'latex', 'fontsize', 12, 'verticalalignment', 'middle');
ylabel('$\phi$', 'Interpreter', 'latex', 'rotation', 0, 'fontsize', 11, 'horizontalalignment', 'right');

% Export figure as a vector graphic
exportgraphics(figs, 'criticalVolumeFractionMonodisperse2.pdf', 'ContentType', 'Vector');
