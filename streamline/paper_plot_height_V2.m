% This script generates Figure 5 from the paper, illustrating streamlines, free surface profiles, and particle volume fractions.

clear
close all

% Create figure window with specific dimensions (3.5'' width, 2'' height) for publication-quality plots
figs = figure('Units','inch','Position',[0 0 3.5 3.5*0.618]);
movegui(figs,'west');

% Set output format for saving figures
outputformat = '.pdf';

% Define axis position parameters for consistent figure layout
leftMargin = 0.085;
BottomMargin = 0.17;
width = 0.75;
hight = 0.8;
axes('Position',[leftMargin, BottomMargin, width, hight]);

% Select test case (1: first column panels (a,c,e), 2: second column panels (b,d,f))
testcase = 2;

switch testcase
    case 1
        % Parameters for the first column of Figure 5 (Panels a, c, e)
        Re = 102.094; % Reynolds number
        Ri = 0.084532; % Richardson number
        alpha = 0.0382; % Channel inclination angle parameter
        R = 0.2; % Dimensionless radius ratio
        Kv = 0.62; % Particle sedimentation coefficient
        Kc = 0.41; % Diffusion coefficient ratio
        phim = 0.61; % Maximum packing fraction of particles
        rhop2 = 3.8; % Density ratio of heavy particles
        rhop1 = 2.5; % Density ratio of light particles
        rhol = 1; % Density of liquid phase
        Rp = 0.213; % Inner radius of particle-laden region
        Rf = 0.502; % Radius where free fluid region begins
        Rl = 0.05; % Initial radius of the flow domain
        height0 = [1]; % Initial height of the free surface
    case 2
        % Parameters for the second column of Figure 5 (Panels b, d, f)
        Re = 233;
        Ri = 0.025;
        alpha = 1;
        R = 0.5;
        Kv = 0.62;
        Kc = 0.41;
        phim = 0.61;
        rhop2 = 3.8;
        rhop1 = 2.5;
        rhol = 1;
        Rp = 0.0965;
        Rf = 0.458;
        Rl = 2;
        height0 = [0.87];
end

%% Compute critical volume fraction for particle-laden flow region
critical_volume(Rl, rhop2, alpha, R, Kc, phim);

%% Compute free surface profile for different regions

% Compute the free surface for the heavy-particle-liquid mixture
rspan = Rl + [0, Rp];
[rp2, heightp2] = ode45(@(r,h) heightEquation(r,h,Re,Ri,rhop2,alpha,R,Kv,Kc,phim), rspan, height0);
height0 = heightp2(end);

% Compute the free surface for the light-particle-liquid mixture
rspan = Rl + [Rp, Rf];
[rp1, heightp1] = ode45(@(r,h) heightEquation(r,h,Re,Ri,rhop1,alpha,R,Kv,Kc,phim), rspan, height0);
height0 = heightp1(end);

% Compute the free surface for the pure fluid region
rspan = Rl + [Rf, 1];
[rf, heightf] = ode45(@(r,h) heightEquation(r,h,Re,Ri,rhol,alpha,R,Kv,Kc,phim), rspan, height0);

% Store computed free surface profiles
height = [heightp2; heightp1; heightf];

%% Compute critical particle volume fraction
critical_phi = @(r, rho) min(phim, 0.5 * (sqrt((2 * r.^2 * R / (9 * alpha * Kc) + 1 / (rho - 1)).^2 + (8 * r.^2 * R) / (9 * alpha * Kc)) - (2 * r.^2 * R / (9 * alpha * Kc) + 1 / (rho - 1))));

% Define radius values for plotting
r = [rp2; rp1; rf];
cphi{1} = critical_phi(rp2, rhop2);
cphi{2} = critical_phi(rp1, rhop1);
cphi{3} = 0;

%% Plot streamlines and free surface

% Plot demarcation lines for different regions
plot((Rl+Rp)*[1,1], [0, max(height)], ':k', 'linewidth', 1);
hold on;
plot((Rl+Rf)*[1,1], [0, max(height)], ':k', 'linewidth', 1);

% Compute and plot stream function contours for each region
% Heavy-particle-liquid mixture
zn = linspace(0,1,51)';
zz = zn * heightp2';
rr = repmat(rp2', size(zn));
hh = repmat(heightp2', size(zn));
psi = streamfunction(zz, rr, hh, Re, Ri, rhop2, alpha, R, Kv, Kc, phim);
contour(rr, zz, psi, linspace(0, max(max(abs(psi))) * 0.8, 5), 'LineWidth', 1);

% Light-particle-liquid mixture
zz = zn * heightp1';
rr = repmat(rp1', size(zn));
hh = repmat(heightp1', size(zn));
psi = streamfunction(zz, rr, hh, Re, Ri, rhop1, alpha, R, Kv, Kc, phim);
contour(rr, zz, psi, linspace(0, max(max(abs(psi))) * 0.8, 5), 'LineWidth', 1);

% Pure fluid region
zz = zn * heightf';
rr = repmat(rf', size(zn));
hh = repmat(heightf', size(zn));
psi = streamfunction(zz, rr, hh, Re, Ri, rhol, alpha, R, Kv, Kc, phim);
contour(rr, zz, psi, 'LineWidth', 1);

% Add colorbar and adjust properties
chb = colorbar;
set(gca, 'ColorScale', 'log');
set(chb, 'position', [0.87, 0.08, 0.05, 0.9]);

% Plot free surface
plot(r, height, '-k', 'linewidth', 1.2);
ylabel('$z$', 'Interpreter', 'latex', 'fontsize', 10, 'rotation', 0);
xlabel('$r$', 'Interpreter', 'latex', 'fontsize', 10);
axis tight;

% Save figure
exportgraphics(figs, ['freeSurface', num2str(testcase), outputformat], 'ContentType', 'Vector');

%% Plot particle volume fraction profiles
figs = figure('Units', 'inch', 'Position', [0 0 3.5 3.5 * 0.618]);
movegui(figs, 'west');
axes('Position', [0.11, 0.17, 0.85, 0.8]);

plot(r, [cphi{1}; 0 * rp1; 0 * rf], '--b', 'linewidth', 1);
hold on;
plot(r, [0 * rp2; cphi{2}; 0 * rf], '-r', 'linewidth', 1);
hold off;
xlabel('$r$', 'Interpreter', 'latex', 'fontsize', 12);
legend('Heavy particle', 'Light particle', 'location', 'northeast');
axis tight;

exportgraphics(figs, ['volumefraction', num2str(testcase), outputformat], 'ContentType', 'Vector');


%% produce figure (e) or (f), the schematic of spiral geometries with respect to dimensionless variables.
figs = figure('Units','inch','Position',[0 0 3.5 3.5*0.618]); % 3.5'' width , 2'' length
movegui(figs,'northwest');

% Define parameters for the spiral
theta = linspace(0, 2*pi, 100); % Angle parameter
r = linspace(Rl, Rl+1, 10); % Height parameter
[tt,rr]=meshgrid(theta,r);
% Parametric equations for a spiral
x = rr .* cos(tt);
y = rr .* sin(tt);
z=alpha*tt;
% Create a 3D surface plot

surf(x, y, z);
axis tight
view([-10,75])
% Customize the plot
% title('Spiral Surface');
xlabel('$x$','Interpreter','latex','fontsize',12)
ylabel('$y$','Interpreter','latex','fontsize',12)
zlabel('$z$','Interpreter','latex','fontsize',12)

exportgraphics(figs,['spiral',num2str(testcase),outputformat],'resolution',300)
