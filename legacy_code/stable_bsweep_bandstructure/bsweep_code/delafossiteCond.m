%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% delafossiteCond.m
%
% Aaron Sharpe
% 1/29/18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0 = 100;
theta = 0;

% Bandstructure inputs
phi = 0*pi/6;
n_theta = 250; % Number of points for generating Fermi surface
n_interp = 250; % Number of points to interpolate Fermi surface


%% Generate closed caustic loop

% For positive fields
delp = delafossiteBandstructure;
delp.phi = phi;
delp.B = 1;
delp.n_theta = n_theta;
delp.n_interp = n_interp;
delp.init();
% Plot orbits to check shape
delp.KtoR(1)
delp.plotOrbits(2)

[Gp, thetaP] = delp.cond(w0, theta);
Gp



deln = delafossiteBandstructure;
deln.phi = phi;
deln.B = -1;
deln.n_theta = n_theta;
deln.n_interp = n_interp;
deln.init();
% Plot orbits to check shape
deln.KtoR(-1)

[Gn, thetaN] = deln.cond(w0, theta);
Gn



