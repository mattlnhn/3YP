clc; clear

n = 50; % no. of nodes = no. of cells per side
L = 1e-4; % length of side, m
theta = 100e-6; % thickness, m

mat.A = [26.982, 51.996, 15.999]; % atomic mass, g mol-1
mat.Z = [13, 24, 8]; % atomic no.
mat.I = [166e-6, 315.8e-6, 95e-6]; % mean ionisation energy, MeV
mat.w = [0.52533, 0.00509, 0.46959]; % mass fraction
mat.rho = 3.85; % density, g cm-3
mat.epsilon = .8; % emissivity coefficient
mat.k = 30; % thermal conductivity, W m-1 K-1
mat.c_p = 900; % specific heat capacity J kg-1 K-1

beam.M = 938.27208816; % particle mass, MeV c-2
beam.gamma = 7460.69; % Lorentz factor
beam.f = 11.2e3;  % bunches per second, Hz
beam.np = 1.15e11; % protons per bunch (nominal 25ns spacing)
beam.sigma = 2.5e-6; % beam width, m
beam.pos = [-5*beam.sigma 0]; % beam centre pos w.r.t. midpoint of left edge, m

dl = L/n;
dt = .25*(mat.rho*1000*mat.c_p*dl*dl)/mat.k; % time increment, s
nt = 100000; % no. of time steps
T = dt*nt;

fprintf('dt = %.10f with %d time steps for total sim time %d seconds', dt, nt, T)

T = heating2(n, L, theta, dt, nt, mat, beam);

figure
contourf(T(:, :, end)-T(:, :, 1))
colorbar

