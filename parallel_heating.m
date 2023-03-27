clear
close all
clc

parpool('Processes', 4)
L = 1e-2;
n = 800;

theta = 500e-6; % thickness, m

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
beam.nb = 2808; % no. of bunches
%beam.sigma = 2.5e-6; % beam width, m
beam.beta_star = 0.15;
beam.epsilon_n = 2.50e-6;
beam.sigma = sqrt(beam.beta_star*beam.epsilon_n/beam.gamma);
%beam.pos = -3*beam.sigma; % beam centre pos w.r.t. midpoint of left edge, m

sigma_range = [3.0, 3.1, 3.2, 3.3];
beampos = -sigma_range.*beam.sigma;

dl = L/(n*8);
tau = .25;
dt = tau*(mat.rho*1000*mat.c_p*dl*dl)/mat.k; % time increment, s
T = .1;
nt = round(T/dt);
T = dt*nt;

dT = zeros(n, 2*n, 100, 4);

total_time = tic;
parfor i = 1:4
    dT(:, :, :, i) = heating2(n, L, theta, dt, nt, mat, beam, beampos(i));
end
toc(total_time)

delete(gcp)