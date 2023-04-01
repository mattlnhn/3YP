clc; clear; close all;

L = 1e-2;
dl = 2.5e-5;
Lf = 1e-4;
dlf = 1e-6;
m = round(dl/dlf);

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
beam.beta_star = 0.15;
beam.epsilon_n = 2.50e-6;
beam.sigma = sqrt(beam.beta_star*beam.epsilon_n/beam.gamma);
beam.pos = -3*beam.sigma; % beam centre pos w.r.t. midpoint of left edge, m

tau = .25;
dt = tau*(mat.rho*1000*mat.c_p*dlf*dlf)/mat.k; % time increment, s
T = 1e-2;
nt = round(T/dt);
T = dt*nt;

fprintf('dt = %d s with %d time steps for total simulation time %d s\nEnter to continue...\n', dt, nt, T)
pause()

total_time = tic;
[dT1, dT2, dT3, dT4] = heating2(dlf, Lf, dl, L, theta, dt, nt, 100, mat, beam);
toc(total_time)

dT2up = kron(dT2, ones(m));
dT3up = kron(dT3, ones(m));
dT4up = kron(dT4, ones(m));

dTL = [dT2up; dT1; dT3up];
dTfinal = [dTL dT4up];

figure()
imagesc(dTfinal)
colorbar
axis square